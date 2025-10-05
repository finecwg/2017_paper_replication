#!/usr/bin/env julia
# scripts/progress_report.jl
# Progress report: spectral D2N diagnostics + short no-contact integration with plots
#
# 실행:
#   julia --project=. scripts/progress_report.jl

# --- 0) 기본 패키지 로드 (+Plots 없으면 설치) -------------------------------------
using BouncingDroplet
using LinearAlgebra, SparseArrays
using SpecialFunctions
using Random
begin
    try
        @eval using Plots
    catch e
        @warn "Plots not found in this environment. Installing Plots..." exception=(e, catch_backtrace())
        import Pkg
        Pkg.add("Plots")
        @eval using Plots
    end
end
Plots.default(titlefont=font(12), legendfontsize=9)

# --- 1) 스펙트럴 D2N 로더 ---------------------------------------------------------
# 프로젝트 루트 기준으로 scripts/ 아래에서 실행되므로, 상대경로 사용
include(joinpath(@__DIR__, "..", "src", "operators", "d2n_spectral.jl"))
using .D2N_Spectral

# --- 2) 연산자 구성 ---------------------------------------------------------------
nr, R = 256, 1.2
mesh = BouncingDroplet.RadialMesh(nr, R)
ΔH   = BouncingDroplet.build_laplacian_operator(mesh)
op   = D2N_Spectral.build_D2N_spectral(ΔH, mesh)

# as_matrix(op)는 밀집행렬(Matrix)을 반환 → 진단/플롯용
N_mat = D2N_Spectral.as_matrix(op)
# time_step_no_contact! 는 CSC 스파스 기대 → 적분용
N_csc = sparse(N_mat)
# 질량 행렬 (대각)
M = D2N_Spectral.mass_matrix(mesh)

# --- 3) 스펙트럴 진단: 부호/정밀도 & 레일리 몫 ----------------------------------
mset = [2, 4, 8, 16, 32, 64]
rq             = zeros(Float64, length(mset))  # Rayleigh quotient: φᵀ N φ / φᵀ φ
rq_relerr      = similar(rq)
spec_residuals = similar(rq)                   # ||Nφ + √λ φ|| / ||φ||
for (i, m) in enumerate(mset)
    φm = @view op.V[:, m]
    Nφ = N_mat * φm
    rq[i] = dot(φm, Nφ) / dot(φm, φm)
    rq_relerr[i] = abs(rq[i] + op.sqrtλ[m]) / (1 + op.sqrtλ[m])   # 기대: rq ≈ -√λ
    spec_residuals[i] = norm(D2N_Spectral.apply_D2N(op, φm) .+ op.sqrtλ[m] .* φm) / norm(φm)
end

# --- 4) 구조적 성질: M-자기수반 & 음정확성 ---------------------------------------
MN = M * N_mat
sym_residual = norm(MN - MN', Inf) / max(norm(MN, Inf), 1.0)  # 기대 ~1e-12

rng = MersenneTwister(2025)
nsamples = 20
qvals = zeros(Float64, nsamples)
for i in 1:nsamples
    x = randn(rng, mesh.nr)
    qvals[i] = dot(x, M * (N_mat * x))          # 기대: ≤ 0 (수치오차 포함)
end
q_max = maximum(qvals)  # 0 위로 뜨는 부분은 수치오차 수준이어야 함
q_min = minimum(qvals)

# --- 5) 짧은 무접촉 적분(run): 매우 보수적인 Δt -----------------------------------
params = BouncingDroplet.DimensionlessNumbers(Fr=1.0, We=5.0, Re=1000.0, M=1.0)
# 초기 자유파 (J0)
k  = 8.0
η0 = [besselj(0, k*r) for r in mesh.r]
state = BouncingDroplet.SimulationState(η0, zeros(mesh.nr), zeros(mesh.nr), 0.0, 0.0, 0, 0.0)

# 작은 외력: A*sin(ω t)
A, ω = 0.02, 40.0
forcing(t) = A * sin(ω*t)

# 안정성 고려한 보수적 Δt (경험적으로 잘 동작)
λmax = maximum(op.sqrtλ)        # ~ 2.7e5
δt   = 0.02 / λmax              # ~ 7e-8
nstep = 2000                    # 총 적분 단계

ts = Vector{Float64}(undef, nstep)
ηnorms = similar(ts)
φnorms = similar(ts)
for i in 1:nstep
    BouncingDroplet.time_step_no_contact!(state, mesh, ΔH, N_csc, params, δt; forcing=forcing)
    ts[i] = state.t
    ηnorms[i] = norm(state.η)
    φnorms[i] = norm(state.φ)
end

# --- 6) 플로팅 --------------------------------------------------------------------
# (a) 고유모드 검증: spec_residuals (semilogy)
p1 = plot(
    mset, spec_residuals;
    seriestype = :scatter, markersize=6,
    yscale = :log10,
    xlabel = "mode index m",
    ylabel = "||Nφ + √λ φ|| / ||φ||",
    title  = "Eigenmode residual (expect ~1e-12–1e-10)",
    label  = "spectral residual",
)
plot!(p1, mset, rq_relerr; seriestype=:scatter, markershape=:diamond, label="Rayleigh quotient rel. err")

# (b) 구조 성질: M-자기수반 & 음정확성
#    - 좌표 1: symmetry residual (log-scale)
#    - 좌표 2: qvals 분포 (dot)
p2a = bar([1.0], [max(sym_residual, 1e-16)];
    yscale=:log10, xticks=([1.0], ["M*N symmetry residual"]),
    ylabel="value (log10)", label=false, title="Structural checks")
hline!(p2a, [1e-12], l=:dash, label="~1e-12 ref")

p2b = scatter(1:nsamples, qvals; xlabel="sample", ylabel="xᵀ M N x", label="quadratic form", legend=:bottomright)
hline!(p2b, [0.0], l=:dash, label="0")

# (c) 시간 적분: ||η||, ||φ|| vs t
p3 = plot(ts, ηnorms; label="‖η‖₂")
plot!(p3, ts, φnorms; label="‖φ‖₂")
xlabel!(p3, "time t")
ylabel!(p3, "norm")
title!(p3, "No-contact integration (δt = $(round(δt, sigdigits=3)), steps = $nstep)")

plt = plot(p1, plot(p2a, p2b, layout=(2,1), size=(600,600)), p3; layout=(2,2), size=(1200,900))

outpath = joinpath(@__DIR__, "progress_report.png")
savefig(plt, outpath)

# --- 7) 콘솔 요약 -----------------------------------------------------------------
println("\n=== Progress summary ===")
println("Mesh: nr=$(mesh.nr), R=$(mesh.R)")
println("λmax ≈ $(maximum(op.sqrtλ))  → δt = $(δt)  (steps = $nstep, t_final ≈ $(ts[end]))")
println("M-symmetry residual (∞-norm, relative): ", sym_residual)
println("Quadratic form xᵀ M N x over $nsamples samples: min=$(q_min), max=$(q_max)")
println("Eigenmode checks for m ∈ $(mset):")
for (i, m) in enumerate(mset)
    println("  m=$(lpad(m, 3)): residual=$(spec_residuals[i]):0.3e,  rq_relerr=$(rq_relerr[i]):0.3e")
end
println("Time series (last):  ‖η‖₂=$(ηnorms[end]),  ‖φ‖₂=$(φnorms[end])")
println("Figure saved to: $outpath")