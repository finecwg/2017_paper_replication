#!/usr/bin/env julia
# scripts/wave_evolution.jl
# Axisymmetric free-wave evolution: GIF (2D disc + 1D section) and r–t heatmap

using BouncingDroplet, LinearAlgebra, SparseArrays, SpecialFunctions
begin
    # plotting deps (자동 설치)
    try
        @eval using Plots
    catch
        import Pkg; Pkg.add("Plots"); @eval using Plots
    end
    try
        @eval using FFMPEG   # GIF/MP4 backend
    catch
        import Pkg; Pkg.add("FFMPEG"); @eval using FFMPEG
    end
end

# --- spectral D2N loader (로컬 파일) ---
include(joinpath(@__DIR__, "..", "src", "operators", "d2n_spectral.jl"))
using .D2N_Spectral

# ---------- 1) 문제 세팅 ----------
nr, Rset = 256, 1.2
mesh = BouncingDroplet.RadialMesh(nr, Rset)
ΔH   = BouncingDroplet.build_laplacian_operator(mesh)
op   = D2N_Spectral.build_D2N_spectral(ΔH, mesh)
N    = sparse(D2N_Spectral.as_matrix(op))      # CSC (time_step에서 기대)
params = BouncingDroplet.DimensionlessNumbers(Fr=1.0, We=5.0, Re=1000.0, M=1.0)

# 초기자유파: J0(kr)
k  = 8.0
η0 = [besselj(0, k*r) for r in mesh.r]
state = BouncingDroplet.SimulationState(η0, zeros(nr), zeros(nr), 0.0, 0.0, 0, 0.0)

# 외력 = 0 (시각화용 안정)
forcing(t) = 0.0

# 매우 보수적 Δt (explicit 안정성 고려)
λmax = maximum(op.sqrtλ)        # ~2.7e5
δt   = 0.01 / λmax              # 더 보수적으로 0.01/λmax
Tfinal = 1.2e-4                 # 짧은 구간만 애니메이션
nsteps = Int(ceil(Tfinal/δt))

# ---------- 2) r→(x,y) 보간 유틸 (mesh.r 만 사용) ----------
function radial_to_grid(η::AbstractVector, mesh; nx::Int=220)
    rvec = mesh.r
    nr   = length(rvec)
    R    = rvec[end]
    dr   = rvec[2] - rvec[1]          # 균일격자 가정 (프로젝트 구조상 OK)

    xs = range(-R, R, length=nx)
    ys = range(-R, R, length=nx)
    Z  = fill(NaN, nx, nx)

    @inbounds for j in 1:nx, i in 1:nx
        r = hypot(xs[i], ys[j])
        if r <= R
            # 선형보간 (경계 처리)
            # r ∈ [0,R], idx ∈ [1, nr-1]
            idxf = clamp(Int(floor(r/dr)) + 1, 1, nr-1)
            r0, r1 = rvec[idxf], rvec[idxf+1]
            α = (r - r0) / (r1 - r0 + eps())
            Z[j,i] = (1-α)*η[idxf] + α*η[idxf+1]
        end
    end
    return xs, ys, Z
end

# ---------- 3) 애니메이션 ----------
frames = 120
steps_per_frame = max(1, div(nsteps, frames))
@info "evolution params" δt nsteps frames steps_per_frame Tfinal

anim = @animate for f in 1:frames
    # 여러 스텝 진전
    for _ in 1:steps_per_frame
        BouncingDroplet.time_step_no_contact!(state, mesh, ΔH, N, params, δt; forcing=forcing)
    end

    # 2D 디스크 컬러맵 + 1D 단면 동시 표시
    xs, ys, Z = radial_to_grid(state.η, mesh; nx=220)

    # NaN-safe amplitude for symmetric color limits
    vals = vec(Z[.!isnan.(Z)])
    amp  = isempty(vals) ? 1.0 : maximum(abs.(vals))
    amp  = amp == 0 ? 1e-6 : amp

    p1 = heatmap(xs, ys, Z; aspect_ratio=1, framestyle=:box, c=:balance,
                 clims = (-amp, amp),
                 xlabel="x", ylabel="y",
                 title="η(x,y,t) at t=$(round(state.t, sigdigits=3))")
    p2 = plot(mesh.r, state.η; xlabel="r", ylabel="η(r,t)", legend=false,
              title="radial section", lw=2)
    plot(p1, p2; layout=(1,2), size=(1100,480), margin=3Plots.mm)
end

gifpath = joinpath(@__DIR__, "wave_evolution.gif")
gif(anim, gifpath, fps=24)

# ---------- 4) 공간-시간(colormap) PNG ----------
# 프레임마다 η 저장 → heatmap(r, t, η)
save_every = steps_per_frame
nshots = fld(nsteps, save_every)     # == nsteps ÷ save_every (정수 나눗셈)

H  = zeros(nr, nshots)
ts = zeros(nshots)

# 초기 상태 되돌리고 기록 루프 (외력=0 유지)
state = BouncingDroplet.SimulationState(η0, zeros(nr), zeros(nr), 0.0, 0.0, 0, 0.0)

for s in 1:nsteps
    BouncingDroplet.time_step_no_contact!(state, mesh, ΔH, N, params, δt; forcing=forcing)
    if s % save_every == 0
        j = s ÷ save_every           # 1,2,3,...,nshots
        @inbounds H[:, j] = state.η
        ts[j] = state.t
    end
end

p_rt = heatmap(ts, mesh.r, H;
    xlabel="time t", ylabel="radius r",
    c=:balance, framestyle=:box, size=(900,600),
    title="Space–time of η(r,t) (free wave)"
)
pngpath = joinpath(@__DIR__, "spacetime_eta.png")
savefig(p_rt, pngpath)

println("\nSaved:")
println("  GIF : $gifpath")
println("  PNG : $pngpath")