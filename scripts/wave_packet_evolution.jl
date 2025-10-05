#!/usr/bin/env julia
using BouncingDroplet, LinearAlgebra, SparseArrays
using Plots, Printf
gr(); default(size=(1200,520), dpi=140)

# --- D2N (로컬 include가 필요하면 주석 해제)
include(joinpath(@__DIR__, "..", "src", "operators", "d2n_spectral.jl"))
using .D2N_Spectral

# ---------- 유틸 ----------
project_coeffs(η, V, M; G=nothing) = (G === nothing ? (V'*(M*V)) : G) \ (V'*(M*η))

"""
밴드제한 초기조건 (고유모드 m0±bw만 사용).
η0 = V*c0,  φ0 = -β * η0  (외향 진행 감각을 주는 간단한 위상)
"""
function bandlimited_bump(op; m0=24, bw=8, amp=1.0, β=0.15)
    V = op.V
    m = size(V,2)
    idx = 2:m                                     # m=1(상수) 제외
    w   = exp.(-0.5*((idx .- m0)./bw).^2)         # 모달 가우시안
    c0  = zeros(m); c0[idx] .= amp .* w
    η0  = V * c0
    φ0  = -β .* η0
    return η0, φ0
end

"""
r=0 근방 2차 보간(η ≈ a + b r^2) + 그 외 구간 선형보간, 컬러 스케일 고정
"""
function radial_to_grid(η::AbstractVector, mesh::BouncingDroplet.RadialMesh; nx=251)
    R = mesh.R
    x  = range(-R, R; length=nx);  y = x
    Z  = zeros(nx, nx)

    # 원점 2차 보간 파라미터 (η ≈ a + b r^2)
    η1, η2, η3 = η[1], η[2], η[3]
    r1, r2, r3 = mesh.r[1], mesh.r[2], mesh.r[3]
    b = (η3-η2)/(r3^2 - r2^2)
    a = η2 - b*r2^2

    @inbounds for j in eachindex(y), i in eachindex(x)
        r = hypot(x[i], y[j])
        if r <= r3
            Z[j,i] = a + b*r^2
        else
            k = searchsortedlast(mesh.r, r)
            k = clamp(k, 2, mesh.nr-1)
            rL, rR = mesh.r[k], mesh.r[k+1]
            t = (r - rL)/(rR - rL)
            Z[j,i] = (1-t)*η[k] + t*η[k+1]
        end
    end
    return x, y, Z
end

"""
아주 약한 모달 지수 필터 (선택). 고주파 누적 방지, 파형 거의 보존.
"""
function exp_filter!(f, V, M; α=6.0, p=8, G=nothing)
    G = (G === nothing) ? (V'*(M*V)) : G
    c = G \ (V'*(M*f))
    m = length(c)
    σ = @. exp(-α * ( (0:m-1)/(m-1) )^p)
    c .*= σ
    f .= V * c
    return nothing
end

# ---------- 설정 ----------
mesh = BouncingDroplet.RadialMesh(256, 1.2)
ΔH   = BouncingDroplet.build_laplacian_operator(mesh)
op   = D2N_Spectral.build_D2N_spectral(ΔH, mesh)
N    = sparse(D2N_Spectral.as_matrix(op))
M    = D2N_Spectral.mass_matrix(mesh)
V    = op.V
G    = V'*(M*V)  # ~I

params = BouncingDroplet.DimensionlessNumbers(Fr=1.0, We=5.0, Re=500.0, M=1.0)
η0, φ0 = bandlimited_bump(op; m0=24, bw=8, amp=1.0, β=0.15)
state  = BouncingDroplet.SimulationState(η0, φ0, zeros(mesh.nr), 0.0, 0.0, 0, 0.0)

λmax = maximum(op.sqrtλ)
δt   = 0.02/λmax                      # 보수적 안정계수(이전 실험에서 안정)
Tfinal = 1.0e-4
nstep  = ceil(Int, Tfinal/δt)

frames = 160
steps_per_frame = max(1, cld(nstep, frames))

@info "evolution params" δt nstep frames steps_per_frame Tfinal

# ---------- 플로팅 ----------
pltL = heatmap(); pltR = plot()
function snapshot!(plt1, plt2, η, mesh, t; nx=301, clim=(-1.0,1.0))
    x, y, Z = radial_to_grid(η, mesh; nx=nx)
    heatmap!(plt1, x, y, Z; color=:RdBu, clim=clim, ratio=1, xlabel="x", ylabel="y", cbar_title="η")
    title!(plt1, @sprintf("η(x,y,t) at t=%.2e", t))

    plot!(plt2, mesh.r, η; xlabel="r", ylabel="η(r,t)", legend=false, title="radial section")
end

# 공통 컬러스케일을 위해 초기 η 범위로 설정
cl = (-1.0, 1.0)

anim = @animate for k in 1:frames
    # 시간 진보
    for _ in 1:steps_per_frame
        BouncingDroplet.time_step_no_contact!(state, mesh, ΔH, N, params, δt; forcing = t->0.0)
        # 선택적으로 약한 필터
        exp_filter!(state.η, V, M; α=6.0, p=8, G=G)
        exp_filter!(state.φ, V, M; α=6.0, p=8, G=G)
    end

    # 모달 고주파 비율 모니터링 (정보 출력만)
    c = project_coeffs(state.η, V, M; G=G)
    Em = dot(state.η, M*state.η)
    m  = length(c); hi = round(Int, 0.8m):m
    frac_hi = sum(abs2, c[hi]) / Em
    @info "frame $(k)/$frames" t=state.t ηnorm=norm(state.η) frac_hi=frac_hi

    pltL = heatmap(); pltR = plot()
    snapshot!(pltL, pltR, state.η, mesh, state.t; nx=281, clim=cl)
    plot(pltL, pltR; layout=(1,2))
end

gif_path = joinpath(@__DIR__, "wave_packet_traveling.gif")
gif(anim, gif_path, fps=24)
@info "Saved animation" gif_path