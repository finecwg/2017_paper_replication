using BouncingDroplet, LinearAlgebra

# 메쉬·라플라시안
mesh = BouncingDroplet.RadialMesh(256, 1.2)
ΔH   = BouncingDroplet.build_laplacian_operator(mesh)

# 스펙트럴 D2N
include(joinpath(@__DIR__, "..", "src", "operators", "d2n_spectral.jl"))
import .D2N_Spectral: build_D2N_spectral, apply_D2N, as_matrix, mass_matrix

op = build_D2N_spectral(ΔH, mesh)

# 고유벡터 하나 꺼내 검증: op 내부 고유기저와 완벽 일치
# (테스트를 위해 행렬 재료화 없이 직접 비교)
V = op.V; M = op.M; sqrtλ = op.sqrtλ
m = 8
φ = @view V[:, m]              # 이산 고유벡터
Nφ = apply_D2N(op, φ)

err = norm(Nφ .+ sqrtλ[m] .* φ) / norm(φ)
@show err
@assert err < 1e-10
println("D2N spectral mode test passed.")