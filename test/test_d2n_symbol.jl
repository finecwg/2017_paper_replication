using BouncingDroplet
using LinearAlgebra
using SpecialFunctions

# 메쉬 (위치 인자!)
mesh = BouncingDroplet.RadialMesh(256, 1.0)

# D2N 만들기
Nop = BouncingDroplet.DirichletNeumann.build_D2N(mesh)

# 시험 함수: φ(r) = J0(k r) 이면 D2N(φ) ≈ -k φ
k  = 12.3
φ  = [besselj(0, k*r) for r in mesh.r]
Nφ = BouncingDroplet.DirichletNeumann.apply_D2N(Nop, φ)

err = norm(Nφ .+ k .* φ) / norm(φ)
@show err
@assert err < 1e-2
println("D2N symbol test passed.")