module DirichletNeumann

using LinearAlgebra
using SpecialFunctions

export DirichletToNeumannOperator, build_D2N, apply_D2N, as_matrix

"""
Bessel-basis Galerkin D2N (axisymmetric, φ(R)=0).

기저   : J0(ζ_n r / R), ζ_n = n번째 J0 영점
계수   : a = G^{-1} Bᵀ W φ,  (G = Bᵀ W B : 이산 Gram)
사상   : Nφ = - B ∘(k)∘ a,   (k_n = ζ_n / R, ∘(k)∘ = Diag(k))
여기서
- B[i,n] = J0(ζ_n r_i / R)
- W = diag(w_i), w_i = r_i Δr (끝점 1/2; 축대칭 정규화에 맞춤)

※ 구현 포인트:
- G를 팩터(F)로 저장하고, apply할 때마다 a를 푼다(F rhs).
- 필요 시에만 아주 작은 ridge λI를 자동 주입해 수치안정 보장.
"""
struct DirichletToNeumannOperator
    # 저장: 행렬 N 자체가 아니라, 적용에 필요한 요소들
    B::Matrix{Float64}                        # nr × Nm
    W::Diagonal{Float64,Vector{Float64}}      # nr × nr
    k::Vector{Float64}                        # Nm
    F::Factorization{Float64}                 # G(=B'WB)의 안정화된 팩터(cholesky 또는 ldlt)
    zeta::Vector{Float64}                     # Nm
    R::Float64
    r::Vector{Float64}
end

"J0 영점: McMahon 초기값 + Newton 보정(J0'=-J1)."
function _j0_zeros(N::Int)
    z = Vector{Float64}(undef, N)
    @inbounds for n in 1:N
        x = (n - 0.25)*π
        for _ in 1:20
            f  = besselj(0, x)
            fp = -besselj(1, x)
            x -= f / fp
        end
        z[n] = x
    end
    return z
end

"Gram이 SPD처럼 보이지 않으면 작은 λI로 ridge를 주입해 factorization."
function _factor_with_ridge(G::Symmetric{Float64,Matrix{Float64}})
    s = maximum(abs, diag(G))
    λ = (s == 0.0) ? 1e-14 : 1e-14*s
    for _ in 1:6
        try
            return cholesky(G + λ*I)
        catch
            λ *= 100.0
        end
    end
    # 마지막 안전장치
    return ldlt(G + λ*I)
end

"메쉬(균일 r, 끝점 포함)에서 D2N 운영자 구성."
function build_D2N(mesh)
    nr, R = mesh.nr, mesh.R
    r, dr = mesh.r, mesh.δr

    Nm = nr - 1                   # r=R 노드는 J0(ζ_n)=0
    ζ  = _j0_zeros(Nm)
    k  = ζ ./ R

    # B: J0(ζ_n r_i / R)
    B = Matrix{Float64}(undef, nr, Nm)
    @inbounds for i in 1:nr, n in 1:Nm
        B[i,n] = besselj(0, ζ[n]*r[i]/R)
    end

    # 축대칭 가중: r_i Δr (끝점 1/2)
    w = r .* dr
    w[1]  *= 0.5
    w[end]*= 0.5
    W = Diagonal(w)

    # Gram & factor
    G = Symmetric(B' * (W * B))
    F = _factor_with_ridge(G)

    return DirichletToNeumannOperator(B, W, k, F, ζ, R, r)
end

"적용: Nφ = -B * ( (k .* (F  (B' * (W * φ)))) )"
function apply_D2N(op::DirichletToNeumannOperator, φ::AbstractVector{<:Real})
    @assert length(φ) == size(op.B,1) "apply_D2N: length(φ) must be nr"
    rhs = op.B' * (op.W * φ)                 # Nm
    a   = op.F \ rhs                          # Nm  (G^{-1} * rhs)
    return - op.B * (op.k .* a)               # nr
end

"행렬이 꼭 필요할 때만 재료화(밀집 행렬)."
function as_matrix(op::DirichletToNeumannOperator)
    # X = G^{-1} Bᵀ W  (Nm×nr)
    X = op.F \ (op.B' * op.W)
    # N = - B Diag(k) X  (nr×nr)
    return - op.B * (Diagonal(op.k) * X)
end

end # module