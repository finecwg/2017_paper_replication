module D2N_Spectral

using LinearAlgebra, SparseArrays

export SpectralD2N, build_D2N_spectral, apply_D2N, as_matrix, mass_matrix

"축대칭 질량가중(∫₀^R • r dr)용 이산 질량행렬"
function mass_matrix(mesh)
    r, dr, R = mesh.r, mesh.δr, mesh.R
    nr = length(r)
    w = Vector{Float64}(undef, nr)
    @inbounds for i in 1:nr
        rL = max(0.0, r[i] - 0.5*dr)
        rR = min(R,   r[i] + 0.5*dr)
        w[i] = 0.5*(rR*rR - rL*rL)   # ∫ r dr over control volume
    end
    return Diagonal(w)
end

"스펙트럴 D2N 운영자: N = - V diag(√λ) Vᵀ M  (일반화 고유분해 기반)"
struct SpectralD2N
    V::Matrix{Float64}                         # nr×q, M-직교 고유벡터
    sqrtλ::Vector{Float64}                     # 길이 q, √λ
    M::Diagonal{Float64,Vector{Float64}}       # nr×nr, 질량(가중)
end

"""
build_D2N_spectral(ΔH, mesh)

- ΔH : 반경 라플라시안 이산 연산자 (Dirichlet at r=R, r=0 정칙성)
- 일반화 고유분해:  A v = λ M v,  (A = -ΔH)
- 결과: Nφ = - V diag(√λ) (Vᵀ M φ)
"""
function build_D2N_spectral(ΔH::SparseMatrixCSC, mesh)
    M  = mass_matrix(mesh)                 # 이제 SPD
    A  = -Matrix(ΔH)

    # 대칭화
    Minv = Diagonal(1.0 ./ diag(M))
    A_sym = 0.5 * (A + Minv * (A' * M))

    # 혹시 모를 가장자리 수치불안정 대비: M에 극미 릿지 자동 주입 시도
    Mmat = Matrix(M)
    ridge = 0.0
    local F
    try
        F = eigen(Symmetric(A_sym), Symmetric(Mmat))
    catch
        ridge = 1e-14 * maximum(diag(M))
        F = eigen(Symmetric(A_sym), Symmetric(Mmat + ridge*I))
    end

    λ = F.values
    V = Matrix(F.vectors)                  # V' M V ≈ I

    # 양의 고유값만 채택
    keep = λ .> 1e-12
    λ = λ[keep]; V = V[:, keep]

    # 엄밀 정규화
    for j in 1:size(V,2)
        nj = sqrt(dot(V[:,j], M * V[:,j]))
        V[:,j] ./= nj
    end

    return SpectralD2N(V, sqrt.(λ), M)
end

"적용: Nφ = - V ( √λ ∘ (Vᵀ M φ) )"
function apply_D2N(op::SpectralD2N, φ::AbstractVector{<:Real})
    c = op.V' * (op.M * φ)
    y = op.V * (op.sqrtλ .* c)
    return -y
end

"행렬이 필요할 때만 재료화: N = - V diag(√λ) Vᵀ M"
function as_matrix(op::SpectralD2N)
    return - op.V * Diagonal(op.sqrtλ) * (op.V' * op.M)
end

end # module