"""
Build discrete Laplacian operator matrix for radial functions
Based on Appendix B, equations (B.3) and (B.4)
"""
function build_laplacian_operator(mesh::RadialMesh)
    nr = mesh.nr
    δr = mesh.δr
    
    # Sparse matrix
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    # First point (i=0): special 5-point stencil (eq B.3)
    # ΔHφ(0) ≈ 4[φ(δr) - φ(0)] / δr²
    push!(rows, 1); push!(cols, 1); push!(vals, -4.0/δr^2)
    push!(rows, 1); push!(cols, 2); push!(vals, 4.0/δr^2)
    
    # Interior points (i=1 to nr-2): equation (B.4)
    # ΔHφ(iδr) ≈ [φ((i+1)δr) - 2φ(iδr) + φ((i-1)δr)]/δr² 
    #           + [φ((i+1)δr) - φ((i-1)δr)]/(2iδr²)
    for i in 2:(nr-1)
        idx = i  # 1-indexed
        
        # Coefficient for φ((i-1)δr)
        coeff_minus = 1.0/δr^2 - 1.0/(2*(i-1)*δr^2)
        # Coefficient for φ(iδr)
        coeff_center = -2.0/δr^2
        # Coefficient for φ((i+1)δr)
        coeff_plus = 1.0/δr^2 + 1.0/(2*(i-1)*δr^2)
        
        push!(rows, idx); push!(cols, idx-1); push!(vals, coeff_minus)
        push!(rows, idx); push!(cols, idx);   push!(vals, coeff_center)
        push!(rows, idx); push!(cols, idx+1); push!(vals, coeff_plus)
    end
    
    # Last point (i=nr-1): boundary condition φ(R) = 0
    push!(rows, nr); push!(cols, nr); push!(vals, 1.0)
    
    return sparse(rows, cols, vals, nr, nr)
end

"""
Apply Laplacian to a radial function
"""
function apply_laplacian(ΔH_matrix::SparseMatrixCSC, φ::Vector{Float64})
    return ΔH_matrix * φ
end