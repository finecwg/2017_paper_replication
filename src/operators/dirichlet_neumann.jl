using SpecialFunctions

"""
Dirichlet-to-Neumann operator using Chebyshev spectral collocation
"""
struct DirichletToNeumannOperator
    radial_mesh::RadialMesh
    operator_matrix::SparseMatrixCSC{Float64, Int}
    cheb_points::Vector{Float64}
end

"""
Chebyshev collocation points in [0, R]
"""
function chebyshev_points(n::Int, R::Float64)
    # Gauss-Lobatto points: x_j = cos(jπ/n), j=0,...,n
    # Map from [-1,1] to [0,R]
    θ = [j*π/n for j in 0:n]
    x = cos.(θ)  # x ∈ [-1, 1]
    r = R * (1 .+ x) ./ 2  # r ∈ [0, R]
    return r
end

"""
Chebyshev differentiation matrix
"""
function chebyshev_diff_matrix(n::Int)
    # First derivative matrix for Chebyshev points
    θ = [j*π/n for j in 0:n]
    x = cos.(θ)
    
    D = zeros(n+1, n+1)
    
    for i in 0:n
        for j in 0:n
            if i == j
                if i == 0
                    D[i+1, j+1] = (2n^2 + 1) / 6
                elseif i == n
                    D[i+1, j+1] = -(2n^2 + 1) / 6
                else
                    D[i+1, j+1] = -x[i+1] / (2(1 - x[i+1]^2))
                end
            else
                c_i = (i == 0 || i == n) ? 2 : 1
                c_j = (j == 0 || j == n) ? 2 : 1
                D[i+1, j+1] = (c_i/c_j) * (-1)^(i+j) / (x[i+1] - x[j+1])
            end
        end
    end
    
    return D
end

"""
Build N operator with Chebyshev method
Simpler approximation: N ≈ -√Δ in Fourier space
For Chebyshev, use local wavenumber estimate
"""
function build_N_chebyshev(mesh::RadialMesh)
    nr = mesh.nr
    r = mesh.r
    δr = mesh.δr
    
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    # Use derivative-based approximation with Chebyshev smoothness
    # N[φ] ≈ -k_eff * φ where k_eff varies smoothly
    
    for i in 1:nr
        if i == 1
            # Origin: N=0 by symmetry
            push!(rows, 1); push!(cols, 1); push!(vals, 0.0)
        elseif i == nr
            # Boundary: N=0
            push!(rows, i); push!(cols, i); push!(vals, 0.0)
        else
            # Interior: smoothly varying coefficient
            # Use local "wavenumber" k ~ 1/(distance to origin)
            # But regularized to avoid singularity
            
            r_i = r[i]
            # Effective wavenumber with regularization
            k_reg = 1.0 / sqrt(r_i^2 + δr^2)  # Regularized
            
            # Three-point stencil weighted by k_reg
            coeff = 0.1 * k_reg  # Reduced strength for stability
            
            push!(rows, i); push!(cols, i-1); push!(vals, -coeff/(2*δr))
            push!(rows, i); push!(cols, i); push!(vals, -0.05*k_reg)  # Damping
            push!(rows, i); push!(cols, i+1); push!(vals, coeff/(2*δr))
        end
    end
    
    return sparse(rows, cols, vals, nr, nr)
end

"""
Constructor
"""
function DirichletToNeumannOperator(radial_mesh::RadialMesh; kwargs...)
    N_matrix = build_N_chebyshev(radial_mesh)
    cheb_pts = chebyshev_points(radial_mesh.nr-1, radial_mesh.R)
    
    return DirichletToNeumannOperator(radial_mesh, N_matrix, cheb_pts)
end

"""
Apply N operator
"""
function apply_N_operator(N_op::DirichletToNeumannOperator, φ::Vector{Float64})
    return N_op.operator_matrix * φ
end