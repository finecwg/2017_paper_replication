"""
Compute curvature κ[η] for radial function
For radial functions: κ = (ηrr(1 + ηr²) + ηr/r) / (1 + ηr²)^(3/2)
Linear approximation: κ ≈ Δ_H η = ηrr + ηr/r
"""
function compute_curvature_linear(mesh::RadialMesh, η::Vector{Float64})
    nr = mesh.nr
    δr = mesh.δr
    κ = zeros(nr)
    
    # At origin (special case)
    κ[1] = 4 * (η[2] - η[1]) / δr^2
    
    # Interior points
    for i in 2:(nr-1)
        r = mesh.r[i]
        η_rr = (η[i+1] - 2*η[i] + η[i-1]) / δr^2
        η_r = (η[i+1] - η[i-1]) / (2*δr)
        
        κ[i] = η_rr + η_r / r
    end
    
    # Boundary
    κ[nr] = 0.0
    
    return κ
end

"""
Compute full nonlinear curvature
"""
function compute_curvature_nonlinear(mesh::RadialMesh, η::Vector{Float64})
    nr = mesh.nr
    δr = mesh.δr
    κ = zeros(nr)
    
    # At origin
    κ[1] = 4 * (η[2] - η[1]) / δr^2
    
    # Interior points
    for i in 2:(nr-1)
        r = mesh.r[i]
        η_rr = (η[i+1] - 2*η[i] + η[i-1]) / δr^2
        η_r = (η[i+1] - η[i-1]) / (2*δr)
        
        numerator = η_rr * (1 + η_r^2) + η_r / r
        denominator = (1 + η_r^2)^1.5
        
        κ[i] = numerator / denominator
    end
    
    κ[nr] = 0.0
    
    return κ
end