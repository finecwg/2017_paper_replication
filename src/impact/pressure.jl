"""
Compute pressure distribution on contact area
From momentum equation: p = -φ_t - η (Bernoulli)
"""
function compute_contact_pressure(mesh::RadialMesh,
                                 φ::Vector{Float64},
                                 φ_prev::Vector{Float64},
                                 η::Vector{Float64},
                                 params::DimensionlessNumbers,
                                 δt::Float64,
                                 k::Int)
    
    nr = mesh.nr
    ps = zeros(nr)
    
    if k == 0
        return ps
    end
    
    # Time derivative of φ
    φ_t = (φ .- φ_prev) ./ δt
    
    # Pressure from Bernoulli (dimensionless)
    # p = -φ_t - (1/Fr)η + viscous terms
    for i in 1:k
        ps[i] = -φ_t[i] - η[i]/params.Fr
    end
    
    return ps
end

"""
Compute total force on sphere from pressure
F = ∫ p dA = 2π ∫₀^{r_c} p(r) r dr
"""
function compute_pressure_force(mesh::RadialMesh, ps::Vector{Float64}, k::Int)
    if k == 0
        return 0.0
    end
    
    # Trapezoidal integration
    F = 0.0
    for i in 1:k-1
        r_i = mesh.r[i]
        r_ip1 = mesh.r[i+1]
        
        # Average pressure and radius
        p_avg = (ps[i] + ps[i+1]) / 2
        r_avg = (r_i + r_ip1) / 2
        
        # Area element: 2π r dr
        dA = 2π * r_avg * mesh.δr
        
        F += p_avg * dA
    end
    
    return F
end