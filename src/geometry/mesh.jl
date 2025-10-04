"""
Radial mesh for axisymmetric problem
Fields:
- nr: number of radial points (including r=0 and r=R)
- δr: radial spacing
- R: domain radius
- r: radial coordinates [0, δr, 2δr, ..., R]
"""
struct RadialMesh
    nr::Int
    δr::Float64
    R::Float64
    r::Vector{Float64}
    
    function RadialMesh(nr::Int, R::Float64)
        δr = R / (nr - 1)
        r = collect(0:δr:R)
        new(nr, δr, R, r)
    end
end

Base.length(mesh::RadialMesh) = mesh.nr
Base.getindex(mesh::RadialMesh, i::Int) = mesh.r[i]

"""
Get index for a given radius (nearest point)
"""
function radius_to_index(mesh::RadialMesh, r_val::Float64)
    return max(1, min(mesh.nr, round(Int, r_val / mesh.δr) + 1))
end

"""
Check if index is valid
"""
function is_valid_index(mesh::RadialMesh, i::Int)
    return 1 <= i <= mesh.nr
end

"""
Create mesh with minimum points per radius
"""
function create_mesh_for_sphere(R₀::Float64, λ_F::Float64; 
                               points_per_radius::Int=40,
                               domain_size_factor::Float64=10.0)
    R_domain = domain_size_factor * λ_F
    δr = R₀ / points_per_radius
    nr = floor(Int, R_domain / δr) + 1
    
    return RadialMesh(nr, R_domain)
end