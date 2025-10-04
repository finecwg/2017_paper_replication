"""
Spherical solid shape
"""
struct SphericalSolid
    R₀::Float64         # Radius
    
    SphericalSolid(R₀::Float64) = new(R₀)
end

"""
Height of sphere surface z_s(r) 
For a sphere centered at origin: z_s(r) = -√(R₀² - r²)
z_s(0) = -R₀ (bottom of sphere)
"""
function surface_height(sphere::SphericalSolid, r::Float64)
    if r >= sphere.R₀
        return NaN  # Outside sphere
    end
    return -sqrt(sphere.R₀^2 - r^2)
end

"""
Derivative dz_s/dr
"""
function surface_slope(sphere::SphericalSolid, r::Float64)
    if r >= sphere.R₀
        return NaN
    end
    if r == 0.0
        return 0.0
    end
    return r / sqrt(sphere.R₀^2 - r^2)
end

"""
Curvature κ = 2/R₀ (constant for sphere)
"""
function surface_curvature(sphere::SphericalSolid, r::Float64)
    return 2.0 / sphere.R₀
end

"""
Evaluate z_s on radial mesh
"""
function evaluate_on_mesh(sphere::SphericalSolid, mesh::RadialMesh)
    z_s = zeros(mesh.nr)
    for i in 1:mesh.nr
        r = mesh.r[i]
        z_s[i] = r < sphere.R₀ ? surface_height(sphere, r) : NaN
    end
    return z_s
end