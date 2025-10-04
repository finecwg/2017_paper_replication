"""
Polar mesh centered at r_center for evaluating N operator
Fields:
- r_center: center point on radial mesh
- n_r: number of radial rings
- n_θ: number of angular points
- Δr: radial spacing (finer than radial mesh)
- Δθ: angular spacing
- points: (r, θ) pairs
"""
struct PolarMesh
    r_center::Float64
    n_r::Int
    n_θ::Int
    Δr::Float64
    Δθ::Float64
    points::Matrix{SVector{2,Float64}}  # [radial_idx, angular_idx] -> (r, θ)
    
    function PolarMesh(r_center::Float64, δr_radial::Float64, R_max::Float64;
                      refinement::Int=10, n_θ::Int=128)
        Δr = δr_radial / refinement
        Δθ = 2π / n_θ
        
        n_r = ceil(Int, R_max / Δr)
        
        points = Matrix{SVector{2,Float64}}(undef, n_r, n_θ)
        for i in 1:n_r
            r = i * Δr
            for j in 1:n_θ
                θ = (j - 1) * Δθ
                points[i, j] = SVector{2,Float64}(r, θ)
            end
        end
        
        new(r_center, n_r, n_θ, Δr, Δθ, points)
    end
end

"""
Convert polar (r, θ) to Cartesian (x, y) relative to mesh center
"""
function polar_to_cartesian(pmesh::PolarMesh, r::Float64, θ::Float64)
    x = pmesh.r_center + r * cos(θ)
    y = r * sin(θ)
    return (x, y)
end

"""
Get distance in xy plane
"""
function get_distance(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end