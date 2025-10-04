using Interpolations

"""
Cubic interpolation on radial mesh
"""
function cubic_interpolate_radial(mesh::RadialMesh, values::Vector{Float64}, r_query::Float64)
    if r_query < 0 || r_query > mesh.R
        return 0.0  # Outside domain
    end
    
    # Create interpolation object
    itp = interpolate(values, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, mesh.r)
    
    return sitp(r_query)
end

"""
Cubic interpolation for multiple query points
"""
function cubic_interpolate_radial(mesh::RadialMesh, values::Vector{Float64}, r_queries::Vector{Float64})
    itp = interpolate(values, BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, mesh.r)
    
    return [r >= 0 && r <= mesh.R ? sitp(r) : 0.0 for r in r_queries]
end