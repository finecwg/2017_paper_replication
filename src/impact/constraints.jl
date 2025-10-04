"""
Check if sphere penetrates free surface
"""
function check_penetration(mesh::RadialMesh, η::Vector{Float64}, 
                          h::Float64, z_s::Vector{Float64})
    gap = h .+ z_s .- η
    valid_gap = filter(!isnan, gap)
    
    if isempty(valid_gap)
        return (penetration = false, min_gap = NaN, location = 0.0)
    end
    
    min_gap = minimum(valid_gap)
    min_idx = 1
    
    for i in 1:length(gap)
        if !isnan(gap[i]) && gap[i] ≈ min_gap
            min_idx = i
            break
        end
    end
    
    return (
        penetration = min_gap < 0,
        min_gap = min_gap,
        location = mesh.r[min_idx]
    )
end

"""
Check tangency condition at contact radius
"""
function check_tangency(mesh::RadialMesh, η::Vector{Float64}, 
                       z_s::Vector{Float64}, k::Int)
    if k <= 1 || k >= mesh.nr
        return Inf
    end
    
    δr = mesh.δr
    η_r = (η[k+1] - η[k-1]) / (2*δr)
    zs_r = (z_s[k+1] - z_s[k-1]) / (2*δr)
    
    return abs(η_r - zs_r)
end

"""
Find optimal contact radius k
"""
function find_contact_radius(mesh::RadialMesh, η::Vector{Float64},
                            h::Float64, z_s::Vector{Float64};
                            k_init::Int=10, max_iter::Int=20)
    
    pen = check_penetration(mesh, η, h, z_s)
    if !pen.penetration
        return 0
    end
    
    k_max = findfirst(isnan.(z_s))
    if k_max === nothing
        k_max = mesh.nr - 1
    else
        k_max = max(2, k_max - 2)
    end
    
    k_min = 2
    
    if k_max <= k_min
        return 0
    end
    
    k_best = min(k_init, k_max)
    error_best = Inf
    
    for k in k_min:max(1, div(k_max-k_min, 10)):k_max
        error = check_tangency(mesh, η, z_s, k)
        if error < error_best
            error_best = error
            k_best = k
        end
    end
    
    for k in max(k_min, k_best-3):min(k_max, k_best+3)
        error = check_tangency(mesh, η, z_s, k)
        if error < error_best
            error_best = error
            k_best = k
        end
    end
    
    return k_best
end

"""
Check if contact configuration is valid
"""
function is_valid_contact(mesh::RadialMesh, η::Vector{Float64},
                         h::Float64, z_s::Vector{Float64}, k::Int)
    if k == 0
        return true
    end
    
    if k <= 1 || k >= mesh.nr
        return false
    end
    
    for i in 1:k
        if isnan(z_s[i])
            return false
        end
        
        diff = abs(η[i] - (h + z_s[i]))
        if diff > 0.01 * mesh.δr
            return false
        end
    end
    
    tang_error = check_tangency(mesh, η, z_s, k)
    if tang_error > 1.0
        return false
    end
    
    return true
end