"""
Solve one contact time step
"""
function solve_contact_step(mesh::RadialMesh,
                            sphere::SphericalSolid,
                            ΔH::SparseMatrixCSC,
                            N_matrix::SparseMatrixCSC,
                            params::DimensionlessNumbers,
                            state::SimulationState,
                            δt::Float64)
    
    nr = mesh.nr
    z_s = evaluate_on_mesh(sphere, mesh)
    
    # Find contact radius
    k = find_contact_radius(mesh, state.η, state.h, z_s, k_init=state.k_contact)
    
    if k == 0
        # No contact - update both fluid AND sphere
        time_step_no_contact!(state, mesh, ΔH, N_matrix, params, δt)
        
        # Sphere falls under gravity
        F_gravity = -1.0 / params.Fr
        state.ht += δt * F_gravity / params.M
        state.h += δt * state.ht
        state.k_contact = 0
        
        return state
    end
    
    # In contact
    η_new = copy(state.η)
    φ_new = copy(state.φ)
    ps_new = zeros(nr)
    
    # Apply contact constraint
    for i in 1:k
        η_new[i] = state.h + z_s[i]
    end
    
    # Update free surface
    for i in k+1:nr
        κ = compute_curvature_linear(mesh, state.η)
        η_new[i] = state.η[i] + δt * (N_matrix[i,:]' * state.φ)
        φ_new[i] = state.φ[i] + δt * (-state.η[i]/params.Fr + κ[i]/params.We)
    end
    
    # Compute pressure
    ps_new = compute_contact_pressure(mesh, φ_new, state.φ, η_new, params, δt, k)
    
    # Force balance
    F_pressure = compute_pressure_force(mesh, ps_new, k)
    F_gravity = -1.0 / params.Fr
    F_total = F_pressure + F_gravity
    
    # Update sphere motion
    ht_new = state.ht + δt * F_total / params.M
    h_new = state.h + δt * ht_new
    
    # Update state
    state.η = η_new
    state.φ = φ_new
    state.ps = ps_new
    state.h = h_new
    state.ht = ht_new
    state.k_contact = k
    state.t += δt
    
    return state
end