"""
State vector for simulation
"""
@with_kw mutable struct SimulationState
    η::Vector{Float64}      # Free surface elevation
    φ::Vector{Float64}      # Velocity potential at surface
    ps::Vector{Float64}     # Pressure on contact area
    h::Float64              # Solid bottom height
    ht::Float64             # Solid vertical velocity
    k_contact::Int          # Number of contact points
    t::Float64              # Current time
end

"""
Initialize state for sphere impact
"""
function initialize_impact_state(mesh::RadialMesh, 
                                sphere::SphericalSolid,
                                h0::Float64,
                                v0::Float64)
    nr = mesh.nr
    
    return SimulationState(
        η = zeros(nr),
        φ = zeros(nr),
        ps = zeros(nr),
        h = h0,
        ht = v0,
        k_contact = 0,
        t = 0.0
    )
end

"""
Perform one time step WITHOUT contact (for testing)
"""
function time_step_no_contact!(state::SimulationState,
                              mesh::RadialMesh,
                              ΔH::SparseMatrixCSC,
                              N_matrix::SparseMatrixCSC,
                              params::DimensionlessNumbers,
                              δt::Float64;
                              forcing::Function = (t) -> 0.0)
    
    nr = mesh.nr
    I_mat = sparse(I, nr, nr)
    
    # Linear system
    block_11 = I_mat - (2*δt/params.Re) * ΔH
    block_12 = -δt * N_matrix
    block_21 = (δt/params.Fr) * I_mat
    block_22 = I_mat - (2*δt/params.Re) * ΔH
    
    Q = [block_11 block_12;
         block_21 block_22]
    
    # Curvature
    κ = compute_curvature_linear(mesh, state.η)
    
    # RHS - forcing을 벡터로 broadcast
    forcing_vec = fill(forcing(state.t), nr)  # 스칼라를 벡터로 변환
    
    F = [state.η; 
         state.φ .+ (δt/params.We) .* κ .- (δt/params.Fr) .* forcing_vec]
    
    # Solve
    W_new = Q \ F
    
    state.η = W_new[1:nr]
    state.φ = W_new[nr+1:2*nr]
    state.t += δt
    
    return state
end
"""
Perform one time step WITH contact (placeholder)
"""
function time_step_with_contact!(state::SimulationState,
                                mesh::RadialMesh,
                                sphere::SphericalSolid,
                                ΔH::SparseMatrixCSC,
                                N_matrix::SparseMatrixCSC,
                                params::DimensionlessNumbers,
                                δt::Float64)
    
    z_s = evaluate_on_mesh(sphere, mesh)
    
    # Check if in contact
    min_gap = minimum(state.h .+ z_s .- state.η)
    
    if min_gap > 1e-6
        # No contact
        time_step_no_contact!(state, mesh, ΔH, N_matrix, params, δt)
    else
        # In contact - placeholder
        state.h += state.ht * δt
        state.ht -= (δt/params.Fr)
    end
    
    return state
end