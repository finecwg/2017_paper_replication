"""
Assemble the system matrix for fluid evolution (without contact)
"""
function assemble_fluid_matrix(mesh::RadialMesh,
                               ΔH::SparseMatrixCSC,
                               N_matrix::SparseMatrixCSC,  # DirichletToNeumannOperator 대신 직접 matrix
                               params::DimensionlessNumbers,
                               δt::Float64)
    nr = mesh.nr
    I_mat = sparse(I, nr, nr)
    
    # System: [η; φ]
    block_11 = I_mat - (2*δt/params.Re) * ΔH
    block_12 = -δt * N_matrix
    block_21 = (δt/params.Fr) * I_mat - (δt/params.We) * ΔH
    block_22 = I_mat - (2*δt/params.Re) * ΔH
    
    # Assemble
    Q = [block_11 block_12;
         block_21 block_22]
    
    return Q
end

"""
Placeholder for full contact assembly
"""
function assemble_system_matrix(mesh::RadialMesh, 
                               ΔH::SparseMatrixCSC,
                               N_matrix::SparseMatrixCSC,
                               params::DimensionlessNumbers,
                               k::Int,
                               δt::Float64,
                               sphere::SphericalSolid,
                               z_s::Vector{Float64})
    # Placeholder - return identity for now
    nr = mesh.nr
    return sparse(I, 2*nr, 2*nr)
end