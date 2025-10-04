# test/test_integration.jl
using Test
using BouncingDroplet

@testset "Integration Tests" begin
    @testset "Free Wave Propagation" begin
        mesh = RadialMesh(101, 5.0)
        ΔH = build_laplacian_operator(mesh)
        N_op = DirichletToNeumannOperator(mesh, polar_refinement=6, n_θ=32)
        
        fluid = water_properties()
        solid = SolidProperties(R₀=0.001, ρ_s=1000.0)
        params = dimensionless_sphere_impact(fluid, solid, 0.1)
        
        state = SimulationState(
            η = 0.01 * exp.(-mesh.r.^2),
            φ = zeros(mesh.nr),
            ps = zeros(mesh.nr),
            h = 1.0,
            ht = 0.0,
            k_contact = 0,
            t = 0.0
        )
        
        # Run a few steps
        δt = 0.001
        for i in 1:10
            time_step_no_contact!(state, mesh, ΔH, N_op, params, δt)
        end
        
        # Check stability
        @test !any(isnan.(state.η))
        @test !any(isinf.(state.η))
        @test maximum(abs.(state.η)) < 1.0  # Reasonable magnitude
        
        println("✓ Free wave propagation test passed")
    end
end