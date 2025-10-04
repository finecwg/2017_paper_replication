using BouncingDroplet
using Plots

function test_free_wave()
    println("="^60)
    println("FREE SURFACE WAVE TEST (SIMPLIFIED)")
    println("="^60)
    
    # Smaller domain, fewer points for stability
    mesh = RadialMesh(101, 5.0)
    println("\nMesh: $(mesh.nr) points, R = $(mesh.R)")
    
    # Build operators
    println("Building operators...")
    ΔH = build_laplacian_operator(mesh)
    N_op = DirichletToNeumannOperator(mesh, polar_refinement=6, n_θ=32)
    N_matrix = N_op.operator_matrix  # Matrix 추출
    
    # Physical parameters (water)
    fluid = water_properties()
    solid = SolidProperties(R₀=0.001, ρ_s=1000.0)
    params = dimensionless_sphere_impact(fluid, solid, 0.1)
    
    println("\nDimensionless parameters:")
    println("  Re = $(params.Re)")
    println("  Fr = $(params.Fr)")
    println("  We = $(params.We)")
    
    # Initial condition: smaller Gaussian bump
    state = SimulationState(
        η = 0.001 * exp.(-(mesh.r.^2)/0.1),  # Smaller amplitude, narrower
        φ = zeros(mesh.nr),
        ps = zeros(mesh.nr),
        h = 1.0,
        ht = 0.0,
        k_contact = 0,
        t = 0.0
    )
    
    # Much smaller time step for stability
    δt = 0.0001  # Reduced from 0.001
    n_steps = 1000
    save_interval = 20
    
    println("\nTime stepping:")
    println("  δt = $δt")
    println("  Total steps = $n_steps")
    
    # Storage
    η_history = [copy(state.η)]
    t_history = [state.t]
    energies = [sum(state.η.^2) * mesh.δr]
    
    # Run simulation
    println("\nRunning simulation...")
    for step in 1:n_steps
        time_step_no_contact!(state, mesh, ΔH, N_matrix, params, δt)
        
        # Check for instability
        if any(isnan.(state.η)) || any(isinf.(state.η))
            println("  ERROR: Simulation became unstable at step $step")
            break
        end
        
        if maximum(abs.(state.η)) > 1.0
            println("  WARNING: Large amplitude at step $step: $(maximum(abs.(state.η)))")
        end
        
        if step % save_interval == 0
            push!(η_history, copy(state.η))
            push!(t_history, state.t)
            push!(energies, sum(state.η.^2) * mesh.δr)
            
            if step % 200 == 0
                println("  Step $step / $n_steps, t = $(round(state.t, digits=5)), max|η| = $(round(maximum(abs.(state.η)), digits=6))")
            end
        end
    end
    
    println("\nSimulation complete!")
    println("  Final max|η| = $(maximum(abs.(state.η)))")
    
    # Plot results
    println("\nGenerating plots...")
    
    # Energy evolution
    plot(t_history, energies,
         xlabel="t", ylabel="Energy",
         title="Energy Evolution (should decay or stay constant)",
         lw=2, legend=false, yscale=:log10)
    savefig("free_wave_energy.png")
    println("  Saved: free_wave_energy.png")
    
    # Final state
    plot(mesh.r, η_history[end],
         xlabel="r", ylabel="η",
         title="Final Surface Profile",
         lw=2, legend=false)
    savefig("free_wave_final.png")
    println("  Saved: free_wave_final.png")
    
    # Space-time (if stable)
    if length(η_history) > 10
        η_matrix = hcat(η_history...)
        heatmap(t_history, mesh.r[1:50], η_matrix[1:50, :],
                xlabel="t", ylabel="r",
                title="Free Surface Evolution",
                color=:RdBu)
        savefig("free_wave_spacetime.png")
        println("  Saved: free_wave_spacetime.png")
    end
    
    println("\n" * "="^60)
    println("TEST COMPLETE")
    println("="^60)
    
    return state, η_history, t_history, energies
end

if abspath(PROGRAM_FILE) == @__FILE__
    test_free_wave()
end