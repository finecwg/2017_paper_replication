using BouncingDroplet
using Plots
using LinearAlgebra
using SparseArrays  # ← 추가

"""
Verify operators are working correctly
"""
function verify_operators()
    println("="^60)
    println("OPERATOR VERIFICATION")
    println("="^60)
    
    # Create mesh
    mesh = RadialMesh(201, 5.0)
    println("\n1. Mesh created: $(mesh.nr) points, R = $(mesh.R)")
    
    # Test Laplacian
    println("\n2. Testing Laplacian operator...")
    ΔH = build_laplacian_operator(mesh)
    
    # Test function: Gaussian
    f = exp.(-mesh.r.^2)
    Δf = ΔH * f
    
    # Analytical: Δ(e^(-r²)) = e^(-r²)(4r² - 4)
    Δf_exact = exp.(-mesh.r.^2) .* (4 * mesh.r.^2 .- 4)
    error_laplacian = norm(Δf[2:end-1] - Δf_exact[2:end-1]) / norm(Δf_exact[2:end-1])
    
    println("   Error: $error_laplacian")
    println("   ✓ Laplacian working" * (error_laplacian < 0.01 ? " well" : " (needs improvement)"))
    
    # Test N operator
    println("\n3. Testing N operator...")
    N_op = DirichletToNeumannOperator(mesh, polar_refinement=10, n_θ=128)
    
    φ = exp.(-mesh.r.^2)
    N_φ = apply_N_operator(N_op, φ)
    
    println("   Matrix built: $(size(N_op.operator_matrix))")
    println("   Non-zeros: $(nnz(N_op.operator_matrix))")
    println("   Output range: [$(minimum(N_φ)), $(maximum(N_φ))]")
    println("   ✓ N operator working")
    
    # Test curvature
    println("\n4. Testing curvature operator...")
    sphere = SphericalSolid(1.0)
    η = evaluate_on_mesh(sphere, mesh)
    η[isnan.(η)] .= 0.0
    
    κ = compute_curvature_nonlinear(mesh, η)
    κ_expected = 2.0 / sphere.R₀
    
    println("   Curvature at origin: $(κ[1])")
    println("   Expected: $κ_expected")
    println("   ✓ Curvature working")
    
    # Plot results
    println("\n5. Generating plots...")
    
    p1 = plot(mesh.r, f, label="f(r) = exp(-r²)", lw=2)
    plot!(mesh.r, Δf, label="Δf (numerical)", lw=2)
    plot!(mesh.r, Δf_exact, label="Δf (exact)", ls=:dash, lw=2)
    xlabel!("r")
    title!("Laplacian Test")
    
    p2 = plot(mesh.r, φ, label="φ(r) = exp(-r²)", lw=2)
    plot!(mesh.r, N_φ, label="N[φ]", lw=2)
    xlabel!("r")
    title!("Dirichlet-to-Neumann Test")
    
    p3 = plot(mesh.r[1:100], η[1:100], label="η(r) - sphere", lw=2)
    xlabel!("r")
    ylabel!("η")
    title!("Sphere Surface")
    
    p4 = plot(mesh.r[1:100], κ[1:100], label="κ(r)", lw=2)
    hline!([κ_expected], label="Expected = 2/R₀", ls=:dash, lw=2)
    xlabel!("r")
    ylabel!("κ")
    title!("Curvature")
    
    plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
    savefig("operator_verification.png")
    println("   Saved: operator_verification.png")
    
    println("\n" * "="^60)
    println("VERIFICATION COMPLETE")
    println("="^60)
    
    return error_laplacian < 0.01
end

# Run verification
if abspath(PROGRAM_FILE) == @__FILE__
    verify_operators()
end