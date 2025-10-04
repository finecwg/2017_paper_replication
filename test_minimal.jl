# test_minimal.jl
println("Step 1: Testing basic Julia...")
println("Julia version: ", VERSION)

println("\nStep 2: Testing file existence...")
files_to_check = [
    "src/BouncingDroplet.jl",
    "src/physics/parameters.jl",
    "src/geometry/mesh.jl",
    "src/geometry/polar_mesh.jl",
    "src/geometry/sphere.jl",
    "src/operators/laplacian.jl",
    "src/operators/curvature.jl",
    "src/utils/interpolation.jl"
]

for file in files_to_check
    if isfile(file)
        println("  ✓ $file")
    else
        println("  ✗ $file MISSING!")
    end
end

println("\nStep 3: Testing individual includes...")
try
    include("src/physics/parameters.jl")
    println("  ✓ parameters.jl loaded")
catch e
    println("  ✗ parameters.jl failed: ", e)
end

try
    include("src/geometry/mesh.jl")
    println("  ✓ mesh.jl loaded")
catch e
    println("  ✗ mesh.jl failed: ", e)
end

println("\nStep 4: Testing struct creation...")
try
    mesh = RadialMesh(51, 1.0)
    println("  ✓ RadialMesh created: $(mesh.nr) points")
catch e
    println("  ✗ RadialMesh failed: ", e)
end

try
    water = water_properties()
    println("  ✓ Water properties: ρ=$(water.ρ) kg/m³")
catch e
    println("  ✗ Water properties failed: ", e)
end