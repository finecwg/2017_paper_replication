module BouncingDroplet

using LinearAlgebra
using SparseArrays
using StaticArrays
using Parameters
using Interpolations

# Physics
include("physics/parameters.jl")
include("physics/fluid_equations.jl")
include("physics/solid_motion.jl")

# Geometry
include("geometry/mesh.jl")
include("geometry/polar_mesh.jl")
include("geometry/sphere.jl")

# Operators
include("operators/laplacian.jl")
include("operators/curvature.jl")
include("operators/dirichlet_neumann.jl")

# Utils
include("utils/interpolation.jl")
include("utils/integration.jl")
include("utils/diagnostics.jl")

# Impact (constraints/pressure/force BEFORE solver)
include("impact/constraints.jl")
include("impact/pressure.jl")
include("impact/force.jl")

# Solver (time_stepping BEFORE contact_solver!)
include("solver/matrix_assembly.jl")
include("solver/time_stepping.jl")        # ← SimulationState 정의됨
include("solver/adaptive_timestep.jl")
include("solver/contact_solver.jl")       # ← 이게 마지막

# Exports
export FluidProperties, SolidProperties, DimensionlessNumbers
export dimensionless_sphere_impact, dimensionless_bouncing_droplet
export water_properties, silicone_oil_20cst

export RadialMesh, PolarMesh, SphericalSolid
export create_mesh_for_sphere, evaluate_on_mesh
export surface_height, surface_slope, surface_curvature

export build_laplacian_operator, apply_laplacian
export DirichletToNeumannOperator, apply_N_operator
export compute_curvature_linear, compute_curvature_nonlinear

export SimulationState, initialize_impact_state
export time_step_no_contact!, time_step_with_contact!
export solve_contact_step

export check_penetration, check_tangency, find_contact_radius, is_valid_contact
export compute_contact_pressure, compute_pressure_force

end # module