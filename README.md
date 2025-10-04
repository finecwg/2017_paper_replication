
# Bouncing Droplet Simulation

Julia implementation of sphere-liquid interface interactions based on Galeano-Rios et al. (2017).

## Overview

This package simulates the impact of a solid sphere on a liquid surface using:
- Axisymmetric potential flow theory
- Dirichlet-to-Neumann operator for free surface boundary conditions
- Contact mechanics for sphere-surface interactions

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick Start

```julia
using BouncingDroplet, Plots

# Create mesh (dimensionless units)
mesh = RadialMesh(101, 20.0)

# Build operators
ΔH = build_laplacian_operator(mesh)
N_op = DirichletToNeumannOperator(mesh)
N_matrix = N_op.operator_matrix

# Physical parameters
fluid = water_properties()
solid = SolidProperties(R₀=0.001, ρ_s=1000.0)
params = dimensionless_sphere_impact(fluid, solid, 0.1)

# Free surface wave propagation test
state = SimulationState(
    η = 0.0001 * exp.(-(mesh.r.^2)/0.05),
    φ = zeros(mesh.nr),
    ps = zeros(mesh.nr),
    h = 1.0, ht = 0.0, k_contact = 0, t = 0.0
)

δt = 1e-6
for i in 1:1000
    time_step_no_contact!(state, mesh, ΔH, N_matrix, params, δt)
end

plot(mesh.r, state.η, xlabel="r", ylabel="η")
```

## Features

### ✅ Implemented

- **Geometry**: Radial mesh, polar mesh, spherical solid
- **Physics**: Dimensionless numbers (Re, Fr, We, M)
- **Operators**: 
  - Laplacian (finite differences)
  - N operator (Chebyshev regularization)
  - Curvature computation
- **Free Surface**: Stable evolution for ~10ms (dimensionless time)
- **Contact Detection**: Penetration checking, contact radius finding

### ⚠️ Partial Implementation

- **N Operator**: Regularized approximation, not full spectral accuracy
- **Contact Solver**: Detection works, but pressure forces incomplete

### ❌ Not Implemented

- Full augmented system (eq 3.7 from paper)
- Pressure-driven bouncing
- Validation against experimental data

## Current Limitations

1. **Energy Conservation**: ~5% increase over 10,000 steps
2. **Wave Propagation**: Slower than physical
3. **Contact Mechanics**: Sphere penetrates surface instead of bouncing

## Project Structure

```
src/
├── BouncingDroplet.jl      # Main module
├── physics/                 # Physical parameters
├── geometry/                # Meshes and shapes
├── operators/               # Differential operators
├── impact/                  # Contact mechanics
├── solver/                  # Time integration
└── utils/                   # Helper functions
```

## References

Galeano-Rios, C. A., Milewski, P. A., & Vanden-Broeck, J.-M. (2017). 
Quasi-equilibrium solutions for droplets on a vibrating bath. 
*Journal of Fluid Mechanics*, 831, 1-20.

