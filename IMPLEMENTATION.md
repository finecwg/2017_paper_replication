
# Implementation Details

## Technical Achievements

### N Operator (Appendix A)

**Challenge**: Singular integral at r=s requires careful treatment.

**Solution**: Regularized wavenumber
```julia
k_reg = 1.0 / sqrt(r² + δr²)  # Avoids 1/r singularity
```

**Result**: Stable for 10,000+ time steps with minimal energy drift.

### Time Stepping

**Method**: Implicit Euler for η and φ equations
```julia
η^{n+1} = η^n + δt * N[φ^{n+1}]
φ^{n+1} = φ^n + δt * (-η^{n+1}/Fr + κ/We)
```

**Stability**: CFL-like condition δt < 1e-6 for mesh spacing δr=0.2

## Key Design Decisions

1. **Chebyshev over Spectral Collocation**: Simpler implementation, sufficient accuracy
2. **Sparse Matrices**: Memory-efficient for 100+ grid points
3. **Modular Structure**: Each operator is independent

## Performance

- **10,000 time steps**: ~30 seconds
- **Memory**: ~50 MB for 101-point mesh
- **N operator construction**: ~2 seconds

## Known Issues

### Issue 1: Contact Solver Incomplete

**Symptom**: Sphere penetrates surface (h < 0)

**Root Cause**: Pressure force calculation doesn't generate sufficient repulsion

**Files**: `src/solver/contact_solver.jl`, `src/impact/pressure.jl`

**Next Steps**: 
1. Verify pressure integral (eq 3.9)
2. Check force balance (M * ht_t = F_pressure + F_gravity)
3. Add constraint enforcement

### Issue 2: N Operator Accuracy

**Symptom**: Wave speed ~100× slower than theory

**Root Cause**: Regularization strength (0.1) too weak

**Attempted Solutions**:
- Polar mesh integration (failed: NaN handling)
- Finite element (not attempted: complexity)
- **Chebyshev** (current: stable but inaccurate)

**Trade-off**: Stability vs Accuracy

## Testing Protocol

### Free Surface Test
```julia
mesh = RadialMesh(51, 1.0)
state = SimulationState(η = 0.0001*exp.(-(mesh.r.^2)/0.1), ...)
for i in 1:1000
    time_step_no_contact!(...)
end
# Check: max|η| < 0.001, no NaN
```

### Contact Detection Test
```julia
sphere = SphericalSolid(1.0)
state.h = 0.5  # Penetrating
k = find_contact_radius(mesh, state.η, state.h, z_s)
# Check: k > 0
```

## Future Work

### Short Term (1-2 days)
1. Fix pressure force magnitude
2. Add augmented system solver
3. Validate bouncing behavior

### Medium Term (1 week)
1. Implement proper spectral N operator
2. Add adaptive time stepping
3. Optimize performance

### Long Term
1. Compare with experimental data
2. Extend to 2D (drop on bath)
3. Add vibrating substrate

## Lessons Learned

1. **Singular integrals are hard**: 4 different methods attempted
2. **Stability ≠ Accuracy**: Regularization helps stability but hurts physics
3. **Modular testing crucial**: Each component tested independently
4. **Documentation essential**: Complex code needs clear explanation