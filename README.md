# 2D-Lid-Driven-Cavity-CFD-Solver
A 2D incompressible Navier-Stokes solver for the lid-driven cavity problem, implemented from scratch in Python using the projection method. Validated against the widely-cited Ghia et al. (1982) benchmark data at Re = 100, 400, and 1000.

# What this is
The lid-driven cavity is a classic CFD benchmark: a square domain where the top wall moves at a constant velocity while the other three walls are stationary. It sounds simple, but the flow develops interesting recirculating structures, and the problem has well-established reference solutions at multiple Reynolds numbers — which makes it ideal for validating a solver.

This notebook builds the solver cell by cell, with all the math shown explicitly. If you want to understand how the projection method actually works (not just run it), this is meant to be readable.

# How it works
The solver uses a **projection (fractional step) method** to decouple the velocity and pressure updates:

1. **Predictor step** — advance velocities using the current pressure field, ignoring the incompressibility constraint
2. **Pressure Poisson equation** — solve for the pressure correction that enforces divergence-free flow
3. **Corrector step** — update velocities using the corrected pressure

Spatial derivatives are discretized on a uniform Cartesian grid using central differences. Time integration is first-order explicit (forward Euler).

# Governing equations
The incompressible Navier-Stokes equations:

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}$$

$$\nabla \cdot \mathbf{u} = 0$$

The pressure Poisson equation derived from the projection step:

$$\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$

# Parameters

| Parameter | Value |
|-----------|-------|
| Grid | 65 × 65 |
| Domain | 1 × 1 |
| Density | 1.0 |
| Lid velocity | 1.0 |
| Poisson iterations | 200 (per time step) |

Reynolds number is controlled via kinematic viscosity `nu`:

| Re | nu | Time steps | dt |
|----|----|------------|----|
| 100 | 0.01 | 11,000 | 0.001 |
| 400 | 0.0025 | 30,000 | 0.001 |
| 1000 | 0.001 | 50,000 | 0.0005 |

# Boundary conditions

**Velocity:**
- Top wall: u = 1 (lid velocity), v = 0
- All other walls: u = 0, v = 0

**Pressure:**
- Top wall: p = 0 (Dirichlet)
- Other walls: dp/dn = 0 (Neumann)

# Error summary (RMSE vs. Ghia et al.)

| Re | RMSE u | RMSE v |
|----|--------|--------|
| 100 | 0.0146 | 0.0117 |
| 400 | 0.0255 | 0.0518 |
| 1000 | 0.0637 | 0.0743 |

The Re = 100 case matches closely. Error grows at higher Reynolds numbers, which is expected — at Re = 1000 the flow is more complex and the 65×65 grid starts to feel the resolution limits. A finer mesh or higher-order scheme would improve this.

# Running the notebook

Requirements:
```
numpy
matplotlib
```

Open in Jupyter or Google Colab and run all cells. The solver will print convergence info every 500 iterations and produce validation plots at the end.

Expected runtime (approximate, on a modern laptop):
- Re = 100: ~1 min
- Re = 400: ~3 min
- Re = 1000: ~8 min

# Code Structure
The notebook is organized as follows:

- **Grid setup** — domain, spacing, mesh arrays
- **Finite difference operators** — forward/backward/central differences, Laplacian
- **Boundary conditions** — velocity and pressure BC functions
- **Predictor step** — `vel_without_pressure()`
- **Pressure Poisson solver** — `pressure_poisson()`
- **Corrector step** — `update_velocity()`
- **Main time-stepping loop** — `simulate_cavity_flow()`
- **Validation** — Ghia et al. data, interpolation, RMSE, plots

---

## Limitations and known issues

- The solver is not optimized for speed. NumPy vectorization helps, but the Poisson solver still uses a fixed iteration count rather than convergence-based stopping.
- At Re = 1000, the RMSE is noticeably higher. The coarse grid is the main culprit; the corner vortices are not well-resolved.
- No adaptive time-stepping. The dt values were chosen conservatively by hand.

---

## Reference

Ghia, U., Ghia, K. N., & Shin, C. T. (1982). High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method. *Journal of Computational Physics*, 48(3), 387–411.
