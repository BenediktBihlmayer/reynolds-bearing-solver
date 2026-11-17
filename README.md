# Reynolds Equation Solver for a Journal Bearing

A C++ implementation of a finite-difference solver for the steady
Reynolds equation in a hydrodynamic journal bearing.

The program calculates the lubricant pressure distribution and integrates
the pressure over the bearing surface to obtain the resultant bearing
force. It was developed in the context of my mechanical engineering
master's thesis to support comparisons with CFD simulations.

## Mathematical Model

The solver uses the steady Reynolds equation for an incompressible,
isoviscous lubricant:

```math
\frac{\partial}{\partial x}\left(h^3 \frac{\partial p}{\partial x}\right) + \frac{\partial}{\partial z}\left(h^3 \frac{\partial p}{\partial z}\right) = 6 \eta U \frac{\partial h}{\partial x}
```

Here, x = Rφ is the circumferential coordinate, z is the axial coordinate,
h is the film thickness, p is the pressure, η is the dynamic viscosity,
and U = ωR is the shaft surface velocity.

The film thickness is prescribed as:

$$
h(\varphi) = C \left[1 + \varepsilon
\cos(\varphi - \delta_\varphi)\right]
$$

### Assumptions and Boundary Conditions

- Steady, laminar lubrication flow
- Incompressible lubricant with constant viscosity
- Prescribed bearing geometry and eccentricity
- Periodic pressure in the circumferential direction
- Zero axial pressure gradient at both axial boundaries
- No cavitation model; negative pressure values are retained
- No thermal effects or elastic deformation

The axial boundary conditions represent sealed bearing ends.
They do not model an open-ended bearing with lubricant leakage.

Although the computational grid has two coordinates, the prescribed
geometry and boundary conditions admit an axially uniform pressure
solution. Bearing width therefore does not introduce the pressure
reduction associated with axial leakage.

## Numerical Method

- Nonuniform circumferential grid
- Uniform axial grid
- Conservative flux-form spatial discretization
- Arithmetic averaging of h³ between circumferential nodes
- Gauss-Seidel iteration
- Stopping criterion based on the maximum pressure change per iteration
- Numerical integration of pressure to obtain Cartesian force components

The stopping criterion measures iterative changes, not discretization
error or the equation residual.

## Default Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Shaft radius R | 0.026979 | m |
| Bearing width B | 0.0182 | m |
| Radial clearance C | 0.000032 | m |
| Eccentricity ratio ε | 0.66 | dimensionless |
| Angular velocity ω | 200.0 | rad/s |
| Dynamic viscosity η | 0.010115 | Pa·s |
| Minimum film-thickness angle | 315.0 | degrees |
| Circumferential nodes | 100 | |
| Axial nodes | 50 | |

Parameters are currently defined directly in the source code.

## Build and Run

Requires a C++17-compatible compiler. No external libraries are required.

Compile with GCC or Clang:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic reynolds_solver.cpp -o reynolds_solver
```

Run on Linux or macOS:

```bash
./reynolds_solver
```

## Output Files

Files are written to the current working directory.

| File | Contents |
|------|----------|
| `gleitlager_druck.csv` | Circumferential angle in radians, axial position in metres, pressure in pascals; no header |
| `gleitlager_kraefte.csv` | Cartesian force components, resultant magnitude, and force direction |
| `gleitlager_druck_verschoben.csv` | Pressure data with a shifted angular coordinate for plotting |

The angular shift is a post-processing operation. Forces are calculated
using the original physical coordinates.

## Current Status and Known Limitations

This source version was reconstructed from a thesis code listing.
Compilation and numerical verification of the reconstructed version
have not yet been documented in this repository.

Known limitations:

- The periodic and Neumann boundary conditions leave the pressure
  undetermined up to an additive constant. An explicit pressure
  reference or mean-pressure constraint is not yet implemented.
- The shifted pressure export appends an endpoint using the original
  zero-angle pressure. This is generally inconsistent with the shifted
  coordinate and requires correction.
- The zero-crossing search uses the first nonpositive grid value
  after the pressure maximum, without interpolation.
- Reaching the iteration limit does not prevent the program from
  exporting results. Check the console output before using them.

Results should be treated as provisional until these issues are
resolved and verification results are documented.

## Author

Benedikt Bihlmayer  
M.Sc. Mechanical Engineering
