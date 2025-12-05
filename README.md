# 2DPOLAR

2DPOLAR is a small 2D polar-coordinate Euler solver that advances a
compressible flow using a finite-volume formulation, HLL Riemann fluxes,
and a stiffened-gas equation of state (SG-EOS). The code is organized into
modular Fortran 2008 units for mesh generation, fluid properties, field
storage, flux evaluation, time stepping, and boundary handling.

## Build and run

```bash
make
./main
```

The `Makefile` builds all headers and source files into the `main`
executable using `gfortran` with Fortran 2008 semantics. Object files and
module interfaces are placed in `build/`.【F:Makefile†L1-L26】 The solver writes
Tecplot-style snapshots under `result/` every 20 steps by default.

## Core data structures

### Mesh configuration (`mesh_types`)
* `PolarGridParameters` bundles grid resolution, guard-cell count, radial
  and angular bounds, equation type (Cartesian transform vs. polar
  Euler), and boundary-condition selectors. It also precomputes index
  ranges for physical and computational domains.【F:header/mesh_types.F90†L24-L77】
* `TimeControlParameters` records CFL settings, time bounds, timestep
  limits, counters, and user-requested output times.【F:header/mesh_types.F90†L79-L97】
* `PolarMesh` stores the geometric arrays (cell centers, face locations,
  spacings, face lengths, cell areas, Cartesian coordinates) allocated on
  index ranges that include ghost layers.【F:header/mesh_types.F90†L104-L152】 The
  helper `allocate_mesh_arrays`/`deallocate_mesh_arrays` routines size all
  geometry arrays based on the configured resolution and guard cells.【F:header/mesh_types.F90†L158-L201】

### Flow fields (`flow_fields`)
* `PrimitiveVariables` holds the cell-centered primitive state `(ρ, p, u,
  v)`. `ConservedVariables` stores `(ρ, ρu, ρv, ρE)` for finite-volume
  updates. Each has paired allocate/deallocate helpers to size arrays on
  provided index bounds.【F:header/flow_fields.F90†L14-L64】

### Flux storage (`flux_types`)
* `FaceFlux2D` collects mass, momentum, and energy flux arrays for one
  face family. `FluxFields` wraps separate radial and angular flux sets.
  Allocation aligns with the mesh’s physical index ranges: radial faces
  are `(Nr+1)×Nθ`, angular faces are `Nr×(Nθ+1)`.【F:header/flux_types.F90†L15-L57】

### Equation-of-state support (`fluid_properties`)
* `StiffenedGas` stores SG-EOS parameters `(γ, p_inf)` plus pressure-floor
  controls. Factory functions return presets for air and water or allow
  custom parameters.【F:header/fluid_properties.F90†L24-L63】
* Utility functions provide sound speed, specific internal/total energy,
  and an optional pressure floor that clips negative pressures for stiff
  fluids.【F:header/fluid_properties.F90†L66-L149】

## Key workflows and routines

### Mesh setup (`mesh_ops`)
`mesh_setup` builds a uniform polar grid, fills geometry arrays, allocates
primitive and conserved fields (including ghost cells), and seeds time
control parameters (CFL, end time, output cadence). It prints a mesh
summary and populates the `time_outputs` schedule.【F:source/mesh_ops.F90†L1-L114】【F:source/mesh_ops.F90†L116-L170】

### Initialization (`initial_ops`)
`initialize_flow_fields` selects an air stiffened-gas model, assigns
initial primitive variables (currently a rotating unit-velocity field with
1 kg/m³ density and 1e5 Pa pressure), converts them to conserved form, and
applies boundary conditions to both views.【F:source/initial_ops.F90†L1-L108】【F:source/initial_ops.F90†L110-L178】

### Primitive–conserved conversion (`field_conversion`)
`update_conserved_from_primitive` computes conserved quantities from the
primitive state using SG-EOS specific total energy. The inverse
`update_primitive_from_conserved` reconstructs primitive variables from
conserved arrays, optionally applying a pressure floor (currently
commented).【F:source/field_conversion.F90†L1-L109】【F:source/field_conversion.F90†L111-L183】

### Boundary conditions (`boundary_update`)
`update_boundary_primitive_conserved` dispatches to radial and angular
handlers based on the mesh’s BC selectors. Supported radial modes include
extrapolation, periodic, and axis symmetry (with sign flips on radial
momentum). Angular boundaries support extrapolation and periodic wrapping.
Both primitive and conserved ghost cells are filled consistently.【F:source/boundary_update.F90†L1-L160】

### Time stepping (`time_step_control`)
`compute_global_timestep` evaluates the CFL-limited timestep using local
radial/azimuthal characteristic speeds `( |u_r|+c, |u_θ|+c )` and effective
length scales `(Δr, rΔθ)`, clamped between user-specified maxima and
minima.【F:source/time_step_control.F90†L10-L102】

### Flux computation (`flux_ops`)
`compute_all_fluxes` selects between Cartesian-mapped and polar Euler
forms. For polar mode, `compute_radial_fluxes_polar` and
`compute_angular_fluxes_polar` call the HLL solver with appropriate normal
vectors and scale by face lengths; a metric-form variant exists for mapped
Cartesian equations.【F:source/flux_ops.F90†L1-L73】【F:source/flux_ops.F90†L119-L193】

### Conserved updates (`field_update_ops`)
`update_conserved_field` routes to metric or polar updates. The metric
path performs pure divergence of fluxes over cell volume. The polar path
adds geometric source terms `(ρ u_θ² + p)/r` and `−ρ u_r u_θ / r` scaled by
cell area to the radial and angular momentum equations, respectively.【F:source/field_update_ops.F90†L1-L85】【F:source/field_update_ops.F90†L87-L167】

### Driver (`main.F90`)
`polar_sim` orchestrates the simulation: it builds the mesh, initializes
fields and flux storage, writes the initial Tecplot output, then advances
until the end time or maximum steps. Each iteration computes the timestep,
fluxes, conserved update, primitive reconstruction, boundary refresh, and
periodic output.【F:source/main.F90†L1-L63】

## Extending the solver

* **Alternative initial conditions:** modify
  `initialize_flow_fields` to set different primitive profiles before the
  conserved conversion call.【F:source/initial_ops.F90†L19-L105】
* **Boundary behaviors:** adjust BC selectors in `mesh_setup` or extend
  `apply_bc_r`/`apply_bc_theta` for new types.【F:source/mesh_ops.F90†L20-L53】【F:source/boundary_update.F90†L35-L160】
* **Output cadence:** change `time%output_step` logic in `main.F90` or
  leverage the precomputed `time%time_outputs` array (currently commented
  in the loop).【F:source/main.F90†L33-L61】【F:source/mesh_ops.F90†L55-L105】
