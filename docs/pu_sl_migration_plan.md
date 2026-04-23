# PU / SL Solver Migration Plan

This repo will migrate the KernelPack PU semi-Lagrangian and coupled
advection-diffusion solvers in dependency order.

## Order of migration

1. `PUSLAdvectionSolver`
   - Start with the 2D smooth-disk case used in
     `examples/SLAdvection/PUBackwardStepperExample.cpp`
   - Support backward transport first
   - Then add forward transport
   - Required validation before moving on:
     - constant-state preservation
     - mass drift check
     - 2D convergence study on the manufactured rotation example

2. `MultiSpeciesPUSLAdvectionSolver`
   - Thin wrapper on top of the scalar PU-SL transport solver
   - Required validation before moving on:
     - species-wise constant preservation
     - agreement with repeated scalar solves on a shared velocity field

3. `PUSLFDAdvectionDiffusionSolver`
   - Couple PU-SL transport with the existing MATLAB `DiffusionSolver`
   - Match the manufactured problem in
     `examples/AdvectionDiffusion/PUSLAdvectionDiffusionCommon.h`
   - Required validation before moving on:
     - 2D convergence study
     - mass sanity checks for homogeneous Neumann data
     - consistency check against the pure-diffusion and pure-advection limits

4. `PUSLPUAdvectionDiffusionSolver`
   - Couple PU-SL transport with a MATLAB `PUDiffusionSolver`
   - Required validation before moving on:
     - 2D convergence study
     - comparison against the FD-diffusion coupled solver on the same test case

## Ground rules

- Do not move to the next solver until the current solver has a convergence
  study.
- If convergence breaks, zoom out:
  - verify exact-solution and forcing consistency
  - check constant-state preservation
  - check patch partition-of-unity normalization
  - check departure-point tracing on a frozen test velocity field
  - check that transport-only and diffusion-only limits reduce correctly

## Shared C++ references

- `includes/KernelPack/solvers/PUSLAdvectionSolver.h`
- `includes/KernelPack/solvers/PULocalizedSLCommon.h`
- `includes/KernelPack/solvers/PUCollocationSLCommon.h`
- `includes/KernelPack/solvers/PUSLFDAdvectionDiffusionSolver.h`
- `includes/KernelPack/solvers/PUSLPUAdvectionDiffusionSolver.h`
- `examples/SLAdvection/PUBackwardStepperExample.cpp`
- `examples/SLAdvection/PUForwardStepperExample.cpp`
- `examples/AdvectionDiffusion/PUSLAdvectionDiffusionCommon.h`
