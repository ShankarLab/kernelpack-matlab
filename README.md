# kernelpack-matlab

`kernelpack-matlab` is a MATLAB companion project to
[KernelPack](https://github.com/VarShankar/kernelpack).

It brings the main geometry, node-generation, polynomial, RBF-FD, and solver
ingredients of KernelPack into a MATLAB codebase that is easier to inspect,
prototype with, and extend.

## What it includes

- Geometry objects for smooth and piecewise-smooth boundaries and surfaces:
  `EmbeddedSurface`, `PiecewiseSmoothEmbeddedSurface`, and `RBFLevelSet`
- Seeded fixed-radius Poisson node generation in axis-aligned boxes, with
  geometry-aware clipping and boundary refinement through `DomainNodeGenerator`
- A compact `DomainDescriptor` for interior, boundary, ghost nodes, normals,
  and nearest-neighbor search structures
- Shared polynomial utilities in `kp.poly`, including Legendre-based
  `PolynomialBasis`
- RBF-FD and weighted least-squares stencil and assembly classes in `kp.rbffd`
- Fixed-domain `PoissonSolver` and `DiffusionSolver`

The main packages live in:

- [`+kp/+geometry`](+kp/+geometry)
- [`+kp/+nodes`](+kp/+nodes)
- [`+kp/+domain`](+kp/+domain)
- [`+kp/+poly`](+kp/+poly)
- [`+kp/+rbffd`](+kp/+rbffd)
- [`+kp/+solvers`](+kp/+solvers)

## Supported workflows

- Smooth closed curves and smooth closed 3D surfaces
- Open curve segments and open 3D surface patches
- Piecewise-smooth planar boundaries and piecewise 3D surfaces
- Fixed-radius Poisson disk sampling in any dimension
- Level-set clipping of box clouds to geometry-defined domains
- Standard and overlapped RBF-FD assembly
- Fixed-domain Poisson solves
- Fixed-domain diffusion stepping with BDF1, BDF2, and BDF3

## Examples and checks

Examples:

- [`examples/geometry_examples.m`](examples/geometry_examples.m)
- [`examples/nodes_examples.m`](examples/nodes_examples.m)
- [`examples/poisson_solver_example.m`](examples/poisson_solver_example.m)
- [`examples/poisson_solver_example_3d.m`](examples/poisson_solver_example_3d.m)
- [`examples/diffusion_solver_example.m`](examples/diffusion_solver_example.m)
- [`examples/poisson_convergence_2d.m`](examples/poisson_convergence_2d.m)
- [`examples/poisson_convergence_2d_neumann.m`](examples/poisson_convergence_2d_neumann.m)
- [`examples/poisson_convergence_3d.m`](examples/poisson_convergence_3d.m)
- [`examples/diffusion_convergence_3d.m`](examples/diffusion_convergence_3d.m)

Checks:

- [`tests/geometry_checks.m`](tests/geometry_checks.m)
- [`tests/nodes_checks.m`](tests/nodes_checks.m)
- [`tests/poly_checks.m`](tests/poly_checks.m)
- [`tests/rbffd_checks.m`](tests/rbffd_checks.m)
- [`tests/poisson_solver_checks.m`](tests/poisson_solver_checks.m)
- [`tests/diffusion_solver_checks.m`](tests/diffusion_solver_checks.m)

## Quick examples

### Smooth 2D geometry

```matlab
% Define a smooth closed planar boundary.
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

% Build the geometric model and level set.
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve,1));
surface.buildLevelSetFromGeometricModel([]);

% Extract boundary samples and normals from the fitted representation.
xb = surface.getSampleSites();
nrmls = surface.getNrmls();
```

### Geometry-clipped interior nodes

```matlab
% Generate an interior-plus-boundary domain from a geometry and target spacing.
generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

% Pull out the packed node sets.
Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
Xg = domain.getGhostNodes();
```

### RBF-FD operator assembly

```matlab
% Ask the code to choose stencil parameters from a target accuracy.
sp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 4, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');

% Record stencil metadata during assembly.
op = kp.rbffd.OpProperties('recordStencils', true);

% Assemble a Laplacian on the domain descriptor.
assembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
assembler.AssembleOp(domain, 'lap', sp, op);
L = assembler.getOp();
```

### End-to-end Poisson solve with pure Neumann data

```matlab
% Build a smooth closed domain.
t = linspace(0, 2*pi, 120).';
t(end) = [];
curve = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve,1));
surface.buildLevelSetFromGeometricModel([]);

% Generate interior, boundary, and ghost nodes.
generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5);

% Set up an RBF-FD Poisson solve on that domain.
solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');

% Use a named target order instead of a magic number.
targetOrder = 4;
solver.init(domain, targetOrder);

% Manufactured pure-Neumann problem on the unit disk.
uExact = @(X) (X(:,1).^2 + X(:,2).^2).^2 - (X(:,1).^2 + X(:,2).^2) + 1/6;
forcing = @(X) 4 - 16*(X(:,1).^2 + X(:,2).^2);
neuCoeff = @(Xb) ones(size(Xb,1), 1);
dirCoeff = @(Xb) zeros(size(Xb,1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) ...
    sum(([4*Xb(:,1).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,1), ...
          4*Xb(:,2).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,2)]).*nr, 2);

% Solve and align the mean for comparison.
result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
Xphys = domain.getIntBdryNodes();
u = result.u;
uTrue = uExact(Xphys);
u = u - mean(u - uTrue);
```

### Diffusion stepping

```matlab
% Set up a fixed-domain diffusion stepper on the same domain.
solver = kp.solvers.DiffusionSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');

% Choose the diffusivity and time step.
nu = 0.25;
dt = 0.02;

% Use a named target order instead of a magic number.
targetOrder = 4;
solver.init(domain, targetOrder, dt, nu);

% Define a manufactured transient problem with Dirichlet data.
uExact = @(time, X) exp(-time) .* (X(:,1).^2 + X(:,2).^2);
forcing = @(nuValue, time, X) -exp(-time) .* (X(:,1).^2 + X(:,2).^2) ...
    - 4 * nuValue * exp(-time);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, time, Xb) uExact(time, Xb);

% Seed the state history and take BDF steps.
solver.setInitialState(uExact(0, domain.getIntBdryNodes()));
u1 = solver.bdf1Step(dt, forcing, neuCoeff, dirCoeff, bc);
u2 = solver.bdf2Step(2 * dt, forcing, neuCoeff, dirCoeff, bc);
u3 = solver.bdf3Step(3 * dt, forcing, neuCoeff, dirCoeff, bc);
```

## Notes

- This is a MATLAB implementation of the main KernelPack ingredients, not a
  full one-to-one port of every C++ path.
- Pure Neumann Poisson problems are solved with the usual nullspace
  augmentation, so comparisons should be made after aligning the constant.
- The repo includes both `RBFStencil` and `WeightedLeastSquaresStencil`, but
  the strongest convergence results so far are on the RBF-FD path.
