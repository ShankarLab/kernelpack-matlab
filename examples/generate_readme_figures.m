function generate_readme_figures(whichFigure)
%GENERATE_README_FIGURES Render the figures shown in the top-level README.

rootDir = fileparts(fileparts(mfilename('fullpath')));
imageDir = fullfile(rootDir, 'docs', 'images');
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end

if nargin < 1
    whichFigure = "all";
end
whichFigure = string(whichFigure);

renderOne(whichFigure, "readme_smooth_2d_geometry", @() renderSmooth2DGeometry(fullfile(imageDir, 'readme_smooth_2d_geometry.png')));
renderOne(whichFigure, "readme_geometry_clipped_nodes", @() renderGeometryClippedNodes(fullfile(imageDir, 'readme_geometry_clipped_nodes.png')));
renderOne(whichFigure, "readme_rbffd_operator", @() renderRbfFdAssembly(fullfile(imageDir, 'readme_rbffd_operator.png')));
renderOne(whichFigure, "readme_poisson_neumann", @() renderPoissonNeumann(fullfile(imageDir, 'readme_poisson_neumann.png')));
renderOne(whichFigure, "readme_diffusion_stepping", @() renderDiffusionStepping(fullfile(imageDir, 'readme_diffusion_stepping.png')));
renderOne(whichFigure, "readme_incompressible_euler_velocity", @() renderIncompressibleEulerVelocity(fullfile(imageDir, 'readme_incompressible_euler_velocity.png')));
end

function renderOne(requested, name, thunk)
if requested == "all" || requested == name
    thunk();
end
end

function renderSmooth2DGeometry(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();

fig = figure('Color', 'w', 'Position', [100 100 1100 360]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(curve(:, 1), curve(:, 2), 'ko', 'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 5);
axis equal;
grid on;
title('Data Sites');

nexttile;
plot(xb(:, 1), xb(:, 2), 'b.', 'MarkerSize', 12);
axis equal;
grid on;
title('Boundary Samples');

nexttile;
plot(xb(:, 1), xb(:, 2), '.', 'MarkerSize', 12);
hold on;
quiver(xb(:, 1), xb(:, 2), 0.05 * nrmls(:, 1), 0.05 * nrmls(:, 2), 0);
axis equal;
grid on;
title('Boundary Normals');

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderGeometryClippedNodes(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
Xg = domain.getGhostNodes();

fig = figure('Color', 'w', 'Position', [100 100 720 560]);
plot(Xi(:, 1), Xi(:, 2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 10);
hold on;
plot(Xb(:, 1), Xb(:, 2), '.', 'Color', [0.85 0.15 0.15], 'MarkerSize', 12);
plot(Xg(:, 1), Xg(:, 2), '.', 'Color', [0.2 0.45 0.9], 'MarkerSize', 10);
axis equal;
grid on;
legend({'Interior', 'Boundary', 'Ghost'}, 'Location', 'best');
title('Geometry-Clipped Domain Nodes');
finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderRbfFdAssembly(outFile)
t = linspace(0, 2*pi, 50).';
t(end) = [];
curve = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.05, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

sp = kp.rbffd.StencilProperties.fromAccuracy( ...
    'Operator', 'lap', ...
    'ConvergenceOrder', 4, ...
    'Dimension', 2, ...
    'Approximation', 'rbf', ...
    'treeMode', 'all', ...
    'pointSet', 'interior_boundary');
op = kp.rbffd.OpProperties('recordStencils', true);

assembler = kp.rbffd.FDDiffOp(@() kp.rbffd.RBFStencil());
assembler.AssembleOp(domain, 'lap', sp, op);
L = assembler.getOp();

fig = figure('Color', 'w', 'Position', [100 100 1000 420]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
Xi = domain.getInteriorNodes();
Xb = domain.getBdryNodes();
plot(Xi(:, 1), Xi(:, 2), '.', 'Color', [0.15 0.15 0.15], 'MarkerSize', 10);
hold on;
plot(Xb(:, 1), Xb(:, 2), '.', 'Color', [0.85 0.15 0.15], 'MarkerSize', 12);
axis equal;
grid on;
title('Domain Nodes');

nexttile;
spy(L);
title(sprintf('Laplacian Sparsity (%d x %d)', size(L, 1), size(L, 2)));

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderPoissonNeumann(outFile)
t = linspace(0, 2*pi, 120).';
t(end) = [];
curve = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.08, ...
    'Seed', 17, ...
    'StripCount', 5);

solver = kp.solvers.PoissonSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
targetOrder = 4;
solver.init(domain, targetOrder);

uExact = @(X) (X(:,1).^2 + X(:,2).^2).^2 - (X(:,1).^2 + X(:,2).^2) + 1/6;
forcing = @(X) 4 - 16*(X(:,1).^2 + X(:,2).^2);
neuCoeff = @(Xb) ones(size(Xb,1), 1);
dirCoeff = @(Xb) zeros(size(Xb,1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, Xb) ...
    sum(([4*Xb(:,1).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,1), ...
          4*Xb(:,2).*(Xb(:,1).^2 + Xb(:,2).^2) - 2*Xb(:,2)]).*nr, 2);

result = solver.solve(forcing, neuCoeff, dirCoeff, bc);
Xphys = domain.getIntBdryNodes();
tri = delaunay(Xphys(:, 1), Xphys(:, 2));
u = result.u;
uTrue = uExact(Xphys);
u = u - mean(u - uTrue);
absErr = abs(u - uTrue);

fig = figure('Color', 'w', 'Position', [100 100 1000 420]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u, u, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Pure-Neumann Poisson Solution');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), absErr, absErr, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Absolute Error');
colorbar;

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderDiffusionStepping(outFile)
t = linspace(0, 2*pi, 80).';
t(end) = [];
curve = [cos(t), 0.8*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(curve);
surface.buildClosedGeometricModelPS(2, 0.06, size(curve, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, 0.1, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.5, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);

solver = kp.solvers.DiffusionSolver( ...
    'LapAssembler', 'fd', ...
    'BCAssembler', 'fd', ...
    'LapStencil', 'rbf', ...
    'BCStencil', 'rbf');
nu = 0.25;
dt = 0.02;
targetOrder = 4;
solver.init(domain, targetOrder, dt, nu);

uExact = @(time, X) exp(-time) .* (X(:,1).^2 + X(:,2).^2);
forcing = @(nuValue, time, X) -exp(-time) .* (X(:,1).^2 + X(:,2).^2) ...
    - 4 * nuValue * exp(-time);
neuCoeff = @(Xb) zeros(size(Xb, 1), 1);
dirCoeff = @(Xb) ones(size(Xb, 1), 1);
bc = @(NeuCoeffs, DirCoeffs, nr, time, Xb) uExact(time, Xb);

tFinal = 0.5;
nSteps = round(tFinal / dt);
Xphys = domain.getIntBdryNodes();
tri = delaunay(Xphys(:, 1), Xphys(:, 2));

solver.setInitialState(uExact(0, Xphys));
states = {solver.currentPhysicalState()};
for step = 1:nSteps
    time = step * dt;
    if step == 1
        uNext = solver.bdf1Step(time, forcing, neuCoeff, dirCoeff, bc);
    elseif step == 2
        uNext = solver.bdf2Step(time, forcing, neuCoeff, dirCoeff, bc);
    else
        uNext = solver.bdf3Step(time, forcing, neuCoeff, dirCoeff, bc);
    end
    states{end + 1, 1} = uNext; %#ok<AGROW>
end

u1 = states{2};
u2 = states{3};
uFinal = states{end};
uTrueFinal = uExact(tFinal, Xphys);
absErr = abs(uFinal - uTrueFinal);

fig = figure('Color', 'w', 'Position', [100 100 1200 380]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u1, u1, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('BDF1 State');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), u2, u2, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('BDF2 State');
colorbar;

nexttile;
trisurf(tri, Xphys(:, 1), Xphys(:, 2), absErr, absErr, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Final Absolute Error');
colorbar;

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function renderIncompressibleEulerVelocity(outFile)
threads = 1;
xi_u = 4;
xi_p = 4;
h = 0.10;
dt = 0.02;
tFinal = 0.04;

velocityStencil = makeEulerStencilProps(xi_u);
pressureStencil = makeEulerStencilProps(xi_p);
dual = makeEulerDomain(h, 0.40, pressureStencil.n, threads);

solver = kp.solvers.PUSLIncompressibleEulerSolver();
solver.init(dual, xi_u, velocityStencil, pressureStencil, dt, threads);
solver.setTangentialFlowBoundary(1.0e-5);

Xu = dual.getVelocityDomain().getIntBdryNodes();
Xp = dual.getPressureDomain().getIntBdryNodes();
problem = kp.solvers.detail.IncompressibleEulerBDFBackend.defaultProblemDefinition();
problem.slip_walls = {kp.solvers.detail.IncompressibleEulerBDFBackend.stationarySlipWall((1:dual.getVelocityDomain().getNumBdryNodes()).')};
problem.gauge_options.mode = "forcepressuremean";

solver.setInitialVelocity(eulerVelocityExact(0.0, Xu), problem);
solver.bdf1Step(dt, @eulerRk4Step, @eulerForcing, problem);
sol = solver.bdf2Step(tFinal, @eulerRk4Step, @eulerForcing, problem);
pExact = eulerPressureExact(Xp);
p = sol.pressure - mean(sol.pressure - pExact);
triP = delaunay(Xp(:, 1), Xp(:, 2));
arrowIdx = 1:4:size(Xu, 1);

fig = figure('Color', 'w', 'Position', [100 100 1120 430]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Xu(:, 1), Xu(:, 2), '.', 'MarkerSize', 8);
hold on;
quiver(Xu(arrowIdx, 1), Xu(arrowIdx, 2), sol.velocity(arrowIdx, 1), sol.velocity(arrowIdx, 2), 0, 'k', 'LineWidth', 0.8);
axis equal;
axis([-1 1 -1 1]);
grid on;
title('PU-SL Euler Velocity Field');

nexttile;
trisurf(triP, Xp(:, 1), Xp(:, 2), p, p, 'EdgeColor', 'none');
shading interp;
view(2);
axis equal tight;
grid on;
title('Pressure Field');
colorbar;

finalizeReadmeFigure(fig);
exportgraphics(fig, outFile, 'Resolution', 180);
close(fig);
end

function sp = makeEulerStencilProps(xi)
sp = kp.rbffd.StencilProperties();
sp.dim = 2;
theta = 2;
sp.ell = max(xi + theta - 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(sp.dim, sp.ell), 1);
sp.n = 2 * sp.npoly + 1;
sp.spline_degree = sp.ell;
if mod(sp.spline_degree, 2) == 0
    sp.spline_degree = sp.spline_degree - 1;
end
sp.spline_degree = max(sp.spline_degree, 5);
sp.treeMode = "interior_boundary";
sp.pointSet = "interior_boundary";
end

function dual = makeEulerDomain(h, pressureFraction, minPressureNodes, ~)
t = linspace(0, 2*pi, 512).';
t(end) = [];
boundary = [cos(t), sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(boundary);
surface.buildClosedGeometricModelPS(2, h, size(boundary, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DualNodeDomainGenerator();
generator.generateSmoothDomainNodesAutoPressure(surface, h, pressureFraction, minPressureNodes, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.75, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
dual = generator.createDualNodeDomainDescriptor();
dual.buildStructs();
end

function a = eulerAmplitude(t)
a = 1.0 + 0.2 * sin(2.0 * t) + 0.15 * cos(7.0 * t);
end

function adot = eulerAmplitudeDt(t)
adot = 0.4 * cos(2.0 * t) - 1.05 * sin(7.0 * t);
end

function U = eulerVelocityShape(X)
x = X(:, 1);
y = X(:, 2);
r2 = x.^2 + y.^2;
s = 1.0 - r2;
omega = (s.^2) .* exp(0.25 * r2);
U = zeros(size(X, 1), 2);
U(:, 1) = -omega .* y;
U(:, 2) = omega .* x;
end

function U = eulerVelocityExact(t, X)
U = eulerAmplitude(t) * eulerVelocityShape(X);
end

function p = eulerPressureExact(X)
p = 0.25 * sin(2.0 * X(:, 1) - X(:, 2)) + 0.1 * X(:, 1) .* X(:, 2);
end

function F = eulerForcing(t, X)
a = eulerAmplitude(t);
F = eulerAmplitudeDt(t) * eulerVelocityShape(X);
V = eulerVelocityShape(X);
epsVal = 1.0e-6;
Xpx = X; Xpx(:, 1) = Xpx(:, 1) + epsVal;
Xmx = X; Xmx(:, 1) = Xmx(:, 1) - epsVal;
Xpy = X; Xpy(:, 2) = Xpy(:, 2) + epsVal;
Xmy = X; Xmy(:, 2) = Xmy(:, 2) - epsVal;
dVdx = (eulerVelocityShape(Xpx) - eulerVelocityShape(Xmx)) / (2.0 * epsVal);
dVdy = (eulerVelocityShape(Xpy) - eulerVelocityShape(Xmy)) / (2.0 * epsVal);
F(:, 1) = F(:, 1) + a * a * (V(:, 1) .* dVdx(:, 1) + V(:, 2) .* dVdy(:, 1));
F(:, 2) = F(:, 2) + a * a * (V(:, 1) .* dVdx(:, 2) + V(:, 2) .* dVdy(:, 2));
phase = 2.0 * X(:, 1) - X(:, 2);
F(:, 1) = F(:, 1) + 0.5 * cos(phase) + 0.1 * X(:, 2);
F(:, 2) = F(:, 2) - 0.25 * cos(phase) + 0.1 * X(:, 1);
end

function Xnext = eulerRk4Step(t, X, dt, velocity)
K1 = velocity(t, X);
X2 = X + 0.5 * dt * K1;
K2 = velocity(t + 0.5 * dt, X2);
X3 = X + 0.5 * dt * K2;
K3 = velocity(t + 0.5 * dt, X3);
X4 = X + dt * K3;
K4 = velocity(t + dt, X4);
Xnext = X + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
end

function finalizeReadmeFigure(fig)
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    disableDefaultInteractivity(axs(k));
end
end
