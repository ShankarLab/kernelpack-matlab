function pusl_incompressible_euler_example()
%PUSL_INCOMPRESSIBLE_EULER_EXAMPLE Small manufactured PU-SL Euler demo on the disk.

threads = 1;
xi_u = 4;
xi_p = 4;
h = 0.10;
dt = 0.02;
final_time = 0.04;

velocity_sp = make_stencil_properties(xi_u);
pressure_sp = make_stencil_properties(xi_p);
dual = make_domain(h, 0.40, pressure_sp.n, threads);

solver = kp.solvers.PUSLIncompressibleEulerSolver();
solver.init(dual, xi_u, velocity_sp, pressure_sp, dt, threads);
solver.setTangentialFlowBoundary(1.0e-5);

Xu = dual.getVelocityDomain().getIntBdryNodes();
problem = kp.solvers.IncompressibleEulerSolver.defaultProblemDefinition();
problem.slip_walls = {kp.solvers.IncompressibleEulerSolver.stationarySlipWall((1:dual.getVelocityDomain().getNumBdryNodes()).')};
problem.gauge_options.mode = "forcepressuremean";

u0 = velocity_exact(0.0, Xu);
solver.setInitialVelocity(u0, problem);
sol = solver.bdf1Step(dt, @rk4_step, @euler_forcing, problem);
sol = solver.bdf2Step(final_time, @rk4_step, @euler_forcing, problem);

u_exact = velocity_exact(final_time, Xu);
Xp = dual.getPressureDomain().getIntBdryNodes();
p_exact = pressure_exact(Xp);
p = sol.pressure - mean(sol.pressure - p_exact);

fprintf('PUSL incompressible Euler example\\n');
fprintf('  velocity nodes: %d\\n', size(Xu, 1));
fprintf('  pressure nodes: %d\\n', size(Xp, 1));
fprintf('  relative velocity error: %.6e\\n', relative_l2(sol.velocity, u_exact));
fprintf('  relative pressure error: %.6e\\n', relative_l2_vec(p, p_exact));
fprintf('  divergence rms/max: %.6e / %.6e\\n', sol.divergence_rms, sol.divergence_max);
fprintf('  wall rms/max: %.6e / %.6e\\n', sol.wall_normal_rms, sol.wall_normal_max);
end

function sp = make_stencil_properties(xi)
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

function dual = make_domain(h, pressure_fraction, min_pressure_nodes, ~)
t = linspace(0, 2*pi, 512).';
t(end) = [];
boundary = [cos(t), sin(t)];
surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(boundary);
surface.buildClosedGeometricModelPS(2, h, size(boundary, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DualNodeDomainGenerator();
generator.generateSmoothDomainNodesAutoPressure(surface, h, pressure_fraction, min_pressure_nodes, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.75, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
dual = generator.createDualNodeDomainDescriptor();
dual.buildStructs();
end

function a = amplitude(t)
a = 1.0 + 0.2 * sin(2.0 * t) + 0.15 * cos(7.0 * t);
end

function adot = amplitude_dt(t)
adot = 0.4 * cos(2.0 * t) - 1.05 * sin(7.0 * t);
end

function U = velocity_shape(X)
x = X(:, 1);
y = X(:, 2);
r2 = x.^2 + y.^2;
s = 1.0 - r2;
omega = (s.^2) .* exp(0.25 * r2);
U = zeros(size(X, 1), 2);
U(:, 1) = -omega .* y;
U(:, 2) = omega .* x;
end

function U = velocity_exact(t, X)
U = amplitude(t) * velocity_shape(X);
end

function p = pressure_exact(X)
p = 0.25 * sin(2.0 * X(:, 1) - X(:, 2)) + 0.1 * X(:, 1) .* X(:, 2);
end

function F = euler_forcing(t, X)
a = amplitude(t);
F = amplitude_dt(t) * velocity_shape(X);
V = velocity_shape(X);
epsVal = 1.0e-6;
Xpx = X; Xpx(:, 1) = Xpx(:, 1) + epsVal;
Xmx = X; Xmx(:, 1) = Xmx(:, 1) - epsVal;
Xpy = X; Xpy(:, 2) = Xpy(:, 2) + epsVal;
Xmy = X; Xmy(:, 2) = Xmy(:, 2) - epsVal;
dVdx = (velocity_shape(Xpx) - velocity_shape(Xmx)) / (2.0 * epsVal);
dVdy = (velocity_shape(Xpy) - velocity_shape(Xmy)) / (2.0 * epsVal);
F(:, 1) = F(:, 1) + a * a * (V(:, 1) .* dVdx(:, 1) + V(:, 2) .* dVdy(:, 1));
F(:, 2) = F(:, 2) + a * a * (V(:, 1) .* dVdx(:, 2) + V(:, 2) .* dVdy(:, 2));
phase = 2.0 * X(:, 1) - X(:, 2);
F(:, 1) = F(:, 1) + 0.5 * cos(phase) + 0.1 * X(:, 2);
F(:, 2) = F(:, 2) - 0.25 * cos(phase) + 0.1 * X(:, 1);
end

function r = relative_l2(U, V)
r = norm(U - V, 'fro') / max(norm(V, 'fro'), 1.0e-14);
end

function r = relative_l2_vec(u, v)
r = norm(u - v) / max(norm(v), 1.0e-14);
end

function Xnext = rk4_step(t, X, dt, velocity)
K1 = velocity(t, X);
X2 = X + 0.5 * dt * K1;
K2 = velocity(t + 0.5 * dt, X2);
X3 = X + 0.5 * dt * K2;
K3 = velocity(t + 0.5 * dt, X3);
X4 = X + dt * K3;
K4 = velocity(t + dt, X4);
Xnext = X + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
end
