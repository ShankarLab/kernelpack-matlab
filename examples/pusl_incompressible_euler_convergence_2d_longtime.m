function results = pusl_incompressible_euler_convergence_2d_longtime(xi_u_vals, xi_p_vals)
%PUSL_INCOMPRESSIBLE_EULER_CONVERGENCE_2D_LONGTIME Longer-time serial convergence study.

threads = 1;
if nargin < 1 || isempty(xi_u_vals)
    xi_u_vals = [2, 4, 6];
end
if nargin < 2 || isempty(xi_p_vals)
    xi_p_vals = xi_u_vals;
end
xi_u_vals = xi_u_vals(:);
xi_p_vals = xi_p_vals(:);
assert(numel(xi_u_vals) == numel(xi_p_vals), ...
    'xi_u_vals and xi_p_vals must have the same number of entries.');

hvals = [0.10, 0.08, 0.07, 0.06];
final_time = 0.20;

results = struct();
for xiIdx = 1:numel(xi_u_vals)
    xi_u = xi_u_vals(xiIdx);
    xi_p = xi_p_vals(xiIdx);
    caseResults = cell(numel(hvals), 1);
    fprintf('xi_u = %d, xi_p = %d\n', xi_u, xi_p);
    for hIdx = 1:numel(hvals)
        h = hvals(hIdx);
        caseResult = run_case(xi_u, xi_p, h, final_time, threads);
        caseResults{hIdx} = caseResult;
        fprintf('  h=%.3f  dt=%.6f  steps=%d  Nu=%d  Np=%d  u_rel=%.6e  p_rel=%.6e  div_rms=%.6e\n', ...
            caseResult.h, caseResult.dt, caseResult.num_steps, ...
            caseResult.Nu, caseResult.Np, caseResult.u_rel, caseResult.p_rel, caseResult.div_rms);
    end
    results.(pairFieldName(xi_u, xi_p)) = vertcat(caseResults{:});
end

imageDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', 'images');
if ~exist(imageDir, 'dir')
    mkdir(imageDir);
end
plotPath = fullfile(imageDir, 'pusl_incompressible_euler_convergence_2d_longtime.png');

fig = figure('Color', 'w', 'Position', [100, 100, 760, 520]);
hold on;
styles = {'-o', '-s', '-d'};
for xiIdx = 1:numel(xi_u_vals)
    xi_u = xi_u_vals(xiIdx);
    xi_p = xi_p_vals(xiIdx);
    R = results.(pairFieldName(xi_u, xi_p));
    plot([R.h], [R.u_rel], styles{xiIdx}, 'LineWidth', 1.4, 'MarkerSize', 7, ...
        'DisplayName', sprintf('\\xi_u = %d, \\xi_p = %d', xi_u, xi_p));
end
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('h');
ylabel('relative velocity L2 error');
title(sprintf('PUSL Incompressible Euler on disk, T = %.2f', final_time));
legend('Location', 'southwest');
exportgraphics(fig, plotPath, 'Resolution', 180);
close(fig);

save(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'pusl_incompressible_euler_convergence_2d_longtime_results.mat'), ...
    'results');

disp('PUSL incompressible Euler longer-time convergence study');
for xiIdx = 1:numel(xi_u_vals)
    xi_u = xi_u_vals(xiIdx);
    xi_p = xi_p_vals(xiIdx);
    R = results.(pairFieldName(xi_u, xi_p));
    fprintf('xi_u = %d, xi_p = %d\n', xi_u, xi_p);
    for k = 1:numel(R)
        if k == 1
            rateText = '-';
        else
            rate = log(R(k-1).u_rel / R(k).u_rel) / log(R(k-1).h / R(k).h);
            rateText = sprintf('%.4f', rate);
        end
        fprintf('  h=%.3f  dt=%.6f  steps=%d  Nu=%d  Np=%d  u_rel=%.6e  p_rel=%.6e  rate=%s\n', ...
            R(k).h, R(k).dt, R(k).num_steps, R(k).Nu, R(k).Np, R(k).u_rel, R(k).p_rel, rateText);
    end
end
fprintf('saved: %s\n', plotPath);
end

function result = run_case(xi_u, xi_p, h, final_time, threads)
velocity_sp = make_stencil_properties(xi_u);
pressure_sp = make_stencil_properties(xi_p);
dual = make_domain(h, max(0.40, min(0.70, 1.15 * pressure_sp.n / 220)), pressure_sp.n, threads);

dt_target = h^(xi_u / 3);
num_steps = max(2, ceil(final_time / dt_target));
dt = final_time / num_steps;

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
if num_steps >= 2
    sol = solver.bdf2Step(2.0 * dt, @rk4_step, @euler_forcing, problem);
end
for step = 3:num_steps
    sol = solver.bdf3Step(step * dt, @rk4_step, @euler_forcing, problem);
end

u_exact = velocity_exact(final_time, Xu);
Xp = dual.getPressureDomain().getIntBdryNodes();
p_exact = pressure_exact(Xp);
pressure = sol.pressure - mean(sol.pressure - p_exact);

u_num = solver.advectionSolver().gatherOutputSamples(sol.velocity);
result = struct();
result.xi_u = xi_u;
result.xi_p = xi_p;
result.h = h;
result.dt = dt;
result.num_steps = num_steps;
result.Nu = size(Xu, 1);
result.Np = size(Xp, 1);
result.u_rel = relative_l2(u_num, u_exact);
result.p_rel = relative_l2_vec(pressure, p_exact);
result.div_rms = sol.divergence_rms;
result.div_max = sol.divergence_max;
result.wall_rms = sol.wall_normal_rms;
result.wall_max = sol.wall_normal_max;
end

function name = pairFieldName(xi_u, xi_p)
name = sprintf('xiu%d_xip%d', xi_u, xi_p);
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
