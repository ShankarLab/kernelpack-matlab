function pusl_multispecies_checks()
%PUSL_MULTISPECIES_CHECKS Basic checks for the multi-species PU-SL wrapper.

domain = create_disk_domain(0.16);

solver = kp.solvers.MultiSpeciesPUSLAdvectionSolver();
solver.init(domain, 4, 0.01, 2);

% Constant preservation per species.
c0 = solver.projectConstants([2.0, -1.0]);
c1 = solver.backwardSLStep(0.0, c0, @velocity_disk, []);
assert(max(abs(c1(:, 1) - 2.0)) < 5.0e-4);
assert(max(abs(c1(:, 2) + 1.0)) < 5.0e-4);

% Joint solve should agree with species-by-species exact transport.
c0 = solver.projectInitial(@initial_values_row);
c1 = solver.backwardSLStep(0.0, c0, @velocity_disk, []);
q = solver.evaluateAtNodes(c1);
qex = exact_solution(0.01, solver.getOutputNodes());
assert(norm(q(:) - qex(:)) / norm(qex(:)) < 3.0e-1, ...
    'Multi-species PU-SL step should track the manufactured solution to modest accuracy.');

disp('pusl multispecies checks passed');
end

function values = initial_values_row(x)
values = [exp(0.3 * x(1) - 0.2 * x(2)) + sin(0.7 * x(1) + 0.4 * x(2)), ...
    cos(0.5 * x(1) - 0.6 * x(2)) + 0.2 * x(1)];
end

function values = initial_values(X)
values = [exp(0.3 * X(:, 1) - 0.2 * X(:, 2)) + sin(0.7 * X(:, 1) + 0.4 * X(:, 2)), ...
    cos(0.5 * X(:, 1) - 0.6 * X(:, 2)) + 0.2 * X(:, 1)];
end

function U = velocity_disk(~, X)
U = zeros(size(X));
U(:, 1) = X(:, 2);
U(:, 2) = -X(:, 1);
end

function q = exact_solution(t, X)
Xd = zeros(size(X));
Xd(:, 1) = X(:, 1) * cos(t) - X(:, 2) * sin(t);
Xd(:, 2) = X(:, 1) * sin(t) + X(:, 2) * cos(t);
q = initial_values(Xd);
end

function domain = create_disk_domain(h)
t = linspace(0, 2 * pi, 101).';
t(end) = [];
data_sites = [cos(t), sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(data_sites);
surface.buildClosedGeometricModelPS(2, h, size(data_sites, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, h, ...
    'Seed', 17, ...
    'StripCount', 5, ...
    'DoOuterRefinement', true, ...
    'OuterFractionOfh', 0.75, ...
    'OuterRefinementZoneSizeAsMultipleOfh', 2.0);
end
