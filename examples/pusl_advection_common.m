function out = pusl_advection_common()
%PUSL_ADVECTION_COMMON Shared manufactured problem for PU-SL advection.

out.smooth_initial_value = @smooth_initial_value;
out.velocity_disk = @velocity_disk;
out.exact_solution = @exact_solution;
out.create_disk_domain = @create_disk_domain;
end

function values = smooth_initial_value(X)
z = zeros(size(X, 1), 1);
if size(X, 2) > 2
    z = X(:, 3);
end
values = exp(0.3 * X(:, 1) - 0.2 * X(:, 2) + 0.1 * z) + ...
    sin(0.7 * X(:, 1) + 0.4 * X(:, 2) - 0.2 * z);
end

function U = velocity_disk(t, X)
r2 = sum(X.^2, 2);
scale = sin(pi * r2) .* sin(pi * t);
U = zeros(size(X));
U(:, 1) = X(:, 2) .* scale;
U(:, 2) = -X(:, 1) .* scale;
end

function q = exact_solution(t, X)
r2 = sum(X.^2, 2);
theta = sin(pi * r2) .* ((1 - cos(pi * t)) / pi);
Xd = zeros(size(X));
Xd(:, 1) = X(:, 1) .* cos(theta) - X(:, 2) .* sin(theta);
Xd(:, 2) = X(:, 1) .* sin(theta) + X(:, 2) .* cos(theta);
if size(X, 2) > 2
    Xd(:, 3:end) = X(:, 3:end);
end
q = smooth_initial_value(Xd);
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
