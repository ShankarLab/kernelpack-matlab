function out = pusl_advection_diffusion_common()
%PUSL_ADVECTION_DIFFUSION_COMMON Shared manufactured problem for upcoming PU-SL ports.
%
% This mirrors the common manufactured problem used by the KernelPack
% advection-diffusion examples. It is intentionally kept independent of any
% solver implementation so it can be reused for convergence studies and math
% sanity checks as the PU-SL solvers are migrated.

out = struct();
out.phi = @phi;
out.phi_x = @phi_x;
out.phi_y = @phi_y;
out.lap_phi = @lap_phi;
out.exact_solution = @exact_solution;
out.boundary_zero_swirl = @boundary_zero_swirl;
out.forcing = @forcing;
out.homogeneous_neumann = @homogeneous_neumann;
out.unit_coeff = @unit_coeff;
out.zero_coeff = @zero_coeff;
out.create_disk_domain = @create_disk_domain;
out.create_sphere_domain = @create_sphere_domain;
end

function values = phi(X)
r2 = sum(X.^2, 2);
values = X(:, 1) .* (1 - r2).^2;
end

function values = phi_x(X)
x = X(:, 1);
r2 = sum(X.^2, 2);
s = 1 - r2;
values = s.^2 - 4 * x.^2 .* s;
end

function values = phi_y(X)
x = X(:, 1);
y = X(:, 2);
r2 = sum(X.^2, 2);
s = 1 - r2;
values = -4 * x .* y .* s;
end

function values = lap_phi(X)
r2 = sum(X.^2, 2);
dim = size(X, 2);
values = 4 * X(:, 1) .* ((dim + 4) * r2 - (dim + 2));
end

function values = exact_solution(t, X)
values = ones(size(X, 1), 1) + exp(-t) .* phi(X);
end

function U = boundary_zero_swirl(~, t, X) %#ok<INUSD>
s = boundary_vanishing_factor(X);
U = zeros(size(X));
U(:, 1) = -4 * X(:, 2) .* s;
U(:, 2) = 4 * X(:, 1) .* s;
end

function values = forcing(~, nu, t, X)
decay = exp(-t);
U = boundary_zero_swirl([], t, X);
adv = U(:, 1) .* phi_x(X) + U(:, 2) .* phi_y(X);
values = decay .* (-phi(X) + adv - nu .* lap_phi(X));
end

function values = homogeneous_neumann(~, ~, nr, t, X) %#ok<INUSD>
values = zeros(size(X, 1), 1);
end

function values = unit_coeff(t, X) %#ok<INUSD>
values = ones(size(X, 1), 1);
end

function values = zero_coeff(t, X) %#ok<INUSD>
values = zeros(size(X, 1), 1);
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
    'StripCount', 5);
end

function domain = create_sphere_domain(h)
data_sites = kp.geometry.fibonacciSphere(800);

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(data_sites);
surface.buildClosedGeometricModelPS(3, h, size(data_sites, 1));
surface.buildLevelSetFromGeometricModel([]);

generator = kp.nodes.DomainNodeGenerator();
domain = generator.buildDomainDescriptorFromGeometry(surface, h, ...
    'Seed', 17, ...
    'StripCount', 6);
end

function s = boundary_vanishing_factor(X)
r2 = sum(X.^2, 2);
s = 1 - r2;
end
