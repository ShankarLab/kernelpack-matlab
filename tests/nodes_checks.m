function nodes_checks()
%NODES_CHECKS Lightweight checks for box Poisson disk sampling.

[x2a, info2a] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 19, 'StripCount', 5);
[x2b, info2b] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 19, 'StripCount', 5);
[x2c, ~] = kp.nodes.generatePoissonNodesInBox(0.075, [0 0], [1 1], ...
    'Seed', 23, 'StripCount', 5);

assert(isequal(x2a, x2b), 'Seeded 2D box sampling should be deterministic.');
assert(~isequal(x2a, x2c), 'Different seeds should usually produce different 2D clouds.');
assert(all(x2a >= 0, 'all') && all(x2a <= 1, 'all'), '2D nodes must stay in the box.');
assert(minPairDistance(x2a) >= 0.075 * (1 - 1e-10), '2D nodes must respect the exclusion radius.');
assert(info2a.deterministic && info2b.deterministic, 'Seeded runs should report deterministic mode.');

[x3, info3] = kp.nodes.generatePoissonNodesInBox(0.16, [0 0 0], [1 1 1], ...
    'Seed', 31, 'StripCount', 5);
assert(size(x3, 2) == 3, '3D nodes should have three coordinates.');
assert(minPairDistance(x3) >= 0.16 * (1 - 1e-10), '3D nodes must respect the exclusion radius.');
assert(info3.strip_count == 5, 'Explicit strip count should be preserved.');

[x4, ~] = kp.nodes.generatePoissonNodesInBox(0.35, zeros(1, 4), ones(1, 4), ...
    'Seed', 5, 'StripCount', 3, 'UseParallel', false);
assert(size(x4, 2) == 4, 'Sampler should work in dimensions beyond 3.');
assert(minPairDistance(x4) >= 0.35 * (1 - 1e-10), '4D nodes must respect the exclusion radius.');

gen = kp.nodes.DomainNodeGenerator();
gen.generatePoissonNodes(0.1, [0 0], [1 1], 'Seed', 13, 'StripCount', 4);
assert(~isempty(gen.getRawPoissonInteriorNodes()), 'DomainNodeGenerator should store the raw Poisson cloud.');

disp('nodes checks passed');
end

function dmin = minPairDistance(X)
    if size(X, 1) < 2
        dmin = inf;
        return;
    end
    D = kp.geometry.distanceMatrix(X, X);
    D(1:size(D, 1)+1:end) = inf;
    dmin = min(D, [], 'all');
end
