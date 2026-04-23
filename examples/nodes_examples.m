function nodes_examples()
%NODES_EXAMPLES Minimal examples for box Poisson disk sampling.

fprintf('Running 2D seeded box Poisson example...\n');
[x2, info2] = kp.nodes.generatePoissonNodesInBox(0.06, [0 0], [1 1], ...
    'Seed', 7, 'StripCount', 5);
fprintf('2D nodes: %d (seed %d, strips %d)\n', size(x2, 1), info2.seed, info2.strip_count);

fprintf('Running 4D seeded box Poisson example...\n');
[x4, info4] = kp.nodes.generatePoissonNodesInBox(0.22, zeros(1, 4), ones(1, 4), ...
    'Seed', 11, 'StripCount', 5, 'UseParallel', false);
fprintf('4D nodes: %d (seed %d, strips %d)\n', size(x4, 1), info4.seed, info4.strip_count);
