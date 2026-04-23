function [X, info] = generatePoissonNodesInBox(radius, x_min, x_max, varargin)
%GENERATEPOISSONNODESINBOX Poisson disk sampling on an axis-aligned box.
%
%   X = kp.nodes.generatePoissonNodesInBox(radius, x_min, x_max)
%   generates a point cloud in the box [x_min, x_max] using a fixed-radius
%   Bridson-style sampler. The box dimension is inferred from the lengths of
%   x_min and x_max and can be any positive integer.
%
%   X = kp.nodes.generatePoissonNodesInBox(..., 'Seed', seed) makes the
%   output deterministic. The same radius, box, strip count, and seed produce
%   the same node set independent of parfor scheduling.
%
%   [X, info] = ... also returns sampler metadata.

    opts = parseInputs(radius, x_min, x_max, varargin{:});

    if any(opts.x_max <= opts.x_min)
        X = zeros(0, numel(opts.x_min));
        info = emptyInfo(opts);
        return;
    end

    stripBoxes = buildStripBoxes(opts.x_min, opts.x_max, opts.radius, opts.strip_count);
    localClouds = cell(opts.strip_count, 1);
    stripSeeds = uint32(mod(double(opts.base_seed) + 104729 * (0:opts.strip_count-1), 2^32 - 1));

    radius = opts.radius;
    attempts = opts.attempts;
    if opts.use_parallel
        parfor k = 1:opts.strip_count
            sample_min = stripBoxes{k}.sample_min;
            sample_max = stripBoxes{k}.sample_max;
            localClouds{k} = bridsonPoissonBox(radius, sample_min, sample_max, attempts, stripSeeds(k));
        end
    else
        for k = 1:opts.strip_count
            localClouds{k} = bridsonPoissonBox(radius, stripBoxes{k}.sample_min, ...
                stripBoxes{k}.sample_max, attempts, stripSeeds(k));
        end
    end

    merged = mergeLocalClouds(localClouds, opts.x_min, opts.x_max, opts.radius);
    X = merged.points;

    info = struct( ...
        'dimension', numel(opts.x_min), ...
        'radius', opts.radius, ...
        'attempts', opts.attempts, ...
        'seed', opts.seed, ...
        'deterministic', opts.deterministic, ...
        'strip_count', opts.strip_count, ...
        'used_parallel', opts.use_parallel, ...
        'num_candidates', merged.num_candidates, ...
        'num_points', size(X, 1));
end

function opts = parseInputs(radius, x_min, x_max, varargin)
    parser = inputParser();
    parser.addRequired('radius', @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}));
    parser.addRequired('x_min', @(x) validateattributes(x, {'numeric'}, {'vector', 'real', 'finite'}));
    parser.addRequired('x_max', @(x) validateattributes(x, {'numeric'}, {'vector', 'real', 'finite'}));
    parser.addParameter('Attempts', 30, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', '>=', 1}));
    parser.addParameter('Seed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x)));
    parser.addParameter('Deterministic', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
    parser.addParameter('StripCount', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x) && isfinite(x) && x >= 1));
    parser.addParameter('UseParallel', true, @(x) islogical(x) || isnumeric(x));
    parser.parse(radius, x_min, x_max, varargin{:});

    opts = parser.Results;
    opts.x_min = opts.x_min(:).';
    opts.x_max = opts.x_max(:).';
    if numel(opts.x_min) ~= numel(opts.x_max)
        error('kp:nodes:BoxSizeMismatch', 'x_min and x_max must have the same length.');
    end

    if isempty(opts.Seed)
        opts.seed = double(randi(2^31 - 1));
    else
        opts.seed = double(opts.Seed);
    end
    opts.base_seed = uint32(mod(opts.seed, 2^32));

    if isempty(opts.Deterministic)
        opts.deterministic = ~isempty(opts.Seed);
    else
        opts.deterministic = logical(opts.Deterministic);
    end

    if isempty(opts.StripCount)
        opts.strip_count = defaultStripCount(logical(opts.UseParallel), opts.deterministic);
    else
        opts.strip_count = max(1, floor(opts.StripCount));
    end

    opts.radius = opts.radius;
    opts.attempts = opts.Attempts;
    opts.use_parallel = logical(opts.UseParallel) && opts.strip_count > 1 && canUseParfor();
end

function info = emptyInfo(opts)
    info = struct( ...
        'dimension', numel(opts.x_min), ...
        'radius', opts.radius, ...
        'attempts', opts.attempts, ...
        'seed', opts.seed, ...
        'deterministic', opts.deterministic, ...
        'strip_count', opts.strip_count, ...
        'used_parallel', false, ...
        'num_candidates', 0, ...
        'num_points', 0);
end

function stripCount = defaultStripCount(useParallel, deterministic)
    if deterministic
        stripCount = 5;
        return;
    end
    if ~useParallel || ~canUseParfor()
        stripCount = 1;
        return;
    end
    pool = gcp('nocreate');
    if isempty(pool)
        stripCount = 1;
    else
        stripCount = max(1, pool.NumWorkers);
    end
end

function tf = canUseParfor()
    tf = license('test', 'Distrib_Computing_Toolbox');
end

function stripBoxes = buildStripBoxes(x_min, x_max, radius, stripCount)
    dim0 = (x_max(1) - x_min(1)) / stripCount;
    overlap = radius;
    stripBoxes = cell(stripCount, 1);
    for k = 1:stripCount
        core_min = x_min;
        core_max = x_max;
        core_min(1) = x_min(1) + (k - 1) * dim0;
        core_max(1) = x_min(1) + k * dim0;
        sample_min = core_min;
        sample_max = core_max;
        if k > 1
            sample_min(1) = max(x_min(1), sample_min(1) - overlap);
        end
        if k < stripCount
            sample_max(1) = min(x_max(1), sample_max(1) + overlap);
        end
        stripBoxes{k} = struct( ...
            'core_min', core_min, ...
            'core_max', core_max, ...
            'sample_min', sample_min, ...
            'sample_max', sample_max);
    end
end

function out = mergeLocalClouds(localClouds, x_min, x_max, radius)
    dim = numel(x_min);
    totalCount = sum(cellfun(@(x) size(x, 1), localClouds));
    if totalCount == 0
        out = struct('points', zeros(0, dim), 'num_candidates', 0);
        return;
    end

    allPoints = zeros(totalCount, dim);
    row0 = 1;
    for k = 1:numel(localClouds)
        pts = localClouds{k};
        m = size(pts, 1);
        if m == 0
            continue;
        end
        rows = row0:(row0 + m - 1);
        allPoints(rows, :) = pts;
        row0 = row0 + m;
    end
    allPoints = allPoints(1:row0-1, :);
    inBox = all(allPoints >= x_min, 2) & all(allPoints <= x_max, 2);
    allPoints = allPoints(inBox, :);
    allPoints = sortrows(allPoints);
    out.points = acceptByGlobalGrid(allPoints, x_min, x_max, radius);
    out.num_candidates = size(allPoints, 1);
end

function X = bridsonPoissonBox(radius, x_min, x_max, attempts, seed)
    dim = numel(x_min);
    cellSize = radius / sqrt(dim);
    gridSize = max(ones(1, dim), ceil((x_max - x_min) / cellSize));
    grid = containers.Map('KeyType', 'char', 'ValueType', 'uint32');
    points = zeros(128, dim);
    active = zeros(128, 1, 'uint32');

    stream = RandStream('mt19937ar', 'Seed', double(seed));
    x0 = x_min + rand(stream, 1, dim) .* (x_max - x_min);
    [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, x0, x_min, cellSize, 0, 0);

    while nActive > 0
        pick = randi(stream, nActive);
        activeIdx = active(pick);
        base = points(activeIdx, :);
        accepted = false;
        for attempt = 1:attempts
            candidate = proposeCandidate(base, radius, dim, stream);
            if any(candidate < x_min) || any(candidate > x_max)
                continue;
            end
            if hasNeighbor(candidate, radius, points, grid, x_min, cellSize, gridSize)
                continue;
            end
            [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, candidate, x_min, cellSize, nPoints, nActive);
            accepted = true;
            break;
        end
        if ~accepted
            active(pick) = active(nActive);
            nActive = nActive - 1;
        end
    end

    X = points(1:nPoints, :);
end

function candidate = proposeCandidate(base, radius, dim, stream)
    direction = randn(stream, 1, dim);
    direction = direction ./ max(norm(direction, 2), eps);
    shellRadius = radius * (1 + rand(stream, 1) * (2^dim - 1))^(1 / dim);
    candidate = base + shellRadius * direction;
end

function [points, active, grid, nPoints, nActive] = addPoint(points, active, grid, point, x_min, cellSize, nPoints, nActive)
    if nPoints == 0
        nPoints = 1;
    else
        nPoints = nPoints + 1;
    end
    if nPoints > size(points, 1)
        points = [points; zeros(size(points, 1), size(points, 2))];
    end
    points(nPoints, :) = point;

    nActive = nActive + 1;
    if nActive > numel(active)
        active = [active; zeros(numel(active), 1, 'uint32')];
    end
    active(nActive) = uint32(nPoints);

    idx = pointToCell(point, x_min, cellSize);
    grid(cellKey(idx)) = uint32(nPoints);
end

function tf = hasNeighbor(point, radius, points, grid, x_min, cellSize, gridSize)
    idx = pointToCell(point, x_min, cellSize);
    reach = max(1, ceil(radius / cellSize));
    ranges = cell(1, numel(idx));
    for d = 1:numel(idx)
        ranges{d} = max(1, idx(d) - reach):min(gridSize(d), idx(d) + reach);
    end
    tf = false;
    neighborCells = enumerateNeighborCells(ranges);
    for k = 1:size(neighborCells, 1)
        key = cellKey(neighborCells(k, :));
        if ~isKey(grid, key)
            continue;
        end
        j = grid(key);
        if norm(point - points(j, :), 2) < radius
            tf = true;
            return;
        end
    end
end

function idx = pointToCell(point, x_min, cellSize)
    idx = floor((point - x_min) ./ cellSize) + 1;
    idx = max(idx, 1);
end

function key = cellKey(idx)
    key = sprintf('%d_', idx);
end

function cells = enumerateNeighborCells(ranges)
    dim = numel(ranges);
    grids = cell(1, dim);
    [grids{:}] = ndgrid(ranges{:});
    n = numel(grids{1});
    cells = zeros(n, dim);
    for d = 1:dim
        cells(:, d) = grids{d}(:);
    end
end

function accepted = acceptByGlobalGrid(points, x_min, x_max, radius)
    dim = size(points, 2);
    if isempty(points)
        accepted = zeros(0, dim);
        return;
    end

    cellSize = radius / sqrt(dim);
    gridSize = max(ones(1, dim), ceil((x_max - x_min) / cellSize));
    grid = containers.Map('KeyType', 'char', 'ValueType', 'uint32');
    keep = false(size(points, 1), 1);
    accepted = zeros(size(points));
    nAccepted = 0;
    for i = 1:size(points, 1)
        pt = points(i, :);
        if hasNeighbor(pt, radius, accepted, grid, x_min, cellSize, gridSize)
            continue;
        end
        nAccepted = nAccepted + 1;
        accepted(nAccepted, :) = pt;
        idx = pointToCell(pt, x_min, cellSize);
        grid(cellKey(idx)) = uint32(nAccepted);
        keep(i) = true;
    end
    accepted = accepted(1:nAccepted, :);
end
