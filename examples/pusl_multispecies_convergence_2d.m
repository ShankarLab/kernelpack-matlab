function results = pusl_multispecies_convergence_2d()
%PUSL_MULTISPECIES_CONVERGENCE_2D Convergence study for the multi-species PU-SL wrapper.

orders = [2, 4, 6];
hValues = [0.20, 0.14, 0.10];
T = 0.25;

results = struct();
results.orders = orders;
results.h = hValues;
results.rows = repmat(struct( ...
    'order', 0, ...
    'h', 0, ...
    'dt', 0, ...
    'nPhysical', 0, ...
    'rel', 0, ...
    'inf', 0), numel(orders), numel(hValues));

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hValues)
        h = hValues(ih);
        dt = 0.1 * h;
        domain = create_disk_domain(h);

        solver = kp.solvers.MultiSpeciesPUSLAdvectionSolver();
        solver.init(domain, xi, dt, 2);

        c = solver.projectInitial(@initial_values_row);
        nsteps = max(ceil(T / max(dt, 1.0e-14)), 1);
        local_dt = T / nsteps;
        solver.setStepSize(local_dt);
        t = 0.0;
        for step = 1:nsteps
            c = solver.backwardSLStep(t, c, @velocity_disk, []);
            t = t + local_dt;
        end

        q = solver.evaluateAtNodes(c);
        Xout = solver.getOutputNodes();
        qex = exact_solution(T, Xout);

        err = q - qex;
        row = struct();
        row.order = xi;
        row.h = h;
        row.dt = local_dt;
        row.nPhysical = size(Xout, 1);
        row.rel = norm(err(:)) / norm(qex(:));
        row.inf = norm(err(:), inf);
        results.rows(io, ih) = row;
    end
end

results.rates = estimateRates(results.rows);
printResults(results);
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

function rates = estimateRates(rows)
rates = struct('order', [], 'rel', [], 'inf', []);
for io = 1:size(rows, 1)
    h = [rows(io, :).h];
    rel = [rows(io, :).rel];
    infv = [rows(io, :).inf];
    rateRow = struct();
    rateRow.order = rows(io, 1).order;
    rateRow.rel = NaN(1, numel(h) - 1);
    rateRow.inf = NaN(1, numel(h) - 1);
    for k = 1:numel(h) - 1
        rateRow.rel(k) = log(rel(k) / rel(k + 1)) / log(h(k) / h(k + 1));
        rateRow.inf(k) = log(infv(k) / infv(k + 1)) / log(h(k) / h(k + 1));
    end
    rates(io) = rateRow; %#ok<AGROW>
end
end

function printResults(results)
fprintf('\n2D rigid-rotation multi-species PU-SL convergence study\n');
for io = 1:size(results.rows, 1)
    fprintf('Order %d\n', results.rows(io, 1).order);
    fprintf('  h        dt        Nphys      rel error       inf error       rel rate    inf rate\n');
    for ih = 1:size(results.rows, 2)
        row = results.rows(io, ih);
        if ih == 1
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10s  %-10s\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, '-', '-');
        else
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10.4f  %-10.4f\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, ...
                results.rates(io).rel(ih - 1), results.rates(io).inf(ih - 1));
        end
    end
    fprintf('\n');
end
end
