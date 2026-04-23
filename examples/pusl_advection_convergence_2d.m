function results = pusl_advection_convergence_2d()
%PUSL_ADVECTION_CONVERGENCE_2D Convergence study for backward PU-SL advection.

common = pusl_advection_common();
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
    'inf', 0, ...
    'mass', 0), numel(orders), numel(hValues));

for io = 1:numel(orders)
    xi = orders(io);
    for ih = 1:numel(hValues)
        h = hValues(ih);
        dt = min(0.01, 0.5 * h);

        domain = common.create_disk_domain(h);
        solver = kp.solvers.PUSLAdvectionSolver();
        solver.init(domain, xi, dt);

        c = solver.projectInitial(@(x) common.smooth_initial_value(x));
        m0 = solver.totalMass(c);

        nsteps = max(ceil(T / max(dt, 1.0e-14)), 1);
        local_dt = T / nsteps;
        solver.setStepSize(local_dt);
        t = 0.0;
        for step = 1:nsteps
            c = solver.backwardSLStep(t, c, ...
                @(time, X) common.velocity_disk(time, X), ...
                @(time, X, dlt, vel) rk4sl(time, X, dlt, vel));
            t = t + local_dt;
        end

        q = solver.evaluateAtNodes(c);
        Xout = solver.getOutputNodes();
        qex = common.exact_solution(T, Xout);
        m1 = solver.totalMass(c);

        row = struct();
        row.order = xi;
        row.h = h;
        row.dt = local_dt;
        row.nPhysical = size(Xout, 1);
        row.rel = norm(q - qex) / norm(qex);
        row.inf = norm(q - qex, inf);
        row.mass = abs(m1 - m0);
        results.rows(io, ih) = row;
    end
end

results.rates = estimateRates(results.rows);
printResults(results);
end

function Xnext = rk4sl(t, X, dt, velocity)
K1 = velocity(t, X);
K2 = velocity(t + 0.5 * dt, X + 0.5 * dt * K1);
K3 = velocity(t + 0.5 * dt, X + 0.5 * dt * K2);
K4 = velocity(t + dt, X + dt * K3);
Xnext = X + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
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
fprintf('\n2D backward PU-SL advection convergence study\n');
for io = 1:size(results.rows, 1)
    fprintf('Order %d\n', results.rows(io, 1).order);
    fprintf('  h        dt        Nphys      rel error       inf error       rel rate    inf rate    mass drift\n');
    for ih = 1:size(results.rows, 2)
        row = results.rows(io, ih);
        if ih == 1
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10s  %-10s  %-10.3e\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, '-', '-', row.mass);
        else
            fprintf('  %-7.3f  %-8.5f  %-9d  %-14.6e  %-14.6e  %-10.4f  %-10.4f  %-10.3e\n', ...
                row.h, row.dt, row.nPhysical, row.rel, row.inf, ...
                results.rates(io).rel(ih - 1), results.rates(io).inf(ih - 1), row.mass);
        end
    end
    fprintf('\n');
end
end
