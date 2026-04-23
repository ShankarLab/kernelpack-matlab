function pusl_advection_checks()
%PUSL_ADVECTION_CHECKS Basic checks for the MATLAB PU-SL advection solver.

common = pusl_advection_common();
domain = common.create_disk_domain(0.16);

solver = kp.solvers.PUSLAdvectionSolver();
solver.init(domain, 4, 0.01);

% Constant preservation sanity check.
c0 = solver.projectConstant(2.0, 1);
c1 = solver.backwardSLStep(0.0, c0, ...
    @(t, X) common.velocity_disk(t, X), ...
    @(t, X, dt, velocity) rk4sl(t, X, dt, velocity));
assert(max(abs(c1 - 2.0)) < 5.0e-4, ...
    'Backward PU-SL step should keep constants nearly invariant for this tangential flow case.');

% One manufactured solve sanity check.
c0 = solver.projectInitial(@(x) common.smooth_initial_value(x));
c1 = solver.backwardSLStep(0.0, c0, ...
    @(t, X) common.velocity_disk(t, X), ...
    @(t, X, dt, velocity) rk4sl(t, X, dt, velocity));
q = solver.evaluateAtNodes(c1);
qex = common.exact_solution(0.01, solver.getOutputNodes());
assert(norm(q - qex) / norm(qex) < 2.5e-1, ...
    'Backward PU-SL step should give a modestly accurate one-step update.');

disp('pusl advection checks passed');
end

function Xnext = rk4sl(t, X, dt, velocity)
K1 = velocity(t, X);
K2 = velocity(t + 0.5 * dt, X + 0.5 * dt * K1);
K3 = velocity(t + 0.5 * dt, X + 0.5 * dt * K2);
K4 = velocity(t + dt, X + dt * K3);
Xnext = X + (dt / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
end
