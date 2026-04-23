function pusl_pu_advection_diffusion_checks()
%PUSL_PU_ADVECTION_DIFFUSION_CHECKS Basic checks for PU-SL + PU diffusion.

common = pusl_advection_diffusion_common();
domain = common.create_disk_domain(0.16);

solver = kp.solvers.PUSLPUAdvectionDiffusionSolver();
solver.init(domain, 4, 4, 0.01, 0.05, "backward");
solver.enableHomogeneousNeumannMassConservation(true);

fdSolver = kp.solvers.PUSLFDAdvectionDiffusionSolver();
fdSolver.init(domain, 4, 4, 0.01, 0.05, "backward");
fdSolver.enableHomogeneousNeumannMassConservation(true);

Xout = solver.getOutputNodes();
u0 = common.exact_solution(0.0, Xout);
solver.setInitialState(u0);
fdSolver.setInitialState(u0);

u1 = solver.bdf1Step(0.01, ...
    @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
    @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
    @(t, X) common.unit_coeff(t, X), ...
    @(t, X) common.zero_coeff(t, X), ...
    @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));

u1fd = fdSolver.bdf1Step(0.01, ...
    @(t, X) common.boundary_zero_swirl(domain, t, X), [], ...
    @(nuValue, t, X) common.forcing(domain, nuValue, t, X), ...
    @(t, X) common.unit_coeff(t, X), ...
    @(t, X) common.zero_coeff(t, X), ...
    @(NeuCoeffs, DirCoeffs, nr, t, X) common.homogeneous_neumann(NeuCoeffs, DirCoeffs, nr, t, X));

uTrue = common.exact_solution(0.01, Xout);
assert(norm(u1 - uTrue) / norm(uTrue) < 5e-1, ...
    'PU-SL + PU diffusion BDF1 should track the manufactured solution to modest accuracy.');
assert(norm(u1 - u1fd) / norm(u1fd) < 3e-1, ...
    'PU-SL + PU diffusion should stay reasonably close to the FD-diffusion coupled path on the same test.');

disp('pusl pu advection diffusion checks passed');
end
