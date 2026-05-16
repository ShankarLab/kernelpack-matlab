classdef PUSLIncompressibleEulerSolver < handle
    %PUSLINCOMPRESSIBLEEULERSOLVER PU-SL wrapper around the MATLAB Euler backend.

    properties
        Domain kp.domain.DualNodeDomainDescriptor = kp.domain.DualNodeDomainDescriptor()
        xi_sl (1,1) double = 0
        dt (1,1) double = 0
        num_omp_threads (1,1) double = 1
        Advection kp.solvers.PUSLAdvectionSolver = kp.solvers.PUSLAdvectionSolver()
        Backend kp.solvers.IncompressibleEulerSolver = kp.solvers.IncompressibleEulerSolver()
        AdvectionDomain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        physical_nodes double = zeros(0, 0)
        global_physical_nodes double = zeros(0, 0)
        interior_nodes double = zeros(0, 0)
        coeff_nm2 double = zeros(0, 0)
        coeff_nm1 double = zeros(0, 0)
        coeff_n double = zeros(0, 0)
        completed_steps_ (1,1) double = 0
    end

    methods
        function clearAdvectionBoundaryCondition(obj)
            obj.Advection.clearAdvectionBoundaryCondition();
        end

        function setTangentialFlowBoundary(obj, normal_velocity_tolerance)
            if nargin < 2 || isempty(normal_velocity_tolerance)
                normal_velocity_tolerance = 1.0e-10;
            end
            obj.Advection.setTangentialFlowBoundary(normal_velocity_tolerance);
        end

        function setInflowDirichletBoundary(obj, inflow_value, normal_velocity_tolerance)
            if nargin < 3 || isempty(normal_velocity_tolerance)
                normal_velocity_tolerance = 1.0e-10;
            end
            obj.Advection.setInflowDirichletBoundary(inflow_value, normal_velocity_tolerance);
        end

        function init(obj, domain, xi_sl, velocity_sp, pressure_sp, dt, num_omp_threads)
            if nargin < 7 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end
            if ~(dt > 0)
                error('kp:solvers:BadDt', 'PUSLIncompressibleEulerSolver requires dt > 0.');
            end
            obj.Domain = domain;
            obj.Domain.buildStructs();
            obj.xi_sl = xi_sl;
            obj.dt = dt;
            obj.num_omp_threads = max(1, num_omp_threads);

            obj.AdvectionDomain = makePhysicalAdvectionDomain(obj.Domain.getVelocityDomain());
            obj.Advection.init(obj.AdvectionDomain, xi_sl, dt);
            obj.Backend.init(obj.Domain, velocity_sp, pressure_sp, dt, obj.num_omp_threads);
            obj.physical_nodes = obj.Advection.getOutputNodes();
            obj.global_physical_nodes = obj.physical_nodes;
            obj.interior_nodes = obj.Domain.getVelocityDomain().getInteriorNodes();
            if size(obj.physical_nodes, 2) ~= size(obj.global_physical_nodes, 2) || ...
                    size(obj.global_physical_nodes, 2) ~= size(obj.interior_nodes, 2)
                error('kp:solvers:BadDimension', ...
                    'PUSLIncompressibleEulerSolver advection/backend dimension mismatch.');
            end

            obj.coeff_nm2 = zeros(0, size(obj.physical_nodes, 2));
            obj.coeff_nm1 = zeros(0, size(obj.physical_nodes, 2));
            obj.coeff_n = zeros(0, size(obj.physical_nodes, 2));
            obj.completed_steps_ = 0;
        end

        function setStepSize(obj, dt)
            if ~(dt > 0)
                error('kp:solvers:BadDt', 'PUSLIncompressibleEulerSolver requires dt > 0.');
            end
            obj.dt = dt;
            obj.Advection.setStepSize(dt);
            obj.Backend.setStepSize(dt);
        end

        function setInitialVelocity(obj, physical_velocity, problem)
            local_velocity = normalizeLocalPhysicalVelocity(obj, physical_velocity, 'setInitialVelocity');
            obj.coeff_nm2 = obj.Advection.projectSamples(local_velocity);
            obj.coeff_nm1 = zeros(0, size(local_velocity, 2));
            obj.coeff_n = zeros(0, size(local_velocity, 2));
            obj.completed_steps_ = 0;
            obj.Backend.setInitialVelocityOwned(local_velocity);
        end

        function setVelocityHistory(obj, varargin)
            problem = varargin{end};
            velocities = varargin(1:end-1);
            %#ok<NASGU>
            for k = 1:numel(velocities)
                velocities{k} = normalizeLocalPhysicalVelocity(obj, velocities{k}, 'setVelocityHistory');
            end
            obj.coeff_nm2 = obj.Advection.projectSamples(velocities{1});
            if numel(velocities) >= 2
                obj.coeff_nm1 = obj.Advection.projectSamples(velocities{2});
                obj.completed_steps_ = 1;
            else
                obj.coeff_nm1 = zeros(0, size(velocities{1}, 2));
                obj.completed_steps_ = 0;
            end
            if numel(velocities) >= 3
                obj.coeff_n = obj.Advection.projectSamples(velocities{3});
                obj.completed_steps_ = 2;
            else
                obj.coeff_n = zeros(0, size(velocities{1}, 2));
            end
            obj.Backend.setVelocityHistoryOwned(velocities{:});
        end

        function out = getInteriorNodes(obj), out = obj.interior_nodes; end
        function out = getOutputNodes(obj), out = obj.physical_nodes; end
        function out = getDomain(obj), out = obj.Domain; end
        function out = advectionSolver(obj), out = obj.Advection; end
        function out = eulerBackend(obj), out = obj.Backend; end

        function sol = bdf1Step(obj, t_next, rk, body_force, problem)
            if isempty(obj.coeff_nm2)
                error('kp:solvers:MissingHistory', ...
                    'PUSLIncompressibleEulerSolver::bdf1Step requires setInitialVelocity() first.');
            end
            transported_nm2 = transportToPhysical(obj, t_next - obj.dt, obj.coeff_nm2, 1, extrapolatedVelocityCallback(obj, t_next, 1), rk);
            obj.Backend.setVelocityHistoryOwned(normalizeLocalPhysicalVelocity(obj, transported_nm2, 'bdf1Step'));
            sol = obj.Backend.bdf1Step(forceAt(obj, t_next, body_force), problem);
            obj.coeff_nm1 = obj.Advection.projectSamples(sol.velocity);
            obj.completed_steps_ = 1;
        end

        function sol = bdf2Step(obj, t_next, rk, body_force, problem)
            if obj.completed_steps_ < 1 || isempty(obj.coeff_nm1)
                error('kp:solvers:MissingHistory', ...
                    'PUSLIncompressibleEulerSolver::bdf2Step requires one prior step.');
            end
            velocity = extrapolatedVelocityCallback(obj, t_next, 2);
            transported_nm2 = transportToPhysical(obj, t_next - 2.0 * obj.dt, obj.coeff_nm2, 2, velocity, rk);
            transported_nm1 = transportToPhysical(obj, t_next - obj.dt, obj.coeff_nm1, 1, velocity, rk);
            obj.Backend.setVelocityHistoryOwned( ...
                normalizeLocalPhysicalVelocity(obj, transported_nm2, 'bdf2Step'), ...
                normalizeLocalPhysicalVelocity(obj, transported_nm1, 'bdf2Step'));
            sol = obj.Backend.bdf2Step(forceAt(obj, t_next, body_force), problem);
            obj.coeff_n = obj.Advection.projectSamples(sol.velocity);
            obj.completed_steps_ = 2;
        end

        function sol = bdf3Step(obj, t_next, rk, body_force, problem)
            if obj.completed_steps_ < 2 || isempty(obj.coeff_n)
                error('kp:solvers:MissingHistory', ...
                    'PUSLIncompressibleEulerSolver::bdf3Step requires two prior steps.');
            end
            velocity = extrapolatedVelocityCallback(obj, t_next, 3);
            transported_nm2 = transportToPhysical(obj, t_next - 3.0 * obj.dt, obj.coeff_nm2, 3, velocity, rk);
            transported_nm1 = transportToPhysical(obj, t_next - 2.0 * obj.dt, obj.coeff_nm1, 2, velocity, rk);
            transported_n = transportToPhysical(obj, t_next - obj.dt, obj.coeff_n, 1, velocity, rk);
            obj.Backend.setVelocityHistoryOwned( ...
                normalizeLocalPhysicalVelocity(obj, transported_nm2, 'bdf3Step'), ...
                normalizeLocalPhysicalVelocity(obj, transported_nm1, 'bdf3Step'), ...
                normalizeLocalPhysicalVelocity(obj, transported_n, 'bdf3Step'));
            sol = obj.Backend.bdf3Step(forceAt(obj, t_next, body_force), problem);
            obj.coeff_nm2 = obj.coeff_nm1;
            obj.coeff_nm1 = obj.coeff_n;
            obj.coeff_n = obj.Advection.projectSamples(sol.velocity);
            obj.completed_steps_ = max(obj.completed_steps_, 3);
        end

        function out = currentInteriorVelocity(obj)
            coeffs = currentCoefficients(obj);
            if isempty(coeffs)
                out = zeros(0, size(obj.interior_nodes, 2));
                return;
            end
            out = obj.Advection.evaluateAtPoints(coeffs, obj.interior_nodes);
        end

        function out = currentPhysicalVelocity(obj)
            coeffs = currentCoefficients(obj);
            if isempty(coeffs)
                out = zeros(0, size(obj.physical_nodes, 2));
                return;
            end
            out = obj.Advection.evaluateAtPoints(coeffs, obj.physical_nodes);
        end

        function out = evaluateCurrentVelocityAtPoints(obj, X)
            coeffs = currentCoefficients(obj);
            if isempty(coeffs)
                out = zeros(size(X, 1), obj.Domain.getDim());
                return;
            end
            out = obj.Advection.evaluateAtPoints(coeffs, X);
        end
    end
end

function local_velocity = normalizeLocalPhysicalVelocity(obj, velocity, where)
dim = obj.Domain.getDim();
local_sized = size(velocity, 1) == size(obj.physical_nodes, 1);
global_sized = size(velocity, 1) == size(obj.global_physical_nodes, 1);
if (~local_sized && ~global_sized) || size(velocity, 2) ~= dim || any(~isfinite(velocity), 'all')
    error('kp:solvers:BadVelocityState', ...
        'PUSLIncompressibleEulerSolver::%s received an invalid physical velocity.', where);
end
local_velocity = velocity;
end

function domain = makePhysicalAdvectionDomain(velocityDomain)
dim = velocityDomain.getDim();
domain = kp.domain.DomainDescriptor();
domain.setNodes(velocityDomain.getInteriorNodes(), velocityDomain.getBdryNodes(), zeros(0, dim));
domain.setNormals(velocityDomain.getNrmls());
domain.setSepRad(velocityDomain.getSepRad());
domain.setOuterLevelSet(velocityDomain.getOuterLevelSet());
domain.setBoundaryLevelSets(velocityDomain.getBoundaryLevelSets());
domain.buildStructs();
end

function out = forceAt(~, t, body_force)
out = @(X) body_force(t, X);
end

function w = extrapolationWeights(theta, order)
if order <= 1
    w = [1.0, 0.0, 0.0];
elseif order == 2
    w = [theta + 1.0, -theta, 0.0];
else
    w = [0.5 * (theta + 1.0) * (theta + 2.0), -theta * (theta + 2.0), 0.5 * theta * (theta + 1.0)];
end
end

function velocity = extrapolatedVelocityCallback(obj, t_next, order)
coeff0 = obj.coeff_nm2;
coeff1 = obj.coeff_nm1;
coeff2 = obj.coeff_n;
dt = obj.dt;
velocity = @(t, X) extrapolatedVelocityEval(obj, t_next, dt, order, coeff0, coeff1, coeff2, t, X);
end

function V = extrapolatedVelocityEval(obj, t_next, dt, order, coeff0, coeff1, coeff2, t, X)
if isempty(coeff0)
    V = zeros(size(X, 1), obj.Domain.getDim());
    return;
end
if order <= 1 || isempty(coeff1)
    V = obj.Advection.evaluateAtPoints(coeff0, X);
    return;
end
theta = (t - (t_next - dt)) / dt;
w = extrapolationWeights(theta, ternary(isempty(coeff2), 2, order));
if isempty(coeff2)
    V = w(1) * obj.Advection.evaluateAtPoints(coeff1, X) + w(2) * obj.Advection.evaluateAtPoints(coeff0, X);
else
    V = w(1) * obj.Advection.evaluateAtPoints(coeff2, X) + ...
        w(2) * obj.Advection.evaluateAtPoints(coeff1, X) + ...
        w(3) * obj.Advection.evaluateAtPoints(coeff0, X);
end
end

function transported = transportToPhysical(obj, t_start, coeff_state, num_steps, velocity, rk)
if num_steps <= 0
    error('kp:solvers:BadSteps', 'PUSLIncompressibleEulerSolver transport requires a positive step count.');
end
transported = backwardTraceEvaluateAtPoints(obj.Advection, t_start, num_steps, coeff_state, obj.physical_nodes, velocity, rk);
end

function coeffs = currentCoefficients(obj)
if obj.completed_steps_ <= 0
    coeffs = obj.coeff_nm2;
elseif obj.completed_steps_ == 1
    coeffs = obj.coeff_nm1;
else
    coeffs = obj.coeff_n;
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
