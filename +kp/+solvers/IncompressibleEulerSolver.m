classdef IncompressibleEulerSolver < handle
    %INCOMPRESSIBLEEULERSOLVER Single-process MATLAB analogue of KernelPack's Euler backend.

    properties
        Domain kp.domain.DualNodeDomainDescriptor = kp.domain.DualNodeDomainDescriptor()
        VelocityStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        PressureStencilProperties kp.rbffd.StencilProperties = kp.rbffd.StencilProperties()
        dt (1,1) double = 0
        num_omp_threads (1,1) double = 1
        VelocityPhysicalDomain kp.domain.DomainDescriptor = kp.domain.DomainDescriptor()
        Xphys double = zeros(0, 0)
        Xb double = zeros(0, 0)
        nr double = zeros(0, 0)
        Xp double = zeros(0, 0)
        GradOps cell = {}
        DivOps cell = {}
        velocity_node_range_ (1,2) double = [1 0]
        pressure_node_range_ (1,2) double = [1 0]
        unm2 double = zeros(0, 0)
        unm1 double = zeros(0, 0)
        un double = zeros(0, 0)
        completed_steps_ (1,1) double = 0
        velocity_history_node_owned_ (1,1) logical = false
    end

    methods
        function init(obj, dualDomain, velocity_sp, pressure_sp, dt, num_omp_threads)
            if nargin < 6 || isempty(num_omp_threads)
                num_omp_threads = 1;
            end
            if ~(dt > 0)
                error('kp:solvers:BadDt', 'IncompressibleEulerSolver requires dt > 0.');
            end
            obj.Domain = dualDomain;
            obj.Domain.buildStructs();
            obj.VelocityStencilProperties = velocity_sp;
            obj.PressureStencilProperties = pressure_sp;
            obj.dt = dt;
            obj.num_omp_threads = max(1, num_omp_threads);

            velDomain = obj.Domain.getVelocityDomain();
            obj.VelocityPhysicalDomain = stripGhosts(velDomain);
            obj.Xphys = obj.VelocityPhysicalDomain.getIntBdryNodes();
            obj.Xb = obj.VelocityPhysicalDomain.getBdryNodes();
            obj.nr = obj.VelocityPhysicalDomain.getNrmls();
            obj.Xp = obj.Domain.getPressureDomain().getIntBdryNodes();
            obj.velocity_node_range_ = [1, size(obj.Xphys, 1)];
            obj.pressure_node_range_ = [1, size(obj.Xp, 1)];

            obj.assembleOperators();
            obj.unm2 = zeros(0, size(obj.Xphys, 2));
            obj.unm1 = zeros(0, size(obj.Xphys, 2));
            obj.un = zeros(0, size(obj.Xphys, 2));
            obj.completed_steps_ = 0;
            obj.velocity_history_node_owned_ = false;
        end

        function setStepSize(obj, dt)
            if ~(dt > 0)
                error('kp:solvers:BadDt', 'IncompressibleEulerSolver requires dt > 0.');
            end
            obj.dt = dt;
        end

        function setInitialVelocity(obj, u0)
            validateVelocityState(obj, u0, 'setInitialVelocity');
            obj.velocity_history_node_owned_ = false;
            obj.unm2 = u0;
            obj.unm1 = zeros(0, size(u0, 2));
            obj.un = zeros(0, size(u0, 2));
            obj.completed_steps_ = 0;
        end

        function setInitialVelocityOwned(obj, u0_local)
            obj.setInitialVelocity(u0_local);
            obj.velocity_history_node_owned_ = true;
        end

        function setVelocityHistory(obj, varargin)
            obj.setVelocityHistoryImpl(false, varargin{:});
        end

        function setVelocityHistoryOwned(obj, varargin)
            obj.setVelocityHistoryImpl(true, varargin{:});
        end

        function out = getOwnedPressureRange(obj)
            out = obj.pressure_node_range_;
        end

        function out = gatherPressureSamples(~, pressure)
            out = pressure;
        end

        function sol = bdf1Step(obj, forcing, problem)
            if isempty(obj.unm2)
                error('kp:solvers:MissingHistory', 'IncompressibleEulerSolver::bdf1Step requires setInitialVelocity() first.');
            end
            sol = obj.solveBdfStep(1.0, [1.0, 0.0, 0.0], forcing, problem);
        end

        function sol = bdf2Step(obj, forcing, problem)
            if obj.completed_steps_ < 1
                error('kp:solvers:MissingHistory', 'IncompressibleEulerSolver::bdf2Step requires one prior step.');
            end
            sol = obj.solveBdfStep(1.5, [2.0, -0.5, 0.0], forcing, problem);
        end

        function sol = bdf3Step(obj, forcing, problem)
            if obj.completed_steps_ < 2
                error('kp:solvers:MissingHistory', 'IncompressibleEulerSolver::bdf3Step requires two prior steps.');
            end
            sol = obj.solveBdfStep(11.0 / 6.0, [3.0, -1.5, 1.0 / 3.0], forcing, problem);
        end
    end

    methods (Static)
        function bc = stationarySlipWall(boundary_indices)
            bc = struct();
            bc.boundary_indices = boundary_indices(:);
            bc.normal_velocity = @(X, nr) zeros(size(X, 1), 1); %#ok<INUSD>
        end

        function problem = defaultProblemDefinition()
            problem = struct();
            problem.slip_walls = {};
            problem.gauge_options = struct( ...
                'mode', "automatic", ...
                'nullspace_tol', 1.0e-7, ...
                'row_scale_detection', true);
        end
    end

    methods (Access = private)
        function setVelocityHistoryImpl(obj, ownedFlag, varargin)
            for k = 1:nargin - 2
                validateVelocityState(obj, varargin{k}, 'setVelocityHistory');
            end
            obj.velocity_history_node_owned_ = ownedFlag;
            obj.unm2 = varargin{1};
            if numel(varargin) >= 2
                obj.unm1 = varargin{2};
                obj.completed_steps_ = 1;
            else
                obj.unm1 = zeros(0, size(obj.unm2, 2));
                obj.completed_steps_ = 0;
            end
            if numel(varargin) >= 3
                obj.un = varargin{3};
                obj.completed_steps_ = 2;
            else
                obj.un = zeros(0, size(obj.unm2, 2));
            end
        end

        function assembleOperators(obj)
            dim = obj.Domain.getDim();
            obj.GradOps = cell(1, dim);
            obj.DivOps = cell(1, dim);
            opProps = kp.rbffd.OpProperties('decompose', false, 'storeWeights', true, 'recordStencils', false);
            for d = 1:dim
                opProps.selectdim = d - 1;
                gradOp = kp.rbffd.CrossNodeDiffOp(@() kp.rbffd.RBFStencil());
                gradOp.AssembleOp(obj.Domain.getPressureDomain(), obj.VelocityPhysicalDomain, 'grad', ...
                    withPointSetTree(obj.PressureStencilProperties, "interior_boundary", "interior_boundary"), opProps);
                obj.GradOps{d} = gradOp.getOp();

                divOp = kp.rbffd.CrossNodeDiffOp(@() kp.rbffd.RBFStencil());
                divOp.AssembleOp(obj.VelocityPhysicalDomain, obj.Domain.getPressureDomain(), 'grad', ...
                    withPointSetTree(obj.VelocityStencilProperties, "interior_boundary", "interior_boundary"), opProps);
                obj.DivOps{d} = divOp.getOp();
            end
        end

        function sol = solveBdfStep(obj, alpha0, historyCoeffs, forcing, problem)
            problem = normalizeProblem(problem);
            validateProblem(obj, problem);

            dim = obj.Domain.getDim();
            Nui = obj.VelocityPhysicalDomain.getNumInteriorNodes();
            Nub = obj.VelocityPhysicalDomain.getNumBdryNodes();
            Nphys = obj.VelocityPhysicalDomain.getNumIntBdryNodes();
            Npib = obj.Domain.getPressureDomain().getNumIntBdryNodes();
            pressureOffset = dim * Nphys;
            ghostOffset = pressureOffset + Npib;

            history_rhs = buildHistoryRhs(obj, historyCoeffs, 1.0 / alpha0);
            forcing_u = forcing(obj.Xphys);
            if ~isequal(size(forcing_u), size(obj.Xphys))
                error('kp:solvers:BadForcing', 'IncompressibleEulerSolver forcing callback returned invalid data.');
            end
            forcing_u = (obj.dt / alpha0) * forcing_u;
            prepared = prepareWallData(obj, problem);
            boundaryNormalVelocity = prepared.boundary_normal_velocity;

            addGauge = pressureConstantIsNull(obj, problem.gauge_options);
            totalRows = dim * Nphys + Npib + Nub + double(addGauge);
            totalCols = totalRows;

            tripletCap = dim * Nphys * (2 + 4 * max(cellfun(@nnz, obj.GradOps))) + ...
                Npib * (dim * max(cellfun(@nnz, obj.DivOps))) + Nub * (dim + 1) + addGauge * (Npib + 1);
            rows = zeros(tripletCap, 1);
            cols = zeros(tripletCap, 1);
            vals = zeros(tripletCap, 1);
            rhs = zeros(totalRows, 1);
            cursor = 1;

            for node = 1:Nphys
                isBoundary = node > Nui;
                bidx = node - Nui;
                for d = 1:dim
                    row = (d - 1) * Nphys + node;
                    rows(cursor) = row; cols(cursor) = row; vals(cursor) = 1; cursor = cursor + 1;
                    [gcols, gvals] = sparseRow(obj.GradOps{d}, node);
                    m = numel(gcols);
                    rows(cursor:cursor+m-1) = row;
                    cols(cursor:cursor+m-1) = pressureOffset + gcols;
                    vals(cursor:cursor+m-1) = (obj.dt / alpha0) * gvals;
                    cursor = cursor + m;
                    if isBoundary
                        rows(cursor) = row;
                        cols(cursor) = ghostOffset + bidx;
                        vals(cursor) = obj.nr(bidx, d);
                        cursor = cursor + 1;
                    end
                    rhs(row) = forcing_u(node, d) + history_rhs(node, d);
                end
            end

            for prow = 1:Npib
                row = pressureOffset + prow;
                for d = 1:dim
                    [dcols, dvals] = sparseRow(obj.DivOps{d}, prow);
                    m = numel(dcols);
                    rows(cursor:cursor+m-1) = row;
                    cols(cursor:cursor+m-1) = (d - 1) * Nphys + dcols;
                    vals(cursor:cursor+m-1) = dvals;
                    cursor = cursor + m;
                end
                if addGauge
                    rows(cursor) = row;
                    cols(cursor) = dim * Nphys + Npib + Nub + 1;
                    vals(cursor) = 1;
                    cursor = cursor + 1;
                end
            end

            for bidx = 1:Nub
                row = ghostOffset + bidx;
                boundaryNode = Nui + bidx;
                for d = 1:dim
                    rows(cursor) = row;
                    cols(cursor) = (d - 1) * Nphys + boundaryNode;
                    vals(cursor) = obj.nr(bidx, d);
                    cursor = cursor + 1;
                end
                rhs(row) = boundaryNormalVelocity(bidx);
            end

            if addGauge
                row = dim * Nphys + Npib + Nub + 1;
                rows(cursor:cursor+Npib-1) = row;
                cols(cursor:cursor+Npib-1) = pressureOffset + (1:Npib);
                vals(cursor:cursor+Npib-1) = 1;
                cursor = cursor + Npib;
            end

            A = sparse(rows(1:cursor-1), cols(1:cursor-1), vals(1:cursor-1), totalRows, totalCols);
            x = gmresWithFallback(A, rhs, zeros(totalCols, 1));
            if any(~isfinite(x))
                error('kp:solvers:BadSolve', 'IncompressibleEulerSolver solve returned non-finite values.');
            end

            velocity = zeros(Nphys, dim);
            for d = 1:dim
                velocity(:, d) = x((d - 1) * Nphys + (1:Nphys));
            end
            pressure = x(pressureOffset + (1:Npib));

            if obj.completed_steps_ <= 0
                obj.unm1 = velocity;
                obj.completed_steps_ = 1;
            elseif obj.completed_steps_ == 1
                obj.un = velocity;
                obj.completed_steps_ = 2;
            else
                obj.unm2 = obj.unm1;
                obj.unm1 = obj.un;
                obj.un = velocity;
                obj.completed_steps_ = max(obj.completed_steps_, 3);
            end

            [divRms, divMax] = divergenceDiagnostics(obj, velocity);
            [wallRms, wallMax] = wallDiagnostics(obj, velocity, boundaryNormalVelocity);

            sol = struct();
            sol.velocity = velocity;
            sol.pressure = pressure;
            sol.pressure_is_local = false;
            sol.pressure_begin = 1;
            sol.pressure_end = Npib;
            sol.divergence_rms = divRms;
            sol.divergence_max = divMax;
            sol.wall_normal_rms = wallRms;
            sol.wall_normal_max = wallMax;
        end
    end
end

function domain = stripGhosts(velocityDomain)
dim = velocityDomain.getDim();
domain = kp.domain.DomainDescriptor();
domain.setNodes(velocityDomain.getInteriorNodes(), velocityDomain.getBdryNodes(), zeros(0, dim));
domain.setNormals(velocityDomain.getNrmls());
domain.setSepRad(velocityDomain.getSepRad());
domain.setOuterLevelSet(velocityDomain.getOuterLevelSet());
domain.setBoundaryLevelSets(velocityDomain.getBoundaryLevelSets());
domain.buildStructs();
end

function sp = withPointSetTree(sp, pointSet, treeMode)
sp.pointSet = pointSet;
sp.treeMode = treeMode;
end

function validateVelocityState(obj, velocity, caller)
expectedRows = obj.VelocityPhysicalDomain.getNumIntBdryNodes();
expectedCols = obj.Domain.getDim();
if size(velocity, 1) ~= expectedRows || size(velocity, 2) ~= expectedCols || any(~isfinite(velocity), 'all')
    error('kp:solvers:BadVelocityState', ...
        'IncompressibleEulerSolver::%s received an invalid velocity state.', caller);
end
end

function problem = normalizeProblem(problem)
if isempty(problem)
    problem = kp.solvers.IncompressibleEulerSolver.defaultProblemDefinition();
end
if ~isfield(problem, 'slip_walls')
    problem.slip_walls = {};
end
if ~isfield(problem, 'gauge_options')
    problem.gauge_options = struct('mode', "automatic", 'nullspace_tol', 1.0e-7, 'row_scale_detection', true);
end
end

function validateProblem(obj, problem)
Nub = obj.VelocityPhysicalDomain.getNumBdryNodes();
covered = zeros(Nub, 1);
for k = 1:numel(problem.slip_walls)
    idx = problem.slip_walls{k}.boundary_indices(:);
    if any(idx < 1 | idx > Nub)
        error('kp:solvers:BadSlipWall', 'IncompressibleEulerSolver boundary index out of range.');
    end
    covered(idx) = covered(idx) + 1;
end
if any(covered ~= 1)
    error('kp:solvers:BadSlipWallCover', ...
        'IncompressibleEulerSolver slip walls must cover each boundary node exactly once.');
end
end

function prepared = prepareWallData(obj, problem)
Xb = obj.VelocityPhysicalDomain.getBdryNodes();
nr = obj.VelocityPhysicalDomain.getNrmls();
Nub = size(Xb, 1);
prepared.boundary_normal_velocity = zeros(Nub, 1);
for wallIdx = 1:numel(problem.slip_walls)
    wall = problem.slip_walls{wallIdx};
    idx = wall.boundary_indices(:);
    Xloc = Xb(idx, :);
    nloc = nr(idx, :);
    v = wall.normal_velocity(Xloc, nloc);
    v = v(:);
    if numel(v) ~= numel(idx) || any(~isfinite(v))
        error('kp:solvers:BadWallVelocity', ...
            'IncompressibleEulerSolver slip wall callback returned invalid normal velocity data.');
    end
    prepared.boundary_normal_velocity(idx) = v;
end
end

function history = buildHistoryRhs(obj, coeffs, scale)
if obj.completed_steps_ <= 0
    state0 = obj.unm2;
    state1 = [];
    state2 = [];
elseif obj.completed_steps_ == 1
    state0 = obj.unm1;
    state1 = obj.unm2;
    state2 = [];
else
    state0 = obj.un;
    state1 = obj.unm1;
    state2 = obj.unm2;
end
history = zeros(size(obj.Xphys, 1), obj.Domain.getDim());
if coeffs(1) ~= 0
    history = history + coeffs(1) * state0;
end
if coeffs(2) ~= 0
    history = history + coeffs(2) * state1;
end
if coeffs(3) ~= 0
    history = history + coeffs(3) * state2;
end
history = scale * history;
end

function tf = pressureConstantIsNull(obj, gaugeOptions)
mode = lower(string(gaugeOptions.mode));
switch mode
    case {"forcepressuremean", "force_pressure_mean"}
        tf = true;
        return;
    case "none"
        tf = false;
        return;
end
dim = obj.Domain.getDim();
Nphys = obj.VelocityPhysicalDomain.getNumIntBdryNodes();
Npib = obj.Domain.getPressureDomain().getNumIntBdryNodes();
y = zeros(dim * Nphys + Npib + obj.VelocityPhysicalDomain.getNumBdryNodes(), 1);
for d = 1:dim
    y((d - 1) * Nphys + (1:Nphys)) = obj.GradOps{d} * ones(Npib, 1);
end
residual = norm(y) / max(1.0, sqrt(Npib));
tf = residual <= gaugeOptions.nullspace_tol;
end

function [cols, vals] = sparseRow(A, row)
[~, cols, vals] = find(A(row, :));
cols = cols(:);
vals = vals(:);
end

function [rmsVal, maxVal] = divergenceDiagnostics(obj, velocity)
Npib = obj.Domain.getPressureDomain().getNumIntBdryNodes();
div = zeros(Npib, 1);
for d = 1:obj.Domain.getDim()
    div = div + obj.DivOps{d} * velocity(:, d);
end
rmsVal = norm(div) / sqrt(max(Npib, 1));
maxVal = max(abs(div));
end

function [rmsVal, maxVal] = wallDiagnostics(obj, velocity, boundaryNormalVelocity)
Nui = obj.VelocityPhysicalDomain.getNumInteriorNodes();
Nub = obj.VelocityPhysicalDomain.getNumBdryNodes();
vals = zeros(Nub, 1);
for bidx = 1:Nub
    vals(bidx) = obj.nr(bidx, :) * velocity(Nui + bidx, :).';
end
resid = vals - boundaryNormalVelocity(:);
rmsVal = norm(resid) / sqrt(max(Nub, 1));
maxVal = max(abs(resid));
end

function x = gmresWithFallback(A, b, x0)
if nargin < 3 || isempty(x0)
    x0 = zeros(size(b));
end
[x, flag] = gmres(A, b, [], 1.0e-10, min(200, size(A, 1)), [], [], x0);
if flag ~= 0 || any(~isfinite(x))
    x = A \ b;
end
end
