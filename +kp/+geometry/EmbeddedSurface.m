classdef EmbeddedSurface < handle
    %EMBEDDEDSURFACE KernelPack-like geometric description of one boundary object.

    properties
        ubox struct = struct('p', [], 'V', [], 'D', [])
        cbox struct = struct('p', [], 'V', [], 'D', [])
        ebox struct = struct('p', [], 'V', [], 'D', [])
        data_sites double = zeros(0, 0)
        data_site_nrmls double = zeros(0, 0)
        uniform_sample_sites double = zeros(0, 0)
        sample_sites double = zeros(0, 0)
        sample_sites_s double = zeros(0, 0)
        thick_copy double = zeros(0, 0)
        uniform_Nrmls double = zeros(0, 0)
        tangents1 double = zeros(0, 0)
        tangents2 double = zeros(0, 0)
        tangents1_s double = zeros(0, 0)
        tangents2_s double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        Nrmls_s double = zeros(0, 0)
        geom_model struct = struct()
        Nd (1,1) double = 0
        N (1,1) double = 0
        surf_dim (1,1) double = 0
        sep_rad (1,1) double = NaN
        uniformTree struct = struct('Points', [])
        tallTree struct = struct('Points', [])
        thick_tree struct = struct('Points', [])
        levelSet kp.geometry.RBFLevelSet = kp.geometry.RBFLevelSet()
    end

    methods
        function setDataSites(obj, dataSites)
            obj.data_sites = dataSites;
            obj.Nd = size(dataSites, 1);
            obj.surf_dim = size(dataSites, 2) - 1;
        end

        function setSampleSites(obj, sampleSites)
            obj.sample_sites = sampleSites;
            obj.N = size(sampleSites, 1);
        end

        function setUniformSampleSites(obj, sampleSites)
            obj.uniform_sample_sites = sampleSites;
        end

        function buildClosedGeometricModelPS(obj, dim, rad, Nb, Ne, ~, method, supersample_fac, ~, ~)
            if dim ~= 2
                error('EmbeddedSurface:NotImplemented', ...
                    'Clean restart currently implements buildClosedGeometricModelPS for dim=2 first.');
            end

            obj.sep_rad = rad;
            obj.Nd = min(Nb, size(obj.data_sites, 1));
            obj.surf_dim = dim - 1;

            dataSites = obj.data_sites(1:obj.Nd, :);
            t = kp.geometry.chordLengthParam(dataSites, true);
            theta = 2 * pi * t;
            [r, ~] = kp.geometry.periodicChordDistance(theta, theta);
            degree = 7;
            K = kp.geometry.phsKernel(r, degree);
            reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
            weights = (K + reg * eye(size(K))) \ dataSites;

            obj.geom_model = struct( ...
                'type', 'closed-curve-sbf', ...
                'degree', degree, ...
                'theta', theta, ...
                'weights', weights);

            if method == 1
                Ns = max(round(supersample_fac * 1.5 * Ne), Ne);
                ts = linspace(0, 1, Ns + 1).';
                ts = ts(1:end-1);
                ptss = obj.evalClosedCurve(ts);
                [tanS, nrmlS] = obj.evalClosedCurveFrame(ts);
                keep = round(linspace(1, Ns, Ne));
                obj.sample_sites_s = ptss;
                obj.sample_sites = ptss(keep, :);
                obj.tangents1_s = tanS;
                obj.Nrmls_s = nrmlS;
                obj.tangents1 = tanS(keep, :);
                obj.Nrmls = nrmlS(keep, :);
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
            else
                tEval = linspace(0, 1, Ne + 1).';
                tEval = tEval(1:end-1);
                obj.sample_sites = obj.evalClosedCurve(tEval);
                [obj.tangents1, obj.Nrmls] = obj.evalClosedCurveFrame(tEval);
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
            end

            obj.N = size(obj.sample_sites, 1);
            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildGeometricModelPS(obj, dim, rad, Nb, Ne, ~, method, supersample_fac, ~, ~)
            if dim ~= 2
                error('EmbeddedSurface:NotImplemented', ...
                    'Clean restart currently implements buildGeometricModelPS for dim=2 first.');
            end

            obj.sep_rad = rad;
            obj.Nd = min(Nb, size(obj.data_sites, 1));
            obj.surf_dim = dim - 1;

            dataSites = obj.data_sites(1:obj.Nd, :);
            u = kp.geometry.chordLengthParam(dataSites, false);
            degree = 7;
            K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(u, u), degree);
            P = [ones(numel(u), 1), u];
            reg = 1e-12 * max(1.0, max(abs(K), [], 'all'));
            A = [K + reg * eye(numel(u)), P; P.', zeros(2, 2)];
            coeffs = A \ [dataSites; zeros(2, 2)];

            obj.geom_model = struct( ...
                'type', 'open-curve-rbf', ...
                'degree', degree, ...
                'u', u, ...
                'rbfWeights', coeffs(1:numel(u), :), ...
                'polyCoeffs', coeffs(numel(u) + 1:end, :));

            if method == 1
                Ns = max(round(supersample_fac * 1.5 * Ne), Ne);
                us = linspace(0, 1, Ns).';
                ptss = obj.evalOpenCurve(us);
                [tanS, nrmlS] = obj.evalOpenCurveFrame(us);
                keep = round(linspace(1, Ns, Ne));
                obj.sample_sites_s = ptss;
                obj.sample_sites = ptss(keep, :);
                obj.tangents1_s = tanS;
                obj.Nrmls_s = nrmlS;
                obj.tangents1 = tanS(keep, :);
                obj.Nrmls = nrmlS(keep, :);
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
            else
                uEval = linspace(0, 1, Ne).';
                obj.sample_sites = obj.evalOpenCurve(uEval);
                [obj.tangents1, obj.Nrmls] = obj.evalOpenCurveFrame(uEval);
                obj.uniform_sample_sites = obj.sample_sites;
                obj.uniform_Nrmls = obj.Nrmls;
            end

            obj.N = size(obj.sample_sites, 1);
            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildLevelSetFromGeometricModel(obj, lambda)
            if nargin < 2 || isempty(lambda)
                pts = obj.uniform_sample_sites;
                nrmls = obj.uniform_Nrmls;
            else
                if strcmp(obj.geom_model.type, 'closed-curve-sbf')
                    pts = obj.evalClosedCurve(lambda);
                    [~, nrmls] = obj.evalClosedCurveFrame(lambda);
                else
                    pts = obj.evalOpenCurve(lambda);
                    [~, nrmls] = obj.evalOpenCurveFrame(lambda);
                end
            end
            obj.levelSet = kp.geometry.RBFLevelSet();
            obj.levelSet.BuildLevelSetFromCFI(pts, nrmls);
        end

        function computeBoundingBox(obj)
            obj.cbox = kp.geometry.pcaOrientedBoundingBox(obj.sample_sites);
        end

        function computeUniformBoundingBox(obj)
            obj.ubox = kp.geometry.pcaOrientedBoundingBox(obj.uniform_sample_sites);
        end

        function computeExtendedBoundingBox(obj, h)
            obj.thick_copy = obj.uniform_sample_sites + h * obj.uniform_Nrmls;
            obj.thick_tree = struct('Points', obj.thick_copy);
            obj.ebox = kp.geometry.pcaOrientedBoundingBox(obj.thick_copy);
        end

        function buildTree(obj)
            obj.tallTree = struct('Points', obj.sample_sites);
        end

        function buildUniformTree(obj)
            obj.uniformTree = struct('Points', obj.uniform_sample_sites);
        end

        function flipNormals(obj)
            obj.Nrmls = -obj.Nrmls;
            obj.uniform_Nrmls = -obj.uniform_Nrmls;
            obj.Nrmls_s = -obj.Nrmls_s;
            obj.data_site_nrmls = -obj.data_site_nrmls;
        end

        function out = getSampleSites(obj), out = obj.sample_sites; end
        function out = getUniformSampleSites(obj), out = obj.uniform_sample_sites; end
        function out = getNrmls(obj), out = obj.Nrmls; end
        function out = getUniformNrmls(obj), out = obj.uniform_Nrmls; end
        function out = getBoundingBox(obj), out = obj.cbox.p; end
        function out = getUniformBoundingBox(obj), out = obj.ubox.p; end
        function out = getExtendedBoundingBox(obj), out = obj.ebox.p; end
        function out = getKDTree(obj), out = obj.tallTree; end
        function out = getUniformKDTree(obj), out = obj.uniformTree; end
        function out = getThickTree(obj), out = obj.thick_tree; end
        function out = getLevelSet(obj), out = obj.levelSet; end
        function out = getN(obj), out = obj.N; end
    end

    methods (Access = private)
        function x = evalClosedCurve(obj, t)
            tw = kp.geometry.wrapPeriodicParameter(t(:));
            theta = 2 * pi * tw;
            [r, ~] = kp.geometry.periodicChordDistance(theta, obj.geom_model.theta);
            K = kp.geometry.phsKernel(r, obj.geom_model.degree);
            x = K * obj.geom_model.weights;
        end

        function [xt, n] = evalClosedCurveFrame(obj, t)
            tw = kp.geometry.wrapPeriodicParameter(t(:));
            theta = 2 * pi * tw;
            [r, delta] = kp.geometry.periodicChordDistance(theta, obj.geom_model.theta);
            degree = obj.geom_model.degree;
            if degree == 1
                dphi = sin(delta) ./ max(r, eps);
                dphi(r == 0) = 0;
            else
                dphi = degree * sin(delta) .* (r .^ (degree - 2));
                dphi(r == 0) = 0;
            end
            xt = (2 * pi) * (dphi * obj.geom_model.weights);
            n = kp.geometry.normalizeRows([-xt(:, 2), xt(:, 1)]);
        end

        function x = evalOpenCurve(obj, u)
            u = min(max(u(:), 0), 1);
            K = kp.geometry.phsKernel(kp.geometry.distanceMatrix(u, obj.geom_model.u), obj.geom_model.degree);
            x = K * obj.geom_model.rbfWeights + [ones(size(u, 1), 1), u] * obj.geom_model.polyCoeffs;
        end

        function [xt, n] = evalOpenCurveFrame(obj, u)
            u = min(max(u(:), 0), 1);
            du = u - obj.geom_model.u.';
            r = abs(du);
            degree = obj.geom_model.degree;
            if degree == 1
                dphi = sign(du);
            else
                dphi = degree * sign(du) .* (r .^ (degree - 1));
            end
            dphi(r == 0 & degree > 1) = 0;
            xt = dphi * obj.geom_model.rbfWeights + repmat(obj.geom_model.polyCoeffs(2, :), size(u, 1), 1);
            n = kp.geometry.normalizeRows([-xt(:, 2), xt(:, 1)]);
        end
    end
end
