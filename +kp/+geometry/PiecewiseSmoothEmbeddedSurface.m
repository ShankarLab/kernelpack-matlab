classdef PiecewiseSmoothEmbeddedSurface < handle
    %PIECEWISESMOOTHEMBEDDEDSURFACE KernelPack-like piecewise boundary object.

    properties (SetAccess = private)
        segments cell = {}
        Xb double = zeros(0, 0)
        Xb_uniform double = zeros(0, 0)
        Nrmls double = zeros(0, 0)
        Nrmls_uniform double = zeros(0, 0)
        Nrmls_raw double = zeros(0, 0)
        Nrmls_uniform_raw double = zeros(0, 0)
        Nrmls_smoothed double = zeros(0, 0)
        Nrmls_uniform_smoothed double = zeros(0, 0)
        thick_copy double = zeros(0, 0)
        N (1,1) double = 0
        segment_map (:,1) double = zeros(0, 1)
        corner_flags (:,1) double = zeros(0, 1)
        cbox struct = struct('p', [], 'V', [], 'D', [])
        ubox struct = struct('p', [], 'V', [], 'D', [])
        ebox struct = struct('p', [], 'V', [], 'D', [])
        bdryTree struct = struct('Points', [])
        uniform_bdryTree struct = struct('Points', [])
        thick_tree struct = struct('Points', [])
        levelSet kp.geometry.RBFLevelSet = kp.geometry.RBFLevelSet()
    end

    methods
        function generatePiecewiseSmoothSurfaceBySegment(obj, bdry_segments, flip_normal, radius, method, supersample_fac, mode)
            if nargin < 5 || isempty(method), method = 1; end
            if nargin < 6 || isempty(supersample_fac), supersample_fac = 2; end
            if nargin < 7 || isempty(mode), mode = 2; end

            nSeg = numel(bdry_segments);
            obj.segments = cell(size(bdry_segments));
            obj.segment_map = zeros(0, 1);
            obj.corner_flags = zeros(0, 1);

            dim = size(bdry_segments{1}, 2);
            if ~isscalar(radius)
                radius = min(radius(:));
            end

            localSegments = cell(size(bdry_segments));
            if obj.shouldUseParfor(nSeg)
                parfor k = 1:nSeg
                    localSegments{k} = kp.geometry.PiecewiseSmoothEmbeddedSurface.buildSegmentModel( ...
                        bdry_segments{k}, flip_normal(k), dim, radius, method, supersample_fac, mode);
                end
            else
                for k = 1:nSeg
                    localSegments{k} = kp.geometry.PiecewiseSmoothEmbeddedSurface.buildSegmentModel( ...
                        bdry_segments{k}, flip_normal(k), dim, radius, method, supersample_fac, mode);
                end
            end
            obj.segments = localSegments;

            counts = zeros(numel(obj.segments), 1);
            for k = 1:numel(obj.segments)
                counts(k) = size(obj.segments{k}.getSampleSites(), 1);
            end
            totalCount = sum(counts);
            Xb_s = zeros(totalCount, dim);
            Nrmls_s = zeros(totalCount, dim);
            obj.segment_map = zeros(totalCount, 1);
            obj.corner_flags = zeros(totalCount, 1);
            row0 = 1;
            for k = 1:numel(obj.segments)
                seg = obj.segments{k};
                pts = seg.getSampleSites();
                nrmls = seg.getNrmls();
                m = size(pts, 1);
                rows = row0:(row0 + m - 1);
                Xb_s(rows, :) = pts;
                Nrmls_s(rows, :) = nrmls;
                obj.segment_map(rows) = k;
                flags = zeros(m, 1);
                flags(1) = 1;
                flags(end) = 1;
                obj.corner_flags(rows) = flags;
                row0 = row0 + m;
            end

            keep = obj.dedupMask(Xb_s, 0.2 * radius);
            obj.Xb = Xb_s(keep, :);
            obj.Nrmls = kp.geometry.normalizeRows(Nrmls_s(keep, :));
            obj.segment_map = obj.segment_map(keep);
            obj.corner_flags = obj.corner_flags(keep);
            obj.Xb_uniform = obj.Xb;
            obj.Nrmls_uniform = obj.Nrmls;
            obj.Nrmls_raw = obj.Nrmls;
            obj.Nrmls_uniform_raw = obj.Nrmls_uniform;
            obj.Nrmls_smoothed = obj.Nrmls;
            obj.Nrmls_uniform_smoothed = obj.Nrmls_uniform;
            obj.N = size(obj.Xb, 1);

            obj.buildTree();
            obj.buildUniformTree();
        end

        function buildLevelSet(obj)
            obj.levelSet = kp.geometry.RBFLevelSet();
            obj.levelSet.BuildLevelSetFromCFI(obj.Xb_uniform, obj.Nrmls_uniform);
        end

        function buildTree(obj)
            obj.bdryTree = struct('Points', obj.Xb);
        end

        function buildUniformTree(obj)
            obj.uniform_bdryTree = struct('Points', obj.Xb_uniform);
        end

        function buildThickTree(obj)
            obj.thick_tree = struct('Points', obj.thick_copy);
        end

        function computeBoundingBox(obj)
            obj.cbox = kp.geometry.pcaOrientedBoundingBox(obj.Xb);
        end

        function computeUniformBoundingBox(obj)
            obj.ubox = kp.geometry.pcaOrientedBoundingBox(obj.Xb_uniform);
        end

        function computeExtendedBoundingBox(obj, h)
            obj.thick_copy = obj.Xb_uniform + h * obj.Nrmls_uniform;
            obj.buildThickTree();
            obj.ebox = kp.geometry.pcaOrientedBoundingBox(obj.thick_copy);
        end

        function out = getSampleSites(obj), out = obj.Xb; end
        function out = getUniformSampleSites(obj), out = obj.Xb_uniform; end
        function out = getNrmls(obj), out = obj.Nrmls; end
        function out = getUniformNrmls(obj), out = obj.Nrmls_uniform; end
        function out = getCornerFlags(obj), out = obj.corner_flags; end
        function out = getKDTree(obj), out = obj.bdryTree; end
        function out = getUniformKDTree(obj), out = obj.uniform_bdryTree; end
        function out = getThickTree(obj), out = obj.thick_tree; end
        function out = getLevelSet(obj), out = obj.levelSet; end
        function out = getBdryNodes(obj), out = obj.Xb; end
        function out = getUniformBdryNodes(obj), out = obj.Xb_uniform; end
        function out = getBdryNrmls(obj), out = obj.Nrmls; end
        function out = getUniformBdryNrmls(obj), out = obj.Nrmls_uniform; end
        function out = getBoundingBox(obj), out = obj.cbox.p; end
        function out = getUniformBoundingBox(obj), out = obj.ubox.p; end
        function out = getExtendedBoundingBox(obj), out = obj.ebox.p; end
    end

    methods (Access = private)
        function tf = shouldUseParfor(~, nSeg)
            tf = false;
            if nSeg < 2
                return;
            end
            if ~license('test', 'Distrib_Computing_Toolbox')
                return;
            end
            tf = true;
        end

        function keep = dedupMask(~, X, tol)
            keep = true(size(X, 1), 1);
            for i = 1:size(X, 1)
                if ~keep(i)
                    continue;
                end
                d = sqrt(sum((X(i+1:end, :) - X(i, :)).^2, 2));
                dup = find(d <= tol) + i;
                keep(dup) = false;
            end
        end
    end

    methods (Static, Access = private)
        function boundary = buildSegmentModel(bdryPts, flipNormal, dim, radius, method, supersample_fac, mode)
            boundary = kp.geometry.EmbeddedSurface();
            boundary.setDataSites(bdryPts);
            boundary.setSampleSites(bdryPts);
            boundary.computeBoundingBox();
            box = boundary.getBoundingBox();
            if dim == 2
                w = abs(max(box(:, 1)) - min(box(:, 1)));
                h = abs(max(box(:, 2)) - min(box(:, 2)));
                bdry_size = 2 * (w + h);
                seg_N = max(8, round(bdry_size / radius));
                boundary.buildGeometricModelPS(dim, radius, size(bdryPts, 1), seg_N, ...
                    'eval + eval_first_ders', method, supersample_fac, mode, 1);
            elseif dim == 3
                w = abs(max(box(:, 1)) - min(box(:, 1)));
                h = abs(max(box(:, 2)) - min(box(:, 2)));
                l = abs(max(box(:, 3)) - min(box(:, 3)));
                bdry_size = 2 * (w * h + h * l + l * w);
                seg_N = max(16, round(bdry_size / max(radius^2, eps)));
                boundary.buildGeometricModelPS(dim, radius, size(bdryPts, 1), seg_N, ...
                    'eval + eval_first_ders', method, supersample_fac, mode, 1);
            else
                error('PiecewiseSmoothEmbeddedSurface:NotImplemented', ...
                    'Unsupported segment dimension.');
            end
            if flipNormal
                boundary.flipNormals();
            end
        end
    end
end
