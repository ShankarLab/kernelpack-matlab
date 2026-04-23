classdef DomainNodeGenerator < handle
    %DOMAINNODEGENERATOR KernelPack-shaped box Poisson node generator.

    properties (SetAccess = private)
        Xi_orig double = zeros(0, 0)
        Xi_pds_raw double = zeros(0, 0)
        s_dim (1,1) double = 0
        last_info struct = struct()
    end

    methods
        function generatePoissonNodes(obj, radius, x_min, x_max, varargin)
            [X, info] = kp.nodes.generatePoissonNodesInBox(radius, x_min, x_max, varargin{:});
            obj.Xi_orig = X;
            obj.Xi_pds_raw = X;
            obj.s_dim = size(X, 2);
            obj.last_info = info;
        end

        function out = getRawPoissonInteriorNodes(obj)
            out = obj.Xi_pds_raw;
        end
    end
end
