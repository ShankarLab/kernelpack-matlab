function out = puPatchGeometry(domain, xi, patch_spacing_factor, patch_radius_factor)
%PUPATCHGEOMETRY Shared PU patch geometry used by the MATLAB PU solvers.

if nargin < 3 || isempty(patch_spacing_factor)
    patch_spacing_factor = 0.0;
end
if nargin < 4 || isempty(patch_radius_factor)
    patch_radius_factor = 0.0;
end

domain.buildStructs();
Xf = domain.getAllNodes();
h = domain.getSepRad();

spacing = choosePatchSpacing(h, patch_spacing_factor);
radius = choosePatchRadius(h, patch_radius_factor);
min_nodes = chooseMinimumPatchNodes(size(Xf, 2), xi);
centers = choosePatchCenters(Xf, spacing);
node_ids = buildPatchNodeIds(domain, centers, radius, min_nodes);

out = struct();
out.centers = centers;
out.radius = radius;
out.spacing = spacing;
out.min_patch_nodes = min_nodes;
out.node_ids = node_ids;
end

function spacing = choosePatchSpacing(h, patch_spacing_factor)
if patch_spacing_factor > 0
    spacing = patch_spacing_factor * h;
else
    spacing = 2.0 * h;
end
end

function radius = choosePatchRadius(h, patch_radius_factor)
if patch_radius_factor > 0
    radius = patch_radius_factor * h;
else
    radius = 3.0 * h;
end
end

function centers = choosePatchCenters(X, spacing)
if isempty(X)
    centers = zeros(0, size(X, 2));
    return;
end

remaining = true(size(X, 1), 1);
centers = zeros(0, size(X, 2));
for i = 1:size(X, 1)
    if ~remaining(i)
        continue;
    end
    xi = X(i, :);
    centers = [centers; xi]; %#ok<AGROW>
    d = sqrt(sum((X - xi).^2, 2));
    remaining(d <= spacing) = false;
end
end

function patch_node_ids = buildPatchNodeIds(domain, centers, radius, min_patch_nodes)
patch_node_ids = cell(size(centers, 1), 1);
num_all = domain.getNumTotalNodes();
for p = 1:size(centers, 1)
    ids = domain.queryBall("all", centers(p, :), radius);
    patch_node_ids{p} = ids{1}(:);
    if numel(patch_node_ids{p}) < min_patch_nodes
        [idsKnn, ~] = domain.queryKnn("all", centers(p, :), min(min_patch_nodes, num_all));
        patch_node_ids{p} = idsKnn(1, :).';
    end
end
end

function min_nodes = chooseMinimumPatchNodes(dim, xi)
ell = max(xi + 1, 2);
npoly = size(kp.poly.total_degree_indices(dim, ell), 1);
min_nodes = 2 * npoly + 1;
end
