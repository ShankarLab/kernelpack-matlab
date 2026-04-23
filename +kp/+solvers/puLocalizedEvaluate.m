function values = puLocalizedEvaluate(domain, patchData, xi, coeffs, Xq)
%PULOCALIZEDEVALUATE Partition-of-unity evaluation on arbitrary query points.

if isempty(Xq)
    values = zeros(0, size(coeffs, 2));
    return;
end

Xnodes = domain.getAllNodes();
radius = patchData.radius;
centers = patchData.centers;
node_ids = patchData.node_ids;
dim = size(Xnodes, 2);
nc = size(coeffs, 2);

sp = kp.rbffd.StencilProperties();
sp.dim = dim;
sp.ell = max(xi + 1, 2);
sp.npoly = size(kp.poly.total_degree_indices(dim, sp.ell), 1);
sp.spline_degree = max(5, sp.ell);
if mod(sp.spline_degree, 2) == 0
    sp.spline_degree = sp.spline_degree - 1;
end

values = zeros(size(Xq, 1), nc);
weight_sum = zeros(size(Xq, 1), 1);
for p = 1:size(centers, 1)
    center = centers(p, :);
    dist = sqrt(sum((Xq - center).^2, 2));
    mask = dist < radius;
    if ~any(mask)
        continue;
    end
    ids = node_ids{p};
    stencil = kp.rbffd.RBFStencil();
    stencil.InitializeGeometry(Xnodes(ids, :), sp);
    wloc = puPatchWeight(dist(mask) ./ radius);
    vloc = stencil.EvalStencil(sp, Xq(mask, :), coeffs(ids, :), false);
    values(mask, :) = values(mask, :) + wloc .* vloc;
    weight_sum(mask) = weight_sum(mask) + wloc;
end

missing = weight_sum <= 1.0e-14;
if any(missing)
    [idx, ~] = domain.queryKnn("all", Xq(missing, :), 1);
values(missing, :) = coeffs(idx(:, 1), :);
    weight_sum(missing) = 1.0;
end

values = values ./ weight_sum;
end

function w = puPatchWeight(r)
w = zeros(size(r));
mask = r < 1;
t = 1 - r(mask);
rm = r(mask);
w(mask) = t.^8 .* (32 * rm.^3 + 25 * rm.^2 + 8 * rm + 1);
end
