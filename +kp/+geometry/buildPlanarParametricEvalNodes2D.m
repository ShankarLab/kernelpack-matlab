function Xe = buildPlanarParametricEvalNodes2D(interpNodes, N)
%BUILDPLANARPARAMETRICEVALNODES2D Grid-like evaluation nodes in a projected hull.

[X, ~, ~] = kp.geometry.projectToBestFitPlane(interpNodes);
x0 = min(X(:, 1));
xm = max(X(:, 1));
y0 = min(X(:, 2));
ym = max(X(:, 2));

m = max(4, ceil(sqrt(max(2 * N, 4))));
[xx, yy] = meshgrid(linspace(x0, xm, m), linspace(y0, ym, m));
P = [xx(:), yy(:)];
inside = kp.geometry.pointsInConvexHull2D(P, X);
Xe = P(inside, :);

if size(Xe, 1) > N
    idx = round(linspace(1, size(Xe, 1), N));
    Xe = Xe(idx, :);
end

