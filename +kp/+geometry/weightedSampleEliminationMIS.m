function keep = weightedSampleEliminationMIS(X, h)
%WEIGHTEDSAMPLEELIMINATIONMIS KernelPack-style weighted sample elimination.

n = size(X, 1);
if n == 0
    keep = false(0, 1);
    return;
end

D = kp.geometry.distanceMatrix(X, X);
weights = zeros(n, 1);
conflictNeighbors = cell(n, 1);
useParfor = localShouldUseParfor(n);

if useParfor
    parfor i = 1:n
        di = D(i, :);
        weightNeighbors = (di <= 2 * h) & (di > 0);
        c = min(di(weightNeighbors), 2 * h);
        weights(i) = sum((1 - c ./ (2 * h)).^8);
        conflictNeighbors{i} = find((di <= h) & (di > 0));
    end
else
    for i = 1:n
        di = D(i, :);
        weightNeighbors = (di <= 2 * h) & (di > 0);
        c = min(di(weightNeighbors), 2 * h);
        weights(i) = sum((1 - c ./ (2 * h)).^8);
        conflictNeighbors{i} = find(di <= h & di > 0);
    end
end

active = true(n, 1);
winners = zeros(n, 1);
nWinners = 0;

while true
    roundWinner = false(n, 1);
    for i = 1:n
        if ~active(i)
            continue;
        end
        dominates = true;
        nbrs = conflictNeighbors{i};
        for jj = 1:numel(nbrs)
            j = nbrs(jj);
            if ~active(j)
                continue;
            end
            if higherPriority(weights(j), j, weights(i), i)
                dominates = false;
                break;
            end
        end
        roundWinner(i) = dominates;
    end

    idx = find(roundWinner);
    if isempty(idx)
        break;
    end

    nw = numel(idx);
    winners(nWinners + (1:nw)) = idx;
    nWinners = nWinners + nw;
    active(idx) = false;

    for ii = 1:nw
        nbrs = conflictNeighbors{idx(ii)};
        active(nbrs) = false;
    end
end

winners = winners(1:nWinners);
if isempty(winners)
    keep = false(n, 1);
    return;
end

[~, order] = sortrows([-weights(winners), winners], [1 2]);
winners = winners(order);
keep = false(n, 1);
keep(winners) = true;

function tf = higherPriority(lhsWeight, lhsIndex, rhsWeight, rhsIndex)
if lhsWeight > rhsWeight
    tf = true;
elseif lhsWeight < rhsWeight
    tf = false;
else
    tf = lhsIndex < rhsIndex;
end

function tf = localShouldUseParfor(nPts)
tf = false;
if nPts < 64
    return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
    return;
end
tf = true;
