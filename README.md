# kernelpack-matlab

`kernelpack-matlab` is a MATLAB companion project to
[KernelPack](https://github.com/VarShankar/kernelpack).

The aim is to bring the core geometry and meshfree building blocks of
KernelPack into a MATLAB codebase that is easier to inspect, prototype with,
and extend.

## Current focus

The repository is currently centered on the geometry layer. The first pieces in
place are modeled after the main geometry objects used in KernelPack:

- `EmbeddedSurface`
- `PiecewiseSmoothEmbeddedSurface`
- `RBFLevelSet`

These classes live in [`+kp/+geometry`](+kp/+geometry).

## What is here now

The current code establishes a KernelPack-shaped starting point for geometry:

- `EmbeddedSurface` stores the data sites, sampled boundary points, normals,
  tangents, bounding boxes, thickened copies, and level-set representation for
  a single boundary object.
- `PiecewiseSmoothEmbeddedSurface` stores the segment list together with the
  assembled boundary cloud, normals, segment map, corner flags, bounding boxes,
  and level-set representation for a piecewise-smooth boundary.
- `RBFLevelSet` provides an implicit boundary representation with evaluation,
  gradient evaluation, inside-outside tests, and Newton projection routines.

At the moment, the implemented geometric-model construction path is the 2D
case:

- smooth closed curves
- open curve segments
- piecewise-smooth planar boundaries assembled from segments

## Basic use

### Smooth closed curve

```matlab
t = linspace(0, 2*pi, 40).';
t(end) = [];
x = [cos(t), 0.7*sin(t)];

surface = kp.geometry.EmbeddedSurface();
surface.setDataSites(x);
surface.buildClosedGeometricModelPS(2, 0.05, size(x,1), 120, ...
    'eval + eval_first_ders', 1, 2, 2, 1);
surface.buildLevelSetFromGeometricModel([]);

xb = surface.getSampleSites();
nrmls = surface.getNrmls();
phi = surface.getLevelSet().Evaluate(xb);
```

### Open curve segment

```matlab
u = linspace(0, 1, 30).';
x = [u, u.^2];

segment = kp.geometry.EmbeddedSurface();
segment.setDataSites(x);
segment.buildGeometricModelPS(2, 0.05, size(x,1), 80, ...
    'eval + eval_first_ders', 1, 2, 2, 1);

xb = segment.getSampleSites();
nrmls = segment.getNrmls();
```

### Piecewise-smooth planar boundary

```matlab
seg1 = [linspace(0,1,20).', zeros(20,1)];
seg2 = [ones(20,1), linspace(0,1,20).'];
seg3 = [linspace(1,0,20).', ones(20,1)];
seg4 = [zeros(20,1), linspace(1,0,20).'];

surface = kp.geometry.PiecewiseSmoothEmbeddedSurface();
surface.generatePiecewiseSmoothSurfaceBySegment( ...
    {seg1, seg2, seg3, seg4}, [false false false false], 0.05, 1, 2, 2);
surface.buildLevelSet();

xb = surface.getBdryNodes();
nrmls = surface.getBdryNrmls();
corners = surface.getCornerFlags();
```

## Project direction

The goal is to keep building outward from the KernelPack geometry contracts
rather than inventing a separate MATLAB abstraction stack. The next major steps
are:

1. strengthen the current geometry objects
2. add the smooth closed 3D surface path
3. extend piecewise geometry handling further in 3D
4. build node-generation and discretization layers on top of that geometry

## Notes

This repository is under active construction. The current contents should be
viewed as an early geometry foundation rather than a full MATLAB port of
KernelPack.
