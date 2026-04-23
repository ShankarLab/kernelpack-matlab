# kernelpack-matlab

Fresh restart of the MATLAB companion to KernelPack.

This pass is being rebuilt around the KernelPack geometry contracts first:

- `EmbeddedSurface`
- `PiecewiseSmoothEmbeddedSurface`
- `RBFLevelSet`

The goal is to mirror the geometry objects and public data they expose before
adding broader solver code.

Current status:

- `EmbeddedSurface` is in place with KernelPack-like stored geometry state.
- `PiecewiseSmoothEmbeddedSurface` is in place with assembled boundary state.
- `RBFLevelSet` is in place with KernelPack-like method names.
- The clean restart currently implements 2D smooth closed and open-segment
  geometric-model construction first.

