```@meta
CurrentModule = ARPES
```

# Dimension-aware analysis utilities

Beyond loading and k-space conversion, `ARPES.jl` provides utilities for reshaping axes,
shifting data, smoothing noisy spectra, and taking numerical derivatives.

## Axis editing

Use these helpers when the axis definition itself must change:

- `add_dim` adds a new scalar axis to a dataset,
- `convert_dim` applies a coordinate transform to an existing axis lookup while
  preserving its label and metadata,
- `rename_dim` renames a dimension label,
- `shift_dim` shifts axis lookup values, and
- `negate_dim` reverses the sign of an axis.

These operations modify the coordinate system while preserving the underlying
array values. They are useful for calibration, unit conversion, offset correction,
and harmonizing imported datasets.

## Data rebinnig

The `rebin` family of functions resample the data values themselves by averaging over
bins along a specified dimension. This is useful for:

- reducing noise by averaging neighboring points,
- matching the resolution of another dataset, and
- preparing data for visualization.

8## Data shifting by interpolation

The `shift` family moves the data values themselves using interpolation.
It supports:

- a uniform shift along one dimension,
- a position-dependent shift controlled by another dimension, and
- shifts specified by a one-dimensional dimensional array.

This is useful for tasks such as:

- energy alignment between scans,
- angle correction varying with another coordinate, and
- compensating for slow drift across a measurement series.

The implementation supports both regularly spaced and monotonic irregular grids.

## Smoothing and denoising

Several 1D filters are applied independently along a chosen dimension:

| Function                 | Typical use                                                  |
| ------------------------ | ------------------------------------------------------------ |
| `kalman_smooth_dim_llpf` | Strong noise suppression with explicit state-space smoothing |
| `boxcar_filter_dim`      | Simple moving average                                        |
| `sg_filter_dim`          | Shape-preserving smoothing with local polynomial fits        |
| `binomial_filter_dim`    | Lightweight repeated smoothing                               |
| `gaussian_filter_dim`    | Gaussian-weighted smoothing                                  |

Window sizes can be expressed either in pixels or in axis units, depending on the function
arguments.

## Numerical derivatives

The derivative utilities operate along physical dimensions instead of raw indices.

- `derivative` computes finite-difference derivatives along one axis,
- `curvature` computes 2D curvature-like enhancement, and
- `minimum_gradient` normalizes by the local gradient modulus over two dimensions.

These routines preserve array shape and metadata, which makes them convenient for
visualization and feature enhancement pipelines.

## Concatenation and metadata consensus

`cat_arpes` concatenates multiple datasets along a chosen dimension while preserving the
dimensional-array structure.

When metadata from multiple inputs must be combined, `merge_consensus` keeps only the
key-value pairs that are consistent across all inputs.

This is a practical way to build combined datasets without silently keeping conflicting
metadata.

See the [API reference](../api.md) for the detailed docstrings of these functions.
