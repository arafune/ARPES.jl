```@meta
CurrentModule = ARPES
```

# API reference

This page collects the main public API that is intended for everyday use.
For conceptual guidance, see the pages in the user guide first.

## Core type

```@docs
ARPESData
```

## I/O

```@docs
load
```

## Axis utilities

```@docs
add_dim
convert_dim
rename_dim
shift_dim
negate_dim
```

## Data transforms

```@docs
rebin
shift
rebuild_with_slice
cat_arpes
```

## Filtering

```@docs
kalman_smooth_dim_llpf
boxcar_filter_dim
sg_filter_dim
binomial_filter_dim
gaussian_filter_dim
```

## Derivative utilities

```@docs
derivative
curvature
minimum_gradient
```

## Dictionary helpers

```@docs
merge_consensus
```

## k-space conversion

```@docs
k_conversion
```
