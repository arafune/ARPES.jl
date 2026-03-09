using Interpolations

"""
    fast_interpolate(p0, p1, [p2], c0, c1, [c2], values)

Perform n-dimensional linear interpolation on a 2D or 3D grid. 

This function automatically dispatches to the most efficient algorithm based on the 
input types and dimensions:
- **Gridded(Linear)**: Used when coordinates are non-uniformly spaced (`AbstractVector`).
- **BSpline(Linear)**: Used when coordinates are uniformly spaced (`AbstractRange`), 
  providing significant performance gains.

# Requirements
- `p0`, `p1`, and (optionally) `p2` **must have the same shape and size**.
- `c0`, `c1`, and `c2` must be monotonic.
- `values` must match the dimensions of the grid defined by `c0`, `c1`, and `c2`.


# Arguments
- `p0, p1, [p2]`: `AbstractArray` of desired positions for each dimension.
- `c0, c1, [c2]`: `AbstractVector` (irregular) or `AbstractRange` (regular) 
  defining the original grid coordinates.
- `values`: `AbstractMatrix` (2D) or `AbstractArray{T, 3}` (3D) containing the 
  data values at grid points.

# Returns
- An array of interpolated values with the same shape as `p0`. 
  Points outside the grid boundaries are returned as `NaN`.
"""
function fast_interpolate end

# --- 2D Bilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function fast_interpolate(
    p0::AbstractArray,
    p1::AbstractArray,
    c0::AbstractRange,
    c1::AbstractRange,
    values::AbstractMatrix,
)
    # Use BSpline interpolation for optimized performance on uniform grids
    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c0, c1)
    return extrapolate(sitp, NaN).(p0, p1)
end


# Irregular grid: coords are general vectors
function fast_interpolate(
    p0::AbstractArray,
    p1::AbstractArray,
    c0::AbstractVector,
    c1::AbstractVector,
    values::AbstractMatrix,
)
    # Use Gridded interpolation for non-uniform spacing
    itp = interpolate((c0, c1), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p0, p1)
end

# --- 3D Trilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function fast_interpolate(
    p0::AbstractArray,
    p1::AbstractArray,
    p2::AbstractArray,
    c0::AbstractRange,
    c1::AbstractRange,
    c2::AbstractRange,
    values::AbstractArray{T,3},
) where {T}
    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c0, c1, c2)
    return extrapolate(sitp, NaN).(p0, p1, p2)
end
# Irregular grid: coords are general vectors
function fast_interpolate(
    p0::AbstractArray,
    p1::AbstractArray,
    p2::AbstractArray,
    c0::AbstractVector,
    c1::AbstractVector,
    c2::AbstractVector,
    values::AbstractArray{T,3},
) where {T}
    itp = interpolate((c0, c1, c2), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p0, p1, p2)
end


