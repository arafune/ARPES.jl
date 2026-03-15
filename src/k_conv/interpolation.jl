using Interpolations

"""
    _interpolate(p0, p1, [p2], c0, c1, [c2], values)

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
function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    c0::AbstractRange{<:Real},
    c1::AbstractRange{<:Real},
    values::AbstractMatrix{<:Real},
)
    # Use BSpline interpolation for optimized performance on uniform grids
    @debug "size(values), length(c0), length(c1)" size(values) length(c0) length(c1)
    values = step(c0) < 0 ? reverse(values, dims = 1) : values
    c0 = step(c0) < 0 ? reverse(c0) : c0
    values = step(c1) < 0 ? reverse(values, dims = 2) : values
    c1 = step(c1) < 0 ? reverse(c1) : c1

    @debug "c0, c1" c0 c1
    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c0, c1)
    return extrapolate(sitp, NaN).(p0, p1)
end


# Irregular grid: coords are general vectors
function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    c0::AbstractVector{<:Real},
    c1::AbstractVector{<:Real},
    values::AbstractMatrix{<:Real},
)
    # Use Gridded interpolation for non-uniform spacing
    @debug "size(values), length(c0), length(c1)" size(values) length(c0) length(c1)
    itp = interpolate((c0, c1), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p0, p1)
end

# --- 3D Trilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    c0::AbstractRange{<:Real},
    c1::AbstractRange{<:Real},
    c2::AbstractRange{<:Real},
    values::AbstractArray{T,3},
) where {T<:Real}
    @debug "size(values), length(c0), length(c1)" size(values) length(c0) length(c1)
    values = step(c0) < 0 ? reverse(values, dims = 1) : values
    c0 = step(c0) < 0 ? reverse(c0) : c0
    values = step(c1) < 0 ? reverse(values, dims = 2) : values
    c1 = step(c1) < 0 ? reverse(c1) : c1
    values = step(c2) < 0 ? reverse(values, dims = 3) : values
    c2 = step(c2) < 0 ? reverse(c2) : c2

    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c0, c1, c2)
    return extrapolate(sitp, NaN).(p0, p1, p2)
end

function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    c0::AbstractVector{<:Real},
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    values::AbstractArray{T,3},
) where {T<:Real}
    @debug "size(values), length(c0), length(c1)" size(values) length(c0) length(c1)
    itp = interpolate((c0, c1, c2), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p0, p1, p2)
end

function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    _,
    c0::AbstractRange{<:Real},
    c1::AbstractRange{<:Real},
    _,
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end

function _interpolate(
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    _,
    c0::AbstractVector{<:Real},
    c1::AbstractVector{<:Real},
    _,
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end

function _interpolate(
    p0::AbstractArray{<:Real},
    _,
    p1::AbstractArray{<:Real},
    c0::AbstractRange{<:Real},
    _,
    c1::AbstractRange{<:Real},
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end

function _interpolate(
    p0::AbstractArray{<:Real},
    _,
    p1::AbstractArray{<:Real},
    c0::AbstractVector{<:Real},
    _,
    c1::AbstractVector{<:Real},
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end


function _interpolate(
    _,
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    _,
    c0::AbstractRange{<:Real},
    c1::AbstractRange{<:Real},
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end

function _interpolate(
    _,
    p0::AbstractArray{<:Real},
    p1::AbstractArray{<:Real},
    _,
    c0::AbstractVector{<:Real},
    c1::AbstractVector{<:Real},
    values::AbstractMatrix{<:Real},
)
    _interpolate(p0, p1, c0, c1, values)
end


