using Interpolations

"""
    _interpolate(p1, p2, [p3], c1, c2, [c3], values)

Perform n-dimensional linear interpolation on a 2D or 3D grid. 

This function automatically dispatches to the most efficient algorithm based on the 
input types and dimensions:
- **Gridded(Linear)**: Used when coordinates are non-uniformly spaced (`AbstractVector`).
- **BSpline(Linear)**: Used when coordinates are uniformly spaced (`AbstractRange`), 
  providing significant performance gains.

# Requirements
- `p1`, `p2`, and (optionally) `p3` **must have the same shape and size**.
- `c1`, `c2`, and `c3` must be monotonic.
- `values` must match the dimensions of the grid defined by `c1`, `c2`, and `c3`.


# Arguments
- `p1, p2, [p3]`: `AbstractArray` of desired positions for each dimension.
- `c1, c2, [c3]`: `AbstractVector` (irregular) or `AbstractRange` (regular) 
  defining the original grid coordinates.
- `values`: `AbstractMatrix` (2D) or `AbstractArray{T, 3}` (3D) containing the 
  data values at grid points.

# Returns
- An array of interpolated values with the same shape as `p1`. 
  Points outside the grid boundaries are returned as `NaN`.
"""
function fast_interpolate end

# --- 2D Bilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    c1::AbstractRange{<:Real},
    c2::AbstractRange{<:Real},
    values::AbstractMatrix{<:Real},
)
    # Use BSpline interpolation for optimized performance on uniform grids
    @debug "size(values), length(c1), length(c2)" size(values) length(c1) length(c2)
    @debug "Input Range Check" extrema(p1) extrema(p2) extrema(c1) extrema(c2) size(values)
    values = step(c1) < 0 ? reverse(values, dims = 1) : values
    c1 = step(c1) < 0 ? reverse(c1) : c1
    values = step(c2) < 0 ? reverse(values, dims = 2) : values
    c2 = step(c2) < 0 ? reverse(c2) : c2
    @debug "Post-Normalization" c1_range=extrema(c1) c2_range=extrema(c2)

    @debug "c1, c2" c1 c2
    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c1, c2)
    res = extrapolate(sitp, NaN).(p1, p2)
    @debug "Result Summary" count_nan=count(isnan, res) total_elements=length(res)
    return res
end


# Irregular grid: coords are general vectors
function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    values::AbstractMatrix{<:Real},
)
    # Use Gridded interpolation for non-uniform spacing
    @debug "size(values), length(c1), length(c2)" size(values) length(c1) length(c2)
    c1 = collect(c1)
    c2 = collect(c2)
    itp = interpolate((c1, c2), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p1, p2)
end

# --- 3D Trilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    p3::AbstractArray{<:Real},
    c1::AbstractRange{<:Real},
    c2::AbstractRange{<:Real},
    c3::AbstractRange{<:Real},
    values::AbstractArray{T,3},
) where {T<:Real}
    @debug "size(values), length(c1), length(c2)" size(values) length(c1) length(c2)
    values = step(c1) < 0 ? reverse(values, dims = 1) : values
    c1 = step(c1) < 0 ? reverse(c1) : c1
    values = step(c2) < 0 ? reverse(values, dims = 2) : values
    c2 = step(c2) < 0 ? reverse(c2) : c2
    values = step(c3) < 0 ? reverse(values, dims = 3) : values
    c3 = step(c3) < 0 ? reverse(c3) : c3

    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c1, c2, c3)
    return extrapolate(sitp, NaN).(p1, p2, p3)
end

function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    p3::AbstractArray{<:Real},
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    c3::AbstractVector{<:Real},
    values::AbstractArray{T,3},
) where {T<:Real}
    @debug "size(values), length(c1), length(c2)" size(values) length(c1) length(c2)
    c1 = collect(c1)
    c2 = collect(c2)
    c3 = collect(c3)
    itp = interpolate((c1, c2, c3), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p1, p2, p3)
end

# --- 4D Quadrilinear Interpolation ---

# Regular grid: coords are ranges (e.g., start:step:stop)
function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    p3::AbstractArray{<:Real},
    p4::AbstractArray{<:Real},
    c1::AbstractRange{<:Real},
    c2::AbstractRange{<:Real},
    c3::AbstractRange{<:Real},
    c4::AbstractRange{<:Real},
    values::AbstractArray{T,4},
) where {T<:Real}
    @debug "size(values), length.(c1,c2,c3,c4)" size(values) length(c1) length(c2) length(c3) length(c4)
    values = step(c1) < 0 ? reverse(values, dims = 1) : values
    c1 = step(c1) < 0 ? reverse(c1) : c1
    values = step(c2) < 0 ? reverse(values, dims = 2) : values
    c2 = step(c2) < 0 ? reverse(c2) : c2
    values = step(c3) < 0 ? reverse(values, dims = 3) : values
    c3 = step(c3) < 0 ? reverse(c3) : c3
    values = step(c4) < 0 ? reverse(values, dims = 4) : values
    c4 = step(c4) < 0 ? reverse(c4) : c4

    itp = interpolate(values, BSpline(Linear()))
    sitp = scale(itp, c1, c2, c3, c4)
    return extrapolate(sitp, NaN).(p1, p2, p3, p4)
end

# Irregular grid: coords are general vectors
function _interpolate(
    p1::AbstractArray{<:Real},
    p2::AbstractArray{<:Real},
    p3::AbstractArray{<:Real},
    p4::AbstractArray{<:Real},
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    c3::AbstractVector{<:Real},
    c4::AbstractVector{<:Real},
    values::AbstractArray{T,4},
) where {T<:Real}
    @debug "size(values), length.(c1,c2,c3,c4)" size(values) length(c1) length(c2) length(c3) length(c4)
    c1 = collect(c1)
    c2 = collect(c2)
    c3 = collect(c3)
    c4 = collect(c4)
    itp = interpolate((c1, c2, c3, c4), values, Gridded(Linear()))
    return extrapolate(itp, NaN).(p1, p2, p3, p4)
end
