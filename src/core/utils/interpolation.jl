"""Bilinear interpolation."""
using Base.Threads

"""
    fast_bilinear_interpolate(
        desired_pos_dim0::AbstractArray,
        desired_pos_dim1::AbstractArray,
        orig_coords_dim0::AbstractVector,
        orig_coords_dim1::AbstractVector,
        orig_values::AbstractMatrix,
    ) -> AbstractArray

Performs bilinear interpolation on a 2D grid for arbitrary (not necessarily rectilinear) axes.

# Arguments
- `desired_pos_dim0::AbstractArray`: Desired positions along the first dimension (x-axis) to interpolate.
- `desired_pos_dim1::AbstractArray`: Desired positions along the second dimension (y-axis) to interpolate.
- `orig_coords_dim0::AbstractVector`: Original grid coordinates along the first dimension (x-axis), must be sorted.
- `orig_coords_dim1::AbstractVector`: Original grid coordinates along the second dimension (y-axis), must be sorted.
- `orig_values::AbstractMatrix`: Matrix of values defined on the original grid.

# Returns
- `AbstractArray`: Interpolated values at the desired positions, with the same shape as the input position arrays.

Returns `NaN` for points outside the original grid.

Threaded for performance.
"""
using Base.Threads

function fast_bilinear_interpolate(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_values::AbstractMatrix,
)
    # Check axis order and set flags
    rev_dim0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev_dim1 = orig_coords_dim1[1] > orig_coords_dim1[end]

    coords_dim0 = rev_dim0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords_dim1 = rev_dim1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    values = orig_values
    if rev_dim0
        values = reverse(values, dims = 1)
    end
    if rev_dim1
        values = reverse(values, dims = 2)
    end

    len_dim0 = length(coords_dim0)
    len_dim1 = length(coords_dim1)
    desired_shape = size(desired_pos_dim0)
    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    n_points = length(xs)
    result = Vector{eltype(orig_values)}(undef, n_points)

    @threads for idx = 1:n_points
        x = xs[idx]
        y = ys[idx]

        x1_idx = searchsortedfirst(coords_dim0, x) - 1
        x2_idx = x1_idx + 1
        y1_idx = searchsortedfirst(coords_dim1, y) - 1
        y2_idx = y1_idx + 1

        if x1_idx < 1 || x2_idx > len_dim0 || y1_idx < 1 || y2_idx > len_dim1
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1 = coords_dim0[x1_idx]
            x2 = coords_dim0[x2_idx]
            y1 = coords_dim1[y1_idx]
            y2 = coords_dim1[y2_idx]

            Q11 = values[x1_idx, y1_idx]
            Q12 = values[x1_idx, y2_idx]
            Q21 = values[x2_idx, y1_idx]
            Q22 = values[x2_idx, y2_idx]

            denom = (x2 - x1) * (y2 - y1)
            result[idx] =
                (
                    Q11 * (x2 - x) * (y2 - y) +
                    Q21 * (x - x1) * (y2 - y) +
                    Q12 * (x2 - x) * (y - y1) +
                    Q22 * (x - x1) * (y - y1)
                ) / denom
        end
    end

    return reshape(result, desired_shape)
end


"""
    fast_bilinear_interpolate_rectilinear(
        desired_pos_dim0::AbstractArray,
        desired_pos_dim1::AbstractArray,
        orig_coords_dim0::AbstractVector,
        orig_coords_dim1::AbstractVector,
        orig_values::AbstractMatrix,
    ) -> AbstractArray

Performs fast bilinear interpolation on a 2D rectilinear grid (uniformly spaced axes).

# Arguments
- `desired_pos_dim0::AbstractArray`: Desired positions along the first dimension (x-axis) to interpolate.
- `desired_pos_dim1::AbstractArray`: Desired positions along the second dimension (y-axis) to interpolate.
- `orig_coords_dim0::AbstractVector`: Original grid coordinates along the first dimension (x-axis), must be sorted and uniformly spaced.
- `orig_coords_dim1::AbstractVector`: Original grid coordinates along the second dimension (y-axis), must be sorted and uniformly spaced.
- `orig_values::AbstractMatrix`: Matrix of values defined on the original grid.

# Returns
- `AbstractArray`: Interpolated values at the desired positions, with the same shape as the input position arrays.

Returns `NaN` for points outside the original grid.

Threaded for performance.
"""
function fast_bilinear_interpolate_rectilinear(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_values::AbstractMatrix,
)
    rev_dim0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev_dim1 = orig_coords_dim1[1] > orig_coords_dim1[end]

    coords_dim0 = rev_dim0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords_dim1 = rev_dim1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    values = orig_values
    if rev_dim0
        values = reverse(values, dims = 1)
    end
    if rev_dim1
        values = reverse(values, dims = 2)
    end

    len_dim0 = length(coords_dim0)
    len_dim1 = length(coords_dim1)
    desired_shape = size(desired_pos_dim0)
    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    n_points = length(xs)
    result = Vector{eltype(orig_values)}(undef, n_points)

    step_dim0 = (coords_dim0[end] - coords_dim0[1]) / (len_dim0 - 1)
    step_dim1 = (coords_dim1[end] - coords_dim1[1]) / (len_dim1 - 1)
    x0 = coords_dim0[1]
    y0 = coords_dim1[1]

    @threads for idx = 1:n_points
        x = xs[idx]
        y = ys[idx]

        x1_idx = floor(Int, (x - x0) / step_dim0) + 1
        x2_idx = x1_idx + 1
        y1_idx = floor(Int, (y - y0) / step_dim1) + 1
        y2_idx = y1_idx + 1

        if x1_idx < 1 || x2_idx > len_dim0 || y1_idx < 1 || y2_idx > len_dim1
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1 = coords_dim0[x1_idx]
            x2 = coords_dim0[x2_idx]
            y1 = coords_dim1[y1_idx]
            y2 = coords_dim1[y2_idx]

            Q11 = values[x1_idx, y1_idx]
            Q12 = values[x1_idx, y2_idx]
            Q21 = values[x2_idx, y1_idx]
            Q22 = values[x2_idx, y2_idx]

            denom = (x2 - x1) * (y2 - y1)
            result[idx] =
                (
                    Q11 * (x2 - x) * (y2 - y) +
                    Q21 * (x - x1) * (y2 - y) +
                    Q12 * (x2 - x) * (y - y1) +
                    Q22 * (x - x1) * (y - y1)
                ) / denom
        end
    end

    return reshape(result, desired_shape)
end


