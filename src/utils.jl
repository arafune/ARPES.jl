using DimensionalData
using DimensionalData: AbstractDimArray, rebuild

_default_rtol(T) = T <: AbstractFloat ? √eps(T) : 0

function _is_equal_spacing(x::AbstractVector; atol = 0.0, rtol = _default_rtol(eltype(x)))
    if length(x) < 2
        return true
    end
    diffs = diff(x)
    return all(isapprox.(diff(diffs), 0; atol = atol, rtol = rtol))
end

# -------------------------------
# Utility: apply 1D function along a dimension
# -------------------------------
"""
    _apply_along_dim(A::AbstractDimArray, dim, f::Function)

Apply `f` to each 1D slice of `A` along dimension `dim` and write the returned
slice back into the corresponding location.

Assumes `f(y)` returns a vector/array with the same length as `y`.
Keeps the DimArray metadata via `rebuild(A, out)`.
"""
function _apply_along_dim(A::AbstractDimArray, dim, f::Function)
    d = DimensionalData.dims(A, dim)
    dim_index = DimensionalData.dimnum(A, d)

    # Keep the same parent-array type/shape as A
    out = similar(parent(A))

    # Iterate over all indices of dimensions other than `dim_index`
    other_axes = ntuple(i -> (i == dim_index ? Base.OneTo(1) : axes(A, i)), ndims(A))

    @views for I in CartesianIndices(other_axes)
        # Build an index like (:, y, z) when dim_index==1
        idx = ntuple(i -> (i == dim_index ? Colon() : I[i]), ndims(A))

        # Take the 1D slice (no copy), then pass a Vector copy to f (safe)
        y = vec(copy(A[idx...]))

        r = f(y)

        # Write back into the slice
        out[idx...] .= r
    end

    return rebuild(A, out)
end
