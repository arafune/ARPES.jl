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
function _apply_along_dim(A::AbstractDimArray, dim, f::Function)
    d = DimensionalData.dims(A, dim)
    dim_index = DimensionalData.dimnum(A, d)

    out = similar(parent(A), eltype(A))

    other_sizes = ntuple(i -> i == dim_index ? 1 : size(A, i), ndims(A))

    for I in CartesianIndices(other_sizes)
        idx = ntuple(i -> i == dim_index ? Colon() : I[i], ndims(A))
        y = Array(A[idx...])
        out[idx...] = f(y)
    end

    return rebuild(A, out)
end


