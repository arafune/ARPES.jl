_default_rtol(T) = T <: AbstractFloat ? √eps(T) : 0

function _is_equal_spacing(x::AbstractVector; atol = 0.0, rtol = _default_rtol(eltype(x)))
    if length(x) < 2
        return true
    end
    diffs = diff(x)
    return all(isapprox.(diff(diffs), 0; atol = atol, rtol = rtol))
end
