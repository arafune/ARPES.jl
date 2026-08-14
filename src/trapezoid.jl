using DimensionalData
using DimensionalData: Dimension
using DimensionalData: basetypeof
using DimensionalData: metadata, name
using DimensionalData.Lookups
using .KConversion: prepare_for_broadcast, _interpolate
export trapezoid

const TrapezoidCorner = NamedTuple{(:eV, :phi),Tuple{Float64,Float64}}
const TrapezoidCorners = NTuple{4,TrapezoidCorner}
const BaseCorners = Tuple{Real,Real}


"""
    trapezoid(A, trapezoid_corners, base_corners)
    trapezoid(A, base_corners, trapezoid_corners)

Applies the trapezoidal correction in angular units by linearly interpolating slices.

This function shares some code with standard coordinate conversion routines, such as
those used for momentum conversion, because it can be viewed as a coordinate
conversion between two angular coordinate systems: the measured angles and the true
angles.



          (UL)_____________ (UR)                 +--------+
        ↑     \\           /                      |        |
        |      \\         /        ⇄              |        |
        eV      \\_______/               (L_Rect) +--------+  (R_Rect)
            (LL)          (LR)

                ----------→ phi

# Argumetns

- `A`: The input `ARPESData` to be transformed.
   - `trapezoid_corners`: The coordinate of the trapezoid corners (UL, UR, LL, LR).
   The tuple of four Named tuple that specifies corners of the trapezoid, the key must be both `:eV` and `:phi`.
   The corner order is inferred automatically from the `(eV, phi)` values, so the input does not need to
   already be ordered as `(UL, UR, LL, LR)`.
- `base_corners`: The tuple of two the phi values of the trensposed rectangle corners. (i.e. L_Rect and R_Rect).
   if not specified (None), use the extrema(A, :phi). Defaults to None.
   (As the `eV` axis does not change and only `phi` axis changes, specifying L_Rect and R_Rect is enough.)

- If `base_corners` comes first, transforms from rectangle to trapezoid.
- If `trapezoid_corners` comes first, transforms from trapezoid to rectangle.

# Returns:

A new array converted.
"""
trapezoid(
    A::ARPESData{T,2} where {T},
    trapezoid_corners::TrapezoidCorners,  # (UL, UR, LL, LR)
    base_corners::Union{BaseCorners,Nothing} = nothing,
) = _trapezoid(A, trapezoid_corners, base_corners, Val(:to_rectangle))

trapezoid(
    A::ARPESData{T,2} where {T},
    base_corners::BaseCorners,
    trapezoid_corners::TrapezoidCorners,  # (UL, UR, LL, LR)
) = _trapezoid(A, trapezoid_corners, base_corners, Val(:to_trapezoid))

function _trapezoid(
    A::ARPESData{T,2} where {T},
    trapezoid_corners::TrapezoidCorners,
    base_corners::Nothing,
    direction::Val{:to_rectangle},
)
    return _trapezoid(A, trapezoid_corners, extrema(lookup(A, :phi)), direction)
end

function _trapezoid(
    A::ARPESData{T,2} where {T},
    trapezoid_corners::TrapezoidCorners,
    base_corners::BaseCorners,
    direction::Union{Val{:to_rectangle},Val{:to_trapezoid}},
)
    _validate_trapezoid_input(A)
    trapezoid_corners = _normalize_trapezoid_corners(trapezoid_corners)
    _validate_base_corners(base_corners)

    phi_lookup = lookup(A, :phi)
    eV_lookup = lookup(A, :eV)
    step_phi = _step(dims(A, :phi))
    min_to_phi, max_to_phi = _output_phi_bounds(
        phi_lookup,
        eV_lookup,
        trapezoid_corners,
        base_corners,
        direction,
    )
    @debug "max_to_phi: $max_to_phi, min_to_phi: $min_to_phi"
    to_phi_range = range(start = min_to_phi, step = abs(step_phi), stop = max_to_phi)
    to_phi_grid, ek_grid = prepare_for_broadcast(to_phi_range, parent(eV_lookup))
    from_phi_grid =
        _source_phi_grid(to_phi_grid, ek_grid, trapezoid_corners, base_corners, direction)
    data_transposed = _interpolate(
        from_phi_grid,
        ek_grid,
        parent(phi_lookup),
        parent(eV_lookup),
        parent(A),
    )
    return _rebuild_with_phi(A, to_phi_range, data_transposed)
end

function _validate_trapezoid_input(A)
    @assert hasdim(A, :eV) && hasdim(A, :phi) " A must have dimensions :eV and :phi"
    @assert _is_equal_spacing(lookup(A, :phi)) "A must be equally spaced along :phi dimension"
end

function _normalize_trapezoid_corners(trapezoid_corners::TrapezoidCorners)
    sorted_corners =
        sort(collect(trapezoid_corners); by = corner -> (corner.eV, corner.phi))
    @assert sorted_corners[1].eV == sorted_corners[2].eV "Two trapezoid corners must share the same upper-edge eV value"
    @assert sorted_corners[3].eV == sorted_corners[4].eV "Two trapezoid corners must share the same lower-edge eV value"
    @assert sorted_corners[1].eV != sorted_corners[3].eV "Upper and lower trapezoid edges must have different eV values"
    return tuple(sorted_corners...)
end

function _validate_base_corners(base_corners::BaseCorners)
    L_rect, R_rect = base_corners
    @assert L_rect < R_rect "L_Rect must be less than R_Rect"
end

function _output_phi_bounds(
    phi_lookup,
    eV_lookup,
    trapezoid_corners::TrapezoidCorners,
    base_corners::BaseCorners,
    direction::Union{Val{:to_rectangle},Val{:to_trapezoid}},
)
    min_input_phi = minimum(phi_lookup)
    max_input_phi = maximum(phi_lookup)
    min_to_phi = minimum(
        _phi_to_phi(min_input_phi, eV_lookup, trapezoid_corners, base_corners, direction),
    )
    max_to_phi = maximum(
        _phi_to_phi(max_input_phi, eV_lookup, trapezoid_corners, base_corners, direction),
    )
    return min_to_phi, max_to_phi
end

_phi_to_phi(p, Ek, trapezoid_corners, base_corners, ::Val{:to_rectangle}) =
    _phi_to_phi(p, Ek, trapezoid_corners, base_corners)

_phi_to_phi(p, Ek, trapezoid_corners, base_corners, ::Val{:to_trapezoid}) =
    _phi_to_phi(p, Ek, base_corners, trapezoid_corners)

_source_phi_grid(
    to_phi_grid,
    ek_grid,
    trapezoid_corners,
    base_corners,
    ::Val{:to_rectangle},
) = _phi_to_phi(to_phi_grid, ek_grid, base_corners, trapezoid_corners)

_source_phi_grid(
    to_phi_grid,
    ek_grid,
    trapezoid_corners,
    base_corners,
    ::Val{:to_trapezoid},
) = _phi_to_phi(to_phi_grid, ek_grid, trapezoid_corners, base_corners)

function _rebuild_with_phi(A, to_phi_range, data_transposed)
    @debug "trapezoid output check" size(data_transposed) length(to_phi_range)
    new_phi_dim = Dim{:phi}(to_phi_range; metadata = metadata(dims(A, :phi)))
    new_dims = Base.setindex(dims(A), new_phi_dim, dimnum(A, :phi))

    return rebuild(A; data = data_transposed, dims = new_dims)
end

function trapezoid(
    A::ARPESData{T,2} where {T},
    base_corners::Tuple{Real,Real},
    trapezoid_corners::NTuple{4,NamedTuple{(:phi, :eV),Tuple{Float64,Float64}}},  # (UL, UR, LL, LR)
)
    reorddered_corners = (
        (eV = trapezoid_corners[1].eV, phi = trapezoid_corners[1].phi),
        (eV = trapezoid_corners[2].eV, phi = trapezoid_corners[2].phi),
        (eV = trapezoid_corners[3].eV, phi = trapezoid_corners[3].phi),
        (eV = trapezoid_corners[4].eV, phi = trapezoid_corners[4].phi),
    )
    return trapezoid(A, base_corners, reorddered_corners)
end

function trapezoid(
    A::ARPESData{T,2} where {T},
    trapezoid_corners::NTuple{4,NamedTuple{(:phi, :eV),Tuple{Float64,Float64}}},  # (UL, UR, LL, LR)
    base_corners::Tuple{Real,Real},
)
    reorddered_corners = (
        (eV = trapezoid_corners[1].eV, phi = trapezoid_corners[1].phi),
        (eV = trapezoid_corners[2].eV, phi = trapezoid_corners[2].phi),
        (eV = trapezoid_corners[3].eV, phi = trapezoid_corners[3].phi),
        (eV = trapezoid_corners[4].eV, phi = trapezoid_corners[4].phi),
    )
    return trapezoid(A, reorddered_corners, base_corners)
end

"""
    _phi_to_phi(p, Ek, trapezoid_corners, base_corners)   # rectangle => trapezoid
    _phi_to_phi(p, Ek, base_corners, trapezoid_corners)   # trapezoid => rectangle

Transforms coordinates between a trapezoidal region (defined by four corners in (eV, phi) space)
and a rectangular region (defined by two base corners). The direction of the transformation is
determined by the argument order:

- If `base_corners` comes first, transforms from rectangle to trapezoid.
- If `trapezoid_corners` comes first, transforms from trapezoid to rectangle.

Arguments:
- `p`: phi coordinate(s) to transform
- `Ek`: energy value(s)
- `trapezoid_corners`: NTuple of 4 NamedTuples (UL, UR, LL, LR) with fields `:eV` and `:phi`
- `base_corners`: Tuple of two real numbers defining the rectangle base

Returns:
- Transformed coordinate(s) in the target region.
"""
_phi_to_phi(
    p,
    Ek,
    trapezoid_corners::TrapezoidCorners,  # (UL, UR, LL, LR)
    base_corners::BaseCorners,
) = _phi_to_phi_impl(
    p,
    Ek,
    _normalize_trapezoid_corners(trapezoid_corners),
    base_corners,
    true,
)

_phi_to_phi(
    p,
    Ek,
    base_corners::BaseCorners,
    trapezoid_corners::TrapezoidCorners,  # (UL, UR, LL, LR)
) = _phi_to_phi_impl(
    p,
    Ek,
    _normalize_trapezoid_corners(trapezoid_corners),
    base_corners,
    false,
)

function _phi_to_phi_impl(
    p,
    Ek,
    trapezoid_corners::TrapezoidCorners,  # (UL, UR, LL, LR)
    base_corners::BaseCorners,
    from_trapezoid::Bool,
)
    @debug "p: typeof: $(typeof(p)), length: $(length(p)), extrema: $(extrema(p))"
    @debug "base_corners: $(base_corners)"
    UL, UR, LL, LR = trapezoid_corners
    left_corner, right_corner = extrema(base_corners)
    slope_left_edge = (UL.phi - LL.phi) / (UL.eV - LL.eV)
    slope_right_edge = (UR.phi - LR.phi) / (UR.eV - LR.eV)
    left_edge = @. slope_left_edge * (Ek - UL.eV) + UL.phi
    @debug "left_edge: length: $(length(left_edge)), extrema: $(extrema(left_edge))"
    right_edge = @. slope_right_edge * (Ek - UR.eV) + UR.phi
    @debug "right_edge: length: $(length(right_edge)), extrema: $(extrema(right_edge))"
    if from_trapezoid # from trapezoid to rectangle
        c = @. (p - left_edge) / (right_edge - left_edge)
        @debug "c: length: $(length(c)), extrema: $(extrema(c))"
        return @. left_corner + c * (right_corner - left_corner)
    else  # from rectangle to trapezoid 
        dac_da = (right_edge - left_edge) / (right_corner - left_corner)
        @debug "dac_da: length: $(length(dac_da)), extrema: $(extrema(dac_da))"
        return @. (p - left_corner) * dac_da + left_edge
    end
end
