using DimensionalData
using DimensionalData: Dimension
using DimensionalData: basetypeof
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

export trapezoid


"""
    trapezoid(A, trapezoid_corners, base_corners; from_trapezoid=true)

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
- `trapezoid_corners`: The coordinate of the trapezoid corners.
   The tuple of four dict that specifies corners of the trapezoid, the key must be both `:eV` and `:phi`.
- `base_corners`: The tuple of two the phi values of the trensposed rectangle corners. (i.e. L_Rect and R_Rect).
   if not specified (None), use the extrema(A, :phi). Defaults to None.
   (As the `eV` axis does not change and only `phi` axis changes, specifying L_Rect and R_Rect is enough.)

# Returns:

A new array converted.
"""

function trapezoid(
    A::ARPESData{T,2} where {T},
  trapezoid_corners::NTuple{4,NamedTuple{(:eV, :phi),Tuple{Float64,Float64}}},  # (UL, UR, LL, LR)
  base_corners::Union{Tuple{Real, Real}, Nothing} = nothing;,
) end


function _rectangle_to_trapezoid(p,
                                 Ek,
                                 trapezoid_corners::NTuple{4,NamedTuple{(:eV, :phi),Tuple{Float64,Float64}}},  # (UL, UR, LL, LR)
                                 base_corners::Tuple{Real, Real})
  UL, UR, LL, LR = trapezoid_corners
  slope_left_edge = UL.phi - LL.phi / (UL.eV - LL.eV)
  slope_right_edge = UR.phi - LR.phi / (UR.eV - LR.eV)
  left_edge = @. slope_left_edge * (Ek - UL.eV) + UL.phi
  right_edge = @. slope_right_edge + (Ek- UR.eV) + UR.phi
  dac_da = (right_edge - left_edge) / (maximim(base_corners) - minimum(base_corners))
  return @. (p - minimum(base_corners)) * dac_da+ left_edge
end

function _traezoid_to_rectangle()
end

