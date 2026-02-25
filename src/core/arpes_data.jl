using DimensionalData
import Base: size, axes, getindex, iterate
using DimensionalData.Dimensions: @dim
using Dates
import DimensionalData: dims, name

#--- Dimension names used for ARPES ---
@dim kx "1/Å"
@dim ky "1/Å"
@dim kz "1/Å"
@dim phi "degrees"
@dim psi "degrees"
@dim eV "eV"
@dim detector_ch
@dim ch2
@dim delay
@dim spin

"""
Abstract type representing the unit of intensity in ARPES data.
"""
abstract type IntensityUnit end

"""
Intensity unit representing raw counts.
"""
struct Counts <: IntensityUnit end

"""
Intensity unit representing counts per second (CPS).
"""
struct CPS <: IntensityUnit end

abstract type AnalyzerConfiguration end

struct TypeI <: AnalyzerConfiguration end
struct TypeII <: AnalyzerConfiguration end
struct TypeIp <: AnalyzerConfiguration end
struct TypeIIp <: AnalyzerConfiguration end

"""
Container for ARPES data.

# Type Parameters
- `T<:AbstractDimArray`: The type of the intensity data array.
- `U<:IntensityUnit`: The unit of intensity (e.g., Counts, CPS).
- `Conf`: The analyzer configuration type (should be a subtype of `AnalyzerConfiguration`).

# Fields
- `intensity::T`: The intensity values, stored as an AbstractDimArray.
- `unit::U`: The unit of intensity.
- `type::Conf`: The analyzer configuration type label.

# Constructor
The inner constructor checks that `type` is a subtype of `AnalyzerConfiguration`.
"""
s
# Fields
- `intensity::T`: The intensity values, stored as an AbstractDimArray.
- `unit::U`: The unit of intensity, subtype of IntensityUnit.
"""
struct ARPESData{T<:AbstractDimArray,U<:IntensityUnit,Conf}
    intensity::T
    unit::U
    type::Conf

    function ARPESData(intensity::T, unit::U, type::Conf) where {T,U,Conf}
        if !(type isa Type && type <: AnalyzerConfiguration)
            throw(
                ArgumentError("type must be a subtype of AnalyzerConfiguration, got $type"),
            )
        end
        new{T,U,Conf}(intensity, unit, type)
    end
end


# delegate methods
# -------------------

size(d::ARPESData) = size(d.intensity)
axes(d::ARPESData) = axes(d.intensity)
iterate(d::ARPESData, args...) = iterate(d.intensity, args...)
getindex(d::ARPESData, I...) = getindex(d.intensity, I...)
dims(d::ARPESData) = dims(d.intensity)
name(d::ARPESData) = name(d.intensity)
