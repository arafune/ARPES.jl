using DimensionalData
import Base: size, axes, getindex, iterate
using DimensionalData.Dimensions: @dim
using Dates
import DimensionalData: dims, name, metadata

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

"""
Container for ARPES data.

# Fields
- `intensity::T`: The intensity values, stored as an AbstractDimArray.
- `unit::U`: The unit of intensity, subtype of IntensityUnit.
"""
struct ARPESData{T<:AbstractDimArray,U<:IntensityUnit}
    intensity::T
    unit::U
end

# delegate methods
# -------------------

size(d::ARPESData) = size(d.intensity)
axes(d::ARPESData) = axes(d.intensity)
iterate(d::ARPESData, args...) = iterate(d.intensity, args...)
getindex(d::ARPESData, I...) = getindex(d.intensity, I...)
dims(d::ARPESData) = dims(d.intensity)
name(d::ARPESData) = name(d.intensity)
metadata(d::ARPESData) = metadata(d.intensity)
