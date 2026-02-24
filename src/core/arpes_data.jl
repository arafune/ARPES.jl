using DimensionalData
using DimensionalData.Dimensions: @dim
using Dates

#--- Dimension names used for ARPES ---
@dim kx "1/Å"
@dim ky "1/Å"
@dim kz "1/Å"
@dim phi "degrees"
@dim psi "degrees"
@dim eV "eV"
@dim ch
@dim ch2
@dim delay
@dim spin

abstract type IntensityUnit end
struct Counts <: IntensityUnit end
struct CPS <: IntensityUnit end

struct ARPESData{T<:AbstractDimArray,U<:IntensityUnit}
    intensity::T
    unit::U
end

