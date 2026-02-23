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
@dim Cycle "cycle"
@dim Ch2 "ch2"
@dim delay

abstract type IntensityUnit end
struct Counts <: IntensityUnit end
struct CPS <: IntensityUnit end

struct ARPESData{T<:AbstractDimArray,U<:IntensityUnit}
    intensity::T
    unit::U
end

