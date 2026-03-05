using DimensionalData
import Base: size, axes, getindex, iterate, parent
using DimensionalData.Dimensions: @dim
using Dates
import DimensionalData: dims, name, metadata, rebuild
using Makie
import Makie: convert_arguments

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
type representing the unit of intensity in ARPES data.
"""
abstract type IntensityUnit end
struct Counts <: IntensityUnit end
struct CPS <: IntensityUnit end



"""
type representing the analyzer configuration defined in Rev. Sci. Instrum. **89**, 043903 (2018).
"""
abstract type AnalyzerConfiguration end
struct TypeI <: AnalyzerConfiguration end
struct TypeII <: AnalyzerConfiguration end
struct TypeIp <: AnalyzerConfiguration end
struct TypeIIp <: AnalyzerConfiguration end

@enum EnergyDefinition begin
    BindingEnergy
    FinalStateEnergy
    KineticEnergy
    IntermediateEnergy
end

"""
Container for ARPES data.

# Type Parameters
- `T<:AbstractDimArray`: The type of the data array (intensity values).
- `D`: The type of the dimension descriptor.
- `N`: The type of the name field.
- `M`: The type of the metadata field.

# Fields
- `data::T`: The intensity values, stored as an AbstractDimArray.
    * The metadata of this array should include the following symbols and values:
        * :hv - the photon energy in eV.
        * :workfunction - the work function of the analyzer in eV.
    * The energy axis must be eV, and the emission angles must be phi and psi.
      (phi is the emission angle of photoelectrons parallel to the slit direction.)
    * If the data is after momentum conversion, the dimensions should be kx, ky, kz, kp.
- `dims::D`: The dimension descriptor, typically matching the dims of the data array.
- `name::N`: The name label for the dataset.
- `metadata::M`: Additional metadata as a dictionary or other structure.
  - metadata should include the following keys:
    * :analyzer_configuration - one of the AnalyzerConfiguration types (TypeI, TypeII, etc.)
    * :energy_definition - one of the EnergyDefinition enum values (BindingEnergy, FinalStateEnergy, etc.)
    * :hv - the photon energy in eV (should match the :hv in the data array metadata)  (Note: in future, consider to use hv instead of hv)
"""
struct ARPESData{T,D,R,N,M}
    data::T
    dims::D
    refdims::R
    name::N
    metadata::M

end



# delegate methods
# -------------------

Base.parent(A::ARPESData) = A.intensity
dims(A::ARPESData) = A.dims
name(A::ARPESData) = A.name
metadata(A::ARPESData) = A.metadata

size(d::ARPESData) = size(d.intensity)
axes(d::ARPESData) = axes(d.intensity)
iterate(d::ARPESData, args...) = iterate(d.intensity, args...)
getindex(d::ARPESData, I...) = getindex(d.intensity, I...)

Base.ndims(A::ARPESData) = ndims(parent(A))
Base.IndexStyle(::Type{<:ARPESData{A}}) where {A} = Base.IndexStyle(A)
Base.Broadcast.broadcastable(A::ARPESData) = Base.Broadcast.broadcastable(parent(A))

function ARPESData(intensity, dims; name = :ARPES, metadata = Dict(), refdims = nothing)
    dims = DimensionalData.format(dims, intensity)
    return ARPESData(intensity, dims, refdims, name, metadata)
end

function ARPESData(intensity::AbstractDimArray)
    return ARPESData(
        parent(intensity),
        intensity.dims,
        refdims = intensity.refdims,
        name = intensity.name,
        metadata = intensity.metadata,
    )
end

function Makie.convert_arguments(P::Type{<:Makie.AbstractPlot}, data::ARPESData)
    return Makie.convert_arguments(P, parent(data))
end
