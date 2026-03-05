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

"""
    @enum EnergyDefinition

Enumeration of possible energy definitions used in ARPES data analysis.

- `BindingEnergy`: Electron binding energy relative to the Fermi level.
- `FinalStateEnergy`: Energy of the electron in the final state after photoemission (Referred to the Fermi level).
- `KineticEnergy`: Kinetic energy of the emitted electron (Referred to the vacuum level of the sample).
- `IntermediateEnergy`: Energy in an intermediate state (e.g., in pump-probe experiments).
"""
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
- `intensity::T`: The intensity values, stored as an AbstractArray.
- `dims::D`: The dimension descriptor, typically matching the dims of the data array.
    * The energy axis must be eV, and the emission angles must be phi and psi.
      (phi is the emission angle of photoelectrons parallel to the slit direction.)
    * If the data is after momentum conversion, the dimensions should be kx, ky, kz, kp.
- `name::N`: The name label for the dataset.
- `metadata::M`: Additional metadata as a dictionary or other structure.
  - metadata should include the following keys:
    * :intensity_unit - one of the IntensityUnit types (Counts, CPS)
    * :analyzer_configuration - one of the AnalyzerConfiguration types (TypeI, TypeII, etc.)
    * :energy_definition - one of the EnergyDefinition enum values (BindingEnergy, FinalStateEnergy, etc.)
    * :hv - the photon energy in eV (should match the :hv in the data array metadata)  (Note: in future, consider to use :hν instead of :hv)
"""
struct ARPESData{T,D,R,N,M}
    intensity::T
    dims::D
    refdims::R
    name::N
    metadata::M
end

"""
  Create ARPESData from a raw intensity array and a dimension descriptor.

# Arguments
- intensity: N-dimensional array of intensity values.
- dims: Dimension descriptor compatible with DimensionalData.jl.

# Keywords
- name = :ARPES: Dataset name (Symbol or String).
- metadata = Dict(): Additional metadata (e.g., :hv, :analyzer_configuration, :energy_definition).
- refdims = nothing: Optional reference-dimension descriptor.

# Behavior
- Normalizes dims via DimensionalData.format(dims, intensity).
- Does not perform metadata or axis validation here.
"""
function ARPESData(intensity, dims; name = :ARPES, metadata = Dict(), refdims = nothing)
    dims = DimensionalData.format(dims, intensity)
    return ARPESData(intensity, dims, refdims, name, metadata)
end

"""
  Construct ARPESData from an AbstractDimArray.

# Arguments
- intensity::AbstractDimArray: Array with dims, refdims, name, and metadata.

# Behavior
- Uses parent(intensity) as the data field.
- Copies dims, refdims, name, and metadata from the given array.
- No additional validation is performed.
"""
function ARPESData(
    intensity::AbstractDimArray;
    intensity_unit::Type{<:IntensityUnit} = Counts,
    analyzer_config::Type{<:AnalyzerConfiguration} = TypeI,
    energy_def::EnergyDefinition = BindingEnergy,
    additional_metadata::AbstractDict = Dict(),
)
    extra_metadata = Dict(
        :analyzer_configuration => analyzer_config,
        :energy_definition => energy_def,
        :intensity_unit => intensity_unit,
    )
    new_metadata =
        merge(Dict(metadata(intensity)), extra_metadata, Dict(additional_metadata))
    return ARPESData(
        parent(intensity),
        intensity.dims,
        refdims = intensity.refdims,
        name = intensity.name,
        metadata = new_metadata,
    )
end

function Makie.convert_arguments(P::Type{<:Makie.AbstractPlot}, data::ARPESData)
    return Makie.convert_arguments(P, parent(data))
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

Base.eltype(A::ARPESData) = eltype(parent(A))
Base.ndims(A::ARPESData) = ndims(parent(A))
Base.IndexStyle(::Type{<:ARPESData{A}}) where {A} = Base.IndexStyle(A)
Base.Broadcast.broadcastable(A::ARPESData) = Base.Broadcast.broadcastable(parent(A))


