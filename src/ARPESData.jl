using DimensionalData
import Base: size, axes, getindex, iterate, parent
using DimensionalData.Dimensions: @dim
import DimensionalData: dims, name, metadata, rebuild

"""
Container for ARPES data.

# Type Parameters
- `A<:AbstractArray{T,N}`: The type of the data array (intensity values).
- `D`: The type of the dimension descriptor.
- `R`: The type of the reference dimension descriptor.
- `Na`: The type of the name field.
- `Me`: The type of the metadata field.

# Fields
- `data::A`: The intensity values, stored as an AbstractArray.
- `dims`: The dimension descriptor, in the ARPESData nomencalture of the substantial axes such as `eV`, `phi` are defined.
    * `eV` is the energy dimension, which can be defined as either binding energy or kinetic energy depending on the `energy_definition` in metadata.
      - Note: the `eV` dim should include the dict as metadata, which includes `:energy_definition` as the key.
        `energy_definition` should be one of the `EnergyDefinition` types (BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy).
    * `phi` is the emission angle of photoelectrons parallel to the slit direction.
    * `psi` is the emission angle of photoelectrons perpendicular to the slit direction.
    * `kx`, `ky`, `kz`, or `kp` are the momentum dimensions, which can be defined based on the angles and photon energy in metadata.
    * `delay` is the pump-probe delay time, which is used for time-resolved ARPES data.
- `name::Na`: The name label for the dataset.
- `metadata::Me`: Additional metadata as a dictionary or other structure.
  - metadata should include the following keys:
    * `:intensity_unit` - one of the `IntensityUnit` types (Counts, CPS)
    * `:analyzer_configuration` - one of the `AnalyzerConfiguration` types (TypeI, TypeII, etc.)
    * `:hv` - the photon energy in eV (should match the :hv in the data array metadata)
      (Note: in future, consider to use :hν instead of :hv)
    * `:history` - a list of transformations applied to the data, for provenance tracking.
    * `:β`: used for momentum conversion. if the dims include psi, angle `β` is not required.
    * `:ξ`: used for momentum conversion.
    * `:χ`: used for momentum conversion.
    * `:δ`: used for momentum conversion.  (About the definition of β, ξ, χ and δ, see Rev. Sci. Instrum. **89**, 043903 (2018).)
"""
struct ARPESData{T,N,D,R,A<:AbstractArray{T,N},Na,Me} <: AbstractDimArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    function ARPESData(
        data::A,
        dims::D,
        refdims::R,
        name::Na,
        metadata::Me,
    ) where {T,N,D,R,A<:AbstractArray{T,N},Na,Me}
        formatted_dims = DimensionalData.format(dims, data)
        return new{T,N,typeof(formatted_dims),R,A,Na,Me}(
            data,
            formatted_dims,
            refdims,
            name,
            metadata,
        )
    end
end

"""
  Create ARPESData from a raw intensity array and a dimension descriptor.

# Arguments
- data: N-dimensional array of intensity values.
- dims: Dimension descriptor compatible with DimensionalData.jl.

# Keywords
- name = :ARPES: Dataset name (Symbol or String).
- metadata = Dict(): Additional metadata (e.g., :hv, :analyzer_configuration, :energy_definition).
- refdims = nothing: Optional reference-dimension descriptor.

# Behavior
- Normalizes dims via DimensionalData.format(dims, data).
- Does not perform metadata or axis validation here.
"""
function ARPESData(data, dims; refdims = (), name = :ARPES, metadata = Dict())
    dims = DimensionalData.format(dims, data)
    ARPESData(data, dims, refdims, name, metadata)
end

"""
  Construct ARPESData from an AbstractDimArray.

# Arguments
- data::AbstractDimArray: Array with dims, refdims, name, and metadata.

# Behavior
- Uses parent(data) as the data field.
- Copies dims, refdims, name, and metadata from the given array.
- No additional validation is performed.
"""
function ARPESData(
    data::AbstractDimArray;
    intensity_unit::Type{<:IntensityUnit} = Counts,
    analyzer_config::Type{<:AnalyzerConfiguration} = TypeI,
    energy_def::EnergyDefinition = BindingEnergy,
    additional_metadata::AbstractDict = Dict(),
)
    extra_metadata =
        Dict(:analyzer_configuration => analyzer_config, :intensity_unit => intensity_unit)
    new_metadata = merge(Dict(metadata(data)), extra_metadata, Dict(additional_metadata))

    # Add :energy_definition to the metadata of the eV dimension
    new_dims = map(dims(data)) do d
        if d isa eV
            curr_dim_meta = Dict(DimensionalData.metadata(d))
            new_dim_meta = merge(curr_dim_meta, Dict(:energy_definition => energy_def))
            curr_lookup = lookup(d)
            new_lookup = rebuild(curr_lookup; metadata = new_dim_meta)
            return rebuild(d, new_lookup)
        else
            return d
        end
    end

    return ARPESData(
        parent(data),
        new_dims,
        refdims = refdims(data),
        name = data.name,
        metadata = new_metadata,
    )
end


# delegate methods
# -------------------
function rebuild(A::ARPESData, data, dims, refdims, name, metadata)
    ARPESData(data, dims, refdims, name, metadata)
end
Base.parent(A::ARPESData) = A.data
dims(A::ARPESData) = A.dims
name(A::ARPESData) = A.name
metadata(A::ARPESData) = A.metadata

Base.eltype(A::ARPESData) = eltype(parent(A))
Base.IndexStyle(::Type{<:ARPESData{T,N,D,R,A}}) where {T,N,D,R,A} = Base.IndexStyle(A)
Base.Broadcast.broadcastable(A::ARPESData) = Base.Broadcast.broadcastable(parent(A))

