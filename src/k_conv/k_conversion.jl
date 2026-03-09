module KConversion
using DimensionalData
using DimensionalData: dims, hasdim, dimnum
using ..ARPES: ARPESData, kx, ky, kz, phi, psi, eV
using ..ARPES: BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy
using ..ARPES: shift_dim

include("mapping.jl")
include("interpolation.jl")

export k_conversion

function k_conversion(
    data::ARPESData;
    kx_range::Union{AbstractArray,Nothing} = nothing,
    ky_range::Union{AbstractArray,Nothing} = nothing,
    kz_range::Union{AbstractArray,Nothing} = nothing,
    eV_range::Union{AbstractArray,Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert _check_arpesdata(data)
    # 1. required variables
    # 1.1. workfunction and kinetic energy

    workfunction = metadata(data)[:workfunction]

    #  (at present, kz conversion from the photon energy depencence is not implemented,
    #   so hv is not used in the conversion. It is included here for future extension.)
    hv = metadata(data)[:hv]

    if metadata(data)[:energy_definition] == FinalStateEnergy
        ke = shift_dim(dims(data)[dimnum(data, :eV)], workfunction)
    elseif metadata(data)[:energy_definition] == BindingEnergy
        ke = shift_dim(dims(data)[dimnum(data, :eV)], hv-workfunction)
    else
        throw(ArgumentError("Not Impremented for $(metadata(data)[:energy_definition])"))
    end

    # 1.2. angles  * α, β_, χ_, δ_, ξ_
    α = _deg2rad(parent(lookup(dims(data)[dimnum(data, :phi)])))
    if hasdim(data, :psi)
        β = parent(lookup(dims(data)[dimnum(data, :psi)]))
    else
        β = metadata(data)[:β]
    end
    β_ = _deg2rad(β .- β0)
    ξ_ = _deg2rad(metadata(data)[:ξ] - ξ0)
    δ_ = _deg2rad(metadata(data)[:δ] - δ0)
    χ_ = _deg2rad(metadata(data)[:χ] - χ0)

    # 1.3. apply necessary transformations to the angles if needed.
    #  * apply offset
    #  * degree -> radian
    # 2. determine k_region if kx, ky are not provided.
    # 3. apply interpolation to get the intensity values on the k grid.
end

function _check_arpesdata(data::ARPESData)
    # check if the required dimensions and metadata are present
    required_dims = [phi, eV]
    for dim in required_dims
        idx = findfirst(d -> name(d) == name(dim), dims(data))
        if idx === nothing
            throw(ArgumentError("Missing required dimension: $(name(dim))"))
        end
    end
    required_metadata = [:workfunction, :hv]
    for meta in required_metadata
        if !haskey(metadata(data), meta)
            throw(ArgumentError("Missing required metadata: $meta"))
        end
    end
    for meta in [:ξ, :δ, :χ]
        if !haskey(metadata(data), meta)
            metadata(data)[meta] = 0.0
            @warn "Missing metadata: $meta. Using default value of 0.0 for conversion."
        end
    end
    if !haskey(metadata(data), :β) && !hasdim(data, :psi)
        throw(
            ArgumentError(
                "Missing required metadata: either :β in metadata or :psi dim must be present",
            ),
        )
    end
    return true
end

"""
    _deg2rad(angle_in_degrees::AbstractRange) :: AbstractRange
    _deg2rad(angle_in_degrees::AbstractArray) :: AbstractArray
    _deg2rad(angle_in_degrees::Real):: Real

Convert an AbstractRange, AbstractArray or Real of angles from degrees to radians.

# Arguments
- `angle_in_degrees`: An AbstractRange of angles in degrees.

# Returns
- An AbstractRange of angles in radians.
"""
function _deg2rad(angle_in_degrees::AbstractRange)
    return angle_in_degrees * (π / 180)
end

function _deg2rad(angle_in_degrees::AbstractArray)
    return angle_in_degrees .* (π / 180)
end

function _deg2rad(angle_in_degrees::Real)::Real
    return angle_in_degrees * (π / 180)
end

end # module
