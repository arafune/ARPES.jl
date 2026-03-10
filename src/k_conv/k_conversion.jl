module KConversion
using DimensionalData
using DimensionalData: dims, hasdim, dimnum, metadata, name, lookup
using ..ARPES: ARPESData, kx, ky, kz, phi, psi, eV
using ..ARPES: BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy
using ..ARPES: AnalyzerConfiguration, TypeI, TypeII, TypeIp, TypeIIp
using ..ARPES: shift_dim, negate_dim

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
    if !(metadata(data)[:energy_definition] in [BindingEnergy, FinalStateEnergy])
        throw(ArgumentError("Not Impremented for $(metadata(data)[:energy_definition])"))
    end
    @assert _check_arpesdata(data)
    # 1. required variables
    # 1.1. workfunction, photon energy and kinetic energy
    workfunction = metadata(data)[:workfunction]

    #  (at present, kz conversion from the photon energy depencence is not implemented,
    #   so hv is not used in the conversion. It is included here for future extension.)
    hv = metadata(data)[:hv]

    if eV_range !== nothing
        if metadata(data)[:energy_definition] == FinalStateEnergy
            ke = eV_range - workfunction
        elseif metadata(data)[:energy_definition] == BindingEnergy
            # negate_eV = ragne(-eV_range[1], step = -step(eV_range), stop = -eV_range[end])
            # ke = negate_eV - workfunction
            ke = hv + eV_range - workfunction
        end
    else
        if metadata(data)[:energy_definition] == FinalStateEnergy
            ke = shift_dim(dims(data)[dimnum(data, :eV)], workfunction)
        elseif metadata(data)[:energy_definition] == BindingEnergy
            # negaate_eV = negate_dim(dims(data)[dimnum(data, :eV)])
            # ke = shift_dim(negate_eV, hv - workfunction)
            ke = shift_dim(dims(data)[dimnum(data, :eV)], hv-workfunction)
        end
    end

    # 1.2. angles  * α, β_, χ_, δ_, ξ_
    #  * apply offset & and convert degree to radian.
    α = _deg2rad(parent(lookup(dims(data)[dimnum(data, :phi)])))
    if hasdim(data, :psi)
        β = parent(lookup(dims(data)[dimnum(data, :psi)]))
    else
        β = range(start = metadata(data)[:β], stop = metadata(data)[:β])
    end
    β_ = _deg2rad(β, β0)
    ξ_ = _deg2rad(metadata(data)[:ξ], ξ0)
    δ_ = _deg2rad(metadata(data)[:δ], δ0)
    χ_ = _deg2rad(metadata(data)[:χ], χ0)

    # 2. determine k_region if kx, ky are not provided.



    # 3. apply interpolation to get the intensity values on the k grid.

end

# --- internal functions

"""
  _determine_kxky_ragion(α, β_, ξ_, δ_, χ_, ke)


"""
function _determine_kxky_ragion(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    α::AbstractArray,
    β_::AbstractArray,
    ξ_,
    δ_,
    χ_,
    ke::AbstractArray,
)
    # 1. Determine the step of kx and ky  
    # 2. Deteremine the minimum and maximum of kx and ky
    #    k min is determined from lowest ke, while k max is determined from the highest ke.
    # 3. Return the kx and ky range  min_k: step_k: max_k + step_k
    return kx_range, ky_range
end


"""
    _check_arpesdata(data::ARPESData)

Check if the given `ARPESData` object contains the required dimensions and metadata for k-space conversion.

# Arguments
- `data::ARPESData`: The ARPES data object to check.

# Throws
- `ArgumentError` if any required dimension (`phi`, `eV`) or metadata (`:workfunction`, `:hv`) is missing.
- `ArgumentError` if neither `:β` metadata nor `:psi` dimension is present.

# Behavior
- If optional metadata (`:ξ`, `:δ`, `:χ`) is missing, sets them to `0.0` and emits a warning.

# Returns
- `true` if all required checks pass.
"""
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
    if !haskey(metadata(data), :β) && !hasdim(data, :psi)
        throw(
            ArgumentError(
                "Missing required metadata: either :β in metadata or :psi dim must be present",
            ),
        )
    end
    for meta in [:ξ, :δ, :χ]
        if !haskey(metadata(data), meta)
            metadata(data)[meta] = 0.0
            @warn "Missing metadata: $meta. Using default value of 0.0 for conversion."
        end
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
function _deg2rad(angle_in_degrees::AbstractRange, offset_deg::Real = 0.0)
    return ((first(angle_in_degrees)-offset_deg)*(π/180)):(step(angle_in_degrees)*(π/180)):((last(
        angle_in_degrees,
    )-offset_deg)*(π/180))
end

function _deg2rad(
    angle_in_degrees::T,
    offset_deg::Real = 0.0,
) where {T<:Union{AbstractArray,Real}}
    return @. (angle_in_degrees - offset_deg) * (π / 180)
end


end # module
