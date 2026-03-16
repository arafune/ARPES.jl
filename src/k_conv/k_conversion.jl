module KConversion
using DimensionalData
using DimensionalData: dims, hasdim, metadata, name, lookup
using ..ARPES: ARPESData, kx, ky, kz, phi, psi, eV
using ..ARPES: EnergyDefinition
using ..ARPES: BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy
using ..ARPES: AnalyzerConfiguration, TypeI, TypeII, TypeIp, TypeIIp
using ..ARPES: shift_dim, negate_dim

include("mapping.jl")
include("interpolation.jl")
include("preprocess.jl")
include("postprocess.jl")
export k_conversion

function k_conversion(
    data::ARPESData;
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    energy_definition = metadata(dims(data, :eV))[:energy_definition]
    if !(energy_definition in [BindingEnergy, FinalStateEnergy])
        throw(ArgumentError("Not Impremented for $energy_definition"))
    end
    @assert _check_arpesdata(data)
    # 1. required variables
    # 1.1. workfunction, photon energy and kinetic energy
    analyzer_conf = metadata(data)[:analyzer_configuration]
    workfunction = metadata(data)[:workfunction]

    #  (at present, kz conversion from the photon energy depencence is not implemented,
    #   so hv is not used in the conversion. It is included here for future extension.)
    hv = metadata(data)[:hv]

    eV_range = isnothing(eV_range) ? parent(lookup(data, :eV)) : eV_range
    ek_original = _ek_range(energy_definition, dims(data, :eV), workfunction, hv)
    ek = _ek_range(energy_definition, eV_range, workfunction, hv)

    @debug "Kinetic energy ref to analyzer" ek
    # 1.2. angles  * α, β_, χ_, δ_, ξ_
    #  * apply offset & and convert degree to radian.
    α =
        haskey(metadata(data), :negate_alpha) && metadata(data)[:negate_alpha] == true ?
        _deg2rad(parent(negate_dim(dims(data, :phi)))) :
        _deg2rad(parent(lookup(data, :phi)))
    β = hasdim(data, :psi) ? parent(dims(data, :psi)) : metadata(data)[:β]
    β_ = _deg2rad(β, β0)
    ξ_ = _deg2rad(metadata(data)[:ξ], ξ0)
    δ_ = _deg2rad(metadata(data)[:δ], δ0)
    χ_ = _deg2rad(metadata(data)[:χ], χ0)

    @debug "Kinetic energy ref to analyzer, and angles" ek α β_ χ_ ξ_ δ_
    # 2. determine k_regions, and use them if kx_range and ky_range are not provided.
    kx_range =
        isnothing(kx_range) ? _kx_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : kx_range
    ky_range =
        isnothing(ky_range) ? _ky_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : ky_range

    @debug "kx_range, ky_range" kx_range ky_range
    # 3. apply interpolation to get the intensity values on the k grid.
    #   3.1 corresponding \alpha and \beta
    kx_grid, ky_grid, ek_grid = reshape_for_nd(kx_range, ky_range, ek)

    @debug "size of kx_grid, ky_grid ek_grid, " size(kx_grid) size(ky_grid) size(ek_grid)
    α_range = mapped_α(analyzer_conf, kx_grid, ky_grid, ek_grid, _deg2rad(β0), χ_, ξ_, δ_)
    β_range = mapped_β(analyzer_conf, kx_grid, ky_grid, ek_grid, _deg2rad(β0), χ_, ξ_, δ_)
    @debug "size of α_range, β_range" size(α_range) size(β_range)
    #   3.2 interpolate
    data_k =
        _interpolate(α_range, β_range, ek_grid, α, β, parent(ek_original), parent(data))
    # 4. construct the output ARPESData object with kx, ky, and eV dimensions.u
    return _build_arpesband(
        data_k,
        kx_range,
        ky_range,
        eV_range,
        metadata(data),
        energy_definition,
    )
end

# --- internal functions

"""
    _deg2rad(angle_deg::AbstractRange, offset_deg::Real = 0.0) :: AbstractRange
    _deg2rad(angle_deg::T, offset_deg::Real = 0.0) where {T<:Union{AbstractArray,Real}} :: T

Convert an AbstractRange, AbstractArray or Real of angles from degrees to radians.

# Arguments
- `angle_deg`: An AbstractRange (AbstractArray and Real) of angles in degrees.
- `offset_deg`: An optional offset in degrees to be subtracted from the input angles
              before conversion. Default is `0.0`.

# Returns
- An AbstractRange (or AbstractArray or Real) of angles in radians.
"""
function _deg2rad(angle_deg::AbstractRange, offset_deg::Real = 0.0)
    return range(
        (first(angle_deg) - offset_deg) * (π/180);
        step = step(angle_deg) * (π/180),
        stop = (last(angle_deg) - offset_deg) * (π/180),
    )
end

function _deg2rad(angle_deg::T, offset_deg::Real = 0.0) where {T<:Union{AbstractArray,Real}}
    return @. (angle_deg - offset_deg) * (π / 180)
end

end # module
