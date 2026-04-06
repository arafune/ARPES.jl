module KConversion
using DimensionalData
using DimensionalData: dims, hasdim, metadata, name, lookup
using ..ARPES: ARPESData, kx, ky, kz, phi, psi, eV
using ..ARPES: EnergyDefinition
using ..ARPES: BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy
using ..ARPES:
    AbstractAnalyzerConfiguration,
    AbstractAnalyzerConfigurationWithoutDeflector,
    AbstractAnalyzerConfigurationWithDeflector,
    TypeI,
    TypeII,
    TypeIp,
    TypeIIp
using ..ARPES: shift_dim, negate_dim, add_dim


include("mapping.jl")
include("interpolation.jl")
include("preprocess.jl")
include("postprocess.jl")
export k_conversion

"""
    k_conversion(data::ARPESData; kx_range=nothing, ky_range=nothing, kz_range=nothing,
                 eV_range=nothing, β0=0.0, ξ0=0.0, δ0=0.0, χ0=0.0)

Convert ARPES data from angle-space to momentum-space coordinates.

Dispatch is determined by the number of dimensions and which named dimensions are present:

| Input dims | Key dims present             | Output         | Internal call       |
|:---------- |:---------------------------- |:-------------- |:------------------- |
| N=2        | `phi`, `eV`                  | `kx × eV`     | (direct)            |
| N=3        | `phi`, `psi`, `eV`           | `kx × ky × eV`| `kxky_conversion`   |
| N=3        | `phi`, `hv`/`hν`, `eV`      | (not impl.)    | `kpkz_conversion`   |
| N=3        | other                        | slices extra dim | `_slice_and_convert`|
| N=4        | `phi`, `psi`, `hv`/`hν`, `eV`| (not impl.)  | `kxkykz_conversion` |
| N=4        | `phi`, `psi`, `eV` + extra   | slices extra dim | `_slice_and_convert`|
| N=4        | `phi`, `hv`/`hν`, `eV` + extra| slices extra dim| `_slice_and_convert`|
| N≥5        | any                          | slices last extra dim | `_slice_and_convert`|

The input must contain `eV` and `phi` dimensions. The dataset metadata is expected to
include `:analyzer_configuration`, `:workfunction`, `:hv`, `:β`, `:ξ`, `:δ`, and `:χ`.
The `eV` dimension metadata must include `:energy_definition`.

# Keywords
- `kx_range`, `ky_range`, `kz_range`: Optional target momentum axes. If omitted, the code
  estimates suitable in-plane momentum ranges from the input geometry.
- `eV_range`: Optional target energy axis. Defaults to the input `eV` lookup.
- `β0`, `ξ0`, `δ0`, `χ0`: Additional angular offsets in degrees applied during mapping.

# Returns
- `ARPESData`: A new dataset sampled on momentum-space axes, preserving the dataset-level
  metadata from the input.

# Notes
- Currently supports `BindingEnergy` and `FinalStateEnergy` energy definitions.
- `kpkz_conversion` and `kxkykz_conversion` are defined but not yet implemented.
"""
function k_conversion(  # kp version
    data::ARPESData{T,2} where {T};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :phi) && hasdim(data, :eV) "Input data must have 'phi' and 'eV' dimensions."
    energy_definition = metadata(dims(data, :eV))[:energy_definition]
    if !(energy_definition in [BindingEnergy, FinalStateEnergy])
        throw(ArgumentError("Not Implemented for $energy_definition"))
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
    β = metadata(data)[:β]
    β_ = _deg2rad(β, β0)
    ξ_ = _deg2rad(metadata(data)[:ξ], ξ0)
    δ_ = _deg2rad(metadata(data)[:δ], δ0)
    χ_ = haskey(metadata(data), :χ) ? _deg2rad(metadata(data)[:χ], χ0) : NaN

    @debug "Kinetic energy ref to analyzer, and angles" ek α β_ χ_ ξ_ δ_
    # 2. determine k_regions, and use them if kx_range and ky_range are not provided.
    kx_range =
        isnothing(kx_range) ? _kx_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : kx_range
    ky_range =
        isnothing(ky_range) ? _ky_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : ky_range

    @debug "kx_range, ky_range" kx_range ky_range
    # 3. apply interpolation to get the intensity values on the k grid.
    #   3.1 corresponding \alpha and \beta
    kx_grid, ky_grid, ek_grid = prepare_for_broadcast(kx_range, ky_range, ek)

    @debug "size of kx_grid, ky_grid ek_grid, " size(kx_grid) size(ky_grid) size(ek_grid)
    α_range, β_range =
        angle_mapping(analyzer_conf, kx_grid, ky_grid, ek_grid, _deg2rad(β0), χ_, ξ_, δ_)
    @debug "size of α_range, β_range" size(α_range) size(β_range)
    #   3.2 interpolate
    data_k = _interpolate(α_range, ek_grid, α, parent(ek_original), parent(data))
    # 4. construct the output ARPESData object with kx, ky, and eV dimensions.
    return _build_arpesband(data_k, kx_range, eV_range, metadata(data), energy_definition)
end


function k_conversion(  # N=3 version version
    data::ARPESData{T,3} where {T};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :eV) && hasdim(data, :phi) "Input data must have 'phi' and 'eV' dimensions."
    if hasdim(data, :psi)  # kxky conversion from ϕ and ψ
        return kxky_conversion(data; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0)
    end
    if hasdim(data, :hv) || hasdim(data, :hν) # kpz conversion from ϕ and hv
        return kpkz_conversion(data; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0)
    end
    # conversion 
    target_dim = only(otherdims(data, (:eV, :phi)))
    return _slice_and_convert(
        data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
    )
end


function k_conversion(
    data::ARPESData{T,4} where {T};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :eV) && hasdim(data, :phi) "Input data must have 'phi' and 'eV' dimensions."
    if hasdim(data, :psi) && (hasdim(data, :hv) || hasdim(data, :hν))
        return kxkykz_conversion(
            data;
            kx_range,
            ky_range,
            kz_range,
            eV_range,
            β0,
            ξ0,
            δ0,
            χ0,
        )
    end
    if hasdim(data, :psi) && !(hasdim(data, :hv) || hasdim(data, :hν))
        target_dim = only(otherdims(data, (:eV, :phi, :psi)))
        return _slice_and_convert(
            data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
        )
    end
    if (hasdim(data, :hv) || hasdim(data, :hν)) && !hasdim(data, :psi)
        target_dim = only(otherdims(data, (:eV, :phi, :hv, :hν)))
        return _slice_and_convert(
            data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
        )
    end
    target_dim = otherdims(data, (:eV, :phi))[end]
    return _slice_and_convert(
        data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
    )
end


function k_conversion(  # N>=5D
    data::ARPESData{T,N} where {T,N};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :eV) && hasdim(data, :phi) "Input data must have 'phi' and 'eV' dimensions."
    # Slice over the last scan dimension (not involved in k-space conversion) and
    # recurse into (N-1)D. Exclude eV, phi, psi, and hv/hν so that those dimensions
    # are handled by the appropriate lower-level dispatch
    # (N=4 routes psi/hv to kxky/kpkz/kxkykz conversions).
    # Using [end] preserves dimension order: inner recursion handles earlier dims first.
    target_dim = otherdims(data, (:eV, :phi, :psi, :hv, :hν))[end]
    return _slice_and_convert(
        data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
    )
end


"""
    kxky_conversion(data::ARPESData{T,3}; kx_range=nothing, ky_range=nothing,
                    kz_range=nothing, eV_range=nothing, β0=0.0, ξ0=0.0, δ0=0.0, χ0=0.0)

Convert 3D `phi × psi × eV` ARPES data to `kx × ky × eV` momentum coordinates.

Called internally by `k_conversion` when the 3D input has both `phi` and `psi` dimensions.
`phi` provides the first in-plane scattering angle (α in the mapping equations) and `psi`
provides the second (β).

The dataset metadata must include `:analyzer_configuration`, `:workfunction`, `:hv`, `:ξ`,
`:δ`, and optionally `:χ`. The `eV` dimension metadata must include `:energy_definition`
set to `BindingEnergy` or `FinalStateEnergy`.

# Keywords
- `kx_range`, `ky_range`, `kz_range`: Optional target momentum axes. If omitted, the code
  estimates suitable ranges from the input geometry.
- `eV_range`: Optional target energy axis. Defaults to the input `eV` lookup.
- `β0`, `ξ0`, `δ0`, `χ0`: Additional angular offsets in degrees applied during mapping.

# Returns
- `ARPESData`: A new 3D dataset on `kx × ky × eV` axes.
"""
function kxky_conversion(
    data::ARPESData{T,3} where {T};
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
        throw(ArgumentError("Not Implemented for $energy_definition"))
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
    χ_ = haskey(metadata(data), :χ) ? _deg2rad(metadata(data)[:χ], χ0) : NaN

    @debug "Kinetic energy ref to analyzer, and angles" ek α β_ χ_ ξ_ δ_
    # 2. determine k_regions, and use them if kx_range and ky_range are not provided.
    kx_range =
        isnothing(kx_range) ? _kx_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : kx_range
    ky_range =
        isnothing(ky_range) ? _ky_range(analyzer_conf, α, β_, ek, χ_, ξ_, δ_) : ky_range

    @debug "kx_range, ky_range" kx_range ky_range
    # 3. apply interpolation to get the intensity values on the k grid.
    #   3.1 corresponding \alpha and \beta
    kx_grid, ky_grid, ek_grid = prepare_for_broadcast(kx_range, ky_range, ek)

    @debug "size of kx_grid, ky_grid ek_grid, " size(kx_grid) size(ky_grid) size(ek_grid)
    α_range, β_range =
        angle_mapping(analyzer_conf, kx_grid, ky_grid, ek_grid, _deg2rad(β0), χ_, ξ_, δ_)
    @debug "size of α_range, β_range" size(α_range) size(β_range)
    #   3.2 interpolate
    data_k =
        _interpolate(α_range, β_range, ek_grid, α, β, parent(ek_original), parent(data))
    # 4. construct the output ARPESData object with kx, ky, and eV dimensions.
    return _build_arpesband(
        data_k,
        kx_range,
        ky_range,
        eV_range,
        metadata(data),
        energy_definition,
    )
end

"""
    kpkz_conversion(data::ARPESData{T,3}; kx_range=nothing, ky_range=nothing,
                    kz_range=nothing, eV_range=nothing, β0=0.0, ξ0=0.0, δ0=0.0, χ0=0.0)

Convert 3D `phi × hv × eV` ARPES data to parallel/perpendicular momentum coordinates.

Called internally by `k_conversion` when the 3D input has `phi` and `hv` (or `hν`)
dimensions but no `psi` dimension.

!!! note "Not yet implemented"
    This function currently raises an error. The kp–kz conversion via photon-energy
    dependence is not yet implemented.
"""
function kpkz_conversion(
    data::ARPESData{T,3} where {T};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :eV) &&
            hasdim(data, :phi) &&
            (hasdim(data, :hv) || hasdim(data, :hν)) "Input data must have 'phi', 'eV', and 'hv' (or 'hν') dimensions."
    error("kpkz conversion for 3D ARPESData is not implemented yet.")
end

"""
    kxkykz_conversion(data::ARPESData{T,4}; kx_range=nothing, ky_range=nothing,
                      kz_range=nothing, eV_range=nothing, β0=0.0, ξ0=0.0, δ0=0.0, χ0=0.0)

Convert 4D `phi × psi × hv × eV` ARPES data to `kx × ky × kz × eV` momentum coordinates.

Called internally by `k_conversion` when the 4D input has `phi`, `psi`, and `hv` (or `hν`)
dimensions.

!!! note "Not yet implemented"
    This function currently raises an error. The full three-dimensional momentum
    conversion is not yet implemented.
"""
function kxkykz_conversion(
    data::ARPESData{T,4} where {T};
    kx_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    ky_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    kz_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    eV_range::Union{AbstractVector{<:Real},Nothing} = nothing,
    β0::Real = 0.0,
    ξ0::Real = 0.0,
    δ0::Real = 0.0,
    χ0::Real = 0.0,
)
    @assert hasdim(data, :eV) &&
            hasdim(data, :phi) &&
            hasdim(data, :psi) &&
            (hasdim(data, :hv) || hasdim(data, :hν)) "Input data must have 'phi', 'eV', 'psi', and 'hv' (or 'hν') dimensions."
    error("kxkykz conversion for 4D ARPESData is not implemented yet.")
end

# --- internal functions

"""
    _slice_and_convert(data, target_dim; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0)

Slice `data` along `target_dim`, apply `k_conversion` to each slice, and concatenate
the results back along the same dimension. Used to reduce N-dimensional data to
(N-1)-dimensional slices for recursive k-space conversion.
"""
function _slice_and_convert(
    data::ARPESData,
    target_dim;
    kx_range,
    ky_range,
    kz_range,
    eV_range,
    β0,
    ξ0,
    δ0,
    χ0,
)
    unit = haskey(metadata(target_dim), :unit) ? metadata(target_dim)[:unit] : nothing
    k_convs = ARPESData[]
    for (val, arpes_slice) in zip(dims(data, target_dim), eachslice(data; dims = target_dim))
        k_converted = k_conversion(
            arpes_slice; kx_range, ky_range, kz_range, eV_range, β0, ξ0, δ0, χ0
        )
        push!(k_convs, add_dim(k_converted, target_dim, val; unit))
    end
    return cat(k_convs...; dims = target_dim)
end

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
