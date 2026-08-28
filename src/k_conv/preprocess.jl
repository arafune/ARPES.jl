using Statistics
using DimensionalData: Dimension
using ..ARPES: EnergyDefinition
using ..ARPES: BindingEnergy, FinalStateEnergy, KineticEnergy, IntermediateEnergy
using ..ARPES: ARPESData, kx, ky, kz, phi, psi, eV

"""
     prepare_for_broadcast(arrs...)

Reshapes each input array in `arrs` for N-dimensional broadcasting.

For each array argument, this function reshapes it so that its length occupies a unique dimension, with all other dimensions set to 1. This enables element-wise broadcasting across all input arrays. Non-array arguments are returned unchanged.

# Arguments
- `arrs...`: A variable number of arguments, typically arrays, but may include non-array values.

# Returns
- A tuple where each array is reshaped for broadcasting, and non-array arguments are unchanged.

# Example
```julia
a = [1, 2, 3]
b = [4, 5]
prepare_for_broadcast(a, b)
# returns (3×1 Array, 1×2 Array)
```
"""
function prepare_for_broadcast(arrs::Union{AbstractArray,Real}...)
    arrays_count = count(a -> a isa AbstractArray, arrs)
    array_indices = cumsum([a isa AbstractArray for a in arrs])
    N = length(arrs)

    ntuple(
        idx -> begin
            a = arrs[idx]
            if a isa AbstractArray
                curr_idx = array_indices[idx]
                v = vec(a)
                shape = ntuple(j -> j == curr_idx ? length(v) : 1, arrays_count)
                reshape(v, shape)
            else
                a
            end
        end,
        N,
    )
end


"""
    _ek_range(energy_def::EnergyDefinition, ek, workfunction::Real, hv::Real)

Return a shifted kinetic energy axis depending on the energy definition.

- If `energy_def` is `FinalStateEnergy`, shift by `-workfunction`.
- If `energy_def` is `BindingEnergy`, shift by `hv - workfunction`.

# Arguments
- `energy_def`: The energy definition type (`FinalStateEnergy` or `BindingEnergy`).
- `ek`: The kinetic energy axis (can be `Dimension`, `AbstractRange`, or `AbstractVector`).
- `workfunction`: The work function value.
- `hv`: The photon energy.

# Returns
A shifted kinetic energy axis of AbstractArray.
"""
function _ek_range(
    energy_def::EnergyDefinition,
    ek::Dimension,
    workfunction::Real,
    hv::Real,
)
    if energy_def == FinalStateEnergy
        return parent(shift_dim(ek, -workfunction))
    elseif energy_def == BindingEnergy
        return parent(shift_dim(ek, hv - workfunction))
    end
end

function _ek_range(
    energy_def::EnergyDefinition,
    ek::AbstractRange,
    workfunction::Real,
    hv::Real,
)
    if energy_def == FinalStateEnergy
        return range(first(ek) - workfunction; step = step(ek), length = length(ek))
    elseif energy_def == BindingEnergy
        return range(first(ek) + hv - workfunction; step = step(ek), length = length(ek))
    end
end

function _ek_range(
    energy_def::EnergyDefinition,
    ek::AbstractVector,
    workfunction::Real,
    hv::Real,
)
    if energy_def == FinalStateEnergy
        return @. ek - workfunction
    elseif energy_def == BindingEnergy
        return @. ek + hv - workfunction
    end
end

"""
    _kx_range(
        analyzer_conf::Type{<:AbstractAnalyzerConfiguration},
        ek::AbstractVector,
        α::AbstractVector,
        β_::Union{AbstractVector{<:Real},Real},
        χ_::Real,
        ξ_::Real,
        δ_::Real,
    )


Determine the kx range for the given analyzer configuration and parameters.

# Arguments
- `analyzer_conf::Type{<:AbstractAnalyzerConfiguration}`: The analyzer configuration type.
- `ek::AbstractVector`: Range/Vector of kinetic energies.
- `α::AbstractVector`:  Range/Vector of alpha angles.
- `β_::AbstractVector`: Range/Vector of beta angles.
- `χ_`: Chi parameter (type depends on context).
- `ξ_`: Xi parameter (type depends on context).
- `δ_`: Delta parameter (type depends on context).

# Returns
- `kx_range`: Range for kx, determined based on the mapping from the provided parameters.
            if the step size cannot be determined (e.g., due to insufficient data
            in most case, :α is a single value),
            returns the minimum kx value.

# Description
1. Determines the step size for kx using the minimum kinetic energy.
2. Computes the minimum and maximum values for kx over the full parameter space.
3. Constructs and returns the kx range using the calculated step size.
"""
function _kx_range(
    analyzer_conf::Type{<:AbstractAnalyzerConfiguration},
    α::AbstractVector{<:Real},
    β_::Union{AbstractVector{<:Real},Real},
    ek::AbstractVector{<:Real},
    χ_::Real,
    ξ_::Real,
    δ_::Real,
)
    ek_min = minimum(ek)

    kx_ = mapped_kx(analyzer_conf, prepare_for_broadcast(α, β_, ek_min)..., χ_, ξ_, δ_)

    d_kx = abs.(diff(kx_, dims = 1))
    step_kx = isempty(d_kx) ? nothing : minimum(d_kx)
    @debug "min_kx, max_kx, step_kx @ minimum ek" minimum(kx_) maximum(kx_) step_kx ek_min
    kx_ = mapped_kx(analyzer_conf, prepare_for_broadcast(α, β_, ek)..., χ_, ξ_, δ_)
    min_kx, max_kx = minimum(kx_), maximum(kx_)
    @debug "kx_ shape" size(kx_)
    @debug "min_kx, max_kx, step_kx" min_kx max_kx step_kx
    if step_kx === nothing || step_kx == 0.0
        return min_kx
    end
    kx_range = range(start = min_kx, stop = max_kx, step = step_kx)
    @debug "kx_range length" length(kx_range)
    return kx_range
end

"""
    _ky_range(
        analyzer_conf::Type{<:AbstractAnalyzerConfiguration},
        ek::AbstractVector,
        α::AbstractVector,
        β_::Union{AbstractVector{<:Real},Real},
        χ_::Real,
        ξ_::Real,
        δ_::Real,
    )

Determine the ky range for the given analyzer configuration and parameters.

# Arguments
- `analyzer_conf::Type{<:AbstractAnalyzerConfiguration}`: The analyzer configuration type.
- `α::AbstractVector`: Range/Vector of alpha angles.
- `β_::AbstractVector`: Range/Vector of beta angles.
- `ek::AbstractVector`: Range/Vector of kinetic energies.
- `χ_`: Chi parameter (type depends on context).
- `ξ_`: Xi parameter (type depends on context).
- `δ_`: Delta parameter (type depends on context).

# Returns
- `ky_range`:  Range for ky, determined based on the mapping from the provided parameters.
             if the step size cannot be determined (e.g., due to insufficient data.
             in most case, :β is a single value), returns the minimum ky value.

# Description
1. Determines the step size for ky using the minimum kinetic energy.
2. Computes the minimum and maximum values for ky over the full parameter space.
3. Constructs and returns the ky range using the calculated step size.
"""
function _ky_range(
    analyzer_conf::Type{<:AbstractAnalyzerConfiguration},
    α::AbstractVector{<:Real},
    β_::AbstractVector{<:Real},
    ek::AbstractVector{<:Real},
    χ_::Real,
    ξ_::Real,
    δ_::Real,
)
    ek_min = minimum(ek)
    ky_ = mapped_ky(analyzer_conf, prepare_for_broadcast(α, β_, ek_min)..., χ_, ξ_, δ_)
    d_ky = abs.(diff(ky_, dims = β_ isa Real ? 1 : 2))
    step_ky = isempty(d_ky) ? nothing : minimum(d_ky)
    ky_ = mapped_ky(analyzer_conf, prepare_for_broadcast(α, β_, ek)..., χ_, ξ_, δ_)
    min_ky, max_ky = minimum(ky_), maximum(ky_)
    if step_ky === nothing || step_ky == 0.0
        return min_ky
    end
    return range(start = min_ky, stop = max_ky, step = step_ky)
end

function _ky_range(
    analyzer_conf::Type{<:AbstractAnalyzerConfiguration},
    α::AbstractVector{<:Real},
    β_::Real,
    ek::AbstractVector{<:Real},
    χ_::Real,
    ξ_::Real,
    δ_::Real,
)
    ek_min = minimum(ek)
    ky_ = mapped_ky(analyzer_conf, prepare_for_broadcast(α, β_, ek_min)..., χ_, ξ_, δ_)
    return mean(filter(!isnan, ky_))
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
    #
    for dim in [:phi, :eV]
        if !hasdim(data, dim)
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

    for meta in _set_supplemental_angles(metadata(data)[:analyzer_configuration])
        if !haskey(metadata(data), meta)
            metadata(data)[meta] = 0.0
            @warn "Missing metadata: $meta. Using default value of 0.0 for conversion."
        end
    end
    return true
end

_set_supplemental_angles(::Type{<:AbstractAnalyzerConfigurationWithoutDeflector}) = (:ξ, :δ)
_set_supplemental_angles(::Type{<:AbstractAnalyzerConfigurationWithDeflector}) =
    (:ξ, :δ, :χ)
