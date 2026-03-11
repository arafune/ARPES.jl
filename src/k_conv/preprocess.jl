
"""
    _kx_range(
        analyzer_conf::Type{<:AnalyzerConfiguration},
        ek::AbstractVector,
        α::AbstractVector,
        β_::AbstractVector,
        χ_,
        ξ_,
        δ_,
    )

Determine the kx range for the given analyzer configuration and parameters.

# Arguments
- `analyzer_conf::Type{<:AnalyzerConfiguration}`: The analyzer configuration type.
- `ek::AbstractVector`: Range/Vector of kinetic energies.
- `α::AbstractVector`:  Range/Vector of alpha angles.
- `β_::AbstractVector`: Range/Vector of beta angles.
- `χ_`: Chi parameter (type depends on context).
- `ξ_`: Xi parameter (type depends on context).
- `δ_`: Delta parameter (type depends on context).

# Returns
- `kx_range`:  Range for kx, determined based on the mapping from the provided parameters.

# Description
1. Determines the step size for kx using the minimum kinetic energy.
2. Computes the minimum and maximum values for kx over the full parameter space.
3. Constructs and returns the kx range using the calculated step size.
"""
function _kx_range(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    ek::AbstractVector{<:Real},
    α::AbstractVector{<:Real},
    β_::AbstractVector{<:Real},
    χ_,
    ξ_,
    δ_,
)
    ek_min = minimum(ek)
    kx_ = mapped_kx(analyzer_conf, ek_min, reshape(α, :, 1), reshape(β_, 1, :), χ_, ξ_, δ_)
    d_kx = abs.(diff(kx_, dims = 1))
    step_kx = isempty(d_kx) ? nothing : minimum(d_kx)
    @debug "min_kx, max_kx, step_kx @ minimum ek" minimum(kx_) maximum(kx_) step_kx ek_min
    kx_ = mapped_kx(
        analyzer_conf,
        reshape(ek, :, 1, 1),
        reshape(α, 1, :, 1),
        reshape(β_, 1, 1, :),
        χ_,
        ξ_,
        δ_,
    )
    min_kx, max_kx = minimum(kx_), maximum(kx_)
    @debug "kx_ shape" size(kx_)
    @debug "min_kx, max_kx, step_kx" min_kx max_kx step_kx
    if step_kx === nothing || step_kx == 0.0
        return min_kx:min_kx
    end
    return range(start = min_kx, stop = max_kx, step = step_kx)

end

"""
    _ky_range(
        analyzer_conf::Type{<:AnalyzerConfiguration},
        ek::AbstractVector,
        α::AbstractVector,
        β_::AbstractVector,
        χ_,
        ξ_,
        δ_,
    )

Determine the ky range for the given analyzer configuration and parameters.

# Arguments
- `analyzer_conf::Type{<:AnalyzerConfiguration}`: The analyzer configuration type.
- `ek::AbstractVector`: Range/Vector of kinetic energies.
- `α::AbstractVector`: Range/Vector of alpha angles.
- `β_::AbstractVector`: Range/Vector of beta angles.
- `χ_`: Chi parameter (type depends on context).
- `ξ_`: Xi parameter (type depends on context).
- `δ_`: Delta parameter (type depends on context).

# Returns
- `ky_range`:  Range for ky, determined based on the mapping from the provided parameters.

# Description
1. Determines the step size for ky using the minimum kinetic energy.
2. Computes the minimum and maximum values for ky over the full parameter space.
3. Constructs and returns the ky range using the calculated step size.
"""
function _ky_range(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    ek::AbstractVector{<:Real},
    α::AbstractVector{<:Real},
    β_::AbstractVector{<:Real},
    χ_,
    ξ_,
    δ_,
)
    ek_min = minimum(ek)
    ky_ = mapped_ky(analyzer_conf, ek_min, reshape(α, :, 1), reshape(β_, 1, :), χ_, ξ_, δ_)
    d_ky = diff(ky_, dims = 2)
    step_ky = isempty(d_ky) ? nothing : minimum(d_ky)
    ky_ = mapped_ky(
        analyzer_conf,
        reshape(ek, :, 1, 1),
        reshape(α, 1, :, 1),
        reshape(β_, 1, 1, :),
        χ_,
        ξ_,
        δ_,
    )
    min_ky, max_ky = minimum(ky_), maximum(ky_)
    if step_ky === nothing || step_ky == 0.0
        return min_ky:min_ky
    end
    return range(start = min_ky, stop = max_ky, step = step_ky)

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
    for dim in [phi, eV]
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


