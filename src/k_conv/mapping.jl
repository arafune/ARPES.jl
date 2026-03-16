using ARPES: TypeI, TypeII
const K_INV = 0.5123167219534328
export momentum_mapping, mapped_kx, mapped_ky, angle_mapping, mapped_α, mapped_β
using ..ARPES: AnalyzerConfiguration

unorm_sinc(x) = x == 0 ? 1 : sinc(x/π)

function momentum_mapping(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    α,
    β_,
    Ek,
    χ_,
    ξ_,
    δ_,
)
    return mapped_kx(analyzer_conf, α, β_, Ek, χ_, ξ_, δ_),
    mapped_ky(analyzer_conf, α, β_, Ek, χ_, ξ_, δ_)
end

function mapped_kx(::Type{TypeI}, α, β_, Ek, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (sin(δ_) * sin(β_) + cos(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        cos(δ_) * cos(ξ_) * sin(α)
    )
end

function mapped_ky(::Type{TypeI}, α, β_, Ek, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (-cos(δ_) * sin(β_) + sin(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        sin(δ_) * cos(ξ_) * sin(α)
    )
end

function mapped_kx(::Type{TypeII}, α, β_, Ek, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (sin(δ_) * sin(ξ_) + cos(δ_) * sin(β_) * cos(ξ_)) * cos(α) -
        (sin(δ_) * cos(ξ_) - cos(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
end

function mapped_ky(::Type{TypeII}, α, β_, Ek, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (-cos(δ_) * sin(ξ_) + sin(δ_) * sin(β_) * cos(ξ_)) * cos(α) +
        (cos(δ_) * cos(ξ_) + sin(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
end

function angle_mapping(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    kx,
    ky,
    Ek,
    β0,
    χ_,
    ξ_,
    δ_,
)
    return mapped_α(analyzer_conf, kx, ky, Ek, β0, χ_, ξ_, δ_),
    mapped_β(analyzer_conf, kx, ky, Ek, β0, χ_, ξ_, δ_)
end

function mapped_α(::Type{TypeI}, kx, ky, Ek, β0, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. asin(
        (sin(ξ_) * sqrt(k^2 - kx^2 - ky^2) - cos(ξ_) * (kx * cos(δ_) + ky * sin(δ_))) / k,
    )
end

function mapped_β(::Type{TypeI}, kx, ky, Ek, β0, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. β0 + atan(
        (
            sin(δ_) * kx - cos(δ_) * ky
        )/(
            sin(ξ_) * cos(δ_) * kx +
            sin(ξ_) * sin(δ_) * ky +
            cos(ξ_) * sqrt(k^2 - kx^2 - ky^2)
        ),
    )
end
