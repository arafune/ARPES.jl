using ARPES: TypeI, TypeII, TypeIp, TypeIIp
const K_INV = 0.5123167219534328
export momentum_mapping, mapped_kx, mapped_kz, angle_mapping, mapped_α, mapped_β

unorm_sinc(x) = x == 0 ? 1 : sinc(x/π)

function momentum_mapping(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    Ek,
    α,
    β_,
    χ_,
    ξ_,
    δ_,
)
    return mapped_kx(analyzer_conf, Ek, α, β_, χ_, ξ_, δ_),
    mapped_ky(analyzer_conf, Ek, α, β_, χ_, ξ_, δ_)
end

function mapped_kx(::Type{TypeIp}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (sin(δ_) * sin(β_) + cos(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        cos(δ_) * cos(ξ_) * sin(α)
    )
end

function mapped_ky(::Type{TypeIp}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (-cos(δ_) * sin(β_) + sin(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        sin(δ_) * cos(ξ_) * sin(α)
    )
end

function mapped_kx(::Type{TypeIIp}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (sin(δ_) * sin(ξ_) + cos(δ_) * sin(β_) * cos(ξ_)) * cos(α) -
        (sin(δ_) * cos(ξ_) - cos(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
end

function mapped_ky(::Type{TypeIIp}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. k * (
        (-cos(δ_) * sin(ξ_) + sin(δ_) * sin(β_) * cos(ξ_)) * cos(α) +
        (cos(δ_) * cos(ξ_) + sin(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
end

function angle_mapping(
    analyzer_conf::Type{<:AnalyzerConfiguration},
    Ek,
    kx,
    ky,
    β0,
    χ_,
    ξ_,
    δ_,
)
    return mapped_α(analyzer_conf, Ek, kx, ky, β0, χ_, ξ_, δ_),
    mapped_β(analyzer_conf, Ek, kx, ky, β0, χ_, ξ_, δ_)
end

function mapped_α(::Type{TypeI}, Ek, kx, ky, β0, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. arcsin(
        (sin(ξ_) * sqrt(k^2 - kx^2 - ky^2) - cos(ξ_) * (kx * cos(δ_) + ky * sin(δ_))) / k,
    )
end

function mapped_β(::Type{TypeI}, Ek, kx, ky, β0, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    return @. β0 + arctan(
        (
            sin(δ_) * kx - cos(δ_) * ky
        )/(
            sin(ξ_) * cos(δ_) * kx +
            sin(ξ_) * sin(δ_) * ky +
            cos(ξ_) * sqrt(k^2 - kx^2 - ky^2)
        ),
    )
end
