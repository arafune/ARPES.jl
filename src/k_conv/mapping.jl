using ARPES: TypeI, TypeII, TypeIp, TypeIIp
const K_INV = 0.5123167219534328

unorm_sinc(x) = x == 0 ? 1 : sinc(x/π)

function mapping(::Type{TypeI}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    kx = @. k * (
        (+sin(δ_) * sin(β_) + cos(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        cos(δ_) * cos(ξ_) * sin(α)
    )
    ky = @. k * (
        (-cos(δ_) * sin(β_) + sin(δ_) * sin(ξ_) * cos(β_)) * cos(α) -
        sin(δ_) * cos(ξ_) * sin(α)
    )
    return kx, ky
end

function mapping(::Type{TypeII}, Ek, α, β_, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    kx = @. k * (
        (sin(δ_) * sin(ξ_) + cos(δ_) * sin(β_) * cos(ξ_)) * cos(α) -
        (sin(δ_) * cos(ξ_) - cos(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
    ky = @. k * (
        (-cos(δ_) * sin(ξ_) + sin(δ_) * sin(β_) * cos(ξ_)) * cos(α) +
        (cos(δ_) * cos(ξ_) + sin(δ_) * sin(β_) * sin(ξ_)) * sin(α)
    )
    return kx, ky
end


function mapping_inverse(::Type{TypeI}, Ek, kx, ky, β0, χ_, ξ_, δ_)
    k = @. K_INV * sqrt(Ek)
    α = @. arcsin(
        (sin(ξ_) * sqrt(k^2 - kx^2 - ky^2) - cos(ξ_) * (kx * cos(δ_) + ky * sin(δ_))) / k,
    )
    β = @. β0 + arctan(
        (
            sin(δ_) * kx - cos(δ_) * ky
        )/(
            sin(ξ_) * cos(δ_) * kx +
            sin(ξ_) * sin(δ_) * ky +
            cos(ξ_) * sqrt(k^2 - kx^2 - ky-2)
        ),
    )
    return α, β
end


