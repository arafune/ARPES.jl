using DimensionalData
using ARPES: ARPESData, phi, eV, Z, FinalStateEnergy

function test_spd_standard()
    file_path = joinpath(pkgdir(ARPES), "testdata", "spd_standard.itx")
    spd_standard = load(file_path, loc = "SPD")
    return spd_standard
end

function test_arpes_ch_resolved()
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    meta_data = Dict(:beta=>0.0)
    arpes_ch_resolved = load(file_path, loc = "SPD", extra_metadata = meta_data)
    return arpes_ch_resolved
end
