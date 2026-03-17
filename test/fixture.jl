using DimensionalData
using ARPES: ARPESData, phi, eV, Z, FinalStateEnergy

function test_DimArray3D()
    x = range(1.0, stop = 3.0, length = 3)
    y = collect(range(5, stop = 8, length = 4))
    z = 10.0:14.0
    data = 1.0:60.0 |> collect |> d->reshape(d, 3, 4, 5)
    da = DimArray(data, (X = x, Y = y, Z = z))
    return da
end

function test_DimArray2D()
    x = range(1.0, stop = 3.0, length = 3)
    y = collect(range(5, stop = 8, length = 4))
    data = 1.0:12.0 |> collect |> d->reshape(d, 3, 4)
    da = DimArray(data, (X = x, Y = y))
    return da
end


function test_ARPESData()
    x = range(-5, stop = 9, length = 10)
    y = collect(range(6, stop = 8, length = 10))
    z = 1:10
    data = rand(10, 10, 10)
    da = DimArray(
        data,
        (
            Dim{:phi}(x; metadata = Dict(:unit => "deg")),
            Dim{:eV}(y; metadata = Dict(:unit => "eV")),
            Dim{:Z}(z; metadata = Dict(:unit => "arb. units")),
        ),
    )
    da_ARPES = ARPESData(
        da;
        energy_def = FinalStateEnergy,
        additional_metadata = Dict(
            :β=>0.0,
            :ξ=>0.0,
            :δ=>0.0,
            :hv => 0.0,
            :workfunction => 4.4,
        ),
    )
    return da_ARPES
end

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
