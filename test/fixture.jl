using ARPES: ARPESData, phi, eV, Z, FinalStateEnergy

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


