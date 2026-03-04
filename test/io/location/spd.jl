using Test
using ARPES.IO: SPDLoader
using ARPES.IO.Location: load_data

@testset "SPDLoader" begin
    # Test loading a file with the SPD loader
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    data = load_data(SPDLoader, file_path)
    @test !isempty(data.intensity)
end
