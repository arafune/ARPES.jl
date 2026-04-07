using Test
using ARPES.IO: SPDLoader
using ARPES.IO.Location: load_data

@testset "SPDLoader" begin
    # Test loading a file with the SPD loader
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    arpes_data = load_data(SPDLoader, file_path)
    @test !isempty(arpes_data.data)
end
