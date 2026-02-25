using Test
using ARPES

@testset "load" begin
    # Test loading a file with the SPD loader
    file_path = joinpath(pkgdir(ARPES), "testdata", "spd_standard.itx")
    data = load(file_path, loc = "SPD")
    @test !isempty(data.intensity)
end
