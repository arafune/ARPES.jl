using Test
using ARPES.IO: select_loader, get_loader_by_name

@testset "select_loader" begin
    # Test with explicit location
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    @test select_loader("dummy.itx", "SPD") == SPDLoader
    @test select_loader(file_path, nothing) == SPDLoader
end

@testset "get_loader_by_name" begin
    @test get_loader_by_name("SPD") == SPDLoader
    @test_throws ErrorException get_loader_by_name("UNKNOWN")
end

