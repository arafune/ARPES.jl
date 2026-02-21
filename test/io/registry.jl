using Test
include("../../src/io/types.jl")
include("../../src/io/registry.jl")

@testset "select_loader" begin
    # Test with explicit location
    file_path = joinpath(pkgdir(ARPES), "testdata", "arpes_ch_resolved.itx")
    @test select_loader("dummy.itx", "SPD") == SPDLoader
    @test select_loader(file_path, nothing) == SPDLoader
end
#    # Test with OPTIONS.loc
#    OPTIONS[] = (loc = "SPD")
#
#    OPTIONS[] = (loc = nothing)
#    # Test with file extension detection
#    @test select_loader("dummy.itx", nothing) == SPDLoader
#    @test select_loader("dummy.sp2", nothing) == SPDLoader
#    @test select_loader("dummy.unknown", nothing) == GenericLoader
#end

@testset "get_loader_by_name" begin
    @test get_loader_by_name("SPD") == SPDLoader
    @test_throws ErrorException get_loader_by_name("UNKNOWN")
end

#@testset "detect_loader" begin
#    @test detect_loader("file.itx") == SPDLoader
#    @test detect_loader("file.sp2") == SPDLoader
#    @test detect_loader("file.unknown") == GenericLoader
#end

