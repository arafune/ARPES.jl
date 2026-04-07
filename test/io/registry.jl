using Test
using ARPES.IO: select_loader, get_loader_by_name

@testset "get_loader_by_name" begin
    @test get_loader_by_name("SPD") == SPDLoader
    @test_throws ErrorException get_loader_by_name("UNKNOWN")
end

