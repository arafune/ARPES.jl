using Test
using ARPES.IO: SPDLoader
using ARPES: phi, eV, detector_ch, ch2
using ARPES.IO.Location: canonical_dim
# canonical_dimがLocationモジュールにある場合はusing Location: canonical_dimも追加

@testset "canonical_dim with SPDLoader" begin
    loader = SPDLoader()

    @test canonical_dim(loader, :non_energy_channel) == phi
    @test canonical_dim(loader, :angle) == phi
    @test canonical_dim(loader, :theta) == phi
    @test canonical_dim(loader, :kinetic_energy) == eV
    @test canonical_dim(loader, :energy) == eV
    @test canonical_dim(loader, :binding_energy) == eV

    @test canonical_dim(loader, :x) == phi
    @test canonical_dim(loader, :y) == eV
    @test canonical_dim(loader, :z) == detector_ch
    @test canonical_dim(loader, :w) == ch2

    @test canonical_dim(loader, :unknown) === nothing
end

