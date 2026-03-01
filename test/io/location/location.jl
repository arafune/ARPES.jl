using Test
using ARPES.IO: SPDLoader
using ARPES: phi, eV, detector_ch, ch2
using ARPES.IO.Location: canonical_dim, pulse_to_theta, negate, convert_Shin_convention

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

@testset "test for helper function in location.jl" begin
    @test pulse_to_theta(0) == -45.0
    @test pulse_to_theta(-6000) == -44.0
    @test pulse_to_theta(-12000) == -43.0
    @test pulse_to_theta(-270000) == 0
    @test pulse_to_theta(-540000) == 45
    @test negate(5) == -5
    @test negate(5.3) == -5.3
    @test negate(-3) == 3
end

# ------ test for convert_Shin_convention ------ 
@testset "test for convert_Shin_convention" begin
    metadata_dict = Dict(:angle => 300, :energy => 10.0, :other => "value", :a => 200)
    conversion_rule = Dict(
        :angle => pulse_to_theta,
        :energy => x -> x * 2,
        :χ => 0.0,
        :β => [Dict(:theta => pulse_to_theta, :a => pulse_to_theta)],
    )

    converted_dict = convert_Shin_convention(metadata_dict, conversion_rule)

    @test converted_dict[:angle] == pulse_to_theta(300)
    @test converted_dict[:energy] == 20.0
    @test converted_dict[:β] == pulse_to_theta(200)
    @test converted_dict[:χ] == 0.0
end

