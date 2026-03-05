using Test
using ARPES.IO: SPDLoader
using ARPES: phi, eV, detector_ch, ch2
using ARPES.IO.Location: canonical_dim, _pulse_to_theta_spd, negate, standardize_metadata

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
    @test _pulse_to_theta_spd(0) == -45.0
    @test _pulse_to_theta_spd(-6000) == -44.0
    @test _pulse_to_theta_spd(-270000) == 0
    @test _pulse_to_theta_spd(-540000) == 45
    @test negate(5.3) == -5.3
    @test negate(-3) == 3
end

# ------ test for standardize_metadata ------ 
@testset "test for standardize_metadata" begin
    metadata_dict = Dict(:angle => 300, :energy => 10.0, :other => "value", :a => 200)
    conversion_rule = Dict(
        :angle => _pulse_to_theta_spd,
        :energy => x -> x * 2,
        :χ => 0.0,
        :β => [Dict(:theta => _pulse_to_theta_spd, :a => _pulse_to_theta_spd)],
    )

    converted_dict = standardize_metadata(metadata_dict, conversion_rule)

    @test converted_dict[:angle] == _pulse_to_theta_spd(300)
    @test converted_dict[:energy] == 20.0
    @test converted_dict[:β] == _pulse_to_theta_spd(200)
    @test converted_dict[:χ] == 0.0

    blank_conversion_rule = Dict()
    converted_dict = standardize_metadata(metadata_dict, blank_conversion_rule)
    @test isempty(converted_dict)

    blank_vector_conversion_rule = Dict(:beta=>[])
    converted_dict = standardize_metadata(metadata_dict, blank_vector_conversion_rule)
    @test isempty(converted_dict)

    conversion_rule_just_replace = Dict(:another_angle => :angle)
    converted_dict = standardize_metadata(metadata_dict, conversion_rule_just_replace)
    @test haskey(converted_dict, :another_angle)

    donversion_rule_with_nonexistent_key = Dict(:nonexistent => :nonexistent)
end

@testset "standardize_metadata with Dict rule" begin
    metadata = Dict(:a => 10, :b => 20)

    conversion_rule1 = Dict(:new1 => Dict(:a => x -> x * 2, :b => x -> x + 1))
    result1 = standardize_metadata(metadata, conversion_rule1)
    @test result1[:new1] == 20
    @test haskey(result1, :new1)

    conversion_rule2 = Dict(:new2 => Dict(:c => x -> x * 100, :b => x -> x + 5))
    result2 = standardize_metadata(metadata, conversion_rule2)
    @test result2[:new2] == 25

    conversion_rule3 = Dict(:new3 => Dict(:x => x -> x * 2, :y => x -> x + 1))
    result3 = standardize_metadata(metadata, conversion_rule3)
    @test !haskey(result3, :new3)
    @test isempty(result3)
end


@testset "standardize_metadata with Vector rule" begin
    metadata = Dict(:a => 10, :b => 20)
    conversion_rule1 = Dict(:new1 => [:a, :b])
    result1 = standardize_metadata(metadata, conversion_rule1)
    @test result1[:new1] == 10

    conversion_rule2 = Dict(:new2 => [:c, :b])
    result2 = standardize_metadata(metadata, conversion_rule2)
    @test result2[:new2] == 20

    conversion_rule3 = Dict(:new3 => [:x, :y])
    result3 = standardize_metadata(metadata, conversion_rule3)
    @test !haskey(result3, :new3)
    @test isempty(result3)
end
