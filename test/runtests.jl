using ARPES
using Coverage
using Test

@testset "ARPES.jl" begin
    @testset "IO" begin
        @testset "Location" begin
            include("io/location/spd.jl")
        end
        @testset "FORMATS" begin
            include("io/formats/itx.jl")
        end
        @testset "registry" begin
            include("io/registry.jl")
        end
    end
end
