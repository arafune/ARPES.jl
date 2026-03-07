using ARPES
using Coverage
using Test

@testset "ARPES.jl" begin
    @testset "Core" begin
        include("ARPESData.jl")
    end
    @testset "IO" begin
        @testset "IO basic" begin
            include("io/io.jl")
        end
        @testset "Location" begin
            include("io/location/location.jl")
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
