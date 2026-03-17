using ARPES
using Test
using Coverage

include("./fixture.jl")

@testset "ARPES.jl" begin
    @testset "function for Dict as metadata" begin
        include("dict.jl")
    end
    @testset "Core" begin
        include("ARPESData.jl")
        include("dims.jl")
        include("data.jl")
    end
    @testset "k_conv" begin
        include("k_conv/preprocess.jl")
        include("k_conv/interpolation.jl")
        include("k_conv/k_conversion.jl")
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
