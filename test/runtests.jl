using ARPES
using Coverage
using Test

@testset "ARPES.jl" begin
    @testset "IO" begin
        @testset "FORMATS" begin
            include("io/formats/itx.jl")
        end
    end
end
