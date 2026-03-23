using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData
using DimensionalData: Dim, hasdim, lookup
using Random
using Statistics

@testset "Test for differential 1st order" begin
  t_rnd = _random_monotonic_vector(200)
  Z_rnd = DimArray(sin.(t_rnd), (Dim{:t}(t_rnd)))
  Z = DimArray(test_data.sine.clean, (Dim{:t}(test_data.t)))
  #
  Z_rnd_diff_1 = differential(Z_rnd, :t)
  Z_diff_1 = differential(Z, :t)
  cos_rnd = cos.(t_rnd)
  cos_clean = cos.(test_data.t)
  @test maximum(parent(Z_rnd_diff_1) - cos_rnd) <  0.005
  @test maximum(parent(Z_diff_1) - cos_clean) <  0.005
end


