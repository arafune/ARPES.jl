using Test
using ARPES
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using DimensionalData: Dim, hasdim, lookup

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end

@testset "Test for vstack_arpes" begin
  data1 = test_ARPESData()
  data2 = test_ARPESData()
  data3 = test_ARPESData()
  data4 = test_spd_standard()

  data_1_1 = add_dim(data1, :w, 0.0)
  data_2_1 = add_dim(data1, :w, 1.0)
  data_3_1 = add_dim(data1, :w, 3.0)
  stack_1_to_3 = vstack_arpes([data_2_1, data_1_1, data_3_1], :w)

  @test size(stack_1_to_3) == (10, 10, 10, 3)
  @test hasdim(stack_1_to_3, :w)
  @tets lookup(stack_1_to_3, :w) == [0.0, 1.0, 3.0]

  stack_1_to_3_alt = vstack_arpes([data_2_1, data_1_1, data_3_1], :w; sort=false)
  @tets lookup(stack_1_to_3_alt, :w) == [1.0, 0.0, 3.0]
end
