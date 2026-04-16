"""
    module ARPES

Main module for the ARPES.jl package, providing core types, data structures, and I/O functionality
for Angle-Resolved Photoemission Spectroscopy (ARPES) data analysis.
"""
module ARPES
include("./types.jl")
include("./dict.jl")
include("./utils.jl")
include("./ARPESData.jl")
export merge_consensus
export kx, ky, kz, phi, psi, eV, detector_ch, ch2, delay, spin
export CPS, Counts
export ARPESData
include("./pumpprobe.jl")
include("./stitch.jl")
include("./transform.jl")

include("io/io.jl")
using .IO
export load

include("./dims.jl")
include("./filter.jl")
include("./derivative.jl")
include("./k_conv/k_conversion.jl")
using .KConversion
export k_conversion

include("./trapezoid.jl")  # As trapezoid uses functions in KConversion, it should be included after KConversion
end
