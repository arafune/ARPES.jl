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

include("./transform.jl")

include("io/io.jl")
using .IO
export load

include("./dims.jl")
include("./filter.jl")
include("./differential.jl")
include("./k_conv/k_conversion.jl")
using .KConversion
export k_conversion
end
