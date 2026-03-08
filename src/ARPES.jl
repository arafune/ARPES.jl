"""
    module ARPES

Main module for the ARPES.jl package, providing core types, data structures, and I/O functionality
for Angle-Resolved Photoemission Spectroscopy (ARPES) data analysis.
"""
module ARPES
include("./ARPESData.jl")
export kx, ky, kz, phi, psi, eV, detector_ch, ch2, delay, spin
export CPS, Counts
export ARPESData

include("io/io.jl")
using .IO
export load

include("./dims.jl")
end
