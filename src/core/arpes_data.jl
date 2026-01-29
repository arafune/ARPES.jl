using DimensionalData
using DimensionalData.Dimensions: @dim
using Dates

# --- Dimension names used for ARPES ---
@dim Kx "kx (1/Å)"
@dim phi "degrees"
@dim Energy "energy (eV)"
@dim Cycle "cycle"
@dim Ch2 "ch2"
@dim Counts "counts"

struct ARPESData
    intensity::DimensionalData.AbstractDimArray
    meta::Dict{String,Any}
end

export ARPESData
