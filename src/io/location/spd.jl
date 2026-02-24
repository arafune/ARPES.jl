using ..Format: read_itx
using DimensionalData
using DimensionalData: DimArray
using ..IO: SPDLoader
using ARPES: ARPESData, phi, eV, ch, ch2
"""
    load_data(::Type{SPDLoader}, fpath:String)

Load and parse an SPD file from the given file path.

# Arguments
- `::Type{SPDLoader}`: Loader type for SPD files.
- `fpath::String`: Path to the SPD file.

# Returns
A standardized DimArray parsed from the file.

# Notes
Supports `.itx` and `.sp2` file extensions.
"""
function load_data(::Type{SPDLoader}, fpath::String)
    # read the SPD file and parse its content
    ext = lowercase(splitext(fpath)[2])
    if ext == ".itx"
        raw_dimarray = read_itx(fpath)
        return spd_to_standard(raw_dimarray)

    elseif ext == ".sp2"
        raw_dimarray = read_sp2(fpath)
        return spd_to_standard(raw_dimarray)
    end
end


"""
Convert raw DimArray from SPD group to the standard ARPESData format, in which the angle notation
follows Rev.Sci.Instrum. 89, 043903 (2018).

"""
function spd_to_standard(raw::DimArray)
    return raw
end
const DIM_ALIAS = Dict(
    phi => [:non_energy_channel, :angle, :theta],
    eV => [:kinetic_energy, :energy, :binding_energy],
)

const DEFAULT_DIM_MAP = Dict(:x => phi, :y => eV, :z => ch, :w=>ch2)

Location.dim_alias(::SPDLoader) = DIM_ALIAS
Location.default_dim_map(::SPDLoader) = DEFAULT_DIM_MAP
