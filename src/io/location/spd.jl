using ..Format: read_itx
using DimensionalData
using DimensionalData: DimArray
using ..IO: SPDLoader
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI
using ..IO.Location: to_standardize

const DIM_ALIAS = Dict(
    phi => [:non_energy_channel, :angle, :theta],
    eV => [:kinetic_energy, :energy, :binding_energy],
    detector_ch => [:energy_channel, :channel, :ch],
)

const DEFAULT_DIM_MAP = Dict(:x => phi, :y => eV, :z => detector_ch, :w=>ch2)

Location.dim_alias(::Type{SPDLoader}) = DIM_ALIAS
Location.default_dim_map(::Type{SPDLoader}) = DEFAULT_DIM_MAP
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
    standard_array = to_standardize(SPDLoader, raw)
    if haskey(metadata(standard_array), :d_scale)  #made from itx
        if startswith(metadata(standard_array)[:d_scale][:unit], "count")
            return ARPESData(standard_array, Counts(), TypeI)
        end
        return ARPESData(standard_array, CPS(), TypeI)
    end
    return ARPESData(standard_array, Counts(), TypeI)
end
