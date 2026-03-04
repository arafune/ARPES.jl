using ..Format: read_itx
using DimensionalData
using DimensionalData: DimArray
using ..IO: SPDLoader
using ARPES: ARPESData, phi, eV, detector_ch, ch2, CPS, Counts, TypeI, FinalStateEnergy
using ..IO.Location: to_standardize

"""
    negate(x)

Returns the negation of the input `x`.

# Arguments
- `x`: A numeric value.

# Returns
- The negated value of `x`.
"""
negate(x) = -x

"""
    pulse_to_theta(pulses::Integer) -> Float64

Converts a given number of pulses from a rotary encoder to the corresponding emission angle (theta) in degrees.
Assumes:
- `normal_emission_angle` is the reference angle (315.0 degrees).
- Each degree corresponds to 6000 pulses.
- One full rotation is 2160000 pulses.

# Arguments
- `pulses::Integer`: The number of pulses to convert.

# Returns
- `Float64`: The calculated emission angle in degrees.
"""
function _pulse_to_theta_spd(pulses::Integer)
    NORMAL_EMISSION_ANGLE = 315.0
    PULSES_ONE_DEGREE = 6000
    PULSE_FULL_ROTATION = 2160000
    angle = mod(pulses, PULSE_FULL_ROTATION)÷PULSES_ONE_DEGREE
    return mod(-(angle - NORMAL_EMISSION_ANGLE) + 90, 180) - 90
end

const DIM_ALIAS = Dict(
    phi => [:non_energy_channel, :angle, :theta, :theta_y],
    eV => [:kinetic_energy, :energy, :binding_energy],
    detector_ch => [:energy_channel, :channel, :ch],
)

const DEFAULT_DIM_MAP = Dict(:x => phi, :y => eV, :z => detector_ch, :w=>ch2)

const STANDARD_ANGLES = Dict(
    :β => [Dict(:theta => _pulse_to_theta_spd, :a => _pulse_to_theta_spd)],
    :ξ => Dict(:beta => negate),
    :δ => 0.0,
)

Location.dim_alias(::Type{SPDLoader}) = DIM_ALIAS
Location.default_dim_map(::Type{SPDLoader}) = DEFAULT_DIM_MAP
Location.angle_Shin_convention(::Type{SPDLoader}) = STANDARD_ANGLES

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
function load_data(
    ::Type{SPDLoader},
    fpath::String;
    extra_metadata::Union{AbstractDict{Symbol,<:Any},Nothing} = nothing,
)
    # read the SPD file and parse its content
    ext = lowercase(splitext(fpath)[2])
    if ext == ".itx"
        raw_dimarray = read_itx(fpath)
    elseif ext == ".sp2"
        raw_dimarray = read_sp2(fpath)
    end
    if extra_metadata !== nothing
        new_metadata = merge(Dict(metadata(raw_dimarray)), extra_metadata)
        raw_dimarray = DimensionalData.rebuild(raw_dimarray; metadata = new_metadata)
    end
    return _spd_to_standard(raw_dimarray)
end


"""
    _spd_to_standard(raw::DimArray) -> ARPESData

Convert a raw `DimArray` from the SPD group into the standard `ARPESData` format,
in which the angle nomenclature follows Rev. Sci. Instrum. **89**, 043903 (2018).

- Standardizes the input array using `to_standardize(SPDLoader, raw)`.
- Determines the appropriate intensity unit (`Counts` or `CPS`) based on the `:d_scale` metadata:
    - If the unit starts with "count", returns `ARPESData` with `Counts()`.
    - Otherwise, returns `ARPESData` with `CPS()`.
- If no `:d_scale` metadata is present, defaults to `Counts()`.

# Arguments
- `raw::DimArray`: The raw SPD data array to be standardized.

# Returns
- `ARPESData`: The standardized ARPES data object with appropriate intensity units.
"""
function _spd_to_standard(raw::DimArray)::ARPESData
    standard_array = to_standardize(SPDLoader, raw)
    # to convert k-space smoothly,
    # * rename beta -> \xi (and flip the sign?)

    if haskey(metadata(standard_array), :d_scale)  #made from itx
        if startswith(metadata(standard_array)[:d_scale][:unit], "count")
            return ARPESData(standard_array, Counts(), TypeI, FinalStateEnergy)
        end
        return ARPESData(standard_array, CPS(), TypeI, FinalStateEnergy)
    end
    return ARPESData(standard_array, Counts(), TypeI, FinalStateEnergy)
end


