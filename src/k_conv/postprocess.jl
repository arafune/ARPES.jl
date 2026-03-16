using DimensionalData: DimArray, metadata, Dim
function _build_arpesband(
    banddata::AbstractArray,
    kx_range::AbstractVector,
    ky_range::AbstractVector,
    ek_range::AbstractVector,
    metadata_original::Dict,
    energy_def::EnergyDefinition,
)

    dimensions = Dim{:kx}(kx_range, metadata = Dict(:unit => "1/Å")),
    Dim{:ky}(ky_range, metadata = Dict(:unit => "1/Å")),
    Dim{:eV}(ek_range, metadata = Dict(:energy_definition => energy_def))
    da = DimArray(
        banddata,
        dimensions,
        name = metadata_original[:name],
        metadata = metadata_original,
    )
    arpes_data = ARPESData(
        da;
        intensity_unit = metadata_original[:intensity_unit],
        analyzer_config = metadata_original[:analyzer_configuration],
        energy_def = energy_def,
        additional_metadata = Dict(),
    )
    return arpes_data
end

function _build_arpesband(
    banddata::AbstractArray,
    kx_range::AbstractVector,
    _,
    ek_range::AbstractVector,
    metadata_original::Dict,
    energy_def::EnergyDefinition,
)
    dimensions = Dim{:kx}(kx_range, metadata = Dict(:unit => "1/Å")),
    Dim{:eV}(ek_range, metadata = Dict(:energy_definition => energy_def))
    da = DimArray(
        banddata,
        dimensions,
        name = metadata_original[:name],
        metadata = metadata_original,
    )
    arpes_data = ARPESData(
        da;
        intensity_unit = metadata_original[:intensity_unit],
        analyzer_config = metadata_original[:analyzer_configuration],
        energy_def = energy_def,
        additional_metadata = Dict(),
    )
    return arpes_data
end


