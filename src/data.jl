using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

function vstack_arpes(datas::Vector{ARPESData})
    # Check that all ARPESData objects have the same dimensions and metadata
    first_dims = dims(datas[1])
    first_metadata = metadata(datas[1])
    
   
    # Stack the data arrays along a new dimension
    stacked_data = vcat([parent(data) for data in datas]...)
    
    # Create a new ARPESData object with the stacked data and original dimensions/metadata
    return ARPESData(stacked_data, dims=dims(datas[1]), metadata=metadata(datas[1]))
end
