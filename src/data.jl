using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

"""
    cat_arpes(args::ARPESData...; dims::Union{Symbol,Dim})

Concatenate multiple `ARPESData` objects along the specified dimension.

# Arguments
- `args::ARPESData...`: One or more `ARPESData` objects to concatenate.
- `dims::Union{Symbol,Dim}`: The dimension (as a `Symbol` or `Dim`) along which to concatenate.

# Throws
- `ArgumentError` if any of the input data objects do not have the specified dimension.

# Returns
- A new `ARPESData` object resulting from concatenation along the specified dimension.

# Example
```julia
cat_arpes(data1, data2; dims=:energy)
```
"""
function cat_arpes(args::ARPESData...; dims::Union{Symbol,Dim})
    missing_indices = findall(d -> !hasdim(d, dims), args)
    if !isempty(missing_indices)
        error_msg = "Dimension :$dims not found in the following argument(s): $(join(missing_indices, ", "))"
        throw(ArgumentError(error_msg))
    end
    lookups = map(d -> lookup(d, dims), args)
    combined_lookup_data = vcat(map(collect, lookups)...)
    proto_lookup = first(lookups)
    target_dims = Dim{dims}(combined_lookup_data, metadata = metadata(proto_lookup))
    return cat(args...; dims = target_dims)
end
