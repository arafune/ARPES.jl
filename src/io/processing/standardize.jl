using DimensionalData
using DimensionalData.Dimensions: name, val
using ..IO: LocationLoader
using ..IO.Location: canonical_dim

import ..IO: to_standardize

function to_standardize(loader::Type{<:LocationLoader}, raw::DimArray)
    new_dims = map(dims(raw)) do d
        raw_name = name(d)
        ctor=canonical_dim(loader, raw_name)
        isnothing(ctor) ? d : ctor(val(d))
    end

    return rebuild(raw; dims = Tuple(new_dims))
end
