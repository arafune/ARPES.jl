using DimensionalData
using DimensionalData: Dimension
using DimensionalData: basetypeof
using DimensionalData.Lookups
using DimensionalData.Dimensions: label

export stitch_along

"""
    stitch_along(A, B, dim; occupation_ratio=nothing, enhance_a=1.0)

Stitch two `AbstractDimArray`s together along `dim`.

- `occupation_ratio`: fraction of the overlap region assigned to the left array.
  `nothing` (default) means no overlap splitting — arrays are concatenated directly.
- `enhance_a`: multiplicative scale factor applied to `A` before concatenating.
"""
function stitch_along(
    A::AbstractDimArray,
    B::AbstractDimArray,
    dim::Union{Dimension,Symbol},
    occupation_ratio::Union{Real,Nothing},
    enhance_a::Real,
)
    @assert name(dims(A)) == name(dims(B))
    otherdims_A = otherdims(A, dim)
    otherdims_B = otherdims(B, dim)
    @assert otherdims_A == otherdims_B
    extrema_A = extrema(lookup(A, dim))
    extrema_B = extrema(lookup(B, dim))
    @debug "extrema_A" extrema_A
    @debug "extrema_B" extrema_B
    if isnothing(occupation_ratio) ||
       (extrema_A[2] < extrema_B[1]) ||
       (extrema_B[2] < extrema_A[1])
        @debug "concatenating"
        return _cat_and_sort(A, B, dim, enhance_a)
    end
    @assert 0.0 <= occupation_ratio <= 1.0 "occupation_ratio must be between 0 and 1"
    if extrema_A[1] < extrema_B[1]
        left, right = A, B
    elseif extrema_B[1] < extrema_A[1]
        left, right = B, A
    else
        error(
            "The arrays have the same extrema, cannot decide which one is left and which one is right",
        )
    end
    seam_dim =
        (maximum(lookup(left, name(dim))) - minimum(lookup(right, name(dim)))) *
        occupation_ratio + minimum(lookup(right, name(dim)))
    @debug "seam_dim" seam_dim
    left_part = left[Dim{name(dim)}(At(<(seam_dim)))]
    right_part = right[Dim{name(dim)}(At(>=(seam_dim)))]
    return _cat_and_sort(left_part, right_part, dim, enhance_a)
end



stitch_along(
    A::AbstractDimArray,
    B::AbstractDimArray,
    dim::Union{Dimension,Symbol};
    occupation_ratio::Union{Real,Nothing} = nothing,
    enhance_a::Real = 1.0,
) = stitch_along(A, B, dim, occupation_ratio, enhance_a)

function _cat_and_sort(
    A::AbstractDimArray,
    B::AbstractDimArray,
    dim::Union{Dimension,Symbol},
    enhance_a,
)
    along_dim = dims(A, dim)
    along_dim_idx = dimnum(A, dim)
    cat_dim = vcat(lookup(A, along_dim), lookup(B, along_dim))
    new_along_dim = rebuild(
        along_dim,
        Sampled(cat_dim, ForwardOrdered(), Irregular(), Points(), metadata(along_dim)),
    )
    @debug "cat_dim" cat_dim
    @debug "new_along_dim" new_along_dim
    cat_data = cat(parent(A) .* enhance_a, parent(B), dims = along_dim_idx)
    @debug "size_of cat_data" size(cat_data)
    stitch_dimArray = rebuild(
        A;
        data = cat_data,
        dims = Base.setindex(dims(A), new_along_dim, along_dim_idx),
    )
    p = sortperm(lookup(stitch_dimArray, new_along_dim))
    return stitch_dimArray[Dim{name(new_along_dim)}(p)]
end
