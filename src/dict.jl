export merge_consensus

"""
    merge_consensus(dicts::AbstractDict{K,V}...) where {K,V}

Merges multiple dictionaries into a single dictionary, keeping only the key-value pairs that are consistent across all input dictionaries.

If a key appears in more than one dictionary with different values, it is excluded from the result.

# Arguments
- `dicts...`: Any number of dictionaries with the same key and value types.

# Returns
A `Dict{K,V}` containing only the key-value pairs that are not in conflict across all input dictionaries.

# Example
```julia
d1 = Dict(:a => 1, :b => 2)
d2 = Dict(:a => 1, :b => 3)
merge_consensus(d1, d2) # returns Dict(:a => 1)
```
"""
function merge_consensus(dicts::AbstractDict{K,V}...) where {K,V}
    result = Dict{K,V}()
    conflicted = Set{K}()

    for d in dicts
        for (k, v) in d
            if k in conflicted
                continue
            elseif haskey(result, k)
                if result[k] != v
                    delete!(result, k)
                    push!(conflicted, k)
                end
            else
                result[k] = v
            end
        end
    end

    return result
end
