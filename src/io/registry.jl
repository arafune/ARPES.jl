"""
Registry mapping location names to their corresponding loader types.

# Example
LOCATION_REGISTRY["SPD"] => SPDLoader
"""
const LOCATION_REGISTRY = Dict{String,Type{<:LocationLoader}}("SPD" => SPDLoader)

"""
Return Loader type from given file path and location name.

# Arguments
 -`fpath::String`: File path.
 -`loc`::Union{String,Nothing}: Location name.

# Returns
 -`Type{<:LocationLoader}`: Loader type.
"""
function select_loader(fpath::String, loc::Union{String,Nothing})
    if !isnothing(loc)
        return get_loader_by_name(loc)
    end

    return detect_loader(fpath)
end

function get_loader_by_name(name::String)
    if haskey(LOCATION_REGISTRY, name)
        return LOCATION_REGISTRY[name]
    else
        error("Unknown location: '$name'. Available: $(keys(LOCATION_REGISTRY))")
    end
end

"""
Judge the loader type from the file content.

# Arguments
- `fpath::String`: File path

# Returns
- `Type{<:LocationLoader}`: Loader 
"""
function detect_loader(fpath::String)
    ext = lowercase(splitext(fpath)[2])

    if ext == ".itx"
        return detect_itx_loader(fpath)
    elseif ext == ".sp2"
        return detect_sp2_loader(fpath)
    else
        @warn "Unknown file extension: $ext, using GenericLoader"
        return GenericLoader
    end
end


"""
Detect the loader type for ITX files by inspecting the file header.

# Arguments
- `fpath::String`: File path

# Returns
- `Type{<:LocationLoader}`: Loader type
"""
function detect_itx_loader(fpath::String)
    #
    header_lines = String[]
    open(fpath, "r") do io
        for _ = 1:50
            eof(io) && break
            push!(header_lines, readline(io))
        end
    end

    header_text = join(header_lines, "\n")

    if occursin("SpecsLab Prodigy", header_text) || occursin("R. Arafune", header_text)
        return SPDLoader
    else
        @info "Could not determine location from ITX file, using GenericLoader"
        return GenericLoader
    end
end

"""
Detect the loader type for SP2 files.

# Arguments
- `fpath::String`: File path

# Returns
- `Type{<:LocationLoader}`: Loader type
"""
function detect_sp2_loader(fpath::String)
    return SPDLoader
end
