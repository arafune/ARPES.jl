# ✅ より良い実装
const LOCATION_REGISTRY = Dict{String,Type{<:LocationLoader}}("SPD" => SPDLoader)

function select_loader(fpath::String, loc::Union{String,Nothing})
    if !isnothing(loc)
        return get_loader_by_name(loc)
    end

    if !isnothing(OPTIONS.loc)
        return get_loader_by_name(OPTIONS.loc)
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
Return Loader type from given file path and location name.

# Arguments
 -`fpath::String`: File path.
 -`loc`::Union{String,Nothing}: Location name.

# Returns
 -`Type{<:LocationLoader}`: Loader type.
"""
function select_loader(fpath::String, loc::Union{String,Nothing})
    if !isnothing(loc)
        return get_loader_by_name(loc)  # loc = "SPD" → SPDLoader
    end
    #
    if !isnothing(OPTIONS.loc)
        return get_loader_by_name(OPTIONS.loc)
    end

    return detect_loader(fpath)
end

function get_loader_by_name(name::String)
    # name = "SPD"
    return LOCATION_REGISTRY[name]  # LOCATION_REGISTRY["SPD"] → SPDLoader
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
    elseif ext == ".h5" || ext == ".hdf5"
        return detect_hdf5_loader(fpath)
    else
        @warn "Unknown file extension: $ext, using GenericLoader"
        return GenericLoader
    end
end




"""ITX file"""
function detect_itx_loader(fpath::String)
    #
    header_lines = String[]
    open(fpath, "r") do io
        for i = 1:50
            eof(io) && break
            push!(header_lines, readline(io))
        end
    end

    header_text = join(header_lines, "\n")

    if occursin("SpecsLab Prodigy", header_text)
        return SPDLoader
    elseif occursin("MAX IV", header_text) || occursin("Bloch", header_text)
        return BlochLoader
    else
        @info "Could not determine location from ITX file, using GenericLoader"
        return GenericLoader
    end
end
