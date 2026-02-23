using DimensionalData
using DimensionalData: set, Dim, DimArray

export read_itx
"""
    read_itx(fpath::String)

Reads an ITX (Igor Text) file and parses its contents into structured data.

- Saves comment lines as raw data.
- Parses comments, wave declarations, and scale information into a dictionary.
- Extracts wave data and converts it into a `DimArray` with appropriate metadata.

# Arguments
- `fpath::String`: Path to the ITX file.

# Returns
- `DimArray`: A dimensional array containing the parsed wave data and metadata.
"""
function read_itx(fpath::String)
    lines = readlines(fpath)

    # Wave data
    waves_info = Dict{Union{String,Symbol},Any}()
    waves_info[:raw_comments] = String[]
    waves_info[:source_file] = fpath
    wave_data = nothing

    in_data_section = false
    data_lines = String[]

    for line in lines
        original_line = line
        line = strip(line)

        # Comment lines starting from X //
        if startswith(line, "X //")
            comment = replace(line, "X //" => "", count = 1)
            comment = strip(comment)

            push!(waves_info[:raw_comments], comment)
            parsed = _parse_comment_line(comment)
            if parsed !== nothing
                merge!(waves_info, parsed)
            end

        elseif startswith(line, "IGOR")
            continue

        elseif startswith(line, "WAVES")
            parsed = _parse_waves_declaration(line)
            if parsed !== nothing
                waves_info = merge(waves_info, parsed)
            end

        elseif startswith(line, "BEGIN")
            in_data_section = true

        elseif startswith(line, "END")
            in_data_section = false
            wave_data = _parse_wave_data(data_lines)

        elseif startswith(line, "X SetScale")
            scale_info = _parse_setscale(line)
            merge!(waves_info, scale_info)

        elseif in_data_section
            push!(data_lines, line)
        end
    end
    return _to_dimarray(wave_data, waves_info)
end

"""
    _parse_comment_line(comment::AbstractString)

Parses a comment line from the ITX file and extracts key-value pairs or metadata.

- Handles user comments, created date, and created by fields.
- Converts keys to symbols and values to appropriate types.

# Arguments
- `comment::AbstractString`: The comment line to parse.

# Returns
- `Dict{Symbol, Any}`: Dictionary of parsed metadata, or `nothing` if not applicable.
"""
function _parse_comment_line(comment::AbstractString)
    if startswith(comment, "User Comment")
        parts = split(comment, "=", limit = 2)
        additional_waves_info = Dict{Symbol,Any}()
        if length(parts) == 2
            key_raw = strip(parts[1])
            val_raw = strip(parts[2])
            additional_waves_info[:user_comment] = val_raw
            info_from_comments = _parse_comment_line_comment(String(val_raw))
            return merge(additional_waves_info, info_from_comments)
        end
    end

    if occursin("=", comment)
        parts = split(comment, "=", limit = 2)
        if length(parts) == 2
            key_raw = strip(parts[1])
            val_raw = strip(parts[2])
            key = Symbol(lowercase(replace(key_raw, " " => "_")))
            val = _parse_value(String(val_raw))
            return Dict(key => val)
        end
    end

    if startswith(comment, "Created Date")
        parts = split(comment, ":", limit = 2)
        if length(parts) == 2
            return Dict(:created_date => strip(parts[2]))
        end
    end

    if startswith(comment, "Created by")
        parts = split(comment, ":", limit = 2)
        if length(parts) == 2
            return Dict(:created_by => strip(parts[2]))
        end
    end
end

"""
    _build_scale(scale_info::Dict, length::Int)

Builds a numerical range for a dimension based on scale information.

# Arguments
- `scale_info::Dict`: Dictionary containing scale parameters (`:min`, `:max`, `:start`, `:step`).
- `length::Int`: Number of points in the range.

# Returns
- `AbstractRange`: The constructed range for the dimension.
"""
function _build_scale(scale_info::Dict, length::Int)
    if haskey(scale_info, :min)
        return range(scale_info[:min], scale_info[:max], length = length)
    elseif haskey(scale_info, :start)
        return range(scale_info[:start], step = scale_info[:step], length = length)
    else
        error("Invalid scale information")
    end
end



"""
    _parse_comment_line_comment(comment::String)

Parses a semicolon-separated comment string into key-value pairs.

Example:
    "beta:0;Temperature:RT;X:13.5;" → Dict(:beta => 0, :temperature => "RT", :x => 13.5)

# Arguments
- `comment::String`: The comment string to parse.

# Returns
- `Dict{Symbol, Any}`: Dictionary of parsed key-value pairs.

"""
function _parse_comment_line_comment(comment::String)
    additional_waves_info = Dict{Symbol,Any}()
    comments = split(comment, ";")
    for info in comments
        info_item = split(info, ":")
        if length(info_item) >= 2
            key = Symbol(lowercase(replace(info_item[1], " " => "_")))
            additional_waves_info[key] = _parse_value(String(info_item[2]))
        end
    end
    return additional_waves_info
end


"""
    _parse_value(val_str::String)

Converts a string value to an appropriate type (Int, Float64, or String).

# Arguments
- `val_str::String`: The string value to convert.

# Returns
- `Int`, `Float64`, or `String`: The parsed value.
"""
function _parse_value(val_str::String)
    val_str = strip(val_str)

    # Try to interpret as number
    try
        # Integer
        if !occursin(".", val_str) && !occursin("e", lowercase(val_str))
            return parse(Int, val_str)
        else
            # Float
            return parse(Float64, val_str)
        end
    catch
        # If not a number, return as String
        return val_str
    end
end

"""
    _parse_wave_data(data_lines::Vector{String})

Parses lines of wave data into a 1D array of Float64 values.

# Arguments
- `data_lines::Vector{String}`: Lines containing space-separated numeric values.

# Returns
- `Vector{Float64}`: Parsed numeric data.
"""
function _parse_wave_data(data_lines::Vector{String})
    all_values = Float64[]
    for line in data_lines
        if isempty(strip(line))
            continue
        end

        # Extract space-saparated numbers
        values = split(line)
        for v in values
            try
                push!(all_values, parse(Float64, v))
            catch
                @warn "Failed to parse value: $v"
            end
        end
    end
    return all_values
end


"""
    _to_dimarray(data::Array, waves_info::Dict)

Converts parsed data and metadata into a `DimArray` with appropriate dimensions and metadata.

# Arguments
- `data::Array`: The wave data array.
- `waves_info::Dict`: Metadata and dimension information.

# Returns
- `DimArray`: The constructed dimensional array.
"""
function _to_dimarray(data::Array, waves_info::Dict)
    # Construct dimension
    dims_list = []
    if haskey(waves_info, :shape)
        data = reshape(data, waves_info[:shape])
    end
    axes_keys = [:x_scale, :y_scale, :z_scale, :w_scale]
    for (i, key) in enumerate(axes_keys)
        if ndims(data) >= i && haskey(waves_info, key)
            info = waves_info[key]
            n = size(data, i)
            vals = _build_scale(info, n)
            default_name = Symbol(string(key)[1])
            dim_name = isempty(info[:label]) ? default_name : Symbol(info[:label])
            d = Dim{dim_name}(vals)
            if haskey(info, :unit) && !isempty(info[:unit])
                d = Dim{dim_name}(
                    vals;
                    metadata = DimensionalData.Dimensions.Metadata(
                        Dict(:unit => info[:unit]),
                    ),
                )
            else
                d = Dim{dim_name}(vals)
            end
            push!(dims_list, d)
        end
    end
    # Create DimArray
    da = DimArray(
        data,
        Tuple(dims_list);
        name = Symbol(get(waves_info, :name, "data")),
        metadata = DimensionalData.Metadata(waves_info),
    )

    return da
end

"""
    _parse_waves_declaration(line::AbstractString)

Parses a WAVES declaration line from the ITX file and extracts wave metadata.

Example:
    WAVES/S/N=(200,2401) 'GrIr111_1'
    → Dict(:type => :single_precision, :shape => (200, 2401), :name => "GrIr111_1")

# Arguments
- `line::AbstractString`: The WAVES declaration line.

# Returns
- `Dict{Symbol, Any}`: Dictionary of parsed wave metadata.
"""
function _parse_waves_declaration(line::AbstractString)
    info = Dict{Symbol,Any}()

    # Extract flag part (/S, /N=...)
    flags = String[]
    if occursin("/", line)
        flag_part = split(line, " ")[1]  # "WAVES/S/N=(200,2401)"
        flags = split(flag_part, "/")[2:end]  # ["S", "N=(200,2401)"]
    end

    # Extract wave name

    name_match = match(r"'([^']+)'", line)
    if !isnothing(name_match)
        info[:name] = name_match.captures[1]
    else
        parts = split(line, " ")
        if length(parts) > 1
            info[:name] = parts[end]
        end
    end
    for flag in flags
        if startswith(flag, "N=")
            shape_str = replace(flag, "N=" => "")
            if startswith(shape_str, "(") && endswith(shape_str, ")")
                shape_str = replace(shape_str, "(" => "")
                shape_str = replace(shape_str, ")" => "")
                dims = parse.(Int, split(shape_str, ","))
                info[:shape] = Tuple(dims)
                info[:n] = dims[1]
            else
                n = parse(Int, shape_str)
                info[:n] = n
            end
        elseif flag == "S"
            info[:type] = :single_precision
        elseif flag == "D"
            info[:type] = :double_precision
        end
    end
    if !haskey(info, :type)
        info[:type] = :double_precision
    end

    return info
end

"""
    _parse_setscale(line::AbstractString)

Parses a SetScale line from the ITX file and extracts axis scaling information.

Examples:
    X SetScale/I x, -11.44, 12.44, "deg (theta_y)", 'GrIr111_1'
    → x-axis: -11.44 to 12.44, unit: "deg", label: "theta_y"

    SetScale/P x, 5.2, 0.002, "eV", "ID_001"
    → x-axis: start=5.2, step=0.002, unit: "eV", label: nothing

# Arguments
- `line::AbstractString`: The SetScale line to parse.

# Returns
- `Dict{Symbol, Any}`: Dictionary containing axis scale information.
"""
function _parse_setscale(line::AbstractString)
    # SetScale/I x, min, max, "unit (label)", 'wave'
    # SetScale/P x, start, step, "unit", "label"
    parts = split(line, ",")
    if length(parts) >= 3
        # Determine mode: /I (interval) or /P (point+step)
        mode = ""
        if occursin("/I", parts[1])
            mode = :inclusive
        elseif occursin("/P", parts[1])
            mode = :perpoints
        else
            mode = :default
        end

        axis_part = strip(split(parts[1], " ")[end])  # "x"
        axis_key = Symbol("$(axis_part)_scale")

        unit, labe = "", ""
        if length(parts) >= 4
            unit_part = strip(parts[4])
            if occursin("(", unit_part) && occursin(")", unit_part)
                unit_match = match(r"\"([^(]+)\s*\(([^)]+)\)\"", unit_part)
                if !isnothing(unit_match)
                    unit = strip(unit_match.captures[1])
                    label = replace(
                        lowercase(strip(unit_match.captures[2])),
                        "-" => "_",
                        " " => "_",
                    )
                end
            else
                unit = replace(unit_part, "\"" => "")
                label = nothing
            end
        end

        name = ""
        if length(parts) >= 5
            name = strip(parts[5])
            name = replace(name, "\"" => "")
            name = replace(name, "'" => "")
        end

        if mode == :inclusive
            min_val = parse(Float64, strip(parts[2]))
            max_val = parse(Float64, strip(parts[3]))
            return Dict(
                axis_key =>
                    Dict(:min => min_val, :max => max_val, :unit => unit, :label => label),
            )
        elseif mode == :perpoints
            start_val = parse(Float64, strip(parts[2]))
            step_val = parse(Float64, strip(parts[3]))
            return Dict(
                axis_key => Dict(
                    :start => start_val,
                    :step => step_val,
                    :unit => unit,
                    :label => label,
                ),
            )
        end
    end
end

