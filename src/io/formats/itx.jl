using DimensionalData

"""
Read itx file

# Memo:

1. Save comment lines as raw data
2. Parse into structured information
"""
function read_itx(fpath::String)
    lines = readlines(fpath)


    # Wave data
    waves_info = Dict{String,Any}()
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

            merge(waves_info, _parse_comment_line(comment))

        elseif startswith(line, "IGOR")
            continue

        elseif startswith(line, "WAVES")
            waves_info = _parse_waves_declaration(line)

        elseif startswith(line, "BEGIN")
            in_data_section = true

        elseif startswith(line, "END")
            in_data_section = false
            wave_data = _parse_wave_data(data_lines)

        elseif startswith(line, "X SetScale")
            scale_info = _parse_setscale!(line)
            merge!(waves_info, scale_info)

        elseif in_data_section
            push!(data_lines, line)
        end
    end

    return _to_dimarray(wave_data, waves_info)
end

"""
Parse comment line

X //Scan Mode         = Fixed Analyzer Transmission
→ Dict{Symbol, Any}(
:scan_mode => "Fixed Analyzer Transmission"
)
"""
function _parse_comment_line(comment::String)
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
"beta:0;Temperature:RT;X:13.5;Y:21.55;Z:+00346000;theta:-00270000;position:187.4655;UV(P);IR(P);+w;P:113mW;+3w;P:6mW;"
->
Dict{Symbol, Any}(
:beta => 0
:p => "6mW"
:temperature => "RT"
:y => 21.55
:position => 187.4655
:z => 346000
:theta => -270000
:x => 13.5)
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
Convert value to appropriate type
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
Parse data as simple 1D array
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
Convert to DimArray
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
            vals = _build_scale(info, length = n)
            default_name = Symbol(string(key)[1])
            dim_name = isempty(info[:label]) ? default_name : Symbol(info[:label])
            d = Dim{dim_name}(vals)
            if haskey(info, :unit) && !isempty(info[:unit])
                d = set(d, MetaData(Dict(:unit => info[:unit])))
            end
            push!(dims_list, d)
        end
    end

    # Create DimArray
    da = DimArray(
        data,
        Tuple(dims_list);
        name = Symbol(get(waves_info, :name, "data")),
        metadata = metadata,
    )

    return da
end

"""
Parse WAVES declaration

WAVES/S/N=(200,2401) 'GrIr111_1'
→ Dict(:type => :single_precision, :shape => (200, 2401), :name => "GrIr111_1")
"""
function _parse_waves_declaration(line::String)
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
Parse SetScale
X SetScale/I x, -11.44, 12.44, "deg (theta_y)", 'GrIr111_1'
→ x-axis: -11.44 to 12.44, unit: "deg", label: "theta_y"

SetScale/P x, 5.2, 0.002, "eV", "ID_001"
→ x-axis: start=5.2, step=0.002, unit: "eV", label: nothing
"""
function _parse_setscale!(line::String)
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

