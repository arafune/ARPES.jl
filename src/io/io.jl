module IO

function load(fpath::String; loc = nothing)
    resolved_path = resolve_path(fpath)

    loader_type = select_loader(resolved_path, loc)
end
end

