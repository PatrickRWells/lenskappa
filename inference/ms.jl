using CSV
using DataFrames
using Base.Filesystem
using JSON
using Statistics
using Plots
using CurveFit

this_file = abspath(@__FILE__)
this_dir = dirname(this_file)
config = JSON.parsefile(joinpath(this_dir, "config.json"))

function strip_units(value:: AbstractString)::Float64
    split_str = split(value, " ")
    return parse(Float64, split_str[1])
end

function clean_column(column:: Vector)::Vector{Float64}
    try
        parsed_column = convert.(Float64, column)
        return parsed_column
    catch
        parsed_column = strip_units.(column)
        return parsed_column        
    end    
end


function normalize_column(df:: DataFrames.DataFrame, col_name:: AbstractString)::DataFrames.DataFrame
    col = df[!, col_name]
    col = clean_column(col)
    median = Statistics.median(col)
    col = col ./ median
    df[!, col_name] = col
    return df
end

function load_single_weightfile(path:: AbstractString)::DataFrames.DataFrame
    if !isfile(path)
        error("File not found: $path")
    end
    df = CSV.read(path, DataFrame)
    df = df[completecases(df), :]
    for column_name = names(df)
        if column_name in config["columns_to_ignore"]
            continue
        end
        df = normalize_column(df, column_name)
    end
    return df
end

function load_all_weightfiles(path:: AbstractString)::Dict{Tuple{Int, Int}, DataFrames.DataFrame}
    if isfile(path)
        error("Folder not found: $path")
    end
    files = readdir(path)
    df_dict = Dict{Tuple{Int, Int}, DataFrames.DataFrame}()
    re = r"[0-9]_[0-9]"
    for file in files
        if !endswith(file, ".csv")
            continue
        end
        m = match(re, file)
        if isnothing(m)
            continue
        end
        field_name = m.match
        field_tuple = Tuple(parse.(Int, split(field_name, "_")))

        file_path = joinpath(path, file)
        df = load_single_weightfile(file_path)
        df_dict[field_tuple] = df
    end
    return df_dict
end

function load_single_field_kappas(path:: AbstractString)::Matrix{Float64}
    if !isfile(path)
        error("File not found: $path")
    end
    if !endswith(path, ".kappa")
        file_name = basename(path)
        error("Expected .kappa binary file, found $file_name")
    end
    data = zeros(Float32, 4096, 4096)
    read!(path, data)
end

function get_kappas_in_folder(path:: AbstractString)::Dict{Tuple{Int, Int}, String}
    files = readdir(path)
    kappas_dict = Dict{Tuple{Int, Int}, String}()
    re = r"8_[0-9]_[0-9]_N"
    for file in files
        if !endswith(file, ".kappa")
            continue
        end
        m = match(re, file)
        if isnothing(m)
            continue
        end
        field_name = m.match
        tuple_elements = split(field_name, "_")[2:3]
        field_tuple = Tuple(parse.(Int, tuple_elements))
        file_path = joinpath(path, file)
        kappas_dict[field_tuple] = file_path
    end
    return kappas_dict
end

function get_indices_from_coordinates(ra:: Vector{Float64}, dec:: Vector{Float64})::Vector{CartesianIndex{2}}
     if any(ra .< 0) || any(ra .> 4)
        println("Reindexing ra...")
        if all(ra .> -2) 
            ra = ra .+ 2.0
        else
            error("RA values are not in the expected range")
        end
    end

    if any(dec .< 0) || any(dec .> 4)
        println("Reindexing dec...")
        if all(dec .> -2) 
            dec = dec .+ 2.0
        else
            error("DEC values are not in the expected range")
        end
    end
    dec = dec .% 360
    ra = ra .% 360
    # We should now have two arrays with everything between 0 and 4 degrees
    l_field = 4.0 # degrees
    n_side = 4096. # pixels
    l_pix = l_field / n_side # degrees
    x_pix = (ra ./ l_pix) .- 0.5
    y_pix = (dec ./ l_pix) .- 0.5

    return CartesianIndex.(Int32.(round.(y_pix)), Int32.(round.(x_pix)))
end

function attach_kappas(weights_path::AbstractString, kappas_path::AbstractString)
    if !isdir(weights_path)
        error("Folder not found: $folder_path")
    end
    if !isdir(kappas_path)
        error("Expected a folder with kappa files, found $kappas_path")
    end
    weights = load_all_weightfiles(weights_path)
    if isempty(weights)
        error("No weight files found in $weights_path")
    end
    kappas = get_kappas_in_folder(kappas_path)
    good_dataframes = []
    for field_tuple in keys(weights)
        if field_tuple in keys(kappas)
            field_kappas = load_single_field_kappas(kappas[field_tuple])
            indices = get_indices_from_coordinates(weights[field_tuple][!, "ra"], weights[field_tuple][!, "dec"])
            kappa_values = field_kappas[indices]
            weights[field_tuple][!, "kappa"] = kappa_values
            push!(good_dataframes, weights[field_tuple]) #Careful, julia and numpy have different indexing conventions
        else
            println("No kappa found for $field_tuple, skipping...")
        end
    end 
    return vcat(good_dataframes...)
end