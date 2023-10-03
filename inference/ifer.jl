include("./ms.jl")
include("./field.jl")
using Base.Filesystem
using DataFrames
using CSV
using Extremes
using Plots

function load_all(field_weight_path::AbstractString, ms_weight_path:: AbstractString, kappa_path::AbstractString)
    if !isdir(ms_weight_path) || !ispath(ms_weight_path)
        error("Expected a folder at $ms_weight_path")
    end
    if !isdir(kappa_path) || !ispath(kappa_path)
        error("Expected a folder at $kappa_path")
    end
    if !isfile(field_weight_path)
        error("File not found: $field_weight_path")
    end
    ms_weights = attach_kappas(ms_weight_path, kappa_path)
    weights = load_field_weights(field_weight_path)
    weight_model = model_weights(ms_weights, "1/r")
    model_params = params(weight_model)
    bins = Vector(0.0:0.01:4.0)
    loc = model_params[1][1]
    scale = model_params[1][2]
    shape = model_params[1][3]
    pdf_values = pdf(GeneralizedExtremeValue(loc, scale, shape), bins)
    Plots.histogram(ms_weights[!, "1/r"], bins = bins, normed=true, label="Field")
    Plots.plot!(bins, pdf_values, label="Gumbel")
    println("loc: $loc, scale: $scale, shape: $shape")
    Plots.savefig("test.png")
    return weight_model
end


kappa_path = "/Volumes/workspace/data/ms/maps/Plane35/kappa"
field_path = "/Users/patrick/Documents/Documents/Work/Research/Environment/Hierarchical/final/field/SL2SJ0212-0555/120_24.csv"
ms_path = "/Users/patrick/Documents/Documents/Work/Research/Environment/Hierarchical/final/ms/SL2SJ0212-0555/120_24"
load_all(field_path, ms_path, kappa_path)