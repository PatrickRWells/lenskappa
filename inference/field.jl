using DataFrames
using CSV
using Base.Filesystem
using Distributions
using StatsBase
using LinearAlgebra
using Extremes

function load_field_weights(path:: AbstractString)
    if !isfile(path)
        error("File not found: $path")
    end
    df = CSV.read(path, DataFrame)
    df = df[completecases(df), :]
    return df
end

function model_weights(df:: DataFrame, weights:: Vector{String})::Dict{String, MaximumLikelihoodAbstractExtremeValueModel}
    dists = Dict{String, Distributions.Gumbel}()
    for weight in weights
        dist = model_weights(df, weight)
        dists[weight] = dist
    end
    return dists
end

function model_weights(df:: DataFrame, weight:: String, nbins:: Int64=100, min:: Float64=0.0, max::Float64=4.0)::MaximumLikelihoodAbstractExtremeValueModel
    if !(weight in names(df))
        error("Weight not found: $weight")
    end
    result = gevfit(df, Symbol(weight))
    return result
end
