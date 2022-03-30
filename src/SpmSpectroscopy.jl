module SpmSpectroscopy

using CSV
using DataFrames
using DataStructures: OrderedDict
using Dates
using Printf
using Statistics
using TOML
using LinearAlgebra: pinv

export SpmSpectrum, load_spectrum,
    correct_background!, no_correction, subtract_minimum, linear_fit
    # deconvolve_sader_jarvis, deconvolve_matrix,
    # savitzky_golay_filter

include("spectrum_functions.jl")
include("AFM_functions.jl")

const VERSION = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "../Project.toml"))["version"])


@enum Direction bwd fwd
@enum Background no_correction subtract_minimum linear_fit

mutable struct SpmSpectrum
    filename::String
    header::OrderedDict{String,String}
    data::DataFrame
    channel_names::Vector{String}
    channel_units::Vector{String}
    
    position::Vector{Float64}

    bias::Float64
    z_feedback::Bool
    
    start_time::DateTime
end
SpmSpectrum(filename::String) = SpmImage(filename, OrderedDict(), DataFrame(), String[], String[], Float64[], 0., false, DateTime(-1))

mutable struct SpmSpectrumChannel
    channel_name::String
    channel_unit::String
    direction::Direction
    data::Array{Float64}
end


"""
    load_spectrum(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
        index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum


Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
`select` can be used to specify which columns to load (see CSV.jl for an explanation of `select`).
If `header_only` is `true`, then only the header is loaded.
If `index_column` is `true`, then an extra column with indices of type `index_column_type` will be added.
If `remove_missing` is `true`, then missing values are dropped.
"""
function load_spectrum(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
    index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum

    ext = rsplit(filename, "."; limit=2)[end]
    if ext == "dat"
        spectrum = load_spectrum_nanonis(filename, select=select, header_only=header_only, remove_missing=remove_missing, index_column=index_column, index_column_type=index_column_type)
    else
        throw(ArgumentError("Unknown file type \"$ext\""))
    end

    return spectrum
end
precompile(load_spectrum, (String, ))


"""
    load_spectrum_nanonis(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
        index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum

Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
`select` can be used to specify which columns to load (see CSV.jl for an explanation of `select`).
If `header_only` is `true`, then only the header is loaded.
If `index_column` is `true`, then an extra column with indices of type `index_column_type` will be added.
If `remove_missing` is `true`, then missing values are dropped.
"""
function load_spectrum_nanonis(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
    index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum

    contents_data = ""
    header = OrderedDict{String,String}()
    data = DataFrame()
    channel_names = Vector{String}()
    channel_units = Vector{String}()
    position = Float64[NaN,NaN,NaN]
    bias = NaN
    z_feedback = false
    start_time = DateTime(-1)
    open(filename) do f
        while !eof(f)
            l = readline(f)
            if l == "[DATA]"
                l = readline(f)
                channels = split(l, "\t")
                channel_names = Vector{String}(undef, length(channels))
                channel_units = Vector{String}(undef, length(channels))
                for i = 1:length(channels)
                    posl = findlast('(', channels[i])
                    posr = findlast(')', channels[i])
                    channel_names[i] = channels[i][1:posl-2]
                    channel_units[i] = channels[i][posl+1:posr-1]
                end
                break
            else
                header_data = split(l, "\t")
                if length(header_data) >= 2
                    header[header_data[1]] = header_data[2]
                end
            end
        end

        if haskey(header, "X (m)")
            position[1] = parse(Float64, header["X (m)"])
        end
        if haskey(header, "Y (m)")
            position[2] = parse(Float64, header["Y (m)"])
        end
        if haskey(header, "Z (m)")
            position[3] = parse(Float64, header["Z (m)"])
        end
        if haskey(header, "Bias (V)")
            bias = parse(Float64, header["Bias (V)"])
        end
        if haskey(header, "Z-Ctrl hold")
            if header["Z-Ctrl hold"] == "FALSE"
                z_feedback = true
            end
        end
        if haskey(header, "Date")
            start_time = DateTime(header["Date"], dateformat"d.m.Y H:M:S")
        end

        if !header_only
            contents_data = read(f, String)
            if length(select) > 0
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names, missingstring="NaN", types=Float64, select=select)
            else
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names, missingstring="NaN", types=Float64)
            end

            if remove_missing
                dropmissing!(data)
            end

            if index_column
                if index_column_type == Int
                    data[!,"Index"] = 1:size(data,1)
                else
                    data[!,"Index"] = convert.(index_column_type, 1:size(data,1))
                end
            end
        end
        if index_column
            push!(channel_names, "Index")
            push!(channel_units, "")
        end
    end

    return SpmSpectrum(filename, header, data, channel_names, channel_units, position, bias, z_feedback, start_time)
end
precompile(load_spectrum_nanonis, (String, ))


"""
    correct_background!(xdata<:Vector{<:Real}, ydata<:Vector{<:AbstractFloat}, type::Background, offset::Bool=true)::Nothing

Background correction of `ydata` vs. `xdata` with using a correction of type `type`.
If `offset` is `true` (default), then `ydata` will be shifted such that its minimum is 0.
"""
function correct_background!(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:AbstractFloat}, type::Background, offset::Bool=true)::Nothing
    if type == no_correction
        return nothing
    end
    if length(ydata) == 0
        return nothing
    end
    if type == linear_fit
        # https://en.wikipedia.org/wiki/Ordinary_least_squares#Simple_linear_regression_model

        varx = var(xdata)
        if varx > 0
            meanx = mean(xdata)
            β = cov(xdata, ydata) / varx
            α = mean(ydata) - β * meanx
            @. ydata = ydata - α - xdata * β
        end
    end
    if type == subtract_minimum || offset
        m = minimum(filter(!isnan, ydata))
        ydata .-= m
    end
    return nothing
end
precompile(correct_background!, (Vector{Float64}, Vector{Float64}, Background, Bool))

end
