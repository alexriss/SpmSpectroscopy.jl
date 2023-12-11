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

@enum Direction bwd fwd
@enum Background no_correction subtract_minimum linear_fit
    
include("spectrum_functions.jl")
include("AFM_functions.jl")

const VERSION = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "../Project.toml"))["version"])

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
SpmSpectrum(filename::String) = SpmImage(
    filename, OrderedDict(), DataFrame(), String[], String[],
    Float64[],
    0., false,
    DateTime(-1)
)


function Base.show(io::IO, s::SpmSpectrum)
    if get(io, :compact, false)
        print(io, "SpmSpectrum(\"", s.filename, "\")")
    else
        print(io, "SpmSpectrum(\"", s.filename, "\", ")
        if haskey(s.header, "Experiment")
            print(io, "Experiment: \"", s.header["Experiment"], "\", ")
        end
        print(io, length(s.channel_names), " channels, ", size(s.data, 1), " points)")
    end
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
    elseif ext == "vpdata"
        spectrum = load_spectrum_gsxm(filename, select=select, header_only=header_only, remove_missing=remove_missing, index_column=index_column, index_column_type=index_column_type)
    else
        throw(ArgumentError("Unknown file type \"$ext\""))
    end

    return spectrum
end
# precompile(load_spectrum, (String, ))
# ompile(Core.kwfunc(load_spectrum), (Vector{Any}, typeof(load_spectrum), String))


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
        experiment = split(readline(f), '\t')
        if length(experiment) >= 2 && experiment[1] == "Experiment"
            header[experiment[1]] = experiment[2]
        else
            error("This is not a Nanonis .dat file.")
        end
        while !eof(f)
            l = readline(f)
            if l == "[DATA]"
                l = readline(f)
                channel_names, channel_units = get_channel_names_units(l)
                break
            else
                header_data = split(l, '\t')
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
            z_feedback = header["Z-Ctrl hold"] == "FALSE" ? true : false
        end
        if haskey(header, "Date")
            try
                start_time = DateTime(header["Date"], dateformat"d.m.Y H:M:S")
            catch e
            end
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
# precompile(load_spectrum_nanonis, (String, ))
# precompile(Core.kwfunc(load_spectrum_nanonis), (Vector{Any}, typeof(load_spectrum_nanonis), String))


"""
    load_spectrum_gsxm(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
        index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum

Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
`select` can be used to specify which columns to load (see CSV.jl for an explanation of `select`).
If `header_only` is `true`, then only the header is loaded.
If `index_column` is `true`, then an extra column with indices of type `index_column_type` will be added.
If `remove_missing` is `true`, then missing values are dropped.
"""
function load_spectrum_gsxm(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, remove_missing::Bool=false,
    index_column::Bool=false, index_column_type::Type=Int)::SpmSpectrum

    contents_data = ""
    header = OrderedDict{String,String}()
    data = DataFrame()
    channel_names = Vector{String}()
    channel_units = Vector{String}()
    position = Float64[NaN,NaN,NaN]
    bias = NaN
    z_feedback = true
    start_time = DateTime(-1)
    open(filename) do f
        while !eof(f)
            l = readline(f)
            comment = false
            if startswith(l, "#C")
                l = l[3:end]
                comment = true
            elseif startswith(l, "#")
                l = l[2:end]
                comment = true
            end

            if comment && endswith(l, "data=")
                l = readline(f)
                l = l[3:end]  # remove the "#C"
                channel_names, channel_units = get_channel_names_units(l)
                break
            else
                header_data = split(l, "::", limit=2)
                length(header_data) == 2 || (header_data = split(l, "\t", limit=2))

                if length(header_data) == 2
                    key = strip(header_data[1])
                    val = strip(header_data[2])

                    # check if key already exists
                    if haskey(header, key)
                        i_key = 1
                        while haskey(header, key * " $(i_key)")
                            i_key += 1
                        end
                        key *= " $(i_key)"
                    end
    
                    header[key] = val
                end
            end
        end

        if haskey(header, "GXSM-Main-Offset")
            h = header["GXSM-Main-Offset"]
            posl = findfirst("X0=", h)
            posr = findnext("Ang", h, posl[end]+1)
            if !isnothing(posl) && !isnothing(posr)
                position[1] = parse(Float64, h[posl[end]+1:posr[1]-1]) * 1e-10  # convert from Angstrom to m
            end
            posl = findnext("Y0=", h, posr[end]+1)
            posr = findnext("Ang", h, posl[end]+1)
            if !isnothing(posl) && !isnothing(posr)
                position[2] = parse(Float64, h[posl[end]+1:posr[1]-1]) * 1e-10  # convert from Angstrom to m
            end
        end
        # todo: bias
        # todo: z-feedback
        if haskey(header, "Date")
            try
                date_strs = split(header["Date"], "=", limit=2)
                if length(date_strs) == 2
                    start_time = DateTime(date_strs[2], dateformat"e u d H:M:S y")
                end
            catch e
            end
        end

        if !header_only
            contents_data = read(f, String)
            if length(select) > 0
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names, missingstring="NaN", types=Float64, comment="#", select=select)
            else
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names, missingstring="NaN", types=Float64, comment="#")
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


"""
    get_channel_names_units(l::String)::Tuple{Vector{String}, Vector{String}}

Extracts the channel names and units from the line `l` of a Nanonis .dat file or GSXM .vpdata file.
"""
function get_channel_names_units(l::String)::Tuple{Vector{String}, Vector{String}}
    channels = split(l, "\t")
    channel_names = Vector{String}(undef, length(channels))
    channel_units = Vector{String}(undef, length(channels))
    for (i, ch) in enumerate(channels)
        posl = findlast('(', ch)
        posr = findlast(')', ch)
        if isnothing(posl) || isnothing(posr)
            ch_name = strip(ch)
            ch_unit = ""
        else
            ch_name = ch[1:posl-2]
            ch_unit = ch[posl+1:posr-1]
        end
        channel_names[i] = strip(ch_name, '"')
        channel_units[i] = ch_unit
    end
    return channel_names, channel_units
end


end