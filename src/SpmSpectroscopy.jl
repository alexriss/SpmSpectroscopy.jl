module SpmSpectroscopy

using CSV
using DataFrames
using DataStructures:OrderedDict
using Dates

export SpmSpectrum, load_spectrum, get_channel


@enum Direction bwd fwd

mutable struct SpmSpectrum
    filename::String
    header::AbstractDict
    data::DataFrame
    channel_names::Vector{String}
    channel_units::Vector{String}
    
    position::Vector{Float64}

    bias::Float64
    z_feedback::Bool
    
    start_time::DateTime
end
SpmSpectrum(filename::String) = SpmImage(filename, OrderedDict(), DataFrame(), String[], String[], Float64[], 0., false, Dates.now())

mutable struct SpmSpectrumChannel
    channel_name::String
    channel_unit::String
    direction::Direction
    data::Array{Float64}
end


"""
    function load_spectrum(filename::AbstractString; header_only::Bool=false)::SpmSpectrum

Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
If `header_only` is `true`, then only the header is loaded.
"""
function load_spectrum(filename::AbstractString; header_only::Bool=false)::SpmSpectrum
    ext = rsplit(filename, "."; limit=2)[end]
    if ext == "dat"
        spectrum = load_spectrum_nanonis(filename, header_only=header_only)
    else
        throw(ArgumentError("Unknown file type \"$ext\""))
    end

    return spectrum
end


"""
    function load_spectrum_nanonis(filename::AbstractString; header_only::Bool=false)::SpmSpectrum

Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
If `header_only` is `true`, then only the header is loaded.
"""
function load_spectrum_nanonis(filename::AbstractString; header_only::Bool=false)::SpmSpectrum
    contents_data = ""
    header = OrderedDict{String,Any}()
    data = DataFrame()
    channel_names = Vector{String}()
    channel_units = Vector{String}()
    position = Float64[NaN,NaN,NaN]
    bias = NaN
    z_feedback = false
    start_time = Dates.now()
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
            data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names)
        end
    end

    return SpmSpectrum(filename, header, data, channel_names, channel_units, position, bias, z_feedback, start_time)
end




end
