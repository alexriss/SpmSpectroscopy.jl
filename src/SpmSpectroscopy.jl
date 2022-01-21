module SpmSpectroscopy

using CSV
using DataFrames
using DataStructures:OrderedDict
using Dates
using Statistics

export SpmSpectrum, load_spectrum, get_channel
export correct_background!, no_correction, subtract_minimum, linear_fit

@enum Direction bwd fwd
@enum Background no_correction subtract_minimum linear_fit

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
    function load_spectrum(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false,
        index_column::Bool=false, index_column_type::Type=Int64)::SpmSpectrum


Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
`select` can be used to specify which columns to load (see CSV.jl for an explanation of `select`).
If `header_only` is `true`, then only the header is loaded.
If `index_column` is `true`, then an extra column with indices of type `index_column_type` will be added.
"""
function load_spectrum(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false,
    index_column::Bool=false, index_column_type::Type=Int64)::SpmSpectrum

    ext = rsplit(filename, "."; limit=2)[end]
    if ext == "dat"
        spectrum = load_spectrum_nanonis(filename, select=select, header_only=header_only, index_column=index_column, index_column_type=index_column_type)
    else
        throw(ArgumentError("Unknown file type \"$ext\""))
    end

    return spectrum
end


"""
    function load_spectrum_nanonis(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false,
        index_column::Bool=false, index_column_type::Type=Int64)::SpmSpectrum

Loads a spectrum from the file `filename`. Currently, only Nanonis .dat files are supported.
`select` can be used to specify which columns to load (see CSV.jl for an explanation of `select`).
If `header_only` is `true`, then only the header is loaded.
If `index_column` is `true`, then an extra column with indices of type `index_column_type` will be added.
"""
function load_spectrum_nanonis(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false,
    index_column::Bool=false, index_column_type::Type=Int64)::SpmSpectrum

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
            if length(select) > 0
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names, select=select)
            else
                data = CSV.read(IOBuffer(contents_data), DataFrame, header=channel_names)
            end
            
            if index_column
                if index_column_type == Int64
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
    function correct_background!(xdata<:Vector{<:Real}, ydata<:Vector{<:AbstractFloat}, type::Background, offset::Bool=true)::Nothing

Background correction of `ydata` vs. `xdata` with using a correction of type `type`.
If `offset` is `true` (default), then `ydata` will be shifted such that its minimum is 0..
"""
function correct_background!(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:AbstractFloat}, type::Background, offset::Bool=true)::Nothing
    if type == no_correction
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


"""
    function trapz(X::AbstractVector{Float64}, Y::AbstractVector{Float64})::Float64

Trapezoidal integration over a function fiven by dicrete points in the arrays `Y` vs `X`,
where the spacing between the points in `X` is not necessarily constant.
"""
function trapz(X::AbstractVector{Float64}, Y::AbstractVector{Float64})::Float64
    @assert length(X) == length(Y)
    res = 0.
    for i = 1:length(X) - 1
        res += (X[i+1] - X[i]) * (Y[i] + Y[i+1])
    end
    return res / 2.
end


"""
    function rolling_mean(arr::AbstractArray{Float64}, n::Int64)::AbstractArray{Float64}

Computes the rolling mean over `n` points. The length of the output array is decreased by `n-1` points.

Adapted from:
https://stackoverflow.com/questions/59562325/moving-average-in-julia
"""
function rolling_mean(arr::AbstractArray{Float64}, n::Int64)::AbstractArray{Float64}
    so_far = sum(@view arr[1:n])
    out = zeros(Float64, length(arr) - n + 1)
    out[1] = so_far / n
    for (i, (start, stop)) in enumerate(zip(arr, @view arr[n+1:end]))
        so_far += stop - start
        out[i+1] = so_far / n
    end
    return out
end


"""
    function deconvolve_sadar_jarvis(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    
AFM force deconvolution using the Sadar-Jarvis method, as described in [Appl. Phys. Lett. 84, 1801 (2004)](https://aip.scitation.org/doi/10.1063/1.1667267).
Values for the tip-height `z` should be in ascending order. The corresponding values for the frequency shift are given in the vector `Δf`, along with
the experimental parameters `f₀` (resonance frequency), `A` (oscillation amplitude), and `k` (cantilever stiffness).

Based on MATLAB code from [Beilstein J. Nanotechnol. 3, 238 (2012)](https://www.beilstein-journals.org/bjnano/articles/3/27).
"""
function deconvolve_sadar_jarvis(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    @assert issorted(z)
    
    ω = Δf / f₀
    dω_dz = diff(ω) ./ diff(z)
    
    F = Vector{Float64}(undef, length(z) - 2)

    # we use `end-1` for both z and ω because dω_dz is shorter by one element
    for j=1:length(z)-3
        
        # adjust length of z, ω and dω_dz
        t = @view z[j+1:end-1]
        ω_ = @view ω[j+1:end-1]
        dω_dz_ = @view dω_dz[j+1:end]     
        
        integral = trapz(
            t,
            @. (1 + sqrt(A) / (8 * sqrt(π * (t - z[j])))) * ω_ - A^(3/2) / sqrt(2 * (t - z[j])) * dω_dz_
        )
                    
        # correction terms
        corr1 = ω[j] * (z[j+1] -z[j])
        corr2 = 2 * (sqrt(A) / (8 * sqrt(pi))) * ω[j] * sqrt(z[j+1] - z[j])
        corr3 = (-2) * (sqrt(A)^3 / sqrt(2)) * dω_dz[j] * sqrt(z[j+1] - z[j])
        F[j] = 2 * k * (corr1 + corr2 + corr3 + integral)
    end
    
    # adjust length of z
    z_ = z[1:length(F)]

    return z_, F
end


"""
AFM force deconvolution using the Matrix method, as described in [Appl. Phys. Lett. 78, 123 (2001)](https://aip.scitation.org/doi/10.1063/1.1335546).
Values for the tip-height `z` should be qually spaced and in ascending order. The corresponding values for the frequency shift are given in the vector `Δf`, along with
the experimental parameters `f₀` (resonance frequency), `A` (oscillation amplitude), and `k` (cantilever stiffness).
    
Based on MATLAB code from [Beilstein J. Nanotechnol. 3, 238 (2012)](https://www.beilstein-journals.org/bjnano/articles/3/27).
"""
function deconvolve_matrix(z::T, Δf::T, f₀::Float64, A::Float64, k::Float64)::Tuple{T,T} where T<:AbstractVector{Float64}
    @assert issorted(z)

    Δf_r = @view Δf[end:-1:1]
    Δz = z[2] - z[1]

    α = round(Int64, A / Δz)

    N = length(z)
    # W = Array{Float64}(undef, N, N)
    W = zeros(Float64, N, N)

    for i = 1:N
        x = max(i - 2α, 1)
        for j = x:i
            W[i,j] = (f₀ / 2 / k) * (2 / π / A) * 2 / (2α + 1) * (sqrt((2α + 1) * (i - j + 1) - (i - j + 1)^2)-sqrt((2α + 1) * (i - j) - (i - j)^2))

            # τ1 = 1 - 2(i - j + 1) / (2α + 1)
            # τ2 = 1 - 2(i - j) / (2α + 1)
            # W[i,j] = (f₀ / 2 / k) * (2 / π / A) * (-sqrt(1 - τ2^2) + sqrt(1 - τ1^2))
        end
    end

    F = W \ Δf_r
    reverse!(F)
    return z[1:end], F  # return copy of z
end




end
