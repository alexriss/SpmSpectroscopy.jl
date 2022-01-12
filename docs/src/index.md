```@meta
CurrentModule = SpmSpectroscopy
```

# SpmSpectroscopy

Documentation for [SpmSpectroscopy](https://github.com/alexriss/SpmSpectroscopy.jl).

## About

A julia library to analyze scanning tunneling and atomic force spectroscopy data.

## Usage

```julia
using SpmSpectroscopy

s = load_spectrum("Bias_spectrocopy_007.dat")

s.position  # get position of the probe
s.channel_names  # get channel names
s.channel_units  # get channel names
s.header  # get raw header data

s.data  # all data in a DataFrame
s.data.Current  # get data for a "Current" channel
s.data[!, "Current"]  # get data for a "Current" channel
```

A more detailed description can be found in the [Reference](@ref)
