<img width="100" height="100" src="docs/src/assets/logo.svg?raw=true" />

# SpmSpectroscopy.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexriss.github.io/SpmSpectroscopy.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexriss.github.io/SpmSpectroscopy.jl/dev)
[![Build Status](https://github.com/alexriss/SpmSpectroscopy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alexriss/SpmSpectroscopy.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alexriss/SpmSpectroscopy.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alexriss/SpmSpectroscopy.jl)

A julia library to analyze scanning tunneling and atomic force spectroscopy data.

*Currently in early stages, with support for Nanonis spectroscopy files. So far, it is only tested for Bias and Z spectroscopy experiments. Reading of other experiments and file formats will be implemented if needed.*

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

## Related projects

- [SpmImageTycoon.jl](https://github.com/alexriss/SpmImageTycoon.jl): App to organize SPM images and spectra.
- [SpmImages.jl](https://github.com/alexriss/SpmImages.jl): Julia library to read and display SPM images.
- [SpmGrids.jl](https://github.com/alexriss/SpmGrids.jl): Julia library to read and analyze SPM grid spectroscopy.
- [imag*ex*](https://github.com/alexriss/imagex): Python scripts to analyze scanning probe images.
- [grid*ex*](https://github.com/alexriss/gridex): Python scripts to analyze 3D grid data.
