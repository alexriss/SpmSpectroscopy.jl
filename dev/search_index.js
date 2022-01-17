var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [SpmSpectroscopy]\nPrivate = false","category":"page"},{"location":"reference/#SpmSpectroscopy.correct_background!","page":"Reference","title":"SpmSpectroscopy.correct_background!","text":"function correct_background!(xdata<:Vector{<:Real}, ydata<:Vector{<:AbstractFloat}, type::Background, offset::Bool=true)::Nothing\n\nBackground correction of ydata vs. xdata with using a correction of type type. If offset is true (default), then ydata will be shifted such that its minimum is 0..\n\n\n\n\n\n","category":"function"},{"location":"reference/#SpmSpectroscopy.load_spectrum-Tuple{AbstractString}","page":"Reference","title":"SpmSpectroscopy.load_spectrum","text":"function load_spectrum(filename::AbstractString; select::AbstractVector=Bool[], header_only::Bool=false, index_column::Bool=false, index_column_type::Type=Int64)::SpmSpectrum\n\nLoads a spectrum from the file filename. Currently, only Nanonis .dat files are supported. select can be used to specify which columns to load (see CSV.jl for an explanation of select). If header_only is true, then only the header is loaded. If index_column is true, then an extra column with indices of type index_column_type will be added.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = SpmSpectroscopy","category":"page"},{"location":"#SpmSpectroscopy","page":"Introduction","title":"SpmSpectroscopy","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Documentation for SpmSpectroscopy.","category":"page"},{"location":"#About","page":"Introduction","title":"About","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"A julia library to analyze scanning tunneling and atomic force spectroscopy data.","category":"page"},{"location":"#Usage","page":"Introduction","title":"Usage","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"using SpmSpectroscopy\n\ns = load_spectrum(\"Bias_spectrocopy_007.dat\")\n\ns.position  # get position of the probe\ns.channel_names  # get channel names\ns.channel_units  # get channel names\ns.header  # get raw header data\n\ns.data  # all data in a DataFrame\ns.data.Current  # get data for a \"Current\" channel\ns.data[!, \"Current\"]  # get data for a \"Current\" channel","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A more detailed description can be found in the Reference","category":"page"}]
}
