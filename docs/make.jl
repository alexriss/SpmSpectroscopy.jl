using SpmSpectroscopy
using Documenter

DocMeta.setdocmeta!(SpmSpectroscopy, :DocTestSetup, :(using SpmSpectroscopy); recursive=true)

makedocs(;
    modules=[SpmSpectroscopy],
    authors="Alex Riss <00alexx@riss.at>",
    repo="https://github.com/alexriss/SpmSpectroscopy.jl/blob/{commit}{path}#{line}",
    sitename="SpmSpectroscopy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexriss.github.io/SpmSpectrpscopy.jl",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
		"Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/alexriss/SpmSpectroscopy.jl",
    devbranch="main",
)
