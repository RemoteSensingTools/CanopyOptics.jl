push!(LOAD_PATH,"src/")

using Documenter, CanopyOptics
using Literate, UnitfulEquivalences, Distributions
ENV["PYTHON"]=""

function build()
    tutorials = ["bilambertian.jl", "specular.jl"] # , 
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials"))
    end
    
    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
            "Home"                  => "index.md",
            "Tutorials"             => tutorials_md
        ]

    mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(),
            ),
        ))
    
    format = Documenter.HTML(
        #assets = [
        #    asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        #    ],
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = mathengine,
        collapselevel = 1,
        )
    makedocs(
            sitename = "Canopy Optics",
            format = format,
            clean = false,
            modules = [CanopyOptics],
            pages = pages)
end
build()

deploydocs(
    repo = "github.com/RemoteSensingTools/CanopyOptics.jl.git",
    target = "build",
    push_preview = true,
)
