pushfirst!(LOAD_PATH, dirname(@__DIR__))

using Documenter, CanopyOptics
using Literate, UnitfulEquivalences, Distributions

@info "Building docs with CanopyOptics source" path=pathof(CanopyOptics)

function build()
    tutorials = ["bilambertian.jl", "specular.jl", "dielectric.jl"]
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
        size_threshold = nothing,
        size_threshold_warn = nothing,
        example_size_threshold = nothing,
        )
    makedocs(
            sitename = "Canopy Optics",
            format = format,
            clean = true,
            checkdocs = :none,
            warnonly = [:cross_references],
            modules = [CanopyOptics],
            pages = pages)
end
build()

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/RemoteSensingTools/CanopyOptics.jl.git",
        target = "build",
        devbranch = "main",
        push_preview = true,
    )
end
