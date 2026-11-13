using Documenter, CanopyOptics
using Literate, UnitfulEquivalences, Distributions

@info "Building docs with CanopyOptics source" path=pathof(CanopyOptics)

function build()
    tutorials = ["bilambertian.jl", "specular.jl", "foursail.jl",
                 "foursail_hotspot.jl", "foursail_brdf.jl", "dielectric.jl",
                 "tree_canopy.jl"]
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials"))
    end
    
    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
            "Home" => "index.md",
            "Leaf Optics & Fluorescence" => [
                "PROSPECT budgets"               => "pages/prospect_budget.md",
                "Leaf surfaces and BSDF"         => "pages/leaf_surfaces_and_bsdf.md",
                "Dorsiventral and multi-pigment" => "pages/dorsiventral_and_multi_pigment.md",
                "Pigment excitation"             => "pages/pigment_excitation.md",
                "Internal leaf emission"         => "pages/plate_emission.md",
                "Emission spectra"               => "pages/emission_spectrum.md",
                "Model evaluations"              => "pages/leaf_model_evaluations.md",
            ],
            "Canopy Radiative Transfer" => [
                "Angular normalization"          => "pages/angular_normalization.md",
            ],
            "Tutorials" => tutorials_md,
            "API Reference" => [
                "Overview & Index"                 => "api/index.md",
                "Leaf Optics & BSDF"               => "api/leaf_optics.md",
                "Excitation & Leaf Fluorescence"   => "api/fluorescence_and_excitation.md",
                "Canopy Structure & Scattering"    => "api/canopy_scattering.md",
                "FourSAIL & Microwave Dielectrics" => "api/foursail_and_dielectrics.md",
            ],
        ]

    format = Documenter.HTML(
        assets = [
            "assets/canopy_plots.css",
            "assets/canopy_plots.js",
        ],
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = Documenter.KaTeX(),
        collapselevel = 1,
        size_threshold = nothing,
        size_threshold_warn = nothing,
        example_size_threshold = nothing,
    )
    DocMeta.setdocmeta!(CanopyOptics, :DocTestSetup, :(using CanopyOptics); recursive = true)
    makedocs(
            sitename = "Canopy Optics",
            format = format,
            clean = true,
            checkdocs = :exports,
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
