push!(LOAD_PATH,"src/")

using Documenter, CanopyOptics, UnitfulEquivalences, Distributions

makedocs(
         sitename = "CanopyOptics.jl",
         modules  = [CanopyOptics],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/RemoteSensingTools/CanopyOptics.jl",
)