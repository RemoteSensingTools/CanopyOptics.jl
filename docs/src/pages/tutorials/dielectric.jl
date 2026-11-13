# # Microwave Dielectric Models

# CanopyOptics ships a small but extensible family of complex permittivity
# models for the materials that appear in canopy radiative transfer at
# microwave frequencies: water, ice, soil, and fresh leaves. They all share
# one entry point — the `dielectric(material, T_kelvin, f_GHz)` function —
# and dispatch is on the material type.
using CanopyOptics
include(joinpath(@__DIR__, "..", "..", "..", "d3_plots.jl")) # hide

# ## The shared contract
#
# Every material type implements `dielectric(model, T, f)` and returns a
# `Complex` permittivity with the loss in the positive imaginary part
# (physics convention `e^{-iωt}`). Temperature is in kelvin, frequency in
# GHz.
materials = (
    "Pure water"   => LiquidPureWater(),
    "Salt water"   => LiquidSaltWater(S = 35.0),
    "Pure ice"     => PureIce(),
    "Moist soil"   => SoilMW(sand_frac = 0.4, clay_frac = 0.2, mᵥ = 0.30, ρ = 1.5),
    "Fresh leaf"   => LeafUlabyElRayes1987(M_g = 0.5),
)

T = 273.0       # K (0 °C — sits inside every model's validity range)
f = 1.4         # GHz (SMOS L-band)

for (name, mat) in materials
    println(rpad(name, 12), " ε(T=$T K, f=$f GHz) = ", dielectric(mat, T, f))
end

# Each method enforces its own validity domain — calls outside the
# documented range throw `ArgumentError` rather than silently returning
# unsupported values:
try
    dielectric(LeafUlabyElRayes1987(M_g = 0.5), 295.0, 100.0)
catch err
    println("rejected: ", err)
end

# ## Frequency sweep — fresh leaf vs liquid water
#
# A fresh leaf at half gravimetric moisture sits well below the bulk water
# response because only the free-water and bound-water volume fractions
# couple to the EM field; the rest is non-dispersive plant matter.
freqs = range(0.2, 20.0, length = 80)
ε_water = [dielectric(LiquidPureWater(),                280.0, f) for f in freqs]
ε_leaf  = [dielectric(LeafUlabyElRayes1987(M_g = 0.5), 280.0, f) for f in freqs]

D3Figure(
    cols = 1,
    rowHeight = 310,
    panels = [
        d3_cartesian(
            title = "Permittivity vs frequency at T = 280 K",
            xlabel = "frequency (GHz)",
            ylabel = "relative permittivity ε",
            series = [
                d3_line(freqs, real.(ε_water); label = "water ε'",  color = "#1f77b4", width = 2.2),
                d3_line(freqs, imag.(ε_water); label = "water ε''", color = "#1f77b4", dash = "dash", width = 2.0),
                d3_line(freqs, real.(ε_leaf);  label = "leaf ε'",   color = "#2ca02c", width = 2.2),
                d3_line(freqs, imag.(ε_leaf);  label = "leaf ε''",  color = "#2ca02c", dash = "dash", width = 2.0),
            ],
        ),
    ],
)

# ## Moisture sweep — leaf
#
# At fixed C-band frequency the real and imaginary parts of the leaf
# permittivity rise monotonically with gravimetric moisture, dominated by
# the free-water contribution at high `M_g`.
mgs   = range(0.0, 0.7, length = 41)
ε_mg  = [dielectric(LeafUlabyElRayes1987(M_g = m), 295.0, 5.0) for m in mgs]

D3Figure(
    cols = 1,
    rowHeight = 300,
    panels = [
        d3_cartesian(
            title = "Leaf permittivity vs moisture at f = 5 GHz, T = 295 K",
            xlabel = "gravimetric moisture M_g",
            ylabel = "relative permittivity ε",
            series = [
                d3_line(mgs, real.(ε_mg); label = "ε' (real)", color = "#1f77b4", width = 2.3),
                d3_line(mgs, imag.(ε_mg); label = "ε'' (loss)", color = "#d62728", dash = "dash", width = 2.1),
            ],
        ),
    ],
)

# ## Adding a new material
#
# Any new material is one struct + one method. The dispatch pattern keeps
# the call site identical:
#
# ```julia
# struct MyMaterial{FT} <: AbstractVegetation        # or AbstractWater / AbstractSoil
#     # … parameters …
# end
#
# function CanopyOptics.dielectric(mod::MyMaterial, T::Real, f::Real)
#     # … returns Complex{FT} …
# end
# ```
#
# Once defined, `MyMaterial` plugs into anything downstream that takes
# an `<: AbstractMaterial` — the planned microwave-canopy expansion will
# build on this same hierarchy.
