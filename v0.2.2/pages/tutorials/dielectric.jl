# # Microwave Dielectric Models

# CanopyOptics ships a small but extensible family of complex permittivity
# models for the materials that appear in canopy radiative transfer at
# microwave frequencies: water, ice, soil, and fresh leaves. They all share
# one entry point — the `dielectric(material, T_kelvin, f_GHz)` function —
# and dispatch is on the material type.
using CanopyOptics
using CairoMakie

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

fig = Figure(size = (640, 360))
ax = Axis(fig[1, 1]; xlabel = "frequency (GHz)", ylabel = "ε",
          title = "Permittivity vs frequency at T = 280 K")
lines!(ax, freqs, real.(ε_water); label = "water  ε'")
lines!(ax, freqs, imag.(ε_water); label = "water  ε''", linestyle = :dash)
lines!(ax, freqs, real.(ε_leaf);  label = "leaf   ε'")
lines!(ax, freqs, imag.(ε_leaf);  label = "leaf   ε''", linestyle = :dash)
axislegend(ax; position = :rt)
fig

# ## Moisture sweep — leaf
#
# At fixed C-band frequency the real and imaginary parts of the leaf
# permittivity rise monotonically with gravimetric moisture, dominated by
# the free-water contribution at high `M_g`.
mgs   = range(0.0, 0.7, length = 41)
ε_mg  = [dielectric(LeafUlabyElRayes1987(M_g = m), 295.0, 5.0) for m in mgs]

fig2 = Figure(size = (640, 360))
ax2 = Axis(fig2[1, 1]; xlabel = "gravimetric moisture M_g",
           ylabel = "ε", title = "Leaf permittivity at f = 5 GHz")
lines!(ax2, mgs, real.(ε_mg); label = "ε'")
lines!(ax2, mgs, imag.(ε_mg); label = "ε''", linestyle = :dash)
axislegend(ax2; position = :lt)
fig2

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
