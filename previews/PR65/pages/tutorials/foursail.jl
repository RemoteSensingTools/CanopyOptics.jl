# # 4SAIL Bidirectional Canopy Reflectance

# `foursail` is the package's implementation of Verhoef's 4SAIL canopy
# bidirectional reflectance model. It returns a [`FourSAILResult`](@ref)
# carrying both the full BRF and the explicit decomposition into
# direct/direct and diffuse contributions, which makes it easy to
# diagnose what each pathway contributes at a given wavelength.
using CanopyOptics
using CairoMakie

# ## Geometry, leaf, soil
#
# Geometry and leaf-angle integration are wavelength-independent; build
# them once with [`FourSAILGeometry`](@ref) and reuse across wavelengths
# or moisture/structure sweeps.
LD   = spherical_leaves()
geom = FourSAILGeometry(LD; sza_deg = 30.0, vza_deg = 20.0, raa_deg = 40.0)

# Use a small synthetic spectrum (red, NIR, SWIR1) so the table fits in
# one screen. In practice you would feed in a full PROSPECT-PRO leaf.
λ      = ["red 670 nm", "NIR 800 nm", "SWIR 1.6 µm"]
leaf_R = [0.06, 0.45, 0.30]      # leaf reflectance ρ
leaf_T = [0.02, 0.45, 0.25]      # leaf transmittance τ
soil_R = [0.10, 0.18, 0.30]
lai    = 3.0

result = foursail(leaf_R, leaf_T, soil_R, geom, lai; hotspot = 0.05)

# `result.rsot` is the full bidirectional reflectance factor — what a
# real instrument would observe at the chosen sun/view geometry.
println("Full BRF (rsot):")
for (name, v) in zip(λ, result.rsot)
    println("  ", rpad(name, 12), round(v, digits = 4))
end

# ## The rsot = rsost + rsodt decomposition
#
# 4SAIL splits the bidirectional reflectance into two physically
# distinct branches:
#
# * `rsost` — the **direct/direct** branch: single scattering from leaves
#   *plus* directly transmitted soil reflectance reaching the sensor
#   without any diffuse scattering. Treat it as a lower bound on what the
#   sensor sees if you were to ignore multiple scattering.
# * `rsodt` — the **diffuse / multiple-scattering remainder** that has to
#   be added to `rsost` to recover the full BRF.
#
# By construction `rsot ≈ rsost + rsodt`. In the NIR band, where leaves
# are highly scattering (`ρ + τ ≈ 0.9`), most of the signal arrives
# through diffuse multiple scattering — `rsost` alone severely
# underestimates the canopy reflectance.
println("\nrsost (direct/direct) vs rsodt (diffuse) vs rsot (full):")
for (name, rost, rodt, rot) in zip(λ, result.rsost, result.rsodt, result.rsot)
    println("  ", rpad(name, 12),
            "  rsost=", round(rost, digits = 4),
            "  rsodt=", round(rodt, digits = 4),
            "  rsot=",  round(rot,  digits = 4))
end

# Confirm the decomposition holds to machine precision.
maximum(abs, result.rsot .- (result.rsost .+ result.rsodt))

# ## Wavelength sweep — visualizing the decomposition
#
# Sweep a synthetic ρ(λ), τ(λ) that interpolates from a red-band
# absorption regime through the NIR plateau. The gap between `rsost`
# (the "if everything were single-scattering" curve) and `rsot` is the
# diffuse contribution `rsodt`.
λs    = range(0.4, 2.5, length = 80)
ρ_λ   = @. 0.04 + 0.42 / (1 + exp(-(λs - 0.72) * 25))  # red dip → NIR plateau
τ_λ   = 0.9 .* ρ_λ                                     # τ tracks ρ in this toy example
soil_λ = fill(0.20, length(λs))

sweep = foursail(ρ_λ, τ_λ, soil_λ, geom, lai; hotspot = 0.05)

fig = Figure(size = (640, 360))
ax = Axis(fig[1, 1]; xlabel = "wavelength (µm)", ylabel = "reflectance",
          title = "4SAIL BRF and its decomposition  (LAI = 3, sza/vza/raa = 30/20/40°)")
lines!(ax, λs, sweep.rsost; label = "rsost (direct/direct)", linestyle = :dash)
lines!(ax, λs, sweep.rsodt; label = "rsodt (diffuse)",       linestyle = :dot)
lines!(ax, λs, sweep.rsot;  label = "rsot  (full BRF)",      linewidth = 2)
axislegend(ax; position = :rt)
fig

# ## Batched geometries — single-pass directional response
#
# For BRDF studies pass a [`FourSAILGeometrySet`](@ref) instead of a
# single geometry: 4SAIL evaluates all (wavelength × angle) combinations
# in one batched pass (CPU or GPU via KernelAbstractions, depending on
# the array backend).
vzas  = collect(0.0:5.0:60.0)
raas  = fill(40.0, length(vzas))
geoms = FourSAILGeometrySet(LD; sza_deg = 30.0, vza_deg = vzas, raa_deg = raas)
brdf  = foursail(leaf_R, leaf_T, soil_R, geoms, lai; hotspot = 0.05)

# `brdf.rsot` is now a (n_wavelength × n_angle) matrix.
size(brdf.rsot)
