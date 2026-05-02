# # Specular Leaf-Surface Reflection

# `SpecularCanopyScattering` models the directional Fresnel-like component of
# leaf surface reflection. This path is still quadrature-based; the closed-form
# Fourier stack currently applies only to the diffuse bi-Lambertian kernel.
using CanopyOptics
using CairoMakie
using Distributions
using Base64

specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2, nQuad = 48)
LD = CanopyOptics.planophile_leaves2()

# Compute a directional reflection value for one incoming and one outgoing
# direction. Directions are represented by signed `μ = cos(θ)` and azimuth
# angle `ϕ`.
incoming = CanopyOptics.dirVector_μ(0.8, 0.0)
outgoing = CanopyOptics.dirVector_μ(0.4, π / 3)
CanopyOptics.compute_reflection(specular, incoming, outgoing, LD)

# The specular model can also be projected onto Fourier Z matrices one moment
# at a time. Unlike the analytic bi-Lambertian routine, this uses numerical
# azimuth quadrature and should be treated as a separate component until a
# mixed diffuse/specular leaf model is introduced.
μ, w = CanopyOptics.gauleg(8, 0.0, 1.0)
Z⁺⁺₀, Z⁻⁺₀ = CanopyOptics.compute_Z_matrices(specular, μ, LD, 0)
size(Z⁻⁺₀)

# Higher moments are available, but sharp specular lobes generally require
# more Fourier moments than diffuse bi-Lambertian scattering.
Z⁺⁺₃, Z⁻⁺₃ = CanopyOptics.compute_Z_matrices(specular, μ, LD, 3)
maximum(abs.(Z⁻⁺₃))

# ## Adding diffuse and specular terms

# Leaf scattering components are additive at the Z-matrix level. The composite
# model dispatches each component through its own implementation: closed-form
# Fourier moments for the diffuse bi-Lambertian term and azimuth quadrature for
# the specular term.
diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2, nQuad = 64)
mixed_leaf = diffuse + specular
Z_mix⁺⁺, Z_mix⁻⁺ = CanopyOptics.compute_Z_matrices_aniso(mixed_leaf, μ, LD, 3)
size(Z_mix⁻⁺)

# ## Lightweight animation

# The animation below keeps the original visual check: the same incoming
# direction is scattered into a grid of outgoing azimuth and zenith angles as
# the leaf-angle distribution changes.
ϕ_a = range(0.0, 2π, length = 48)
μ_a, _ = CanopyOptics.gauleg(28, 0.0, 1.0)
dirs_a = [CanopyOptics.dirVector_μ(a, b) for a in μ_a, b in ϕ_a]
θ_a = sort(acos.(μ_a))
x = collect(0:0.01:1)
αβ = [(0.8, 5.0), (1.5, 4.0), (2.5, 2.5), (4.0, 1.5), (5.0, 0.8)]

LD_obs = Observable(LD)
pdf_obs = @lift pdf.($LD_obs.LD, x)
R_obs = @lift reverse(
    CanopyOptics.compute_reflection.([specular], [incoming], dirs_a, [$LD_obs]),
    dims = 1)'

fig = Figure(size = (760, 340))
ax0 = Axis(fig[1, 1]; title = "Leaf angle distribution",
           xlabel = "inclination θ (degrees)", limits = (0, 90, 0, 3))
ax1 = Axis(fig[1, 2]; title = "Specular reflection",
           xlabel = "azimuth ϕ (degrees)", ylabel = "zenith θ (degrees)")

lines!(ax0, rad2deg.(π .* x ./ 2), pdf_obs)
hm = heatmap!(ax1, rad2deg.(collect(ϕ_a)), rad2deg.(θ_a), R_obs;
              colorrange = (0, 0.02), colormap = :viridis)
Colorbar(fig[1, 3], hm; label = "R")

path = record(fig, "anim_specular.gif", eachindex(αβ); framerate = 3) do i
    LD_obs[] = CanopyOptics.LeafDistribution(Beta(αβ[i]...), 2 / π)
end
HTML("<img src=\"data:image/gif;base64,$(base64encode(read(path)))\" />")
