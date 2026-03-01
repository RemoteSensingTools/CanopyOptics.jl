# # Specular Canopy Scattering

# Using packages:
using CairoMakie, Distributions, Base64
using CanopyOptics

# Create a specular model with refractive index and "roughness" κ
specularMod = CanopyOptics.SpecularCanopyScattering(nᵣ=1.5, κ=0.2)

# Use standard planophile leaf distribution:
LD = CanopyOptics.planophile_leaves2()

# Create some discretized angles:
# Azimuth (in radians)
ϕ = range(0.0, 2π, length=200);
# Polar angle as cos(Θ)
μ, w = CanopyOptics.gauleg(180, 0.0, 1.0);
# Create directional vectors:
dirs = [CanopyOptics.dirVector_μ(a, b) for a in μ, b in ϕ];

# Compute specular reflectance for a fixed incoming direction:
R = CanopyOptics.compute_reflection.([specularMod], [dirs[50,1]], dirs, [LD]);

# ### Polar plot of reflectance using PolarAxis
# Zenith angle θ is the radial axis (0 = nadir, π/2 = horizon).
# acos.(μ) is descending (μ ascending); sort to ascending and reverse R rows to match.
# PolarAxis surface! expects matrix shape (n_angle × n_radius).
θ       = sort(acos.(μ))            # ascending, 0 → π/2
R_polar = reverse(R, dims=1)'       # (nϕ=200 × nθ=180)

fig = Figure(size=(520, 450))
ax  = PolarAxis(fig[1, 1], title="Reflectance R")
hm  = surface!(ax, collect(ϕ), θ, zeros(size(R_polar)),
               color=R_polar, colorrange=(0, 0.02), shading=NoShading)
rlims!(ax, 0, π/2)
Colorbar(fig[1, 2], hm)
fig

# ### Animation over different Beta leaf distributions
# Use a coarser angular grid for the animation to keep the GIF size small.
ϕ_a     = range(0.0, 2π, length=60)
μ_a, _  = CanopyOptics.gauleg(40, 0.0, 1.0)
dirs_a  = [CanopyOptics.dirVector_μ(a, b) for a in μ_a, b in ϕ_a]
θ_a     = sort(acos.(μ_a))

steps  = 0.1:0.2:5
α_vals = [collect(steps); 5 * ones(length(steps))]
β_vals = [5 * ones(length(steps)); reverse(collect(steps))]
x      = collect(0:0.01:1)

# Observables: updating LD_obs propagates to all derived plots reactively.
LD_obs  = Observable(LD)
pdf_obs = @lift pdf.($LD_obs.LD, x)
R_obs   = @lift reverse(
    CanopyOptics.compute_reflection.([specularMod], [dirs_a[20,1]], dirs_a, [$LD_obs]),
    dims=1)'
# Use map to properly capture m in each iteration (avoids closure-in-loop)
Zpm_obs = map(0:3) do m
    @lift CanopyOptics.compute_Z_matrices(specularMod, μ_a, $LD_obs, m)[2]
end

fig2 = Figure(size=(700, 700))
ax0  = Axis(fig2[1,1], title="Leaf angle distribution",
            xlabel="Θ (degrees)", limits=(0, 90, 0, 3))
ax1  = Axis(fig2[1,2], title="Reflectance R",
            xlabel="Azimuth ϕ (°)", ylabel="Zenith θ (°)")
ax_Z = [Axis(fig2[i+1, j], title="Z⁻⁺ m=$(2(i-1)+(j-1))",
             xlabel="μꜜ", ylabel="μꜛ") for i in 1:2, j in 1:2]

lines!(ax0, rad2deg.(π * x / 2), pdf_obs)
heatmap!(ax1, rad2deg.(collect(ϕ_a)), rad2deg.(θ_a), R_obs, colorrange=(0, 0.02))
for (k, ax) in enumerate(ax_Z)
    heatmap!(ax, μ_a, μ_a, Zpm_obs[k], colorrange=(0, 0.3))
end

path = record(fig2, "anim_specular.gif", eachindex(α_vals); framerate=10) do i
    LD_obs[] = CanopyOptics.LeafDistribution(Beta(α_vals[i], β_vals[i]), 2/π)
end
HTML("<img src=\"data:image/gif;base64,$(base64encode(read(path)))\" />")
