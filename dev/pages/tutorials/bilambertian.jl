# # Bilambertian Canopy Scattering

# Using packages:
using CairoMakie, Distributions
using CanopyOptics

# Compute quadrature points:
μ, w = CanopyOptics.gauleg(20, 0.0, 1.0)

# Use standard planophile leaf distribution:
LD = CanopyOptics.planophile_leaves2()

# Use a Bilambertian model and specify Reflectance R and Transmission T
BiLambMod = CanopyOptics.BiLambertianCanopyScattering(R=0.4, T=0.2)

# Compute Scattering Matrices 𝐙⁺⁺, 𝐙⁻⁺ for the 0th Fourier moment:
𝐙⁺⁺, 𝐙⁻⁺ = CanopyOptics.compute_Z_matrices(BiLambMod, μ, LD, 0)

# Plots matrices (shows single scattering contributions for incoming and outgoing directions)
fig = Figure(size=(820, 350))
ax1 = Axis(fig[1,1], title="Z⁻⁺ (Reflection)",  xlabel="μꜜ", ylabel="μꜛ")
ax2 = Axis(fig[1,2], title="Z⁺⁺ (Transmission)", xlabel="μꜛ", ylabel="μꜛ")
hm1 = heatmap!(ax1, μ, μ, 𝐙⁻⁺)
hm2 = heatmap!(ax2, μ, μ, 𝐙⁺⁺)
Colorbar(fig[1,3], hm1)
Colorbar(fig[1,4], hm2)
fig

# ### Animation over different Beta leaf distributions
steps   = 0.1:0.2:5
α_vals  = [collect(steps); 5 * ones(length(steps))]
β_vals  = [5 * ones(length(steps)); reverse(collect(steps))]
x       = collect(0:0.01:1)

# Use Observables so the figure updates reactively each animation frame
LD_obs  = Observable(LD)
pdf_obs = @lift pdf.($LD_obs.LD, x)
Z_obs   = @lift CanopyOptics.compute_Z_matrices(BiLambMod, μ, $LD_obs, 0)
Zpp_obs = @lift $Z_obs[1]
Zpm_obs = @lift $Z_obs[2]

fig2 = Figure(size=(1100, 320))
ax0 = Axis(fig2[1,1], title="Leaf angle distribution",
           xlabel="Θ (degrees)", limits=(0, 90, 0, 3))
ax1 = Axis(fig2[1,2], title="Z⁻⁺ (Reflection)",  xlabel="μꜜ", ylabel="μꜛ")
ax2 = Axis(fig2[1,3], title="Z⁺⁺ (Transmission)", xlabel="μꜛ", ylabel="μꜛ")
lines!(ax0, rad2deg.(π * x / 2), pdf_obs)
heatmap!(ax1, μ, μ, Zpm_obs, colorrange=(0.5, 2))
heatmap!(ax2, μ, μ, Zpp_obs, colorrange=(0.5, 2))

record(fig2, "anim_fps10.gif", eachindex(α_vals); framerate=10) do i
    LD_obs[] = CanopyOptics.LeafDistribution(Beta(α_vals[i], β_vals[i]), 2/π)
end
fig2
