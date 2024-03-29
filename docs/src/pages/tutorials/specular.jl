# # Specular Canopy Scattering

# Using packages:
using Plots, Distributions
using CanopyOptics
pyplot()
# theme(:ggplot2)

# Create a specular model with refractive index and "roughness" κ
specularMod = CanopyOptics.SpecularCanopyScattering(nᵣ=1.5, κ=0.2)

# Use standard planophile leaf distribution:
LD = CanopyOptics.planophile_leaves2()

# Create some discretized angles:
# Azimuth (in radians)
ϕ = range(0.0, 2π,  length=200);
# Polar angle as cos(Θ)
μ,w = CanopyOptics.gauleg(180,0.0,1.0);
# Create directional vectors:
dirs = [CanopyOptics.dirVector_μ(a,b) for a in μ, b in ϕ];

# Compute specular reflectance:
R = CanopyOptics.compute_reflection.([specularMod], [dirs[50,1]], dirs, [LD]);

# Polar plot of reflectance
contourf(ϕ, acos.(μ), R,proj=:polar, ylim=(0,π/2), alpha=0.8)

# ### Animation over different Beta leaf distributions
steps = 0.1:0.1:5 
α = [collect(steps); 5*ones(length(steps))] 
β = [ 5*ones(length(steps)); reverse(collect(steps));] 
x = 0:0.01:1
l = @layout [a  b; c d; e f]

anim = @animate for i ∈ eachindex(α)
    LD = CanopyOptics.LeafDistribution(Beta(α[i],β[i]), 2/π)
    R = CanopyOptics.compute_reflection.([specularMod], [dirs[120,1]], dirs, [LD])
    pZ = []
    for m=0:3
        Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(specularMod, μ, LD, m)
        push!(pZ, contourf(μ, μ, Z⁻⁺, title="Z⁻⁺ Fourier moment $m", xlabel="μꜜ", aspect_ratio = :equal, ylabel="μꜛ"))
    end
    p0 = plot(rad2deg.(π * x/2), pdf.(LD.LD,x), legend=false, ylim=(0,3), title="Leaf angle distribution", xlabel="Θ (degrees)")
    p1 = contourf(ϕ, acos.(μ), R, proj=:polar, ylim=(0,π/2), clims=(0,0.02), alpha=0.8,legend = :none, title="Polar plot R")
    plot(p0, p1, pZ..., layout = l, margin=5Plots.mm)
    plot!(size=(700,700))
end

gif(anim, "anim_fps10.gif", fps = 10)
