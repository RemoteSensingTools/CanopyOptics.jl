```@meta
EditURL = "bilambertian.jl"
```

# Bilambertian Canopy Scattering

Using packages:

````@example bilambertian
using Plots, Distributions
using CanopyOptics
theme(:ggplot2)
````

Compute quadrature points:

````@example bilambertian
μ,w = CanopyOptics.gauleg(20,0.0,1.0)
````

Use standard planophile leaf distribution:

````@example bilambertian
LD = CanopyOptics.planophile_leaves2()
````

Use a Bilambertian model and specify Reflectance R and Transmission T

````@example bilambertian
BiLambMod = CanopyOptics.BiLambertianCanopyScattering(R=0.4,T=0.2)
````

Compute Scattering Matrices 𝐙⁺⁺, 𝐙⁻⁺ for 0th Fourier moment:

````@example bilambertian
𝐙⁺⁺, 𝐙⁻⁺ = CanopyOptics.compute_Z_matrices(BiLambMod, μ, LD, 0)
````

Plots matrices (shows single scattering contributions for incoming and outgoing directions)

````@example bilambertian
l = @layout [a  b]
p1 = contourf(μ, μ, 𝐙⁻⁺, title="Z⁻⁺ (Reflection)", xlabel="μꜜ", ylabel="μꜛ")
p2 = contourf(μ, μ, 𝐙⁺⁺, title="Z⁺⁺ (Transmission)", xlabel="μꜛ", ylabel="μꜛ")
plot(p1, p2,  layout = l, margin=5Plots.mm)
plot!(size=(900,350))
````

### Animation over different Beta leaf distributions

````@example bilambertian
steps = 0.1:0.1:5
α = [collect(steps); 5*ones(length(steps))]
β = [ 5*ones(length(steps)); reverse(collect(steps));]
x = 0:0.01:1

anim = @animate for i ∈ eachindex(α)
    LD = CanopyOptics.LeafDistribution(Beta(α[i],β[i]), 2/π)
    Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(BiLambMod, μ, LD, 0)
    l = @layout [a  b  c]
    p0 = plot(rad2deg.(π * x/2), pdf.(LD.LD,x), legend=false, ylim=(0,3), title="Leaf angle distribution", xlabel="Θ (degrees)")
    p1 = contourf(μ, μ, Z⁻⁺, title="Z⁻⁺ (Reflection)", xlabel="μꜜ", ylabel="μꜛ",clims=(0.5,2))
    p2 = contourf(μ, μ, Z⁺⁺, title="Z⁺⁺ (Transmission)", xlabel="μꜛ", ylabel="μꜛ",clims=(0.5,2))
    plot(p0, p1, p2,  layout = l, margin=5Plots.mm)
    plot!(size=(1100,300))
end

gif(anim, "anim_fps10.gif", fps = 10)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

