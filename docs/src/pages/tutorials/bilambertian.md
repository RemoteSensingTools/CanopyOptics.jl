```@meta
EditURL = "bilambertian.jl"
```

# Bi-Lambertian Canopy Scattering

CanopyOptics represents the diffuse part of a leaf as a bi-Lambertian
scatterer with leaf reflectance `R` and transmittance `T`.

````@example bilambertian
using CanopyOptics
using CairoMakie
using Distributions
using Base64
````

Use a positive-μ quadrature grid and a standard planophile leaf-angle
distribution.

````@example bilambertian
μ, w = CanopyOptics.gauleg(12, 0.0, 1.0)
LD = CanopyOptics.planophile_leaves2()
````

Build the leaf scattering model. The single-scattering albedo is
`ϖ = R + T`; CanopyOptics divides that factor out of the returned Z matrices
so downstream layer solvers can apply `ϖ` explicitly.

````@example bilambertian
leaf = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2, nQuad = 64)
````

## Closed-form Fourier stack

`compute_Z_matrices_aniso_analytic` returns all cosine Fourier moments from
`m = 0:m_max` in one call. The array layout is `Z[i_out, j_in, m+1]`.

````@example bilambertian
m_max = 8
Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices_aniso_analytic(leaf, μ, LD, m_max)
size(Z⁺⁺)
````

`Z⁺⁺` is same-sign transmission and `Z⁻⁺` is sign-change reflection. For
conservative leaves, the `m = 0` quadrature-weighted column integral of
`Z⁺⁺ + Z⁻⁺` is approximately 2 in this normalization.

````@example bilambertian
column_flux = vec(sum(w .* (Z⁺⁺[:, :, 1] .+ Z⁻⁺[:, :, 1]), dims = 1))
extrema(column_flux)
````

The legacy single-moment entry point delegates to the same analytic closure.

````@example bilambertian
Z⁺⁺₂, Z⁻⁺₂ = CanopyOptics.compute_Z_matrices_aniso(leaf, μ, LD, 2)
maximum(abs.(Z⁺⁺₂ .- Z⁺⁺[:, :, 3]))
````

The older `compute_Z_matrices` function is the azimuthally averaged `m = 0`
Shultis-Myneni assembly. It is retained for compatibility and should match
the analytic stack's first moment.

````@example bilambertian
Z⁺⁺₀, Z⁻⁺₀ = CanopyOptics.compute_Z_matrices(leaf, μ, LD, 0)
maximum(abs.(Z⁻⁺₀ .- Z⁻⁺[:, :, 1]))
````

## Changing the leaf-angle distribution

Leaf-angle distributions can be supplied directly as distributions on
`[0, 1]`, with the standard `2 / π` scaling from inclination angle to the
normalized interval.

````@example bilambertian
erectophile = CanopyOptics.LeafDistribution(Beta(5, 2), 2 / π)
Z_erect⁺⁺, Z_erect⁻⁺ =
    CanopyOptics.compute_Z_matrices_aniso_analytic(leaf, μ, erectophile, m_max)

maximum(abs.(Z_erect⁻⁺[:, :, 1] .- Z⁻⁺[:, :, 1]))
````

## Lightweight animation

CairoMakie renders without a Python backend. This animation sweeps a family
of beta leaf-angle distributions and updates the `m = 0` reflection matrix.

````@example bilambertian
function reflection_matrix_for_beta(a, b)
    LD = CanopyOptics.LeafDistribution(Beta(a, b), 2 / π)
    _, Z⁻⁺ = CanopyOptics.compute_Z_matrices_aniso_analytic(leaf, μ, LD, 0)
    return Z⁻⁺[:, :, 1]
end

αβ = [(0.8, 5.0), (1.5, 4.0), (2.5, 2.5), (4.0, 1.5), (5.0, 0.8)]
Zobs = Observable(reflection_matrix_for_beta(αβ[1]...))

fig = Figure(size = (560, 420))
ax = Axis(fig[1, 1]; xlabel = "outgoing μ", ylabel = "incoming μ",
          title = "Z⁻⁺, m = 0")
hm = heatmap!(ax, μ, μ, Zobs; colormap = :viridis)
Colorbar(fig[1, 2], hm; label = "Z⁻⁺")

path = record(fig, "anim_bilambertian.gif", eachindex(αβ); framerate = 3) do i
    Zobs[] = reflection_matrix_for_beta(αβ[i]...)
    ax.title = "Z⁻⁺, m = 0, Beta$(αβ[i])"
end
HTML("<img src=\"data:image/gif;base64,$(base64encode(read(path)))\" />")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

