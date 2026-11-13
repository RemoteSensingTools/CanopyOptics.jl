# # Bi-Lambertian Canopy Scattering

# CanopyOptics represents the diffuse part of a leaf as a bi-Lambertian
# scatterer with leaf reflectance `R` and transmittance `T`.
using CanopyOptics
using Distributions
include(joinpath(@__DIR__, "..", "..", "..", "d3_plots.jl")) # hide

# Use a positive-μ quadrature grid and a standard planophile leaf-angle
# distribution.
μ, w = CanopyOptics.gauleg(12, 0.0, 1.0)
LD = CanopyOptics.planophile_leaves2()

# Build the leaf scattering model. The single-scattering albedo is
# `ϖ = R + T`; CanopyOptics divides that factor out of the returned Z matrices
# so downstream layer solvers can apply `ϖ` explicitly.
leaf = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2)
quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64)

# ## Closed-form Fourier stack

# Passing a moment range to `compute_Z_matrices` returns all cosine Fourier
# moments in one call. The array layout is `Z[i_out, j_in, m+1]` for the
# common `0:m_max` range.
m_max = 8
# Use the same leaf quadrature for extinction and scattering. Both integrate
# the same normalized leaf-area measure, including for PDFs with endpoint peaks.
# Do not normalize the outgoing Z columns to force energy closure: refine the
# outgoing angular grid to check its remaining integration error.
G_leaf = CanopyOptics.G(μ, LD; quadrature)
Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(leaf, μ, LD, 0:m_max; quadrature)
size(Z⁺⁺)

# `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is sign-change reflection. For
# conservative leaves, the `m = 0` quadrature-weighted column integral of
# `Z⁺⁺ + Z⁻⁺` is approximately 2 in this normalization.
column_flux = vec(sum(w .* (Z⁺⁺[:, :, 1] .+ Z⁻⁺[:, :, 1]), dims = 1))
extrema(column_flux)

# Single-moment calls slice the same analytic closure.
Z⁺⁺₂, Z⁻⁺₂ = CanopyOptics.compute_Z_matrices(leaf, μ, LD, 2; quadrature)
maximum(abs.(Z⁺⁺₂ .- Z⁺⁺[:, :, 3]))

# A single `m = 0` call slices the same analytic stack.
Z⁺⁺₀, Z⁻⁺₀ = CanopyOptics.compute_Z_matrices(leaf, μ, LD, 0; quadrature)
maximum(abs.(Z⁻⁺₀ .- Z⁻⁺[:, :, 1]))

# G now defaults to `CanopyQuadrature().n_leaf`, matching scattering. The legacy
# `nLeg` keyword remains an explicit resolution override for standalone callers;
# when coupling transport, passing the same `quadrature` is preferred.

# ## Changing the leaf-angle distribution

# Leaf-angle distributions can be supplied directly as distributions on
# `[0, 1]`, with the standard `2 / π` scaling from inclination angle to the
# normalized interval.
erectophile = CanopyOptics.LeafDistribution(Beta(5, 2), 2 / π)
Z_erect⁺⁺, Z_erect⁻⁺ =
    CanopyOptics.compute_Z_matrices(leaf, μ, erectophile, 0:m_max; quadrature)

maximum(abs.(Z_erect⁻⁺[:, :, 1] .- Z⁻⁺[:, :, 1]))

# ## Interactive leaf-angle sweep

# Sweep a family of `Beta(α, β)` leaf-angle distributions from planophile to
# erectophile and inspect how the `m = 0` canopy reflection matrix `Z⁻⁺` shifts
# with incoming and outgoing stream cosine `μ`.
function reflection_matrix_for_beta(a, b)
    ld = CanopyOptics.LeafDistribution(Beta(a, b), 2 / π)
    _, Z_mp = CanopyOptics.compute_Z_matrices(leaf, μ, ld, 0:0; quadrature)
    return permutedims(Z_mp[:, :, 1])
end

αβ = [(0.8, 5.0), (1.2, 4.2), (1.8, 3.2), (2.5, 2.5), (3.2, 1.8), (4.2, 1.2), (5.0, 0.8)]
Z_frames = [reflection_matrix_for_beta(a, b) for (a, b) in αβ]
z_hi = ceil(maximum(maximum, Z_frames); digits = 2)

D3Figure(
    title = "Bi-Lambertian Canopy Reflection Matrix Z⁻⁺ (m = 0)",
    subtitle = "Drag the slider or press Play to sweep Beta(α, β) from planophile to erectophile",
    cols = 1,
    rowHeight = 340,
    frames = [
        (;
            label = "Beta(α = $(a), β = $(b))",
            panels = [
                d3_heatmap(
                    μ, μ, Z_frames[i];
                    title = "Z⁻⁺(μ_out, μ_in), m = 0 — Beta($(a), $(b))",
                    xlabel = "outgoing μ",
                    ylabel = "incoming μ",
                    cbarLabel = "Z⁻⁺",
                    clim = [0.0, z_hi],
                    colormap = "viridis",
                ),
            ],
        ) for (i, (a, b)) in enumerate(αβ)
    ],
)
