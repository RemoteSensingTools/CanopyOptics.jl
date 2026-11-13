# # Specular Leaf-Surface Reflection

# `SpecularCanopyScattering` models the directional Fresnel-like component of
# leaf surface reflection. This path is still quadrature-based; the closed-form
# Fourier stack currently applies only to the diffuse bi-Lambertian kernel.
using CanopyOptics
using Distributions
include(joinpath(@__DIR__, "..", "..", "..", "d3_plots.jl")) # hide

specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2)
quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64, n_azimuth = 48)
LD = CanopyOptics.planophile_leaves2()

# Compute a directional reflection value for one incoming and one outgoing
# direction. Directions are represented by signed `μ = cos(θ)` and azimuth
# angle `ϕ`.
incoming = CanopyOptics.DirectionCosine(0.8, 0.0)
outgoing = CanopyOptics.DirectionCosine(0.4, π / 3)
CanopyOptics.compute_reflection(specular, incoming, outgoing, LD)

# The specular model can also be projected onto Fourier Z matrices one moment
# at a time. Unlike the analytic bi-Lambertian routine, this uses numerical
# azimuth quadrature.
μ, w = CanopyOptics.gauleg(8, 0.0, 1.0)
Z⁺⁺₀, Z⁻⁺₀ = CanopyOptics.compute_Z_matrices(specular, μ, LD, 0; quadrature)
size(Z⁻⁺₀)

# Higher moments are available, but sharp specular lobes generally require
# more Fourier moments than diffuse bi-Lambertian scattering.
Z⁺⁺₃, Z⁻⁺₃ = CanopyOptics.compute_Z_matrices(specular, μ, LD, 3; quadrature)
maximum(abs.(Z⁻⁺₃))

# ## Adding diffuse and specular terms

# Leaf scattering components are additive at the Z-matrix level. The composite
# model dispatches each component through its own implementation: closed-form
# Fourier moments for the diffuse bi-Lambertian term and azimuth quadrature for
# the specular term.
diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2)
mixed_leaf = diffuse + specular
Z_mix⁺⁺, Z_mix⁻⁺ = CanopyOptics.compute_Z_matrices(mixed_leaf, μ, LD, 0:3; quadrature)
size(Z_mix⁻⁺)

# ## Interactive leaf-angle sweep

# Sweep `Beta(α, β)` leaf inclination distributions and observe how a fixed
# incoming beam (`μ = 0.8, ϕ = 0`) scatters specularly across outgoing
# `(ϕ, θ)` directions.
ϕ_deg = collect(range(0.0, 360.0, length = 49))
θ_deg = collect(range(2.0, 88.0, length = 29))
u_grid = collect(0.01:0.01:0.99)
θ_lad_deg = 90.0 .* u_grid
αβ = [(0.8, 5.0), (1.5, 4.0), (2.5, 2.5), (4.0, 1.5), (5.0, 0.8)]

D3Figure(
    title = "Specular Leaf-Surface Reflection Across Leaf-Angle Distributions",
    subtitle = "Incoming source direction (μ = 0.8, ϕ = 0°); nᵣ = 1.5, κ = 0.2",
    cols = 2,
    rowHeight = 310,
    frames = [
        begin
            ld_i = CanopyOptics.LeafDistribution(Beta(a, b), 2 / π)
            pdf_y = pdf.(Ref(ld_i.LD), u_grid)
            R_mat = [
                CanopyOptics.compute_reflection(
                    specular,
                    incoming,
                    CanopyOptics.DirectionCosine(cosd(θv), deg2rad(ϕv)),
                    ld_i,
                ) for θv in θ_deg, ϕv in ϕ_deg
            ]
            (;
                label = "Beta(α = $(a), β = $(b))",
                panels = [
                    d3_cartesian(
                        title = "Leaf inclination distribution",
                        xlabel = "inclination θ (degrees)",
                        ylabel = "density on [0, 1]",
                        xlim = [0.0, 90.0],
                        ylim = [0.0, 3.2],
                        series = [
                            d3_line(θ_lad_deg, pdf_y; label = "Beta($(a), $(b))", color = "#1f77b4", width = 2.3),
                        ],
                    ),
                    d3_heatmap(
                        ϕ_deg, θ_deg, R_mat;
                        title = "Specular reflection lobe",
                        xlabel = "azimuth ϕ (degrees)",
                        ylabel = "zenith θ (degrees)",
                        cbarLabel = "f_s",
                        clim = [0.0, 0.02],
                        colormap = "viridis",
                    ),
                ],
            )
        end for (a, b) in αβ
    ],
)
