# # Composing tree canopies: leaves + stems + branches
#
# CanopyOptics' `CanopyComponent` + `MixedCanopy` machinery already
# supports area-index-weighted compositions of multiple scattering
# populations. This tutorial walks through the ergonomic constructors
# `LeafComponent`, `StemComponent`, `BranchComponent`, and the
# one-liner wrapper `TreeCanopy(...)` that produces a `MixedCanopy`
# with sensible per-population defaults.
#
# The motivation: in real forests the scattering populations differ
# along three axes simultaneously — area index, angle distribution,
# and scattering model. Leaves have spherical-ish angle distributions
# and bi-Lambertian scattering (reflectance + transmittance). Trunks
# (stems) are nearly vertical (erectophile) and opaque (Lambertian
# reflection, no transmission). Branches sit between the two, both in
# tilt (plagiophile) and in spectral reflectance.

using CanopyOptics
using CanopyOptics: bulk_G, G, LeafComponent, StemComponent, BranchComponent,
                    TreeCanopy, spherical_leaves, erectophile_leaves,
                    plagiophile_leaves

# ## A leaf-only canopy (the default)

# Passing only `LAI` produces a single-component canopy equivalent to
# the existing `LeafComponent`-only construction.

leaves_only = TreeCanopy(LAI = 4.0)

#-

length(leaves_only.components)

# `bulk_G` over a μ grid then equals `LAI · G(μ, spherical)`:

μ_test = [0.1, 0.3, 0.5, 0.7, 1.0]
bulk_G_leaves = bulk_G(leaves_only, μ_test)

#-

analytical = 4.0 .* vec(G(μ_test, spherical_leaves()))
bulk_G_leaves ≈ analytical

# ## Leaves + stems

# Adding `SAI > 0` (stem area index, m² stem / m² ground) appends a
# `StemComponent` to the canopy.  Stems default to:
#
#   - opaque Lambertian scattering (no transmission)
#   - erectophile angle distribution (mostly vertical)
#   - broadband bark reflectance R = 0.25 (Bonan 2019, Table 14.1)

leaves_and_stems = TreeCanopy(LAI = 4.0, SAI = 0.9)

#-

length(leaves_and_stems.components)

# At low sun angles (small μ), the erectophile stems intercept more
# light per unit area than spherical leaves, so the canopy `bulk_G`
# grows. At overhead sun (μ = 1) stems are seen edge-on and their
# `G` is small.

bulk_G_stems = bulk_G(leaves_and_stems, μ_test)

#-

# Difference from leaves-only:

bulk_G_stems .- bulk_G_leaves

# ## Full tree canopy: leaves + stems + branches

# Branches use plagiophile (~45° tilt) leaves by default, with a
# slightly higher broadband reflectance (R = 0.30) than mature
# stems — younger / smoother bark.

full_canopy = TreeCanopy(LAI = 4.0, SAI = 0.9, BAI = 0.1)

#-

length(full_canopy.components)

# Each component carries its own scatterer, LAD, and area index.

[(typeof(c.scatterer).name.name, c.area_index) for c in full_canopy.components]

# ## Overriding defaults

# Each per-population constructor accepts overrides on every axis.
# `TreeCanopy` forwards `leaf_kwargs`, `stem_kwargs`, `branch_kwargs`
# NamedTuples to the underlying calls. For example, band-resolved
# bark reflectance via a lookup table:

bark_LUT = CanopyOptics.LUTWoodReflectance(
    grid      = [400.0, 800.0, 2500.0],
    R         = [0.10,  0.20,  0.45],
    grid_unit = :nm,
)

#-

canopy_LUT_bark = TreeCanopy(LAI = 4.0, SAI = 0.9,
                             stem_kwargs = (; R = bark_LUT))
stem_component  = canopy_LUT_bark.components[2]
CanopyOptics.wood_reflectance(stem_component.scatterer, 800.0)

# ## What this plugs into

# The returned `MixedCanopy` is a regular CanopyOptics scattering
# object: `bulk_G`, `compute_Z_matrices`, and the higher-level
# radiative-transfer paths in RRTMGP / ClimaLand / CanopyColumn all
# consume it directly. The convenience constructors only set defaults;
# the underlying machinery is unchanged.
