# Canopy Structure & Scattering

Canopy radiative-transfer geometry in `CanopyOptics` combines continuous leaf-angle
distributions (`AbstractLeafDistribution`), clumping and hotspot gap corrections,
multi-component leaf/wood mixtures (`MixedCanopy`), and cosine Fourier scattering
matrices `Z[i_out, j_in]` across bi-Lambertian, dorsiventral, specular, and coupled
`4×4` leaf-BSDF scattering laws.

## `Z`-Matrix Single-Scattering Albedo Conventions

- **Symmetric bi-Lambertian (`BiLambertianCanopyScattering`, `LambertianWoodCanopyScattering`)**: Scalar single-scattering albedo `ϖ = R + T` is divided out of `Z` so downstream layer solvers multiply by `ϖ` separately; every non-black (`R + T > 0`) leaf satisfies `(1/2) ∑_i w_i (Z⁺⁺_{ij,0} + Z⁻⁺_{ij,0}) ≃ 1` (`= 1` under [`DiscreteAngularNormalization`](@ref)).
- **Dorsiventral & coupled BSDF (`compute_Z_matrices_directional`, `compute_Z_matrices_coupled_bsdf`)**: By default (`normalize_by_ω = false`), `Z` is normalized by `G(μ_in)` only and retains directional albedo inside `Z` (`(1/2) ∑_i w_i (Z⁺⁺_{ij,0} + Z⁻⁺_{ij,0}) ≃ ω_down[j]`), composing directly without `ϖ`-mixing ambiguities.
- **Pure surface specular (`SpecularCanopyScattering`)**: Retains Fresnel/microfacet reflectance strength directly inside `Z` in the right-handed `(e_s, e_p, k)` meridian Stokes basis.

```@contents
Pages = ["canopy_scattering.md"]
Depth = 2
```

## Directions & Angular Quadrature

```@docs
DirectionCosine
direction_from_angles
CanopyOptics.dot_product
CanopyQuadrature
CanopyOptics.gauleg
```

## Leaf Angle Distributions & Projection (`G`) Functions

```@docs
AbstractLeafDistribution
LeafDistribution
planophile_leaves
planophile_leaves2
erectophile_leaves
plagiophile_leaves
uniform_leaves
spherical_leaves
flat_leaves
beta_leaves
βparameters
G
G2
bfG
CanopyOptics.H
CanopyOptics.compute_Ψ
```

## Clumping & Hotspot Gap Geometry

```@docs
AbstractClumping
NoClumping
ConstantClumping
EmpiricalDirectionalClumping
ChenLeblancClumping
clumping_index
effective_G
AbstractHotSpot
NoHotSpot
KuuskHotSpot
canopy_extinction
hotspot_separation
hotspot_correction
joint_gap_probability
```

## Multi-Component Canopies & Wood Reflectance

```@docs
CanopyComponent
MixedCanopy
LeafComponent
StemComponent
BranchComponent
TreeCanopy
component_G
bulk_G
AbstractWoodReflectance
ConstantWoodReflectance
AbstractLUTWoodReflectance
LUTWoodReflectance
PolynomialWoodReflectance
wood_reflectance
```

## Canopy Scattering Models & Fourier `Z`-Matrices

```@docs
AbstractCanopyScatteringType
BiLambertianCanopyScattering
DorsiventralBiLambertianCanopyScattering
SpecularCanopyScattering
CompositeCanopyScattering
LambertianWoodCanopyScattering
LambertianWood
compute_Z_matrices
compute_Z_matrices_aniso_analytic
CanopyOptics.compute_Z_matrices_aniso
compute_Z_matrices_directional
compute_Z_matrices_coupled_bsdf
CanopyOptics.compute_Γ
CanopyOptics.compute_Γ_isotropic
compute_Γ_directional
CanopyOptics.compute_lambertian_Γ
compute_reflection
compute_reflection_mueller
CanopyOptics.getSpecularΩ
CanopyOptics.A
CanopyOptics.Asm
```

## Discrete Angular Normalization & Leaf Source Moments

```@docs
AbstractAngularNormalization
ContinuousAngularNormalization
DiscreteAngularNormalization
leaf_angular_scale
leaf_projection_moments
```
