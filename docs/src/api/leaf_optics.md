# Leaf Optics, Surfaces & Coupled BSDF

Leaf-level elastic radiative transfer in `CanopyOptics` separates tabulated spectral
optical constants from leaf structural/biochemical traits, dielectric surface
boundaries, and directional BSDF assembly:

1. **Optical tables** ([`LeafOpticalProperties`](@ref), [`MultiPigmentOpticalProperties`](@ref), [`load_wax_refractive_index`](@ref)) resample refractive indices and specific absorption spectra onto a caller-chosen wavelength grid.
2. **Hemispherical mesophyll & dorsiventral budgets** ([`prospect_budget`](@ref), [`leaf_sip_budget`](@ref), [`dorsiventral_budget`](@ref), [`single_slab_dorsiventral_budget`](@ref)) compute directional-hemispherical or bi-hemispherical reflectance ``R``, transmittance ``T``, and directly integrated nonnegative material absorptance ``A`` (``R + T + A = 1``).
3. **Dielectric boundary models** ([`ConeSurface`](@ref), [`RoughDielectricSurface`](@ref)) evaluate smooth Stern–Allen cone transmittance or Beckmann microfacet reflection/transmission with energy-conserving multiple-scattering compensation.
4. **Coupled ``4\times 4`` polarized leaf BSDF** ([`CoupledLeafBSDF`](@ref), [`prepare_leaf_bsdf`](@ref), [`leaf_bsdf_spectrum`](@ref)) couples specular surface reflection (``4\times 4`` Fresnel Mueller matrix) with interface-transmitted diffuse interior transport while preserving Helmholtz reciprocity and directional energy conservation ``R_{\mathrm{dir}}(\mu_{\mathrm{in}}) + T_{\mathrm{dir}}(\mu_{\mathrm{in}}) + A_{\mathrm{dir}}(\mu_{\mathrm{in}}) = 1``.

```@contents
Pages = ["leaf_optics.md"]
Depth = 2
```

## Spectral Optical Tables & Wax Dispersion

```@docs
LeafOpticalProperties
MultiPigmentOpticalProperties
load_wax_refractive_index
```

## Mesophyll Plate & Leaf-SIP Budgets

```@docs
LeafProspectProProperties
prospect
prospect_budget
prospect_constituent_coefficients
constituent_optical_depths
LeafMultiPigmentProperties
LeafSIPProperties
leaf_sip_budget
CanopyOptics.prospect_primal
```

## Dorsiventral & Asymmetric-Surface Leaf Stacks

```@docs
DorsiventralLeafProperties
AsymmetricSurfaceLeafProperties
dorsiventral_from_traits
dorsiventral_budget
single_slab_dorsiventral_budget
```

## Dielectric Surface Interfaces & Microfacet Optics

```@docs
AbstractLeafSurface
ConeSurface
RoughDielectricSurface
surface_interface_optics
calctav_smooth
CanopyOptics.calctav
CanopyOptics.Fᵣ
CanopyOptics.fresnel_components
CanopyOptics.K
beckmann_D
smith_beckmann_G1
microfacet_directional_reflectance
microfacet_hemispherical_reflectance
```

## Coupled `4×4` Polarized Leaf BSDF

```@docs
CoupledLeafBSDF
PreparedLeafBSDFSpectrum
prepare_leaf_bsdf
surface_brdf_mueller
leaf_directional_budget
leaf_bsdf_spectrum
```
