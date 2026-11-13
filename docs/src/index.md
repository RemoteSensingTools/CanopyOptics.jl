# CanopyOptics.jl
*A package to compute canopy scattering properties*
## Package Features
- Use leaf angle distributions to compute bi-Lambertian and dorsiventral canopy scattering matrices
- Compute closed-form cosine Fourier moments for symmetric and dorsiventral bi-Lambertian canopy kernels
- Choose an opt-in [conservative discrete angular measure](pages/angular_normalization.md) shared by bi-Lambertian scattering and leaf emission geometry
- Couple dielectric leaf surfaces (Stern–Allen smooth cones, Beckmann microfacets, epicuticular wax dispersion) with internal mesophyll budgets into an energy-conserving, Helmholtz-reciprocal 4×4 leaf BSDF
- Compose diffuse and specular leaf-surface scattering models additively or via coupled BSDF quadrature
- Request scalar or Stokes-expanded canopy Z matrices with `npol`
- Combine leaf and wood canopy components with separate angle distributions
- Differentiate leaf budgets and canopy Z matrices with `ForwardDiff`
- Compute leaf reflectance, transmittance, and material absorptance with PROSPECT-PRO, separated chlorophyll *a*/*b* multi-pigment bases, Leaf-SIP, and continuous `N_tot >= 1` dorsiventral stacks
- Compute single-zone and two-zone dorsiventral internal fluorescence escape and photosystem excitation partitioning on CPU and CUDA
- Compute complex microwave dielectric ε(T, f) of water, ice, soil and leaves

## Z-matrix convention

Canopy scattering matrices use `Z[i_out, j_in]`: rows are outgoing streams and
columns are incoming streams. `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is
sign-change reflection. The leaf single-scattering albedo `ϖ = R + T` is not
folded into `Z`; layer solvers multiply by `ϖ` separately. For conservative
bi-Lambertian leaves, the `m = 0` column integral of `Z⁺⁺ + Z⁻⁺` is therefore
approximately `2` under the quadrature weights.

Specular components include their Fresnel/roughness strength in the returned Z
contribution. If a downstream RT solver multiplies canopy Z by a separate
single-scattering albedo, it must provide a matching effective albedo when
using specular or diffuse+specular composite leaves.

## Microwave dielectric models

A growing set of materials implements the [`dielectric`](@ref) function,
returning the complex relative permittivity `ε(T, f)` (loss as positive
imaginary part):

| Material | Type | Domain |
| --- | --- | --- |
| Pure liquid water | [`LiquidPureWater`](@ref) | 0.2–1000 GHz, 265–310 K |
| Salt water | [`LiquidSaltWater`](@ref) | + salinity 0–45 PSU |
| Pure ice | [`PureIce`](@ref) | 0.2–1000 GHz, 233–273 K |
| Moist soil | [`SoilMW`](@ref) | Dobson model, 0.3–18 GHz |
| Fresh leaf | [`LeafUlabyElRayes1987`](@ref) | 0.2–20 GHz, gravimetric moisture 0–0.7 |

All five share the contract `dielectric(model, T_kelvin, f_GHz)`. New
materials are added by extending the `dielectric` function on a fresh
`<: AbstractMaterial` subtype — see `src/utils/dielectric.jl` for the
existing methods. The vegetation hierarchy (`AbstractVegetation`) is
the entry point for a planned microwave-canopy expansion.

## Installation

The latest release of CanopyOptics can be installed from the Julia REPL prompt with

```julia
julia> ]add https://github.com/RemoteSensingTools/CanopyOptics.jl
```

## Documentation & API Reference

| Topic | Conceptual Guides | Structured API Reference |
| --- | --- | --- |
| **Leaf Optics & BSDF** | [PROSPECT budgets](pages/prospect_budget.md), [Leaf surfaces and BSDF](pages/leaf_surfaces_and_bsdf.md), [Dorsiventral and multi-pigment](pages/dorsiventral_and_multi_pigment.md), [Model evaluations](pages/leaf_model_evaluations.md) | [Leaf Optics & BSDF](api/leaf_optics.md) |
| **Excitation & Leaf Fluorescence** | [Pigment excitation](pages/pigment_excitation.md), [Internal leaf emission](pages/plate_emission.md), [Emission spectra](pages/emission_spectrum.md) | [Excitation & Leaf Fluorescence](api/fluorescence_and_excitation.md) |
| **Canopy Structure & Scattering** | [Angular normalization](pages/angular_normalization.md), Tutorials | [Canopy Structure & Scattering](api/canopy_scattering.md) |
| **FourSAIL & Microwave Dielectrics** | FourSAIL & Dielectric Tutorials | [FourSAIL & Microwave Dielectrics](api/foursail_and_dielectrics.md) |
| **Full Module & Symbol Index** | — | [Overview & Symbol Index](api/index.md) |
