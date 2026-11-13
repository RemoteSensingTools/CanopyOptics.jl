# API Reference

This section documents all exported and `public` types, constants, and functions in
`CanopyOptics`, grouped by physical domain rather than alphabetical order.

```@docs
CanopyOptics
```

## Subsystem Map

| Physical Domain | API Reference Page | Guides & Evaluations |
| --- | --- | --- |
| Leaf optical tables, PROSPECT / Multi-Pigment / Leaf-SIP / dorsiventral budgets, smooth & microfacet dielectric surfaces, wax dispersion, and coupled `4×4` polarized leaf BSDF | [Leaf Optics & BSDF](leaf_optics.md) | [PROSPECT budgets](../pages/prospect_budget.md), [Leaf surfaces and BSDF](../pages/leaf_surfaces_and_bsdf.md), [Dorsiventral and multi-pigment](../pages/dorsiventral_and_multi_pigment.md), [Model evaluations](../pages/leaf_model_evaluations.md) |
| Constituent absorption partitioning, antenna selection, PSII/PSI excitation, 1-zone & 2-zone internal photon escape, and conditional emission spectra | [Excitation & Leaf Fluorescence](fluorescence_and_excitation.md) | [Pigment excitation](../pages/pigment_excitation.md), [Internal leaf emission](../pages/plate_emission.md), [Emission spectra](../pages/emission_spectrum.md) |
| Spherical directions, leaf angle distributions, Ross `G` functions, clumping, hotspot, wood/tree components, Fourier `Z`-matrices, and conservative angular normalization | [Canopy Structure & Scattering](canopy_scattering.md) | [Angular normalization](../pages/angular_normalization.md) |
| `FourSAIL` canopy radiative transfer, microwave dielectric permittivity `ε(T, f)`, and deprecated aliases | [FourSAIL & Microwave Dielectrics](foursail_and_dielectrics.md) | Tutorials in sidebar |

## Conventions

- **Exported vs. `public` bindings**: Symbols shown unqualified (`prospect_budget`, `compute_Z_matrices`) are exported by `using CanopyOptics`. Lower-level optical and geometric primitives declared `public` on Julia 1.11+ (`CanopyOptics.calctav`, `CanopyOptics.Fᵣ`, `CanopyOptics.H`, `CanopyOptics.gauleg`, etc.) are part of the supported API when qualified with `CanopyOptics.` or explicitly imported.
- **Precision & automatic differentiation**: Leaf budgets, surface optics, and canopy scattering kernels preserve caller floating-point precision (`Float32`, `Float64`) and propagate `ForwardDiff.Dual` tangents with respect to optical and trait parameters.

## Section Contents

```@contents
Pages = [
    "leaf_optics.md",
    "fluorescence_and_excitation.md",
    "canopy_scattering.md",
    "foursail_and_dielectrics.md",
]
Depth = 2
```

## Alphabetical Symbol Index

```@index
Pages = [
    "index.md",
    "leaf_optics.md",
    "fluorescence_and_excitation.md",
    "canopy_scattering.md",
    "foursail_and_dielectrics.md",
]
```
