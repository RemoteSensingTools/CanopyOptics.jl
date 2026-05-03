# CanopyOptics.jl
<p align="center">
  
  <a href="https://RemoteSensingTools.github.io/CanopyOptics.jl/dev/">
    <img src="https://img.shields.io/badge/docs-latest-blue.svg"
         alt="Docs">
  </a>
  <a href="https://github.com/RemoteSensingTools/CanopyOptics.jl/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/RemoteSensingTools/CanopyOptics.jl"
         alt="License">
  </a>
  <a href="https://github.com/RemoteSensingTools/CanopyOptics.jl/commits/main">
    <img src="https://img.shields.io/github/commit-activity/y/RemoteSensingTools/CanopyOptics.jl"
         alt="Github Commit Frequency">
  </a>
</p>

CanopyOptics.jl computes canopy radiative-transfer building blocks —
leaf and wood optical properties, leaf-angle distributions,
clumping/hotspot corrections, 4SAIL bidirectional reflectance, and
Fourier-decomposed canopy Z matrices — for downstream radiative-transfer
solvers. It is the canopy companion to
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl) and
shares its Z-matrix / Stokes / Fourier conventions so canopy outputs
plug directly into vSmartMOM's layer assembly. Today the package
covers the solar/SWIR spectrum (PROSPECT-PRO leaves, 4SAIL canopy);
microwave dielectric models for water, ice, soil and leaves are in
place as the seed of a planned MW canopy expansion.

## Installation

```julia
julia> ]add CanopyOptics
```

To track the development branch:

```julia
julia> ]add https://github.com/RemoteSensingTools/CanopyOptics.jl
```

## Quickstart

```julia
using CanopyOptics

# Leaf reflectance / transmittance from PROSPECT-PRO
opti = createLeafOpticalStruct(400.0:1.0:2500.0)             # nm grid
leaf = LeafProspectProProperties{Float64}(Ccab = 40.0)
T, R = prospect(leaf, opti)                                  # (transmittance, reflectance)

# Bi-Lambertian canopy Z matrices on a Gauss-Legendre μ grid
LD       = spherical_leaves()
μ, w     = CanopyOptics.gauleg(10, 0.0, 1.0)
canopy   = BiLambertianCanopyScattering(R = 0.45, T = 0.05)
Z⁺⁺, Z⁻⁺ = compute_Z_matrices(canopy, μ, LD, 0:4)            # Fourier stack
```

See the [tutorials](https://RemoteSensingTools.github.io/CanopyOptics.jl/dev/)
for 4SAIL, specular reflection, mixed leaf/wood canopies, clumping,
and microwave dielectric models.

## Conventions

Canopy Z matrices use `Z[i_out, j_in]` and normalize BiLambertian
leaves by `G(μ_in) * (R + T)`, matching the vSmartMOM phase-matrix
convention. Use `compute_Z_matrices(model, μ, LD, 0:m_max)` to build a
full Fourier stack. Numerical integration settings live in
`CanopyQuadrature` and pass via `compute_Z_matrices(...; quadrature)`.

The leaf single-scattering albedo `ϖ = R + T` is **not** folded into
`Z`; layer solvers multiply by `ϖ` separately. Specular kernels do
include their Fresnel/roughness strength in the returned Z, so
consumers that multiply canopy Z by a separate `ϖ` (e.g. vSmartMOM
layer assembly) need a matching effective albedo when using specular
or composite components.

Canopy scattering parameters are `Real`-generic, so `ForwardDiff.Dual`
leaf reflectance/transmittance and specular parameters can flow
through Z assembly.

## v0.2 highlights

- Closed-form Fourier Z matrices for bi-Lambertian canopy scattering.
- Stokes-aware Z expansion via `npol = 1`, `3`, or `4`.
- Additive `CompositeCanopyScattering` for diffuse + specular leaf
  components.
- Mixed leaf/wood canopy components with separate angle distributions
  and area indices.
- Clumping utilities for propagation/extinction and a Kuusk-style
  hotspot utility for direct-beam source-term integration.
- Constant, LUT, and polynomial wood reflectance models.
- 4SAIL bidirectional reflectance kernel with full
  `rsot = rsost + rsodt` decomposition (CPU + GPU via
  KernelAbstractions).
- Microwave dielectric ε(T, f) for water, ice, soil, and fresh leaves
  (Ulaby & El-Rayes 1987).

## Roadmap

- Joint design note pinning the CanopyOptics ↔ vSmartMOM Z / ϖ /
  Stokes / Fourier contract in one place across both repos.
- Microwave canopy expansion: full Stokes phase matrices for oriented
  scatterers (Rayleigh-Gans disk, finite cylinder), ensemble averaging
  into per-layer (τ, ϖ, Z) for vSmartMOM coupling.
- Effective-albedo interface for specular/composite canopy layers in
  vSmartMOM.
- SIF source-term hooks.

## Citation

A DOI is forthcoming. In the meantime, please cite as:

> Frankenberg, C. and contributors. *CanopyOptics.jl: canopy radiative-transfer
> building blocks for vSmartMOM.* RemoteSensingTools, 2025.
> https://github.com/RemoteSensingTools/CanopyOptics.jl
