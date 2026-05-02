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
  <a href="https://github.com/RemoteSensingTools/vSmartMOM.jl/commits/master">
    <img src="https://img.shields.io/github/commit-activity/y/RemoteSensingTools/CanopyOptics.jl"
         alt="Github Commit Frequency">
  </a>
</p>

Tools for the computation of canopy optical parameters, see docs for details
(works in concert with vSmartMOM.jl).

Convention note: canopy Z matrices use `Z[i_out, j_in]` and normalize
BiLambertian leaves by `G(μ_in) * (R + T)`, matching the vSmartMOM
phase-matrix convention. Use `compute_Z_matrices(model, μ, LD, 0:m_max)` to
build a full Fourier stack. Numerical integration settings live in
`CanopyQuadrature` and can be passed as `compute_Z_matrices(...; quadrature)`.

Canopy scattering parameters are `Real`-generic, so `ForwardDiff.Dual` leaf
reflectance/transmittance and specular parameters can flow through Z assembly.

## v0.2 highlights

- Closed-form Fourier Z matrices for bi-Lambertian canopy scattering.
- Stokes-aware Z expansion via `npol = 1`, `3`, or `4`.
- Additive `CompositeCanopyScattering` for diffuse + specular leaf components.
- Mixed leaf/wood canopy components with separate angle distributions and area
  indices.
- Clumping utilities for propagation/extinction and a Kuusk-style hotspot
  utility for direct-beam source-term integration.
- Constant, LUT, and polynomial wood reflectance models.

Specular kernels include the Fresnel/roughness strength in the returned Z
contribution. Consumers that multiply canopy Z by a separate single-scattering
albedo, such as vSmartMOM layer assembly, need to provide a matching effective
albedo when using specular components in full RT.

## Roadmap / ToDo for next releases

- Wire Kuusk hotspot source scaling into vSmartMOM's canopy source-term path.
- Add default literature-backed wood spectra and YAML access.
- Add an explicit effective-albedo interface for specular/composite canopy
  layers in vSmartMOM.
- Add SIF source-term hooks.
