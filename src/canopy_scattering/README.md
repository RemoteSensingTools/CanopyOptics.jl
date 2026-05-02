# Canopy Scattering Source Map

This folder contains canopy-to-Z-matrix assembly code.

- `projection_geometry.jl`: Ross `G` projection geometry and reference checks.
- `quadrature_controls.jl`: `CanopyQuadrature` compatibility helpers.
- `../wood_optics/reflectance.jl`: wood reflectance models and the opaque
  Lambertian wood scatterer used by Z assembly.
- `specular.jl`: specular leaf-surface reflection and its Z moments.
- `bilambertian_fourier.jl`: analytic bi-Lambertian Fourier moments. The
  internal one-sided projection helpers split `Ω ⋅ Ω_L` into the positive and
  negative leaf faces before taking Fourier moments.
- `z_matrices.jl`: public `compute_Z_matrices` dispatch and composite summation.
- `legacy_azimuth_reference.jl`: old direct-azimuth reference and compatibility paths.

Production callers should prefer `compute_Z_matrices(model, μ, LD, m)` or
`compute_Z_matrices(model, μ, LD, 0:m_max)`.
