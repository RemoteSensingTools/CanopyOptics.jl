# Changelog

## v0.2.0

- Added closed-form Fourier Z matrices for bi-Lambertian canopy scattering.
- Added Stokes-aware canopy Z assembly with `npol = 1`, `3`, or `4`.
- Added additive composite canopy scattering, including diffuse + specular
  leaf components.
- Added canopy quadrature controls, clumping utilities, and a Kuusk-style
  hotspot gap-probability utility.
- Added mixed leaf/wood canopy components and wood reflectance model
  infrastructure.
- Exported the public leaf-angle distribution constructors and `G` utilities.

Specular Z matrices include Fresnel/roughness strength in the returned kernel.
Downstream RT solvers that multiply by a separate canopy single-scattering
albedo need a matching effective-albedo path before using specular components
in full layer assembly.
