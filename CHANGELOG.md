# Changelog

## Unreleased

### Added
- Microwave dielectric model for fresh leaves
  (`LeafUlabyElRayes1987` + `dielectric` method) and exported
  `AbstractMaterial`, `AbstractWater`, `AbstractSoil`,
  `AbstractVegetation`, `SoilMW` traits and types.
- Canonical docstring on `foursail!` documenting both the in-place
  vector and batched-matrix dispatch shapes and citing Verhoef
  (1984); resolves the unresolved `@ref foursail!` Documenter
  warning.
- README quickstart, abstract, install block, citation block, and
  microwave-direction note. Fixed badge link (was pointing at
  `vSmartMOM.jl`).

### Changed
- Public `dielectric(...)` methods now throw `ArgumentError` instead
  of `AssertionError` for out-of-range inputs (water, ice, soil,
  leaf). Consumers that explicitly catch `AssertionError` need to
  switch.

### Fixed
- `foursail` no longer returns `NaN` for conservative leaves
  (`ρ + τ = 1`, e.g. `(0.5, 0.5)`). The `m → 0`, `rinf → 1`
  singularity is patched with a `cbrt(eps)` floor on `m` (same
  trick as canonical PROSAIL Fortran). The resulting answer
  matches the well-defined `ρ + τ → 1⁻` limit to ~6 decimal
  places.
- `LeafUlabyElRayes1987` now enforces its documented frequency
  domain (0.2–20 GHz) — calls outside that range throw
  `ArgumentError` instead of returning unsupported values.

### Removed (breaking)
- The unfinished `forest_prototyping/` prototype is gone, along with
  its public exports `wood_forward`, `wood_backward`, `afsal`,
  `asal`, `abs_components`. These functions referenced undefined
  module-scope names and threw `UndefVarError` on first call, so
  they had no working callers; they are removed rather than
  refactored. The Karam/Fung/Antar (1988) algorithms they sketched
  will be reintroduced as part of the planned microwave-canopy
  expansion.
- Dropped `QuadGK` and `YAML` dependencies (only `forest_prototyping`
  used them).

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
