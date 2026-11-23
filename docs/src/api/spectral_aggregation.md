# Spectral Band Aggregation

Plain-array wavenumber band integration, wavelength–wavenumber grid conversion, and
composable leaf spectrum evaluation across UV, shortwave, and thermal bands. See the
[Spectral Band Aggregation guide](../pages/band_aggregation.md) for mathematical
conventions, fallback rules, and pigment-share aggregation.

```@contents
Pages = ["spectral_aggregation.md"]
Depth = 2
```

## Band Averaging & Sub-Band Partition Shares

```@docs
band_average
band_average!
subband_weights
subband_weights!
```

## Spectral Weights & Wavelength–Wavenumber Conversion

```@docs
AbstractSpectralWeight
FlatWeight
PlanckWeight
TabulatedWeight
ReferenceSolar
wavelength_to_wavenumber
wavenumber_to_wavelength
```

## Leaf Spectrum Models & Pigment Shares

```@docs
leaf_spectrum
TabulatedLeaf
ConstantLeaf
CompositeLeaf
leaf_pigment_shares
```

## Qualified Integrals & Extension Hooks

The following symbols are not exported by `using CanopyOptics`; access or extend them
via qualified `CanopyOptics.*` bindings.

```@docs
CanopyOptics.band_weight
CanopyOptics.band_weight!
CanopyOptics.weight_interval_integral
CanopyOptics.spectral_domain
CanopyOptics.SpectralCoverageError
```
