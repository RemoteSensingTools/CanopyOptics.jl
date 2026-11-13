# Excitation & Leaf Fluorescence

Inelastic chlorophyll fluorescence optics are factored into three composable stages
that operate on CPU scalars/arrays and `KernelAbstractions` GPU backends without
embedding photosynthetic physiology:

1. **Constituent absorption partitioning** ([`partition_absorption`](@ref)) splits total leaf material absorptance at each excitation wavelength across pigments, water, and dry matter in proportion to their specific optical depths.
2. **Antenna selection & photosystem allocation** ([`antenna_excitation`](@ref), [`photosystem_excitation`](@ref)) converts absorbed spectral photons into photosynthetic antenna absorption (`ChlorophyllOnly`, `CarotenoidAntenna`, `MultiChlorophyllAntenna`) and PSII/PSI excitation rates ready for physiological models such as `LeafPhotosynthesis.jl`.
3. **Internal photon escape & emission spectra** ([`diffuse_plate_optics`](@ref), [`dorsiventral_plate_optics`](@ref), [`leaf_emission_fates`](@ref), [`emission_probabilities`](@ref)) propagate internally emitted fluorescence photons through 1-zone (`RegularPlateEscape`) or 2-zone (`DorsiventralPlateEscape`) mesophyll stacks with exact photon conservation (`front + back + reabsorbed = 1` for all `t ∈ (0, 1]`) and Stokes-valid emission support.

```@contents
Pages = ["fluorescence_and_excitation.md"]
Depth = 2
```

## Constituent Absorption Partitioning

```@docs
AbstractAbsorptionPartition
HomogeneousAbsorption
AbsorptionPartitionResult
AbsorptionPartitionStatus
AbsorptionPartitionValid
AbsorptionNoMaterial
AbsorptionInvalidInput
AbsorptionIncompatibleInput
partition_absorption
absorption_valid
```

## Antenna Selection & Photosystem Excitation

```@docs
AbstractAntennaSelection
ChlorophyllOnly
CarotenoidAntenna
MultiChlorophyllAntenna
antenna_excitation
select_absorbed_photons
photosystem_excitation
```

## Internal Leaf Photon Escape (1-Zone & 2-Zone Plates)

```@docs
AbstractLeafBirthProfile
UniformBirthProfile
ExponentialBirthProfile
DorsiventralBirthProfile
RegularPlateEscape
DorsiventralPlateEscape
diffuse_plate_optics
dorsiventral_plate_optics
leaf_emission_fates
```

## Conditional Fluorescence Emission Spectra

```@docs
TabulatedPhotonSpectrum
AbstractEmissionSupport
StrictStokes
emission_allowed
emission_probabilities
```
