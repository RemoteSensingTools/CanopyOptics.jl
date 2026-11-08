# Conditional photon emission spectra

An intrinsic spectrum distributes **already produced internal photons** over
wavelength. It is separate from leaf absorption, pigment/PSII allocation,
physiological fluorescence yield, internal escape and canopy transport. The core
supplies no default empirical spectrum.

```julia
# Inputs and their units come from the caller's explicit experiment parameters.
mass = photon_density_per_nm .* quadrature_width_nm
source = TabulatedPhotonSpectrum(wavelength_nm, mass)
probability = emission_probabilities(source, StrictStokes(), excitation_nm)
```

`TabulatedPhotonSpectrum` stores a discrete photon measure, not a density. Multiply
by the spectral integration measure once before construction. Wavelengths and
masses must share a floating-point type; Float32 and Float64 are tested. Inputs
are copied, but stored vectors remain mutable; preparation revalidates them.
Host preparation retains raw relative masses and filters support **before**
max-scaled normalization, so an extremely small retained tail survives even when
excluded values are large. A retained mass smaller than the floating-point range
relative to other retained masses can still round to zero; this is not arbitrary
precision. Probabilities stay on the original grid with exact excluded zeros.

`StrictStokes()` permits only emission wavelengths longer than excitation. It is
an explicit approximation excluding thermal anti-Stokes excitation. Conditioning
changes the reference shape while preserving the separately prescribed photon
yield. Removed reference weights do not establish extra heat or quenching.
Empty positive support is an error. A physically dark source should be skipped
by the caller rather than requesting an undefined conditional distribution.

## Combining with optical escape

For excitation node `e` and emission node `m`, prepare photon probabilities
`p[m,e]`. Evaluate the independent optical fates at each emitted wavelength and
the appropriate excitation-dependent birth profile. The front, back and
reabsorbed operators are

```math
M_{F,m e}=p_{m e}F_{m e},\quad
M_{B,m e}=p_{m e}B_{m e},\quad
M_{A,m e}=p_{m e}A_{m e}.
```

Each excitation column sums over wavelengths and fates to one, within roundoff.
Multiply these operators by already-produced photons in each excitation bin.
Do not multiply absorption or yield again. Share operators across equivalent
optical/source classes instead of storing one spectral tensor per leaf population.
The matrices can be transferred to CPU/GPU backends for batched application;
the vector-owning source object and validation stay on the host.

Convert all three photon outcomes to energy using their **emission** wavelengths.
Strict Stokes support gives a nonnegative excitation-to-emission energy difference
per photon. This source/escape accounting is not a complete biochemical or thermal
budget. Reabsorption terminates the photon; secondary fluorescence is separate.
The regular-plate fates are hemispheric face totals. Directional canopy SIF needs
its own face-radiance closure and canopy transport.

## Replacing a component

Subtype `AbstractEmissionSupport` and implement the scalar Boolean
`emission_allowed(policy, emission_nm, excitation_nm)` for a different support
rule. Inputs provided by preparation are finite positive wavelengths of the same
type. The tests demonstrate a finite Stokes band without changing the source,
escape or physiology code. A future analytic or state-dependent spectrum can
dispatch on the first argument of `emission_probabilities`; no storage registry
or canopy-grid dependency is required.

Scalar validation compares normalization with BigFloat, including support at
equality, rare tails, unequal caller-owned measures, mutable inputs and invalid
precision. Prepared CPU/CUDA matrix tests cover both precisions, dark columns,
multiple excitation distributions, front/rear birth profiles, reordered and empty
batches, photon closure and emitted-energy bounds. These are numerical and
interface checks, not empirical calibration of a spectral shape.

API reference: [`TabulatedPhotonSpectrum`](@ref),
[`AbstractEmissionSupport`](@ref), [`StrictStokes`](@ref),
[`emission_allowed`](@ref), and [`emission_probabilities`](@ref).
