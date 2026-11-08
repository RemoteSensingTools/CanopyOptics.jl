# Internal leaf photon escape

`RegularPlateEscape` computes the fates of a photon **already produced inside**
one symmetric regular leaf plate. It returns front escape, back escape and
terminal material reabsorption. The birth-depth model is a separate argument:

```julia
escape = RegularPlateEscape(emission_depth, interface_transmission)
birth = ExponentialBirthProfile(excitation_depth, front_amplitude, rear_amplitude)
fates = leaf_emission_fates(escape, birth)
uniform_fates = leaf_emission_fates(escape, UniformBirthProfile())
```

All numerical inputs use the same `Float32` or `Float64`. Construct validated
parameters on the host, then apply the scalar function or broadcast it over
arrays of prepared parameters on CPU/CUDA. Both built-in birth choices are
isbits objects; application allocates no host memory. The module is
`src/leaf_optics/plate_emission.jl`.

## What the parameters mean

Depth `z=0` is the fixed anatomical front and `z=1` the back. Effective excitation
depth `a` and emission depth `b` are nonnegative dimensionless attenuation depths.
They describe traversal survival `exp(-a)` and `exp(-b)`. Interface transmission
`t` is the internal-to-outside probability at **each boundary**, distinct from
the traversal survival. Reflection is `r=1-t`. Require `0<t≤1`; a perfectly
trapping lossless plate needs a fourth, retained-photon fate and is excluded.

The exponential birth density is proportional to

```math
f\exp(-az)+g\exp[-a(1-z)].
```

The host constructor normalizes nonnegative relative amplitudes with max scaling;
both normalized components are retained when representable. An all-zero pair is
an error because it cannot define a conditional birth density. Zero excitation
depth gives the uniform profile even when the amplitudes differ. It does not
assert that a leaf with no absorbed light produces fluorescence.

For front illumination of a homogeneous regular plate, relative amplitudes are
`(1, r_exc*exp(-a))`. The excitation depth uses **total material attenuation**.
Homogeneous pigment allocation and depth-independent physiological yield cancel
when conditioning on an internal birth. This effective exponential internal-depth
closure is an approximation, as in [Platespect](https://doi.org/10.1016/j.rse.2025.114990),
not an exact angular solution of light inside a leaf. This component does not
invert fractional PROSPECT stacks into possibly nonpassive effective interfaces.

## Photon paths and numerical limits

Emission initially travels toward either face with probability one half. Let
`F₀,B₀` be first-arrival probabilities after within-plate absorption. Reflected
traversals have survival `q=r*exp(-b)`, giving

```math
F=t(F_0+qB_0)/(1-q^2),\qquad B=t(B_0+qF_0)/(1-q^2).
```

Implementation factors the denominator using the retained transmission and
direct traversal loss `-expm1(-b)`. Reabsorption is computed as first-flight loss
plus absorption during subsequent reflected traversals. It is not reconstructed
solely as a small residual of rounded escaped probabilities. No result is clipped
or rescaled. Lossless emission gives `(1/2,1/2,0)` for any supported birth profile.
The initial reabsorption derivative is `(1+r)/(2t)`, independently of birth depth.

Normalized exponential overlap integrals handle equal depths, finite depths whose
sum overflows and sharply localized sources. Below `b=eps(FT)^(1/4)`, a four-term
moment expansion preserves weak absorption. Positive moment series and a
large-depth recurrence avoid cancellation at small excitation depth. These
thresholds and coefficients are numerical algorithm choices, not physical constants.
Relative reabsorption accuracy near the Taylor/direct-integral switch need not
equal machine epsilon; the independent test envelope is `3e-5` in Float32 and
`2e-11` in Float64, with a few-subnormal-ULP absolute allowance in batch checks.

Tests include independent 2048-bit boundary balances, 90-digit numerical source
depth integration with explicit repeated-reflection paths, numerical switch
neighborhoods, face reversal, zero/weak/opaque absorption, normalized extreme
amplitudes and smallest-positive interface transmission. CPU and CUDA batches
exercise both birth choices, reordered and empty arrays. ForwardDiff is not yet
a supported contract for these preparation types.

The [validation record](../data/plate_emission_validation.toml) includes source and
reference hashes, test logs and 1,230 scalar accuracy probes per precision. Against
references rounded to the tested FT, the measured maxima were:

| Metric | Float32 | Float64 |
| --- | ---: | ---: |
| Absolute photon closure error | 2.38e-7 | 2.22e-16 |
| Relative front/back error, normal outputs | 4.18e-7 | 1.43e-15 |
| Relative reabsorption error, normal outputs | 9.01e-6 | 2.41e-12 |

These are measurements on the archived probe grid, not universal error bounds.

## Composition and extension

Neither the escape model nor the birth profile owns an intrinsic emission
spectrum, absorptance, pigment selection, PSII allocation or physiological yield.
Those quantities enter the surrounding calculation once. A spectral operator
must supply explicit emission measures and wavelength support; emission energy
uses emission wavelength. Photon conservation here is not a complete leaf energy
budget. Reabsorbed photons terminate; secondary fluorescence is a separate model.

These are hemispheric face totals. Angular emission into a canopy requires a
separate face-radiance closure. Two-sided excitation retains distinct source
contributions; anatomical face labels do not change with illumination. A custom
birth type can subtype `AbstractLeafBirthProfile` and implement
`leaf_emission_fates(escape, birth)`; a different escape structure can dispatch on
the first argument. Neither extension changes a photosynthesis model or host grid.

Actual optical-database preparation, multiple-plate assembly, intrinsic PSII
spectra and canopy coupling follow separately. No fractional-stack, directional
SIF or complete coupled energy capability is implied by this scalar component.
