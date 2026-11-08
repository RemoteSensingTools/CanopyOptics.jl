# PROSPECT lossless-limit validation, 2026-11-08

This repair preserves the PROSPECT-PRO plate model and optical database. It removes
artificial epsilon absorption and adds direct material absorptance while preserving
`prospect`'s `(T,R)` interface. The new `prospect_budget` result can support later
pigment attribution without embedding any physiology or canopy-state dependency.

Independent code/bugs and physics reviewers cleared the implementation and
reference construction. The full package suite passes **90,301 checks**, including
**3,178 new checks** for this repair. The new references use 90-digit angular
integration, repeated integer plate adding and general matrix eigendecomposition;
750-digit stress fixtures cover extreme primal values. First derivatives have
independent finite-difference references, including zero absorption and N=1.

The retained legacy Fresnel calculation is tested against independent angular
quadrature over indices 1.25–1.6. Its nearly index-matched limit remains unsupported
by the accuracy claim; invalid surface transmissions are rejected without clipping.
The Float32 E3 precision limitation, derivative domains and separate tolerances
are described in the [numerical documentation](src/pages/prospect_budget.md).

A separate local integration audit used the same three physical leaf profiles,
all 2,101 native wavelength samples and both floating-point types, plus zero-material
controls. For ordinary leaves, the largest changes across those profiles are:

| Quantity | Float32 maximum absolute change | Float64 maximum absolute change |
| --- | ---: | ---: |
| Transmission | 3.88e-7 | 7.22e-16 |
| Reflection | 3.58e-7 | 5.56e-16 |
| Direct absorption versus previous 1−R−T | 6.26e-7 | 1.03e-15 |

For every zero-material control, direct absorptance is exactly zero. The former
rounded residual reached 1.60e-5 in Float32 and 6.83e-10 in Float64. The repaired
R+T+A closure error is at most one machine epsilon for these profiles in either
type. Rounded `1−R−T` may still contain roundoff and should not replace the direct
absorptance when attributing very small material absorption.

These measurements concern the CPU/ForwardDiff optical calculation. They do not
claim CUDA compilation of the exponential integral, a calibrated SIF model,
chlorophyll-only absorption, or a canopy energy-balance solution.

Exact source/reference/environment hashes and validation commands are recorded in
`prospect_lossless_validation.json`. The public tests and reference generator are
self-contained; no private physiology dependency is introduced.
