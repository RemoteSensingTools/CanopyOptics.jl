# Leaf reflection, transmission and material absorption

`prospect(leaf, optics)` still returns `(T, R)`. The additive
[`prospect_budget`](@ref) interface returns spectral arrays named `transmission`,
`reflectance` and `absorptance`. All three describe monochromatic fractions of
incident radiation. At one wavelength, energy and photon fractions are identical.
Absorption includes all material; it is neither chlorophyll-only APAR nor net
energy retained after fluorescence.

Absorptance is evaluated from absorption within the plates and their multiple
reflections. It is not inferred by subtracting rounded reflectance and transmission,
and none of the three outputs is clipped or normalized. This matters for weakly
absorbing material: absorption can be much smaller than the spacing of floating
point numbers near the separately reported reflection and transmission.

The calculation retains the PROSPECT-PRO plate model, the optical coefficient
database, and its conventional first-surface 40° and internal diffuse 90° angular
integrations. It changes the numerical evaluation of the internal attenuation and
fractional stack. In particular, the old epsilon-adjusted exponential integral and
Stokes discriminant are removed; an artificial loss is not needed to solve a
nonabsorbing stack.

## Numerical limits and derivatives

Internal slab transmission is ``\tau=2E_3(k)``. For small nonnegative material
optical depth, its absorptance is evaluated as

```math
1-\tau=-\operatorname{expm1}(-k)+k\exp(-k)-k^2 E_1(k).
```

At zero optical depth, ``\tau=1`` and ``d\tau/dk=-2``. First derivatives therefore
remain available at the nonabsorbing endpoint; the physical second derivative
contains a logarithmic singularity. No finite Hessian at that endpoint is claimed.
For strong absorption, transmission is evaluated independently using ``E_3`` so
that subtracting absorption from one does not erase a small transmitted flux.

For ``m=N-1`` remaining plates, the code evaluates the analytic continuation of
an integer plate stack. A series handles its conservative limit; scaled hyperbolic
functions handle strong absorption or large thickness. The switch accounts for
both optical depth and plate count. A `log1p` expression retains weak attenuation
when a very large plate count makes the full stack absorbing.

The exactly opaque, zero-thickness endpoint is defined as an empty remainder.
Opaque and zero-thickness limits do not commute, so no general finite thickness
derivative is claimed at their joint endpoint. Finite-transmission derivatives at
``N=1`` are retained.

## Input domain and precision

Both entry points now enforce finite ``N\ge1`` and finite, nonnegative material
contents and absorption coefficients. Fractional plate counts are supported.
Inputs are not repaired silently. Float32 and Float64 are tested, including
first-order ForwardDiff calculations; optical coefficients should use the desired
precision as well as leaf parameters. Mixed input precisions retain Julia's
promotion rules.

The existing analytic Fresnel routine is retained and checked against independent
angular quadrature across the optical database's index range. Its nearly
index-matched limit (refractive index approaching one) remains ill-conditioned.
Nonfinite or out-of-range surface transmissions are rejected; this check does not
establish accuracy for all indices outside the tested range. Such applications
need a separately validated surface treatment. No surface coefficient is clipped.

The installed special-function library evaluates ``E_3(5)`` about 7.8 parts per
million low in Float32. Raising the plate transfer to a larger thickness can
amplify that relative transmission error: the independent ``N=8,k=5`` reference
differs by about 65 parts per million. The Float32 reference tests allow ``10^{-4}``
relative error for transmission, with tighter ``8\times10^{-6}`` bounds for
reflection and directly computed absorption. This is an explicitly measured
precision limit, not a conservation correction. Float64 uses ``2\times10^{-12}``
relative bounds for all three values. Derivatives have separate absolute and
relative tolerances because strongly absorbing leaves can have derivatives much
smaller than the spacing of their almost saturated absorptance.

This host-side PROSPECT calculation does not claim CUDA support for the special
functions. Its optical tables can feed separate CPU/GPU canopy and pigment
calculations. No particular canopy grid, physiology type or time integrator is
introduced by the budget interface.

## Independent verification

`test/reference/generate_prospect_reference.py` uses mpmath at 90 decimal digits
(and 750 digits for extreme stress cases). Angular integration checks attenuation;
integer plate counts use repeated flux adding, and fractional counts use a general
matrix eigendecomposition. The committed CSV references need no Python dependency
when running Julia tests. Separate tests compare values and first derivatives,
check both numerical switches, and exercise the actual spectral coefficient data.
The astronomical-thickness and subnormal-transmission stress cases check values;
they do not establish derivative accuracy at every representable extreme. Derivative
references cover ordinary leaf plate counts and additional switch checks up to
1,000 remaining plates.

At very weak absorption the meaningful checks are relative accuracy of the direct
absorptance and absolute closure of ``R+T+A``. Requiring a relatively accurate
``1-R-T`` would instead demand information the rounded fluxes cannot represent.
