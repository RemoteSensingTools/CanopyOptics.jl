# Excitation allocation and several emission sources

A leaf can contain several fluorescence sources, for example photosystem II and
photosystem I. They share the leaf's absorbed photons but differ in emission
spectrum, photon yield and possibly birth-depth profile. Three small functions
compose the existing absorption, escape and spectrum components for any number of
named sources; none of them fixes the number of photosystems:

```julia
allocated = allocate_excitation(antenna_photons, fractions)   # per excitation bin
ops = emission_operators(probabilities, fates)                # once per source
signal = band_signal(response, ops.front, produced_photons)   # once per source and face
```

The caller supplies spectra, fractions, yields, birth profiles and detector weights.
The core supplies no default for any of them and contains no physiology.

## Joint excitation allocation

`allocate_excitation(photons, fractions)` divides photons in one wavelength bin
among named sinks, `sink_i = f_i * photons`, and returns
`(sinks = NamedTuple, unallocated = FT)` with the keys and order of `fractions`.
`photons` and the fractions share one floating-point type; mixed precision is a
`MethodError`, as for the other scalar kernels. Valid input has finite
`photons ≥ 0`, finite fractions in ``[0,1]`` and

```math
\sum_{i=1}^{N} f_i \le 1 + N\,\varepsilon_{\mathrm{FT}} ,
```

where ``\varepsilon_{\mathrm{FT}}`` is `eps(FT)`. The allowance covers the
left-to-right summation of ``N`` caller-rounded fractions, so fractions normalized
to one in the caller's precision are accepted. Otherwise every sink and the
remainder are typed NaN, as in [`photosystem_excitation`](@ref).

The remainder is `unallocated = max(1 - Σf, 0) * photons`. It is never negative and
never exceeds `photons`; within the allowance a sum slightly above one leaves zero.
Closure `Σ sinks + unallocated = photons` holds to about ``(N+1)`` rounding errors.
When the fractions nearly exhaust the photons, the small remainder has an absolute,
not relative, accuracy of that order. With no sinks, valid photons are returned
unallocated. A one-sink call reproduces `photosystem_excitation(photons, f)`
exactly, so `photosystem_excitation` remains the one-sink shorthand.

The scalar call is type stable and allocation free for any number of sinks. It
broadcasts over arrays of photons with shared (`Ref`) or per-bin fractions, on the
CPU and on GPU backends.

### Two allocation conventions

`allocate_excitation` does not know what `photons` are; the caller's fractions
state the convention. Published allocations use one of two denominators:

- fractions ``f_i`` of **antenna excitation**, applied to the output of
  [`antenna_excitation`](@ref);
- fractions ``a_i`` of **leaf-absorbed photons**, applied to total absorption.

With ``h = Q_{\rm ant}/Q_{\rm abs}`` the antenna share of absorbed photons in the
same bin, both describe the same sink when

```math
a_i = h\,f_i ,\qquad a_i\,Q_{\rm abs} = f_i\,Q_{\rm ant}.
```

Applying ``a_i`` after antenna selection counts the transfer loss twice (the sink
is too small by ``h``). Applying ``f_i`` to total absorption omits it (too large by
``1/h``). For compatible allocations ``\sum_i a_i \le h``, with equality only when
the antenna excitation is allocated exhaustively (``\sum_i f_i = 1``).

The remainders differ between the routes:

```math
U_{\rm ant} = \Bigl(1-\sum_i f_i\Bigr) Q_{\rm ant},\qquad
U_{\rm abs} = \Bigl(1-\sum_i a_i\Bigr) Q_{\rm abs} = (1-h)\,Q_{\rm abs} + U_{\rm ant}.
```

The absorbed-photon remainder also contains the photons that never reached the
antenna. Neither remainder by itself specifies nonphotochemical quenching, heat or
fluorescence; it is only the part of the input not assigned to a named sink.

To check published ``a_i`` against a pigment description, compute ``h`` in each
bin from the partition and a caller-chosen antenna selection,
`antenna_excitation(selection, partition) / partition.total_absorbed`, and compare
it with ``\sum_i a_i`` (or with ``a_i/f_i`` where both are given). The comparison
is meaningful only where the pigment data and the physiological state match the
published measurement. In-vitro absorption shares and in-vivo transfer
efficiencies are different quantities and should not be mixed in one fraction.

## Emission operators, one per source

For excitation bin `e` and emission bin `m`, let `p[m, e]` be the conditional
emission probabilities of one source (columns from
[`emission_probabilities`](@ref), any spectrum and support rule) and `F, B, A` the
front, back and reabsorbed fates of a photon born by that source. Then
`emission_operators(p, fates)` returns

```math
M_{F,me}=p_{me}F_{me},\qquad M_{B,me}=p_{me}B_{me},\qquad M_{A,me}=p_{me}A_{me}.
```

`fates` is either a NamedTuple of three matrices `(front = F, back = B,
reabsorbed = A)`, with keys in any order, or a matrix of `(front, back,
reabsorbed)` NamedTuples as produced by broadcasting
`leaf_emission_fates.(reshape(escape, :, 1), reshape(birth, 1, :))`. All arrays
share the size and floating-point type of `p`; mismatches throw on the host. The
arrays may be CPU or GPU arrays. Nothing is clipped or renormalized.

[`emission_closure_residual`](@ref) returns, per excitation column,

```math
r_e=\sum_m\bigl(M_{F,me}+M_{B,me}+M_{A,me}-p_{me}\bigr),
```

which is at roundoff level when the fates close (`F + B + A = 1`): each column of the
three operators sums to one for a conditional spectrum and to zero for a dark
column. The caller chooses the tolerance.

### Dark columns: contract

An excitation bin in which the source has no supported emission on the grid gives
an all-zero (dark) column of `p`; [`emission_probabilities`](@ref) rejects that case,
so only the caller can build such a column. **A dark column is permitted only in
excitation bins where the source produces no photons.** Photons produced in a dark
column have no front, back or reabsorbed outcome: they are silently lost, and the
closure residual is still zero because it does not see production. For an
excitation bin with production and empty emission support, the caller must either
expand the emission grid until the support is represented or treat the bin as an
error. Nothing is clipped or renormalized to hide the loss.

The caller can detect a violation with existing functions, before applying the
operators:

```julia
dark = iszero.(vec(sum(p; dims = 1)))   # dark excitation columns of this source
all(iszero, produced[dark]) || error("photons produced in a dark excitation column")
```

With that check satisfied, the three faces summed over emission and excitation
equal the total produced photons of the source to rounding.

Call `emission_operators` once per source. Sources may differ in probabilities and
in fates: a source with a different birth-depth profile (for example because its
yield or stoichiometry varies with depth) gets its own fates through the
birth-profile argument of [`leaf_emission_fates`](@ref); see
[State-dependent birth profiles](@ref). No source container type is
introduced; a NamedTuple of operators keyed by source name is sufficient.

## Band observation

`band_signal(response, operator, produced_photons)` evaluates

```math
S = \sum_m r_m \sum_e M_{me}\,q_e = r^{\mathsf T} M q ,
```

where `produced_photons[e]` ``= q_e`` are photons already produced by the source in
excitation bin `e` (absorption, allocation and yield applied once, by the caller;
zero in dark columns) and `response[m]` ``= r_m`` are dimensionless
**photon-basis** weights on the emission grid. Emission-bin measures are already
inside the probabilities and are not applied again. The result is in photons, in
the units of `produced_photons`.

A photon-basis weight is, for example, a detector quantum efficiency (electrons
per incident photon) times a filter transmission. A power responsivity in A/W is
not a photon weight. If the caller folds the photon energy ``hc/\lambda_m`` of the
**emission** wavelength into ``r_m`` (for example responsivity times
``hc/\lambda_m``), the result is in energy-based units instead of photons; that
conversion is the caller's and always uses emission, not excitation, wavelengths.

`S` is a linear functional of hemispheric face totals. It knows nothing of
instrument geometry, collection aperture or fluorometers, and it gives no
directional radiance. Dimension or element-type mismatches throw.

## Combining several sources

With a yield ``\phi_s`` per source ``s`` (the caller's value) and the allocated
excitation ``Q_{s,e}``, the produced photons are ``q_{s,e} = \phi_s Q_{s,e}``, and a
band observation of the front face is

```math
S_F = \sum_s r^{\mathsf T} M_{F,s}\,q_s .
```

Compute one operator per source, multiply the yields in the caller, and sum the
band signals. Adding a third source adds one key to the fractions, one operator and
one term; no package code changes.

## State-dependent birth profiles

### Depth coordinate

``z = x/d`` is the normalized geometric depth of a homogeneous plate of thickness
``d``, from the front face (``z = 0``) to the back (``z = 1``). It is the same at every
wavelength and stays defined when ``b = 0``. Because the plate is homogeneous, the
effective excitation and emission optical depths from the front face to depth ``z``
are ``a\,z`` and ``b\,z``. Equal layers in ``z`` are therefore equal in thickness and
in optical depth at every wavelength. Nonuniform pigment changes the local
excitation and emission attenuation; that cannot be represented by changing the
birth weights alone (see the limits below).

### Production: allocation and yield are separate factors

For excitation bin ``e`` let ``Q_e`` be the antenna excitation of the leaf (photons,
for example from [`antenna_excitation`](@ref)) and ``P_{e,k}`` the fraction of it
in layer ``k``, from [`birth_layer_fractions`](@ref) applied to the excitation
profile of that bin. For source ``s`` let ``\alpha_{s,e,k}`` be the fraction of the
local antenna excitation **allocated** to the source and ``\phi_{s,k}`` its
**yield**: emitted photons per excitation photon allocated to the source. Then the
photons produced by source ``s`` in layer ``k`` are

```math
W_{s,e,k} = Q_e\,P_{e,k}\,\alpha_{s,e,k}\,\phi_{s,k},
\qquad q_{s,e} = \sum_k W_{s,e,k},
```

``q_{s,e}`` is the source's production in bin ``e`` (the `produced_photons` of
[`band_signal`](@ref)), and `PiecewiseUniformBirthProfile(W[s, e, :])` normalizes
the same products into its birth profile. Quenching or the open-centre fraction
belong to ``\phi``; photosystem stoichiometry and antenna sharing belong to
``\alpha``.

An equivalent convention first conditions on the allocation:
``Q_{s,e} = Q_e \sum_k P_{e,k}\alpha_{s,e,k}`` is the excitation allocated to the
source and ``P^{(s)}_{e,k} = P_{e,k}\alpha_{s,e,k} / \sum_j P_{e,j}\alpha_{s,e,j}``
its source-conditioned layer fractions; then ``W_{s,e,k} = Q_{s,e}\,P^{(s)}_{e,k}\,\phi_{s,k}``
with the yield alone. With a depth-independent allocation ``f_s``,
``Q_{s,e} = f_s Q_e`` is the output of [`allocate_excitation`](@ref) and
``P^{(s)} = P``. Use one convention: putting the allocation into ``P^{(s)}`` or
``Q_{s,e}`` **and** into the yield counts it twice. A depth-varying allocation makes
the source-conditioned excitation profile non-exponential even at constant yield,
so such a source needs a layered profile even if its yield does not vary.

If neither allocation nor yield varies with depth, both cancel on normalization:
the birth profile is the excitation profile itself, which is exactly what
[`ExponentialBirthProfile`](@ref) describes. With depth-independent allocation
**and** yield the layered fates converge to the exponential ones with second-order
error in ``1/N``
(observed orders 1.95 to 2.00 from ``N = 8``; for ``a = 3``, ``b = 2``, ``t = 0.6``
the largest fate difference is ``1.1\times10^{-3}`` at ``N = 8`` and
``4.3\times10^{-6}`` at ``N = 128``). When both allocation and yield are
depth-independent, use the exponential profile directly; if either varies with
depth, the source needs a layered profile.

With layer values of a smoothly varying allocation or yield (for example at layer
midpoints), ``W_{s,e,k}`` approximates the layer integral of the production density
to second order in ``1/N``.

### Several illuminations and excitation columns

- Keep separate excitation columns where possible; each gets its own ``P_{e,k}``,
  production and profile. Keep them separate whenever the conditional emission
  spectra of the columns differ.
- If columns must be collapsed, combine **production**,
  ``W_{s,k} = \sum_e Q_{s,e}\,P^{(s)}_{e,k}\,\phi_{s,e,k}``, and normalize afterwards.
  Averaging normalized profiles discards the production weights.
- A yield that depends on the physiological state must be computed by the caller
  from the combined illumination of all columns first; the state is nonlinear in
  the light and is not a per-column property.
- For back illumination, swap the forward and backward excitation amplitudes of the
  `ExponentialBirthProfile` (``(r_{\rm exc}e^{-a},\,1)`` instead of
  ``(1,\,r_{\rm exc}e^{-a})``). The layers keep their anatomical order from the front
  face and the yield field keeps its depth indexing; nothing is mirrored.

### Fates and accuracy

The regular-plate fates of a layered profile are closed form. A photon born at depth
``z`` starts toward either face with probability ``1/2``; its first-flight survival
is ``e^{-bz}`` to the front and ``e^{-b(1-z)}`` to the back, averaged exactly over
each layer with ``L(x) = -\operatorname{expm1}(-x)/x``, ``x = b/N``. The first-flight
reabsorption is a sum of nonnegative terms (``-\operatorname{expm1}`` of the
boundary depths and ``1 - L(x)``, the latter by a series below ``x = 1``), so no
near-unit escape probabilities are subtracted when ``b`` is small; for
``b \le \varepsilon_{\mathrm{FT}}`` it is ``b/2``, exact to relative ``O(b)`` for any
profile, which also avoids rounding ``x = b/N`` on the subnormal grid. The layer
sums, and the normalization of the weights, use compensated (Neumaier) summation,
so the relative accuracy does not degrade with ``N``. The multiple reflections
between the faces use the same closure as the exponential profile. The cost is one
exponential per layer.

**Accuracy.** With ``u`` = `nextfloat(zero(FT))`, the smallest subnormal number,
every fate ``\hat F`` satisfies

```math
|\hat F - F| \le 8(1+b)\,\varepsilon_{\mathrm{FT}}\,|F| + (N + 3)\,u ,
```

and the reabsorbed fate has a further ``u/(2t)`` when ``b`` is below
`floatmin(FT)` and ``1 - t \ne 1`` in `FT`.

- The relative term includes the rounding of the exponent ``b\,j/N``, which is why
  it grows with ``b``.
- ``(N + 3)\,u`` is a derived worst case for underflow. Each of the at most ``N``
  layer products per face that falls below `floatmin` is rounded to a multiple of
  ``u``, the face sums are multiplied by factors of at most ``1/2``, and the
  reflection closure adds a few more roundings. It matters only for results below
  about ``(N+3)\,u/\varepsilon_{\mathrm{FT}}``: subnormal results, or results just
  above `floatmin`. Results that are not subnormal carry essentially the relative
  bound.
- ``u/(2t)`` comes from the shared reflection closure: the subnormal product
  ``(1 - t)(1 - e^{-b})`` is divided by about ``t``. The exponential profile has the
  same term. It is below `floatmin(FT)` and arises only for subnormal ``b``.

The tests assert exactly this bound in `Float32` and `Float64` against an
independent closed form evaluated with at least 2048 bits, and with more bits for
subnormal ``b``. They cover ``b`` from ``2u`` to 700 (to 800 in `Float64`), ``N``
from 1 to 1000, ``t`` from 1 down to ``u``, sparse deep sources whose escapes
underflow, and the CPU and CUDA batches. A scan of 3510 `Float32` and 3330 `Float64` cases
(``N`` up to 1000, deep, shallow, sparse and dense sources, ``b`` from ``2u`` to 1500,
``t`` from 1 to ``3u``) gave a largest ratio of error to this bound of 0.33 and 0.84.
Apart from the closure term, no error exceeded the relative term by more than
0.7 ``u``, and results that are normal numbers stayed below
``4.3(1+b)\,\varepsilon_{\mathrm{FT}}``. The underflow term is a worst case and
grows with ``N``. Closure ``F + B + A = 1`` holds within four
roundings, and the stored weights sum to one within two, for any tested ``N``.

`N` is a type parameter. The kernel compiles in about 10 ms for any `N`, runs
without allocation and costs about 10 ns per layer on one CPU core. The constructor
and `birth_layer_fractions` unroll over `N`: their first call compiles in well under
a second up to a few hundred layers, about 1 s at ``N = 1000`` and 2 to 3 s at
``N = 2000``. Tens to a few hundred layers are practical; on a GPU each thread
holds the `N` weights.

```julia
using CanopyOptics
FT = Float64
N = 8
depth = ntuple(k -> (k - FT(0.5)) / N, N)             # layer midpoints, 0 = front face

# Caller's layer yields (synthetic): quenched near the front, one dark layer.
# The allocation is depth-independent here, so the weights are P times the yield.
yield = map(z -> FT(0.012) + FT(0.02) * z, depth)
yield = Base.setindex(yield, zero(FT), 6)

escape = RegularPlateEscape(FT(0.8), FT(0.6))          # one emission bin
excitation = [ExponentialBirthProfile(a, one(FT), FT(0.2)) for a in FT[3, 1, 2]]

profiles = map(excitation) do p
    P = birth_layer_fractions(p, Val(N))               # excitation per layer, sums to 1
    w = map(*, P, yield)                               # production per layer
    (total = sum(w), birth = PiecewiseUniformBirthProfile(w))
end
totals = [x.total for x in profiles]                   # photons per allocated excitation photon
fates = leaf_emission_fates.(Ref(escape), [x.birth for x in profiles])

# With depth-independent allocation and yield the layered profile approaches the
# exponential one.
flat = PiecewiseUniformBirthProfile(birth_layer_fractions(excitation[1], Val(64)))
leaf_emission_fates(escape, flat).front - leaf_emission_fates(escape, excitation[1]).front  # ≈ -7.9e-6
```

All numbers are synthetic. Compared with the depth-independent case, the
front-quenched yields lower the front escape of the first excitation bin from
0.274 to 0.265; a yield falling with depth raises it. Multiplying `totals` by the
allocated excitation of each bin gives the produced photons ``q_{s,e}``.

### What this does and does not represent

- Equal-width layers with uniform production inside each layer are a
  **discretization of the source**, not a model of tissue layers. All layers share
  the plate's excitation and emission attenuation.
- **Dark layers.** A layer with zero yield or zero allocation gets ``W = 0``: it
  produces nothing but still attenuates and reabsorbs emission, because ``b`` is a
  property of the plate. The production ``q_{s,e}`` excludes it.
- **Dark sources.** If every layer of a source is dark in bin ``e``,
  ``q_{s,e} = 0`` and no profile can be formed (the constructor rejects all-zero
  weights). The caller skips the source in that bin or uses any placeholder profile;
  with zero production the placeholder does not affect any signal. This is the
  production side of the dark-column contract above: zero production is what makes
  a dark emission column admissible.
- **Unequal boundaries** (for example a palisade and a spongy layer of different
  thickness) need another representation; refine ``N`` until the boundaries fall on
  layer edges or use a different profile type.
- **Heterogeneous tissue transport** (layers with their own absorption, scattering
  or interfaces) needs layer optical properties and interface accounting, not birth
  weights. The two-zone [`DorsiventralPlateEscape`](@ref) does that for two zones;
  it has no method for this type, and its [`DorsiventralBirthProfile`](@ref) zone
  weights already carry zone-resolved production (zone excitation times zone
  allocation times zone yield, with depth-independent yield inside each zone).
- **Depth-independent emission spectrum.** Production of a source separates as
  ``q_{s,e}(z)\,p_{s,m,e}``: its conditional emission spectrum does not vary with
  depth. If it does (for example through different photosystem compositions in
  different tissues), split the source into components, each with its own spectrum,
  production and profile.
- The excitation density is the package's effective exponential closure;
  `birth_layer_fractions` integrates it exactly but does not make it an angular
  solution of light inside the leaf.
- No allocation, yield, rate constant or light response is supplied. If they depend
  on the local excitation, the caller computes that dependence (for example from
  ``Q_e P_{e,k}``) and any iteration it needs.
- Reabsorbed photons terminate; they are not re-emitted.

## Synthetic worked example

The following block is illustrative: all numbers are synthetic and none is a
recommended value. It runs as written.

```julia
using CanopyOptics
FT = Float64
emission = FT[670, 690, 710, 730, 750, 770]   # emission bin centres, nm
excitation = FT[450, 550, 650]               # excitation bin centres, nm
absorbed = FT[300, 120, 260]                 # leaf-absorbed photons per bin
coefficients = (chlorophyll = FT(0.8), carotenoids = FT(0.2))

# 1. Antenna excitation per bin, then fractions f of antenna excitation.
selection = CarotenoidAntenna(FT(0.7))
antenna = [antenna_excitation(selection,
               partition_absorption(HomogeneousAbsorption(), q, coefficients))
           for q in absorbed]
f = (psii = FT(0.55), psi = FT(0.4))
allocated = [allocate_excitation(q, f) for q in antenna]
sources = NamedTuple{keys(f)}(keys(f))

# 2. The caller applies its own photon yields.
yields = (psii = FT(0.02), psi = FT(0.004))
produced = map(n -> yields[n] .* [a.sinks[n] for a in allocated], sources)

# 3. One operator per source: its own spectrum and its own birth profile.
spectra = (psii = TabulatedPhotonSpectrum(emission, FT[9, 6, 4, 3, 2, 1]),
           psi  = TabulatedPhotonSpectrum(emission, FT[1, 2, 4, 6, 5, 3]))
escape = [RegularPlateEscape(b, FT(0.6)) for b in FT[2, 1.2, 0.6, 0.3, 0.2, 0.1]]
births = (psii = [ExponentialBirthProfile(a, one(FT), FT(0.2)) for a in FT[3, 1, 2]],
          psi  = [ExponentialBirthProfile(a, one(FT), FT(0.4)) for a in FT[3, 1, 2]])
operators = map(sources) do n
    p = reduce(hcat, [emission_probabilities(spectra[n], StrictStokes(), ex)
                      for ex in excitation])
    fates = leaf_emission_fates.(reshape(escape, :, 1), reshape(births[n], 1, :))
    ops = emission_operators(p, fates)
    @assert maximum(abs, emission_closure_residual(ops, p)) <= 8eps(FT)
    ops
end

# 4. A band with unit weights from 680 to 700 nm; sum the sources.
response = FT.(abs.(emission .- 690) .<= 10)
signals = map((op, q) -> band_signal(response, op.front, q), operators, produced)
front_band_photons = sum(values(signals))
```

To use fractions of absorbed photons instead, call
`allocate_excitation(absorbed[e], map(x -> h[e] * x, f))` with `h = antenna ./ absorbed`
from step 1. Passing `antenna` together with these fractions would count the
transfer loss twice.

## What is not provided

- No photon yields, rate constants or quenching states; produced photons are the
  caller's.
- No depth-dependent yield. A caller with depth-varying yield builds a
  [`PiecewiseUniformBirthProfile`](@ref) per source and excitation bin (regular
  plate only), as described above.
- No directional observation, instrument geometry or angular closure; all
  results are hemispheric face totals.
- No default emission spectrum, action spectrum, pigment content, antenna transfer
  fraction or allocation fraction.
- No conversion to energy; that uses emission wavelengths and is the caller's.
- No secondary fluorescence. Reabsorption terminates a tracked photon; re-emission
  of reabsorbed photons would need its own wavelength redistribution and yield
  accounting and is not part of these operators.

## Tests and extension

Tests cover Float32 and Float64; zero, one, three and forty sinks; reordered keys;
fraction sums at one, within and beyond the allowance, negative and nonfinite
values; invalid photons; closure; no promotion; zero allocations and inferred
types; agreement with `photosystem_excitation`; and agreement of the two allocation
conventions after conversion, using the package's partition and antenna selection.
Operator tests use two sources with different spectra, support rules
(`StrictStokes` and a finite Stokes band) and fates (regular and dorsiventral plates),
dark columns, per-column closure (including a perturbed operator entry), closure
against all produced photons, a demonstration of the photon loss when production
sits in a dark column, and an independent BigFloat accumulation. Band
signals are checked against a hand-computed value, for linearity in each argument,
for zero response and against the single-source path composed from existing
functions. CPU and CUDA batches cover shared and per-bin fractions, reordered and
empty batches, both accepted fate forms, matched permutations of excitation
columns and zero excitation columns. A third named source is added in the
tests without code changes.

Layered birth profiles are tested in both precisions for closure, inference, zero
allocation and `isbits`; against a 2048-bit closed form for ``10^{-12} \le b \le 700``,
up to 1000 layers and with dominant-weight blocks, and against Simpson quadrature for
a single lit layer among forty; for reduction to the uniform profile (one layer, or
equal weights), mirror symmetry under reversed weights, the limits ``b = 0``, small
``b`` (reabsorption slope ``(2-t)/(2t)``), large ``b``, ``t = 1`` and transmissions
down to the smallest subnormal (with ``b`` comparable to ``t``); against the accuracy
bound above for subnormal ``b``, subnormal and underflowing results and sparse deep
sources up to 1000 layers; for second-order
convergence to the exponential profile; for layer fractions against quadrature and
the exact integral; for constant, linearly falling and rising yields and a dark
layer; for invalid constructors; and in CPU and CUDA batches. An end-to-end test
runs two sources with depth-varying allocation and yield, two lit excitation bins
and one dark bin through `emission_operators` and `band_signal` and compares the
band signals with an explicit production-weighted sum over layers built from
independently derived spectral probabilities (photon density times bin width) and
excitation layer fractions (closed-form integral in BigFloat); it also checks
that both production conventions agree and that counting the allocation twice
does not.

API reference: [`allocate_excitation`](@ref), [`emission_operators`](@ref),
[`emission_closure_residual`](@ref), [`band_signal`](@ref),
[`PiecewiseUniformBirthProfile`](@ref) and [`birth_layer_fractions`](@ref).
