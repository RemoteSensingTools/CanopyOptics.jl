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
yield or stoichiometry varies with depth) gets its own fates through the existing
birth-profile argument of [`leaf_emission_fates`](@ref). No source container type is
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
- No depth-dependent source profiles beyond the existing birth-profile types; a
  depth-varying yield is represented by the caller's birth profile per source.
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

API reference: [`allocate_excitation`](@ref), [`emission_operators`](@ref),
[`emission_closure_residual`](@ref) and [`band_signal`](@ref).
