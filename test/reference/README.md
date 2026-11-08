# PROSPECT reference fixtures

Regenerate the three CSV files with Python and `mpmath==1.3.0`:

```sh
python test/reference/generate_prospect_reference.py
```

The generator sets precision explicitly and imports no production CanopyOptics
code. Julia tests read the committed text, so Python is not a test dependency.

- `prospect_reference.csv`: 78 combinations of optical depth and integer/fractional
  plate count, with transmission/reflection/material absorption and first
  derivatives. Integer stacks use repeated flux adding; fractional stacks use
  the eigendecomposition of a general 2×2 transfer matrix. Derivatives use
  independent high-precision finite differences, one-sided at zero depth and N=1.
  The exact nonabsorbing identity sets A=0 instead of retaining an arbitrary
  precision subtraction residual. E3 is independently checked by angular
  integration of the attenuation before generating the table.
- `prospect_fresnel_reference.csv`: direct angular quadrature of unpolarized
  Fresnel transmission at five indices spanning the native database range.
- `prospect_stack_reference.csv`: 750-digit stress calculations, including huge
  lossless stacks, large-N weak absorption, small simultaneous transmission and
  absorption, and the smallest positive Float32/Float64 transmission. These
  extreme fixtures validate values, not all possible extreme derivatives.

Reference values use 90 digits (750 for extreme stacks) internally and are written
with 45 significant digits. Julia converts through BigFloat so an unrepresentable
small reference value rounds to zero instead of raising a decimal-parser error.
Reference inputs outside an FT's representable domain are tested in the other FT;
no production result is clipped to make a reference comparison pass.

The E3 library's measured Float32 error and separate transmission/absorption
acceptance criteria are documented in `docs/src/pages/prospect_budget.md`.

## Internal regular-plate emission

`generate_plate_emission_reference.py` requires mpmath 1.3.0 and writes
`plate_emission_reference.csv`. It freezes inputs to Float32/Float64, integrates
the normalized source density numerically at 90 digits, and explicitly sums
successive reflection/absorption paths until the unprocessed probability is less
than 1e-85. The residual is reported, never normalized into an escape fate. The 28
cases cover equal excitation/emission depths, zero and weak loss, opaque material,
uniform and asymmetric birth, and mirrored amplitudes. This generator imports no
production code and does not use the production analytic overlap/escape formulas.

The wider Julia grid independently evaluates integrated exponentials and a two
boundary balance system at 2048 bits. It includes extreme finite depths and
smallest-positive interface transmission beyond the practical quadrature grid.

`generate_plate_preparation_reference.py` also uses mpmath1.3.0, importing the
independent Fresnel quadrature from `generate_prospect_reference.py`. It writes
140 precision-frozen optical cases to `plate_preparation_reference.csv`, with
280 source-depth integrals for distinct 90-degree and 40-degree entry conventions.
Traversal survival and direct loss are independent angular integrals; opaque
survival scales out `exp(-k)` before quadrature to retain relative resolution.
Record both Python sources when binding reference provenance. The Julia helper
rejects E3-library underflow even when the exact traversal probability would round
to a positive subnormal; no unresolved effective depth is silently manufactured.
The [preparation validation record](../../docs/src/data/plate_preparation_validation.json)
binds both generators, the fixture, production sources and completed package gates.
