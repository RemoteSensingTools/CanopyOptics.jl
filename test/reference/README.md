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
