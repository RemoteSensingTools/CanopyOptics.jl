# Pigment allocation validation, 2026-11-08

The optical component now partitions absorption among arbitrary named constituents,
selects chlorophyll-only or explicitly weighted carotenoid antenna excitation, and
applies an independent photosystem fraction. It introduces no physiology, canopy
state, timestepper or database dependency into the scalar partition/selection rules.
A separate host helper shares PROSPECT composition products with the optical solver.

Both independent code/bugs and physics reviewers cleared the final implementation.
The package suite passes **102,469 checks**, including **12,168 new checks**. Explicit
NVIDIA CUDA acceptance passes **6,736 checks** with scalar indexing disabled, in
both Float32 and Float64. The complete documentation build passes (the existing
aggregate search-index size warning remains).

The independent reference integrates species absorption rates along three weighted
optical paths using BigFloat power-series integration. Other checks cover invalid
or absent material, missing pigment keys, mixed precision, extreme coefficient
magnitudes, endpoint transfer fractions, inference/allocations, custom absorber and
selector extension, empty/reordered device batches and matching spectral reductions.

During review, photon-only diagnostics proved insufficient at the smallest positive
photon input: several constituent photon amounts can all round to zero. The result
therefore retains unitless fractions as well. Selection uses these probabilities
before photon multiplication, preserving zero-selected-pigment, half-selected and
all-selected limits without clipping or normalization. No-material fractions are
undefined (NaN), while actual material absorption is zero. Extreme coefficient
ratios can still underflow the fractions; no relative accuracy beyond the selected
FT's resolution is claimed.

The real-data fixture uses three explicit PROSPECT-PRO compositions, native spectral
coefficients, total direct absorptance and matching constituent products. It checks
wavelength-specific photon and energy accounting over the declared PAR bins, with
independently configured carotenoid and PSII fractions. The same fixture runs on
CPU and CUDA after optical tables are prepared on the host. All 12 original optical
profile/control spectra remain bitwise unchanged by the shared composition helper.

Parameters and units are in `test/parameters/pigment_absorption.toml`; equations and
API boundaries are in [the component documentation](src/pages/pigment_excitation.md).
Exact source, test, parameter and environment hashes are recorded in
`pigment_excitation_validation.json`.

This is the optical half of the pigment-to-physiology handoff. Wiring explicit Q2
through LeafPhotosynthesis and adding a combined teaching example remain the next
increment. No calibrated SIF escape matrix, mesophyll segregation model, Metal
execution, or AD of the new partition/selection kernels is claimed.
