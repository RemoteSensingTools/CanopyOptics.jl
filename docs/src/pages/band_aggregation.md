# Spectral Band Aggregation

Atmospheric radiative-transfer solvers (such as correlated-`k` shortwave and longwave
models) carry leaf reflectance and transmittance on coarse wavenumber bands in `cm⁻¹`,
whereas leaf optical models (`PROSPECT-PRO`, measured reflectance libraries) are
defined on fine wavelength grids. `CanopyOptics` provides a plain-array wavenumber
aggregation utility that conserves first-order reflected and transmitted band flux.

## Exact-Overlap Band Averaging

Given `N` piecewise-constant bin means `values` on `N + 1` strictly ascending
wavenumber bin edges `edges` (`cm⁻¹`) and a `(2, nbnd)` band-limit matrix `bands`
(`[\nu_{\mathrm{lo}}; \nu_{\mathrm{hi}}]` in `cm⁻¹`), [`band_average`](@ref) (and its
mutating counterpart [`band_average!`](@ref), which runs with zero heap allocation when
`dest`, input arrays, and built-in `weight` share a matching `Float32` or `Float64`
element type; `dest` must not alias `values`, `edges`, `bands`, `parent_band`, or
`weight`'s storage, while disjoint views of a shared buffer are permitted) evaluates

```math
\bar{v}_b = \frac{\int_{\nu_{\mathrm{lo},b}}^{\nu_{\mathrm{hi},b}} v(\nu)\,w(\nu)\,d\nu}{\int_{\nu_{\mathrm{lo},b}}^{\nu_{\mathrm{hi},b}} w(\nu)\,d\nu}
```

by exact intersection of every fine bin `[edges[i], edges[i + 1]]` with each band
`[\nu_{\mathrm{lo}}, \nu_{\mathrm{hi}}]`, rather than testing whether bin centers fall
inside the band.

- **Coverage contract**: If `[first(edges), last(edges)]` does not cover every band in
  `bands`, `band_average` throws an `ArgumentError` naming the uncovered band.
- **Partition invariance (including mixed coverage)**: Splitting a parent band ``B``
  into contiguous sub-bands ``B_1, \dots, B_m`` preserves the parent band value across
  every weight and coverage regime when both sub-band means and sub-band shares resolve
  the parent band's effective weight:

```julia
v_sub = band_average(v, edges, sub_bands, weight; parent_band = parent_band)
shares = subband_weights(sub_bands, parent_band, weight)
# sum(shares .* v_sub) ≈ only(band_average(v, edges, parent_band, weight))
```

  Passing `parent_band` to `band_average` / `band_average!` ensures that if a
  [`TabulatedWeight`](@ref) covers only part of `parent_band`, the parent's fallback to
  [`FlatWeight`](@ref) is resolved once and applied consistently to every child sub-band
  as well as to [`subband_weights`](@ref).

## Spectral Weights, `ReferenceSolar`, and Dark-Band Rules

Three built-in [`AbstractSpectralWeight`](@ref) types supply ``w(\nu)`` (per `cm⁻¹`):

- [`FlatWeight()`](@ref): uniform per-`cm⁻¹` weighting (``w(\nu) = 1``).
- [`PlanckWeight(T)`](@ref): hemispheric blackbody spectral exitance
  ``M_\nu(\nu, T) = \pi B_\nu(\nu, T)`` at temperature `T` (`1 <= T <= 10^8 K`),
  integrated in relative coordinates across each fine-bin overlap via composite 8-point
  Gauss–Legendre quadrature.
- [`TabulatedWeight(edges, w)`](@ref) and [`ReferenceSolar`](@ref): piecewise-constant
  per-`cm⁻¹` spectral irradiance tables.

### `ReferenceSolar` CSV Format and UV Zero-Irradiance Convention

`ReferenceSolar(grid, irradiance; grid_unit = :nm)` constructs a `TabulatedWeight` from
in-memory arrays. The zero-argument `ReferenceSolar()` (or `ReferenceSolar(path)`) loads
`data/reference_solar_global.csv` if installed, and throws an `ArgumentError` when no
default table is present.

A solar CSV file uses UTF-8 comma-separated values with optional `#` comment lines, an
optional single header row with non-numeric text in both columns
(`wavelength_nm,irradiance_W_m2_nm`), and two numeric columns:
1. `wavelength_nm`: strictly ascending wavelengths in `nm`.
2. `irradiance_W_m2_nm`: clear-sky global (direct + diffuse) surface spectral irradiance
   in ``\mathrm{W\,m^{-2}\,nm^{-1}}`` (`>= 0`).

Include an explicit zero-irradiance row below `280 nm` (for example `195.0,0.0` paired
with `280.0,0.0` at the stratospheric ozone cutoff) so UV bands extending below `280 nm`
are covered by the table rather than treated as out-of-domain:
- A partially illuminated UV band such as `263–345 nm` (`28985.5–38022.8 cm⁻¹`) is
  covered with positive solar weight and integrates only its illuminated `> 280 nm`
  portion (`0` weight below `280 nm`).
- A fully dark UV band such as `200–263 nm` (`38022.8–50000.0 cm⁻¹`) has zero physical
  solar flux (`CanopyOptics.band_weight` returns `0.0`, and `subband_weights` assigns
  share `0.0` inside an illuminated parent band). Because the solar-weighted mean `0/0`
  is undefined when total band solar weight is zero, `band_average` explicitly returns
  the `FlatWeight()` wavenumber mean over that dark band as a finite placeholder optical
  property (and `subband_weights` likewise returns `FlatWeight()` width shares only when
  the entire `parent_band` has zero total weight).

## Wavelength-to-Wavenumber Conversion, `|dλ/dν|` Jacobian, and Native-Grid Weighting

[`wavelength_to_wavenumber`](@ref) maps point samples on a regular strictly ascending
wavelength grid in `nm` (such as `400.0:1.0:2500.0`) to ascending wavenumber bin edges
`10^7 / (λ_k ∓ Δλ/2)` in `cm⁻¹` ([`wavenumber_to_wavelength`](@ref) inverts both modes):

- Intensive quantities (`density = false`, e.g. leaf ``\rho(\lambda), \tau(\lambda)``)
  reverse bin order without rescaling.
- Per-`nm` spectral densities (`density = true`, e.g. solar irradiance ``E_\lambda`` in
  ``\mathrm{W\,m^{-2}\,nm^{-1}}``) multiply each bin by the exact bin-averaged Jacobian
  ``|d\lambda / d\nu|`` over ``[\lambda_{\mathrm{lo}}, \lambda_{\mathrm{hi}}]``:

```math
\bar{E}_{\nu,i} = E_{\lambda,k}\,\frac{\lambda_{\mathrm{lo}}\,\lambda_{\mathrm{hi}}}{10^7},
\qquad \left|\frac{d\lambda}{d\nu}\right| = \frac{10^7}{\nu^2}.
```

- **Whole-bin vs. sub-bin flux conservation**: The Jacobian
  ``\lambda_{\mathrm{lo}}\lambda_{\mathrm{hi}} / 10^7`` preserves **whole-bin** integrated
  flux ``\bar{E}_{\nu,i}\,\Delta\nu_i = E_{\lambda,k}\,\Delta\lambda_k`` exactly (including
  on irregular midpoint-defined `ReferenceSolar` grids). Within each converted bin,
  however, ``\bar{E}_{\nu,i}`` is stored as a constant density per `cm⁻¹`, so a band
  boundary that cuts **through the interior** of a bin integrates that constant per-`cm⁻¹`
  density linearly in ``\nu``. For example, a single uniform `400–800 nm`
  (`12500–25000 cm⁻¹`) bin at `1 W m⁻² nm⁻¹` has total flux `400 W m⁻²` and converted
  density `0.032 W m⁻² (cm⁻¹)⁻¹`; integrating over `600–800 nm` (`12500–16666.67 cm⁻¹`,
  one-third of the wavenumber width) yields `133.33 W m⁻²` rather than the `200 W m⁻²` of
  a constant per-`nm` representation.
- **Weight on the native grid before coarsening**: Each `1 nm` PROSPECT sample is treated
  as piecewise constant across its bin, and resampling via `leaf_spectrum` preserves the
  flat-in-wavenumber (`FlatWeight`) integral. Always apply the physical spectral weight
  (`ReferenceSolar`, `PlanckWeight`, or flat-in-wavelength
  `ReferenceSolar(λ, ones(length(λ)))`) on the native fine (`~1 nm`) grid before
  coarsening. Because ``|d\lambda / d\nu| = 10^7 / \nu^2`` weights longer wavelengths more
  heavily than flat wavenumber, on the default `LeafProspectProProperties()` leaf across
  the `625–778 nm` (`12850–16000 cm⁻¹`) red-edge band, `FlatWeight()` gives
  ``\rho / \tau = 0.1786 / 0.1970`` whereas flat-in-wavelength weighting on the native
  `1 nm` grid gives ``\rho / \tau = 0.1968 / 0.2195``.

## Leaf Spectrum Composition (`CompositeLeaf`)

[`leaf_spectrum(model, edges)`](@ref) evaluates bin-mean `(rho, tau)` on `edges` and
throws an `ArgumentError` if `model` does not cover every bin of `edges`:

- [`LeafProspectProProperties`](@ref) evaluates [`prospect_budget`](@ref) on the
  bundled `1 nm` optical table (`399.5–2500.5 nm`, `3999.20–25031.29 cm⁻¹`) and
  resamples onto `edges` by exact wavenumber overlap.
- [`TabulatedLeaf(edges, rho, tau)`](@ref) resamples measured wavenumber spectra.
- [`ConstantLeaf(rho, tau; ν_min, ν_max)`](@ref) supplies constant regional values.
- [`CompositeLeaf(models...)`](@ref) selects **the first model covering each entire bin**.
  Custom leaf models can specialize `CanopyOptics.spectral_domain(model)` for interval
  domains or throw `CanopyOptics.SpectralCoverageError` (or return `nothing`) for
  uncovered single-bin probes without masking genuine `ArgumentError` exceptions.

When splicing UV, PROSPECT VNIR/SWIR, mid-IR (`2.5–4.0 µm`), and thermal (`> 4.0 µm`)
regions, bound each regional `ConstantLeaf` to its physical wavenumber interval and
insert the regional boundary wavenumbers into the fine grid so no bin straddles a
boundary:

```julia
leaf_pro = LeafProspectProProperties()
ν_pro_lo, ν_pro_hi = CanopyOptics.spectral_domain(leaf_pro) # 3999.20, 25031.29 cm⁻¹

# 400 nm PROSPECT optics extended across the UV (λ < 399.5 nm, ν >= ν_pro_hi)
r_400, t_400 = leaf_spectrum(leaf_pro, wavelength_to_wavenumber(400.0:1.0:401.0))
uv_leaf   = ConstantLeaf(r_400[end], t_400[end]; ν_min = ν_pro_hi)
mwir_leaf = ConstantLeaf(0.04, 0.0; ν_min = 2500.0, ν_max = ν_pro_lo) # 2.5–4.0 µm
tir_leaf  = ConstantLeaf(0.03, 0.0; ν_min = 0.0, ν_max = 2500.0)      # > 4.0 µm

broadband_leaf = CompositeLeaf(leaf_pro, uv_leaf, mwir_leaf, tir_leaf)
fine_edges = sort!(unique!(vcat(
    collect(10.0:10.0:2500.0),
    [2500.0, ν_pro_lo],
    wavelength_to_wavenumber(400.0:1.0:2500.0),
    [ν_pro_hi],
    collect(25050.0:50.0:50000.0),
)))
rho, tau = leaf_spectrum(broadband_leaf, fine_edges)
```

## PROSPECT Pigment Absorption Shares

[`leaf_pigment_shares(leaf, edges)`](@ref) returns `(; chlorophyll, carotenoids)` shares
of total leaf material absorption per bin for `LeafProspectProProperties` via
[`HomogeneousAbsorption`](@ref). At each native `1 nm` bin, the constituent weight is
concentration × specific absorption, with chlorophyll, carotenoids, anthocyanins, brown
material, water, bulk dry matter, protein, and carbon constituents all included in the
denominator ``\sum_j C_j K_j(\lambda)``.

- **Zero-absorption convention**: Native bins where all constituent products vanish
  (`AbsorptionNoMaterial`) return `0.0` for both pigment shares (rather than `NaN`). When
  resampling onto a multi-bin interval whose integrated material absorptance is zero while
  constituent coefficients are nonzero (subnormal optical-depth underflow), the resampled
  share falls back to the flat-in-wavenumber mean of ``C_i K_i / \sum_j C_j K_j``.
- **Use directly computed material absorptance on the native grid**: Obtain fine-grid
  absorptance from `prospect_budget(leaf, LeafOpticalProperties()).absorptance` rather
  than `1 - rho - tau`. In weakly absorbing bins, `1 - rho - tau` can round to `0.0`
  below the spacing of `rho` and `tau` while direct plate absorptance retains full
  relative precision (for example `3.13 × 10⁻¹⁷` at `675 nm` at trace pigment content).
- **Coarse-band pigment shares (energy vs. photon-count)**: A coarse-band pigment *share*
  is the weighted mean of pigment absorptance ``\alpha(\nu)\,f_{\mathrm{pig}}(\nu)``
  **divided by** the weighted mean of total leaf absorptance ``\alpha(\nu)``:

```julia
opt = LeafOpticalProperties()
edges, α = wavelength_to_wavenumber(400.0:1.0:2500.0, prospect_budget(leaf, opt).absorptance)
shares = leaf_pigment_shares(leaf, edges)

# Energy-weighted coarse pigment shares under weight w (per cm⁻¹):
α_band     = band_average(α, edges, par_bands, w)
α_chl_band = band_average(α .* shares.chlorophyll, edges, par_bands, w)
α_car_band = band_average(α .* shares.carotenoids, edges, par_bands, w)
f_chl_band = ifelse.(α_band .> 0, α_chl_band ./ α_band, zero.(α_band))
f_car_band = ifelse.(α_band .> 0, α_car_band ./ α_band, zero.(α_band))
```

  For **photon-count** shares (APAR / antenna excitation), convert the piecewise-constant
  energy weight density ``w(\nu)`` (`W m⁻² (cm⁻¹)⁻¹`) to a photon-flux density by
  dividing by photon energy ``100\,h c\,\nu`` (`J`, with
  ``h = 6.62607015\times 10^{-34}\;\mathrm{J\,s}`` and
  ``c = 299792458\;\mathrm{m\,s^{-1}}``; at `550 nm`, `1 J` equals `4.598 µmol` photons).
  Over any sub-interval ``[\nu_{\mathrm{lo}}, \nu_{\mathrm{hi}}]`` (`\nu_{\mathrm{lo}} > 0`)
  on which ``w(\nu) = w_k`` is constant, integrating ``w_k / (100\,h c\,\nu)``
  analytically yields the exact bin-mean photon-flux density
  ``\bar{q}_{\nu,k} = (w_k / (100\,h c))\,\ln(\nu_{\mathrm{hi}} / \nu_{\mathrm{lo}}) / (\nu_{\mathrm{hi}} - \nu_{\mathrm{lo}})``.
  Evaluating ``\bar{q}_{\nu,k}`` on the merged knots of `w.edges`, `edges`, and
  `par_bands` integrates ``1/\nu`` across every fine leaf bin even when `w` has coarse
  bins (for example, on a single uniform `12500–25000 cm⁻¹` weight bin split at
  `18750 cm⁻¹`, exact ``1/\nu`` integration gives share ``\ln(1.5)/\ln(2) \approx 0.585``
  rather than `0.500`), and guarding `α_photon .> 0` preserves the zero-absorption
  convention:

```julia
h_planck, c_light = 6.62607015e-34, 299792458.0
w_lo, w_hi = first(w.edges), last(w.edges)
ν_u = sort!(unique!(vcat(
    w.edges,
    filter(ν -> w_lo < ν < w_hi, edges),
    filter(ν -> w_lo < ν < w_hi, vec(par_bands)),
)))
w_u = band_average(w.w, w.edges, vcat(ν_u[1:end-1]', ν_u[2:end]'), FlatWeight())
q_u = [
    iszero(wk) ? zero(wk) :
    (wk / (100.0 * h_planck * c_light)) * log(ν_u[k + 1] / ν_u[k]) / (ν_u[k + 1] - ν_u[k])
    for (k, wk) in enumerate(w_u)
]
w_photon     = TabulatedWeight(ν_u, q_u)
α_photon     = band_average(α, edges, par_bands, w_photon)
α_chl_photon = band_average(α .* shares.chlorophyll, edges, par_bands, w_photon)
α_car_photon = band_average(α .* shares.carotenoids, edges, par_bands, w_photon)
f_chl_photon = ifelse.(α_photon .> 0, α_chl_photon ./ α_photon, zero.(α_photon))
f_car_photon = ifelse.(α_photon .> 0, α_car_photon ./ α_photon, zero.(α_photon))
```
