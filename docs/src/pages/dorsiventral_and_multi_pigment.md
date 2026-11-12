# Dorsiventral anatomy, multi-pigment bases and two-zone fluorescence

Homogeneous plate models represent a leaf as a single symmetric stack with one
pooled chlorophyll ``a{+}b`` spectrum and uniform pigment depth distribution.
Dicotyledonous leaves violate these assumptions in three distinct ways:

1. **Dorsiventral mesophyll & cuticle asymmetry**: Dense columnar palisade
   parenchyma beneath a smooth, waxy adaxial epidermis concentrates ``65``–``85\%``
   of leaf chlorophyll, whereas lacunose spongy mesophyll above a rougher
   stomatal abaxial epidermis increases internal air–cell-wall scattering,
   making ``R_{\mathrm{ab}} > R_{\mathrm{ad}}`` in the visible while preserving
   ``T_{\mathrm{ad}} = T_{\mathrm{ab}}`` under hemispherical illumination.
2. **Distinct chlorophyll ``a`` and ``b`` absorption bands**: Chlorophyll ``a``
   (Soret near ``434\text{ nm}``, ``Q_y`` near ``668``–``676\text{ nm}``) and
   chlorophyll ``b`` (Soret near ``466\text{ nm}``, ``Q_y`` near ``650\text{ nm}``)
   differ across blue-green and red-edge wavelengths, affecting both pigment
   retrievals and photosystem excitation partitioning.
3. **Depth-asymmetric fluorescence reabsorption**: Because palisade tissue is
   optically dense at the ``685\text{ nm}`` PSII emission peak and nearly
   transparent at ``740\text{ nm}``, adaxial vs abaxial excitation produces
   strongly wavelength-dependent front/back escape ratios.

## Continuous `N_tot >= 1` dorsiventral and asymmetric-surface leaves

[`DorsiventralLeafProperties`](@ref) couples a palisade sublayer and a spongy
sublayer between independent adaxial and abaxial [`AbstractLeafSurface`](@ref)
boundaries via exact Redheffer star-product adding, returning per-face
reflectance, transmittance, total absorptance, and separate zone absorptances
(`absorptance_palisade`, `absorptance_spongy`) through [`dorsiventral_budget`](@ref).

Because empirical leaf inversions frequently encounter compact sun leaves or
monocots with total effective structure parameter ``1 \le N_{\mathrm{tot}} < 2``,
[`dorsiventral_from_traits`](@ref) and [`DorsiventralLeafProperties`](@ref)
bridge smoothly across ``N_{\mathrm{tot}} = 2``:

- For ``N_{\mathrm{tot}} \ge 2``: ``N_{\mathrm{pal}} = 1 + (N_{\mathrm{tot}} - 2) f_m``
  and ``N_{\mathrm{spo}} = 1 + (N_{\mathrm{tot}} - 2)(1 - f_m)``, evaluated via
  the two-zone fractional Allen stack.
- For ``1 \le N_{\mathrm{tot}} < 2``: convex interpolation
  ``(2 - N_{\mathrm{tot}}) \mathcal{B}_{\mathrm{slab}} + (N_{\mathrm{tot}} - 1) \mathcal{B}_{1,1}``
  between the unbroken two-sublayer slab ([`single_slab_dorsiventral_budget`](@ref))
  at ``N_{\mathrm{tot}} = 1`` and the ``(N_{\mathrm{pal}}, N_{\mathrm{spo}}) = (1, 1)``
  two-plate stack at ``N_{\mathrm{tot}} = 2``.

This bridge is ``C^0`` continuous at ``N_{\mathrm{tot}} = 2`` (``C^\infty`` on
``(1, 2)`` and ``(2, \infty)``), strictly passive (``R + T + A = 1``,
``A_{\mathrm{pal}} + A_{\mathrm{spo}} = A``), `ForwardDiff`-compatible, and
preserves exact bi-hemispherical reciprocity ``T_{\mathrm{ad}} = T_{\mathrm{ab}}``
to machine precision. [`AsymmetricSurfaceLeafProperties`](@ref) applies the same
continuous ``N \ge 1`` bridge when only the adaxial and abaxial cuticles differ.

## Multi-pigment optical bases and dry-matter conventions

[`LeafMultiPigmentProperties`](@ref) and [`MultiPigmentOpticalProperties`](@ref)
separate ``C_{\mathrm{chla}}`` and ``C_{\mathrm{chlb}}`` while retaining
carotenoids, anthocyanins, senescent brown pigments, water, and either bulk dry
matter (`Cm`) or PROSPECT-PRO protein + carbon-based constituents (`Cprot`, `Ccbc`).

`MultiPigmentOpticalProperties(λ_bnds; basis = ...)` constructs [`MultiPigmentOpticalProperties`](@ref)
on any wavelength grid with four selectable basis options:

| `basis` | Pigment spectra (`400–800 nm`) | Units |
| --- | --- | --- |
| `:calibrated_in_vivo` | Empirically separated in-vivo `Kchla`, `Kchlb` (`ANGERS` species-grouped CV) + `OPTI_2021` `Kcar`, `Kant` | `cm² µg⁻¹` |
| `:prospect_mp_2017` | Zhang et al. (2017) `PROSPECT-MP` pseudo-Voigt `Kchla`, `Kchlb`, `Kcar` + `OPTI_2021` `Kant` | `cm² µg⁻¹` |
| `:prospect_mp_plus_2022` | Zhang et al. (2022) `PROSPECT-MP+` pseudo-Voigt `Kchla`, `Kchlb`, `Kcar`, `Kant` | `cm² µg⁻¹` |
| `:prospect_pro_pooled` | Baseline `Kchla == Kchlb == Kcab` from `OPTI_2021` | `cm² µg⁻¹` |

Above ``800\text{ nm}``, all four bases splice cleanly onto `OPTI_2021` so that
unbounded pseudo-Voigt Lorentzian wings cannot leak artificial pigment
absorption into the ``800``–``2500\text{ nm}`` NIR/SWIR plateau. For excitation
partitioning, [`MultiChlorophyllAntenna`](@ref) weights ``C_{\mathrm{chla}}``,
``C_{\mathrm{chlb}}`` and ``C_{\mathrm{car}}`` separately in
[`photosystem_excitation`](@ref).

### Dry-matter mutual exclusivity and anthocyanin mass units

All three leaf property types ([`LeafProspectProProperties`](@ref),
[`LeafMultiPigmentProperties`](@ref), [`LeafSIPProperties`](@ref)) enforce two
input-consistency safeguards on keyword construction:

1. **No silent `Cm` + `Cprot`/`Ccbc` double-counting**: Passing `Cm` explicitly
   (PROSPECT-D bulk dry-matter mode, e.g. `LeafProspectProProperties(Cm = 0.006)`)
   automatically defaults omitted `Cprot` and `Ccbc` to `0` instead of adding the
   default `0.001 + 0.009 = 0.010 g cm⁻²` PROSPECT-PRO dry-matter background on
   top of `Cm`. Passing positive `Cm` together with positive `Cprot` or `Ccbc`
   throws an `ArgumentError`.
2. **Consistent `[µg cm⁻²]` anthocyanin units**: In `OPTI_2021` (`dataSpec_PRO.csv`,
   Féret et al., 2017/2021), `Kant` peaks at ``0.0725\text{ cm}^2\,\mu\text{g}^{-1}``
   near ``535\text{ nm}``; `Canth` is therefore specified in mass units
   `[µg cm⁻²]` across all `CanopyOptics` leaf models.

## Two-zone internal fluorescence escape

[`DorsiventralPlateEscape`](@ref) and [`DorsiventralBirthProfile`](@ref) extend
[`leaf_emission_fates`](@ref) from a single homogeneous plate
([`RegularPlateEscape`](@ref)) to two coupled mesophyll zones (`adaxial_depth`
``b_1``, `abaxial_depth` ``b_2``, outer diffuse exit transmissions ``t_1, t_2 \in (0, 1]``,
and palisade–spongy internal boundary transmission ``t_{\mathrm{int}} \in (0, 1]``).

```julia
using CanopyOptics

dpo_exc = dorsiventral_plate_optics(2.4, 0.8, 1.42; internal_air_gap = true)
dpo_685 = dorsiventral_plate_optics(1.8, 0.6, 1.41; internal_air_gap = true)
birth   = DorsiventralBirthProfile(
    0.85, 0.15,
    dpo_exc.adaxial_effective_depth, 1.0, 0.1,
    dpo_exc.abaxial_effective_depth, 1.0, 0.2,
)
fates = leaf_emission_fates(dpo_685.escape, birth)
```

Local dimensionless depth ``\zeta_i \in [0, 1]`` runs anatomically **adaxial to
abaxial** in both zones (``\zeta_1 = 0`` at the adaxial epidermis to ``\zeta_1 = 1``
at the internal boundary; ``\zeta_2 = 0`` at the internal boundary to ``\zeta_2 = 1``
at the abaxial epidermis). Reversing leaf orientation therefore swaps
``(w_1, a_1, f_1, g_1) \leftrightarrow (w_2, a_2, g_2, f_2)`` and
`front` ``\leftrightarrow`` `back` while leaving `reabsorbed` invariant.

When ``t_{\mathrm{int}} = 1`` and ``t_1 = t_2 = t``, splitting a homogeneous
plate ``b = \eta b + (1-\eta)b`` at any ``\eta \in (0, 1)`` reproduces
[`RegularPlateEscape`](@ref)`(b, t)` to floating-point roundoff. Setting
`internal_air_gap = true` in [`dorsiventral_plate_optics`](@ref) uses the
multiple-reflection diffuse mesophyll–air–mesophyll transmittance
``t_{\mathrm{gap}} = t_{90} / [n_r^2 (2 - t_{90})]``, capturing partial internal
back-reflection at the spongy–palisade intercellular air interface. Both
[`DorsiventralPlateEscape`](@ref) and [`DorsiventralBirthProfile`](@ref) are
`isbits`, allocation-free, subnormal-safe across ``(0, 1]^3``, and execute
identically on CPU and CUDA (`Float32`/`Float64`).
