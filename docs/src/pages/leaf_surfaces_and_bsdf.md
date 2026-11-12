# Dielectric leaf surfaces and energy-conserving BSDF

Standard plate models such as PROSPECT-5, PROSPECT-D and PROSPECT-PRO evaluate
external surface entry with a fixed ``40^\circ`` solid-angle cone
(`calctav(40, nr)`) and assume bi-Lambertian outgoing reflectance and
transmittance. That simplification introduces two physical limitations when
comparing against directional or sphere-integrated measurements:

1. **Observation-geometry mismatch**: Integrating spheres use ``0^\circ`` or
   ``8^\circ`` collimated entry, whereas the legacy ``40^\circ`` cone lowers
   smooth-surface Fresnel reflectance by roughly ``0.004``–``0.006`` across the
   visible spectrum.
2. **Non-conserving directional specular additions**: Adding a microfacet BRDF
   term on top of an unmodified ``40^\circ`` PROSPECT diffuse budget (as in
   uncoupled diffuse+specular overlays or `PROSPECULAR` at oblique solar zenith
   angles) double-counts boundary reflection and violates energy conservation
   (``R(\theta_i) + T(\theta_i) > 1`` in weakly absorbing NIR bands at
   ``\theta_i \gtrsim 55^\circ``).

`CanopyOptics` separates the dielectric boundary operator ([`AbstractLeafSurface`](@ref))
from the internal mesophyll solver and couples both into an energy-conserving,
Helmholtz-reciprocal 4×4 polarized BSDF ([`CoupledLeafBSDF`](@ref)).

## Dielectric boundary models

Two boundary parameterizations subtype [`AbstractLeafSurface`](@ref):

- [`ConeSurface`](@ref)`(α_deg, n_scale)`: analytical Stern–Allen cone-averaged
  smooth Fresnel interface ([`calctav_smooth`](@ref)), continuous down to exact
  normal incidence ``\alpha = 0^\circ`` where
  ``t_0(n) = 4n/(n+1)^2``.
- [`RoughDielectricSurface`](@ref)`(σ, n_scale)`: isotropic Beckmann microfacet
  slope distribution ([`beckmann_D`](@ref)) with separable Smith shadowing–masking
  ([`smith_beckmann_G1`](@ref)), integrated over ``(\mu_o, \Delta\phi)`` via
  Gauss–Legendre quadrature ([`microfacet_directional_reflectance`](@ref),
  [`microfacet_hemispherical_reflectance`](@ref)).

For any external boundary model with hemispherical external transmittance
``t_{\mathrm{hem,ext}}(n_{\mathrm{eff}}, \sigma)`` into a lossless thin cuticle
of effective index ``n_{\mathrm{eff}} = n_{\mathrm{scale}} \, n_{\mathrm{surface}}(\lambda)``
overlying an aqueous mesophyll of index ``n_r(\lambda)``, isothermal detailed
balance (étendue conservation between air ``n=1`` and mesophyll ``n_r``) fixes
the diffuse internal exit transmittance without free parameters:

```math
t_{\mathrm{int,diffuse}} = \frac{t_{\mathrm{hem,ext}}(n_{\mathrm{eff}}, \sigma)}{n_r^2} \in (0, 1),
\qquad
r_{\mathrm{int,diffuse}} = 1 - t_{\mathrm{int,diffuse}}.
```

[`surface_interface_optics`](@ref) evaluates external entry and diffuse internal
exit consistently for [`prospect_budget`](@ref), [`leaf_sip_budget`](@ref) and
[`dorsiventral_budget`](@ref).

### Epicuticular wax refractive-index dispersion

Aqueous mesophyll refractive index ``n_r(\lambda)`` in `OPTI_2021` drops to
``1.27``–``1.34`` across ``1900``–``2500\text{ nm}``, whereas epicuticular wax
cuticles remain near ``1.40``–``1.48`` across ``400``–``2500\text{ nm}``
(Stuckens et al., 2009). [`load_wax_refractive_index`](@ref) averages the
tabulated Stuckens et al. (2009) normal dispersion curve onto any wavelength
grid so it can be supplied via `n_surface` (or per-face `n_surface_adaxial` and
`n_surface_abaxial`):

```julia
using CanopyOptics, Unitful

λ_bnds = (399.5:1.0:2500.5)u"nm"
opti   = LeafOpticalProperties(λ_bnds)
n_wax  = load_wax_refractive_index(λ_bnds)

leaf = LeafProspectProProperties(N = 1.6, Ccab = 42.0, Cw = 0.013, Cm = 0.006)
surf = RoughDielectricSurface(0.16, 1.0)
bd   = prospect_budget(leaf, opti, surf; n_surface = n_wax)
```

In this thin lossless cuticle approximation, smooth internal wax–mesophyll
Fresnel reflection ``(n_{\mathrm{eff}} - n_r)^2 / (n_{\mathrm{eff}} + n_r)^2``
is neglected (``\le 0.05\%`` in the VNIR and ``\le 0.5\%`` at ``2450\text{ nm}``,
with zero external-to-mesophyll total internal reflection whenever ``n_r > 1``),
while ``t_{\mathrm{int,diffuse}} = t_{\mathrm{hem,ext}}(n_{\mathrm{eff}}, \sigma) / n_r^2``
preserves exact passive energy conservation and bi-hemispherical reciprocity
``T_{\mathrm{ad}} = T_{\mathrm{ab}}`` even when only one face carries wax.

## Coupled energy-conserving leaf BSDF

[`CoupledLeafBSDF`](@ref) and [`prepare_leaf_bsdf`](@ref) couple single-bounce
polarized microfacet reflection ([`surface_brdf_mueller`](@ref)) directly to the
angle-dependent entry transmittance
``t_{\mathrm{ext}}(\mu_i) = 1 - R_{\mathrm{surf}}(\mu_i)`` of the underlying
single-zone or two-zone [`DorsiventralLeafProperties`](@ref) mesophyll stack.

For an incident ray with cosine ``\mu_i`` on face ``s \in \{\text{adaxial}, \text{abaxial}\}``:

1. Single-scattering surface reflection carries directional-hemispherical weight
   ``R_{\mathrm{surf},s}(\mu_i)`` and 4×4 Mueller BRDF ``\mathbf{M}_{\mathrm{surf},s}(\Omega_i, \Omega_o)``
   rotated smoothly into the incoming and outgoing meridian-plane Stokes bases.
2. Only the transmitted fraction ``t_{\mathrm{ext},s}(\mu_i) = 1 - R_{\mathrm{surf},s}(\mu_i)``
   enters the interior stack, where multiple scattering between the mesophyll
   and the diffuse inner boundaries ``r_{\mathrm{int},s}`` produces escaping
   interior reflection ``R_{\mathrm{vol},s}(\mu_i)``, transmission
   ``T_{\mathrm{vol},s}(\mu_i)``, and material absorptance ``A_s(\mu_i)``
   satisfying ``R_{\mathrm{surf},s}(\mu_i) + R_{\mathrm{vol},s}(\mu_i) + T_{\mathrm{vol},s}(\mu_i) + A_s(\mu_i) = 1``
   identically for all ``\mu_i \in (0, 1]`` ([`leaf_directional_budget`](@ref)).
3. Interior reflection escapes across face ``s`` with dimensionless normalized
   escape profile ``E_s(\mu) = t_{\mathrm{ext},s}(\mu) / t_{\mathrm{hem,ext},s}``,
   giving interior BRDF contribution
   ``f_{r,\mathrm{vol},s}(\Omega_i \to \Omega_o) = \frac{R_{\mathrm{vol,hem},s}}{\pi} E_s(\mu_i) E_s(\mu_o)``
   and guaranteeing exact directional BRDF reciprocity
   ``f_r(\Omega_i \to \Omega_o) = f_r(\Omega_o \to \Omega_i)`` and 4×4 Hovenier
   reciprocity
   ``\mathbf{M}(\Omega_2, \Omega_1) = \mathbf{P}\,\mathbf{M}(\Omega_1, \Omega_2)^\top \mathbf{P}``
   with ``\mathbf{P} = \operatorname{diag}(1, 1, -1, 1)``.
4. Transmitted radiation splits between the reciprocal diffuse profile
   ``E_{s'}(\mu_o) / \pi`` (weight ``1 - \eta_t``) and a normalized Henyey–Greenstein-like
   forward-peaked mesophyll lobe (weight ``\eta_t \in [0, 1)``, sharpness
   ``\gamma_t > 0``) symmetric in ``\mathbf{k}_{\mathrm{in}} \cdot \mathbf{k}_{\mathrm{out}}``,
   preserving exact cross-face BTDF reciprocity
   ``f_{t,\mathrm{ad}\to\mathrm{ab}}(\Omega_1 \to \Omega_2) = f_{t,\mathrm{ab}\to\mathrm{ad}}(-\Omega_2 \to -\Omega_1)``
   and exact hemispherical integration ([`leaf_bsdf_spectrum`](@ref)).

## Canopy scattering integration

Canopy scattering kernels consume these directional and dorsiventral leaf models
at two levels of fidelity:

- [`SpecularCanopyScattering`](@ref) and [`compute_reflection_mueller`](@ref)
  evaluate leaf-angle-averaged microfacet reflection in the right-handed canopy
  meridian-plane Stokes basis ``(\hat{\mathbf{e}}_s, \hat{\mathbf{e}}_p, \mathbf{k})``
  (with ``+Q`` along ``\hat{\mathbf{e}}_s`` perpendicular to the meridian plane)
  using singularity-free Cartesian specular bisectors (`getSpecularΩ`), preserving
  Hovenier 4×4 reciprocity and ``m \ge 1`` ``Q \leftrightarrow U`` Fourier couplings.
- [`DorsiventralBiLambertianCanopyScattering`](@ref) and
  [`compute_Z_matrices_directional`](@ref) provide closed-form Fourier moments
  for dorsiventral bi-Lambertian leaves (``R_{\mathrm{ad}} \neq R_{\mathrm{ab}}``),
  returning all four directional blocks ``(Z^{++}, Z^{-+}, Z^{+-}, Z^{--})`` and
  separate downward/upward stream single-scattering albedos ``(\omega_{\mathrm{down}}, \omega_{\mathrm{up}})``.
- [`compute_Z_matrices_coupled_bsdf`](@ref) projects a prepared
  [`PreparedLeafBSDFSpectrum`](@ref) over leaf inclination and azimuth onto the
  four canopy Fourier phase-matrix blocks while preserving exact passive stream
  albedo bounds ``0 \le \omega_{\mathrm{down}}, \omega_{\mathrm{up}} \le 1`` and
  4-block canopy Helmholtz reciprocity.
