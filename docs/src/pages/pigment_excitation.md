# Pigment absorption and antenna excitation

The optical solver supplies total absorbed radiation. Pigment allocation then
answers which material absorbed those photons. Antenna selection answers which
absorbed photons transfer excitation into the photosynthetic antenna. Photosystem
allocation is a third, separate choice. Keeping these steps separate allows a
future leaf optical model or physiological model to replace one component without
owning the others.

The core workflow is scalar and applies separately to each wavelength bin:

```julia
partition = partition_absorption(HomogeneousAbsorption(), absorbed_photons, coefficients)
antenna = antenna_excitation(selection, partition)
q2 = photosystem_excitation(antenna, psii_fraction)
```

`coefficients` is a NamedTuple of concentration × specific-absorption products in
a common optical-depth convention. Its values and the photon input use the same
floating-point type. The caller supplies the optical model, spectrum, composition,
transfer probabilities and photosystem fractions. No species profile, light field
or empirical transfer fraction is selected by the core.

## Homogeneous material assumption

If every optical path encounters the same constituent proportions, then at one
wavelength the fraction of material absorption due to constituent ``i`` is

```math
Q_i = Q_{\rm abs}\frac{k_i}{\sum_j k_j}.
```

This follows by integrating the species-specific absorption rate
``k_i\exp[-(\sum_j k_j)s]`` along any path, then summing paths. It remains valid
with multiple reflections when composition is unchanged along those paths. It is
not a model of spatially segregated pigments in heterogeneous mesophyll. A different
partition rule can implement `AbstractAbsorptionPartition` when that physics is
needed.

`HomogeneousAbsorption` first scales coefficients by their maximum to avoid
summing overflowing products. Normalizing these constituent probabilities is the
physical partition definition; it does not alter or normalize total optical
absorption. The optical solver must independently close its radiation budget.
Do not attribute broadband absorption using a broadband-averaged coefficient
ratio: partition the spectral bins before integration.

For PROSPECT-PRO, the host helper `prospect_constituent_coefficients(leaf, optics)`
returns the exact products also used by `prospect_budget`. It names chlorophyll,
carotenoids, anthocyanins, brown material, water, bulk dry matter, protein and carbon
constituents. Their common division by plate count cancels in the fractions.
Set bulk dry matter `Cm=0` when protein/carbon constituents already account for the
same mass. For a different optical model, supply a different named tuple directly;
no PROSPECT fields or database are required by the scalar partition.

Use `prospect_budget(...).absorptance` for weak material absorption. The rounded
residual `1−R−T` cannot retain absorption smaller than the spacing of the separately
reported reflection and transmission. See [the numerical optical contract](prospect_budget.md).

## Antenna and photosystem choices

`ChlorophyllOnly()` selects chlorophyll absorption. `CarotenoidAntenna(eta)` selects

```math
Q_{\rm ant}=Q_{\rm chl}+\eta_{\rm car}Q_{\rm car},\qquad 0\le\eta_{\rm car}\le1.
```

The carotenoid transfer fraction is a **photon-transfer probability**. It is not
an energy efficiency or the physiological NPQ rate ``K_N``. Untransferred carotenoid
absorption stays in the material absorption/energy ledger; it is not automatically
assigned to a particular heat or quenching pathway. These conventions are consistent
with the constituent allocation described by [Platespect (2025)](https://doi.org/10.1016/j.rse.2025.114990).

The partition retains unitless constituent fractions as well as photon diagnostics.
Built-in selection accumulates the smaller selected or unselected probability and
multiplies by total photons, using the complementary probability when appropriate.
This equivalent probability identity avoids a one-ulp excess over total photons near complete
selection. It neither clips photons nor changes constituent diagnostics. This
preserves pigment identity even when constituent photon amounts underflow. At
subnormal input magnitudes, individual constituent values may still round to zero;
extreme coefficient ratios can also underflow the constituent fractions themselves.
No relative-accuracy claim below the selected FT's probability/flux resolution is made.

The separate photosystem allocation computes ``Q_2(\lambda)=\alpha_2(\lambda)
Q_{\rm ant}(\lambda)`` with an explicit finite fraction in ``[0,1]``. The caller
then integrates Q2 over a declared spectral domain. A leaf model receiving this
already allocated Q2 must not apply another PSII fraction. The same excitation
should drive electron transport and the physiological fluorescence source.
This optics package introduces no photosynthesis or fluorescence dependency.

Total material photon absorption, antenna excitation, PSII excitation and absorbed
energy remain different quantities. Energy requires wavelength-specific photon
energy before summation. PAR energy is not the full shortwave energy budget.
Use consistent spectral extents, units and reduction conventions for total PAR
and Q2. The all-selected ``\eta_{\rm car}=\alpha_2=1`` endpoint is tested on both
CPU and CUDA, including exact equality of the resulting integrated photon arrays.

## Status and extension contracts

| Status | Meaning |
| --- | --- |
| `AbsorptionPartitionValid` | Valid material attribution, including zero illumination of absorbing material |
| `AbsorptionNoMaterial` | All coefficient products are zero (or the tuple is empty), and supplied absorption is zero |
| `AbsorptionInvalidInput` | A photon input or coefficient product is negative or nonfinite |
| `AbsorptionIncompatibleInput` | Nonzero absorption was supplied for material with no absorption coefficients |

Invalid results contain NaN numerical payloads. For absent absorption coefficients,
constituent fractions are undefined (NaN) although absorbed photons are zero.
`absorption_valid(result)` checks
status; the antenna wrapper checks it before calling a selection. Built-in policies
return typed NaN for a missing required pigment name, even for an empty partition.
`ChlorophyllOnly` does not require a carotenoid key. `photosystem_excitation` returns
NaN for invalid photon inputs or fractions. Mixed-precision coefficient tuples
are not silently converted.

Add a new absorber by adding a name/value to the coefficient NamedTuple. A custom
`AbstractAntennaSelection` implements `select_absorbed_photons(selection, partition)`
and can use that name without changing a central list. It returns a same-FT photon
amount, or NaN when required inputs are unavailable. Normally call the public
`antenna_excitation` wrapper, which checks result status and the returned photon
bound. Custom selectors remain responsible for the physical meaning and validation
of their own parameters.

The implemented kernels are tested in Float32/Float64 on CPU and NVIDIA CUDA,
including status propagation, empty/reordered batches and a custom selector.
Optical tables and PROSPECT composition are prepared on the host. This does not
claim GPU PROSPECT special functions, Metal execution or automatic differentiation
of these new partition/selection methods. No canopy grid or persistent individual
leaf identity is introduced.
