# CanopyOptics.jl
*A package to compute canopy scattering properties*
## Package Features
- Use leaf angle distributions to compute bi-Lambertian canopy scattering matrices
- Compute closed-form cosine Fourier moments for the bi-Lambertian canopy kernel
- Compose diffuse and specular leaf-surface scattering models additively
- Request scalar or Stokes-expanded canopy Z matrices with `npol`
- Combine leaf and wood canopy components with separate angle distributions
- Differentiate canopy Z matrices with ForwardDiff leaf optical parameters
- Compute polarized specular leaf-surface reflection terms
- Compute leaf reflectance and transmittance based on Prospect-PRO
- Compute complex microwave dielectric ε(T, f) of water, ice, soil and leaves

## Z-matrix convention

Canopy scattering matrices use `Z[i_out, j_in]`: rows are outgoing streams and
columns are incoming streams. `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is
sign-change reflection. The leaf single-scattering albedo `ϖ = R + T` is not
folded into `Z`; layer solvers multiply by `ϖ` separately. For conservative
bi-Lambertian leaves, the `m = 0` column integral of `Z⁺⁺ + Z⁻⁺` is therefore
approximately `2` under the quadrature weights.

Specular components include their Fresnel/roughness strength in the returned Z
contribution. If a downstream RT solver multiplies canopy Z by a separate
single-scattering albedo, it must provide a matching effective albedo when
using specular or diffuse+specular composite leaves.

## Microwave dielectric models

A growing set of materials implements the [`dielectric`](@ref) function,
returning the complex relative permittivity `ε(T, f)` (loss as positive
imaginary part):

| Material | Type | Domain |
| --- | --- | --- |
| Pure liquid water | [`LiquidPureWater`](@ref) | 0.2–1000 GHz, 265–310 K |
| Salt water | [`LiquidSaltWater`](@ref) | + salinity 0–45 PSU |
| Pure ice | [`PureIce`](@ref) | 0.2–1000 GHz, 233–273 K |
| Moist soil | [`SoilMW`](@ref) | Dobson model, 0.3–18 GHz |
| Fresh leaf | [`LeafUlabyElRayes1987`](@ref) | 0.2–20 GHz, gravimetric moisture 0–0.7 |

All five share the contract `dielectric(model, T_kelvin, f_GHz)`. New
materials are added by extending the `dielectric` function on a fresh
`<: AbstractMaterial` subtype — see `src/utils/dielectric.jl` for the
existing methods. The vegetation hierarchy (`AbstractVegetation`) is
the entry point for a planned microwave-canopy expansion.

## Installation

The latest release of CanopyOptics can be installed from the Julia REPL prompt with

```julia
julia> ]add https://github.com/RemoteSensingTools/CanopyOptics.jl
```

## Code docs:

### Types 

```@autodocs
Modules = [CanopyOptics]
Private = false
Order = [:type]
```
### Functions 

```@autodocs
Modules = [CanopyOptics]
Private = false
Order = [:function]
```
