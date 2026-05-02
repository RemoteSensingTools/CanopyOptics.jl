# CanopyOptics.jl
*A package to compute canopy scattering properties*
## Package Features
- Use leaf angle distributions to compute bi-Lambertian canopy scattering matrices
- Compute closed-form cosine Fourier moments for the bi-Lambertian canopy kernel
- Compose diffuse and specular leaf-surface scattering models additively
- Differentiate canopy Z matrices with ForwardDiff leaf optical parameters
- Compute experimental specular leaf-surface reflection terms
- Compute leaf reflectance and transmittance based on Prospect-PRO

## Z-matrix convention

Canopy scattering matrices use `Z[i_out, j_in]`: rows are outgoing streams and
columns are incoming streams. `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is
sign-change reflection. The leaf single-scattering albedo `ϖ = R + T` is not
folded into `Z`; layer solvers multiply by `ϖ` separately. For conservative
bi-Lambertian leaves, the `m = 0` column integral of `Z⁺⁺ + Z⁻⁺` is therefore
approximately `2` under the quadrature weights.

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
