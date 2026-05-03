
"""
$(FUNCTIONNAME)(mod::LiquidSaltWater, T::FT,f::FT)

Computes the complex dieletric of liquid salty walter (Ulaby & Long book)
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct, provide Salinity in PSU 
- `T`  Temperature in `[⁰K]`
- `f`  Frequency in `[GHz]`
- `S`  Salinity in `[PSU]` comes from the `mod` [`LiquidSaltWater`](@ref) struct

# Examples
```julia-repl
julia> w = CanopyOptics.LiquidSaltWater(S=10.0)     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,283.0,10.0)
51.79254073931596 + 38.32304044382495im
```
"""
function dielectric(mod::LiquidSaltWater, T::FT,f::FT) where FT<:Real
    (;S)  = mod
    0 ≤ S ≤ 45     || throw(ArgumentError("Salinity must be ∈ [0, 45] PSU (got $S)"))
    265 ≤ T ≤ 310  || throw(ArgumentError("Temperature must be ∈ [265, 310] K (got $T)"))
    # Equations were tuned for ⁰C but want K input!
    T = T - FT(273)        
    # Compute σ (Eqs 4.21 - 4.21 in Ulaby and Long)
    σ₃₅ = TempPoly(T)
    P   = S * SalPoly1(S) / SalPoly2(S)
    Q   = FT(1) + (α₀PolNom(S)/α₀PolDenom(S) * (T - FT(15)))/(T + α₁(T));
        
    σ = σ₃₅ * P * Q
 
    ϵS   = FT(87.85306) * exp(FT(-0.00456992)*T - aU[1]*S - aU[2]*S^2 - aU[3]*S*T);
    ϵOne =  aU[4] * exp(-aU[5]*T - aU[6]*S - aU[7]*S*T);
    # Relaxation time constants
    τ1   = (aU[8]  + aU[9] *S) * exp(aU[10] / (T+aU[11]));
    τ2   = (aU[12] + aU[13]*S) * exp(aU[14] / (T+aU[15]));
    ϵInf =  aU[16] + aU[17]*T  + aU[18]*S;

    # Complex Permittivity Calculation
    ϵ = ((ϵS-ϵOne) / (FT(1) - 1im * 2π * f * τ1)) + ((ϵOne-ϵInf) / (1-1im * 2π * f * τ2)) + (ϵInf) + 1im*((FT(17.9751)*σ) /f);
end

"""
$(FUNCTIONNAME)(mod::LiquidPureWater, T::FT,f::FT)

Computes the complex dieletric of liquid pure walter (Ulaby & Long book)
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct
- `T`  Temperature in `[⁰K]`
- `f`  Frequency in `[GHz]`

# Examples
```julia-repl
julia> w = CanopyOptics.LiquidPureWater()     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,283.0,10.0)
53.45306674215052 + 38.068642430090044im
```
"""
function dielectric(mod::LiquidPureWater, T::FT,f::FT) where FT<:Real
    265 ≤ T ≤ 310 || throw(ArgumentError("Temperature must be ∈ [265, 310] K (got $T)"))
    # Equations were tuned for ⁰C but want K input!
    T = T - FT(273)

    ϵS   =  FT(87.85306) * exp(FT(-0.00456992)*T);
    ϵOne =  aU[4]  * exp(-aU[5]*T);
    τ1   =  aU[8]  * exp(aU[10] / (T+aU[11]));
    τ2   =  aU[12] * exp(aU[14] / (T+aU[15]));
    ϵInf =  aU[16] + aU[17]*T;
    # Complex Permittivity Calculation
    ϵ = ((ϵS-ϵOne) / (FT(1) - 1im * 2π * f * τ1)) + ((ϵOne-ϵInf)./(1-1im * 2π * f * τ2)) + (ϵInf);
end

"""
$(FUNCTIONNAME)(mod::PureIce, T::FT,f::FT)

Computes the complex dielectric constant of pure ice (Ulaby & Long book, Eqs. 4.22–4.25).
Valid temperature range: 233–273.15 K.
# Arguments
- `mod` a [`PureIce`](@ref) type struct
- `T`  Temperature in `[K]`, must be ∈ [233, 273.15]
- `f`  Frequency in `[GHz]`

# Examples
```julia-repl
julia> w = CanopyOptics.PureIce()
julia> CanopyOptics.dielectric(w, 253.0, 10.0)   # T = -20°C, f = 10 GHz
# Real part ≈ 3.17 (ice has low imaginary dielectric compared to liquid water)
```
"""
function dielectric(mod::PureIce, T::FT, f::FT) where FT
    233 ≤ T ≤ 273.15 || throw(ArgumentError("Temperature must be ∈ [233, 273.15] K (got $T)"))
    # Some constants first
    θ  = (FT(300)/T) - FT(1);
    B1 = FT(0.0207);
    B2 = FT(1.16e-11);
    b  = FT(335);
    # Eqs 4.22-4.25
    α = (FT(0.00504) + FT(0.0062) * θ) * exp(-FT(22.1) * θ);
    βₘ  = (B1/T) * exp(b/T) / (exp(b/T)-FT(1))^2  + B2 * f^2;
    δβ  = exp(-FT(9.963) + FT(0.0372) * (T-FT(273.16)));
    
    β = βₘ + δβ;
    ϵ = FT(3.1884) + FT(9.1e-4) *(T-FT(273)) + (α/f + β*f)im;
end

"""
$(FUNCTIONNAME)(mod::SoilMW, T::FT,f::FT)

Computes the complex dielectric constant of moist soil (Ulaby & Long book, Eqs. 4.66–4.70).
# Arguments
- `mod` a `SoilMW` type struct (fields: `sand_frac`, `clay_frac`, `mᵥ`, `ρ`)
- `T`  Temperature in `[K]`
- `f`  Frequency in `[GHz]`

# Examples
```julia-repl
julia> w = SoilMW(sand_frac=0.2, clay_frac=0.2, mᵥ=0.3, ρ=1.7)
julia> dielectric(w, 283.0, 10.0)
```
"""
function dielectric(mod::SoilMW{FT}, T::FT,f::FT) where FT
    (;sand_frac, clay_frac, mᵥ, ρ) = mod 
    # T = Tk - FT(273)
    f_hz = f * FT(1e9);    #  transform from GHz to Hz

    α  = FT(0.65);   # eq: 4.68a
    β₁ = FT(1.27) - FT(0.519) * sand_frac - FT(0.152) * clay_frac; # eq: 4.68b
    β₂ = FT(2.06) - FT(0.928) * sand_frac - FT(0.255) * clay_frac; # eq: 4.68c 
    
    ϵ₀ = FT(8.854e-12); 

    if f > 1.3
        σ = FT(-1.645) + FT(1.939) * ρ - FT(2.256) * sand_frac + FT(1.594) * clay_frac;  # eq: 4.68d
    elseif 0.3 ≤ f ≤ 1.3
        σ = FT(0.0467) + FT(0.22)  * ρ - FT(0.411) * sand_frac + FT(0.661) * clay_frac;  # eq: 4.70
    else
        σ = FT(0)
    end

    # Dielectric Constant of Pure Water
    ϵw =  dielectric(LiquidPureWater(),T,f)
    # Add conductivity term to ϵ'' (eq. 4.67b)
    ϵw += im*((FT(2.65)-ρ) / FT(2.65 * mᵥ) * σ / (2π * ϵ₀ * f_hz))

    # calculating dielectric constant of soil using eq. 4.66a and 4.66b
    epsr = (FT(1) + FT(0.66) * ρ + mᵥ^β₁ * real(ϵw).^α - mᵥ)^(1/α);
    epsi = mᵥ^β₂ .* imag(ϵw);
    epsr + epsi*im
end

# === Vegetation dielectric models ============================================

"""
$(FUNCTIONNAME)(mod::LeafUlabyElRayes1987, T, f)

Complex relative permittivity ε(f) of a fresh leaf with gravimetric
moisture `mod.M_g`, at frequency `f` [GHz]. Implements the dual-Debye
mixing model of Ulaby & El-Rayes (1987):

```math
\\varepsilon = \\varepsilon_r + v_{fw}\\,\\varepsilon_{fw}
             + v_b\\,\\varepsilon_b
```

where the residual, free-water and bound-water volume fractions are
fitted from `M_g`, and ``\\varepsilon_{fw}``, ``\\varepsilon_b`` are
the free- and bound-water dispersions. The temperature `T` [K] is
accepted for signature consistency with the other [`dielectric`](@ref)
methods but is not used by this model.

Sign convention: `+i ε''` (loss positive in imaginary part), matching
[`LiquidPureWater`](@ref) and [`SoilMW`](@ref) in this package.

# Examples
```julia-repl
julia> leaf = LeafUlabyElRayes1987(M_g = 0.5);

julia> dielectric(leaf, 295.0, 5.0)         # C-band, mid-moisture
14.400774… + 4.689454…im
```

# Reference
Ulaby & El-Rayes 1987, IEEE TGRS 25(5), 550–557 (Eqs. 11–14);
Ulaby & Long 2014, §11-9.
"""
function dielectric(mod::LeafUlabyElRayes1987, T::Real, f::Real)
    (; M_g, σ) = mod
    0 ≤ M_g ≤ 0.7   || throw(ArgumentError("Gravimetric moisture must be ∈ [0, 0.7] (got $M_g)"))
    0.2 ≤ f ≤ 20    || throw(ArgumentError("Frequency must be ∈ [0.2, 20] GHz for Ulaby & El-Rayes 1987 (got $f)"))

    FT = promote_type(typeof(M_g), typeof(σ), typeof(f))

    # Volume fractions of free and bound water (Eqs. 12, 13)
    v_fw = M_g * (FT(0.55) * M_g - FT(0.076))
    v_b  = FT(4.64) * M_g^2 / (FT(1) + FT(7.36) * M_g^2)

    # Non-dispersive residual (Eq. 11)
    ε_r  = FT(1.7) - FT(0.74) * M_g + FT(6.16) * M_g^2

    # Free-water Debye + ionic conductivity loss (Eq. 14a, physics convention)
    # Conductivity coefficient 17.9751 = 1/(2π ε₀) with f in GHz, σ in S/m
    # — same constant used by LiquidSaltWater above.
    ε_fw = FT(4.9) + FT(75.0) / (FT(1) - 1im * f / FT(18.0)) +
           1im * FT(17.9751) * σ / f

    # Bound-water dispersion (Eq. 14b)
    ε_b  = FT(2.9) + FT(55.0) / (FT(1) + sqrt(-1im * f / FT(0.18)))

    return ε_r + v_fw * ε_fw + v_b * ε_b
end
