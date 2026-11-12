"""
    dielectric(mod::LiquidSaltWater, T, f)

Compute the complex relative permittivity `ε(T, f)` of saline liquid water
following Ulaby & Long (2014, Chapter 4, Eqs. 4.19–4.21), using the `+i ε''`
loss sign convention.

# Arguments
- `mod::LiquidSaltWater`: saline water model supplying salinity `S` in `[0, 45]` `PSU`.
- `T::Real`: temperature in `K` (`265 <= T <= 310`).
- `f::Real`: frequency in `GHz`.

# Examples
```jldoctest
julia> w = LiquidSaltWater(S = 10.0);

julia> round(dielectric(w, 283.0, 10.0); digits = 4)
51.7925 + 38.323im
```
"""
function dielectric(mod::LiquidSaltWater{FTS}, T_in::FTT, f_in::FTF) where {FTS<:Real,FTT<:Real,FTF<:Real}
    FT = float(promote_type(FTS, FTT, FTF))
    BFT = _leaf_base_float(FT)
    S = FT(mod.S)
    T_k = FT(T_in)
    f = FT(f_in)
    S_p = prospect_primal(S)
    T_p = prospect_primal(T_k)
    zero(BFT) <= S_p <= BFT(45) ||
        throw(ArgumentError("Salinity must be ∈ [0, 45] PSU (got $S_p)"))
    BFT(265) <= T_p <= BFT(310) ||
        throw(ArgumentError("Temperature must be ∈ [265, 310] K (got $T_p)"))

    T = T_k - BFT(273)
    σ35 = evalpoly(T, BFT.(_TEMP_POLY_COEFFS))
    P   = S * evalpoly(S, BFT.(_SAL_POLY1_COEFFS)) / evalpoly(S, BFT.(_SAL_POLY2_COEFFS))
    α0  = evalpoly(S, BFT.(_ALPHA0_NOM_COEFFS)) / evalpoly(S, BFT.(_ALPHA0_DEN_COEFFS))
    Q   = one(FT) + (α0 * (T - BFT(15))) / (T + evalpoly(T, BFT.(_ALPHA1_COEFFS)))
    σ   = σ35 * P * Q

    au = BFT.(_ULABY_WATER_AU)
    ϵS   = BFT(87.85306) * exp(BFT(-0.00456992) * T - au[1] * S - au[2] * S^2 - au[3] * S * T)
    ϵOne = au[4] * exp(-au[5] * T - au[6] * S - au[7] * S * T)
    τ1   = (au[8]  + au[9]  * S) * exp(au[10] / (T + au[11]))
    τ2   = (au[12] + au[13] * S) * exp(au[14] / (T + au[15]))
    ϵInf = au[16] + au[17] * T + au[18] * S

    two_pi_f = (BFT(2) * FT(π)) * f
    return (ϵS - ϵOne) / Complex(one(FT), -two_pi_f * τ1) +
           (ϵOne - ϵInf) / Complex(one(FT), -two_pi_f * τ2) +
           ϵInf + Complex(zero(FT), (BFT(17.9751) * σ) / f)
end

"""
    dielectric(mod::LiquidPureWater, T, f)

Compute the complex relative permittivity `ε(T, f)` of pure liquid water
following Ulaby & Long (2014, Chapter 4), using the `+i ε''` loss sign
convention.

# Arguments
- `mod::LiquidPureWater`: pure liquid water tag.
- `T::Real`: temperature in `K` (`265 <= T <= 310`).
- `f::Real`: frequency in `GHz`.

# Examples
```jldoctest
julia> round(dielectric(LiquidPureWater(), 283.0, 10.0); digits = 4)
53.4531 + 38.0686im
```
"""
function dielectric(::LiquidPureWater, T_in::FTT, f_in::FTF) where {FTT<:Real,FTF<:Real}
    FT = float(promote_type(FTT, FTF))
    BFT = _leaf_base_float(FT)
    T_k = FT(T_in)
    f = FT(f_in)
    T_p = prospect_primal(T_k)
    BFT(265) <= T_p <= BFT(310) ||
        throw(ArgumentError("Temperature must be ∈ [265, 310] K (got $T_p)"))

    T = T_k - BFT(273)
    au = BFT.(_ULABY_WATER_AU)
    ϵS   = BFT(87.85306) * exp(BFT(-0.00456992) * T)
    ϵOne = au[4]  * exp(-au[5] * T)
    τ1   = au[8]  * exp(au[10] / (T + au[11]))
    τ2   = au[12] * exp(au[14] / (T + au[15]))
    ϵInf = au[16] + au[17] * T

    two_pi_f = (BFT(2) * FT(π)) * f
    return (ϵS - ϵOne) / Complex(one(FT), -two_pi_f * τ1) +
           (ϵOne - ϵInf) / Complex(one(FT), -two_pi_f * τ2) +
           ϵInf
end

"""
    dielectric(mod::PureIce, T, f)

Compute the complex relative permittivity `ε(T, f)` of pure water ice following
Ulaby & Long (2014, Chapter 4, Eqs. 4.22–4.25), valid for `233 <= T <= 273.15` `K`.

# Arguments
- `mod::PureIce`: pure water ice tag.
- `T::Real`: temperature in `K` (`233 <= T <= 273.15`).
- `f::Real`: frequency in `GHz`.

# Examples
```jldoctest
julia> round(real(dielectric(PureIce(), 253.0, 10.0)); digits = 4)
3.1702
```
"""
function dielectric(::PureIce, T_in::FTT, f_in::FTF) where {FTT<:Real,FTF<:Real}
    FT = float(promote_type(FTT, FTF))
    BFT = _leaf_base_float(FT)
    T = FT(T_in)
    f = FT(f_in)
    T_p = prospect_primal(T)
    BFT(233) <= T_p <= BFT(273.15) ||
        throw(ArgumentError("Temperature must be ∈ [233, 273.15] K (got $T_p)"))

    θ  = (BFT(300) / T) - one(FT)
    B1 = BFT(0.0207)
    B2 = BFT(1.16e-11)
    b  = BFT(335)

    α  = (BFT(0.00504) + BFT(0.0062) * θ) * exp(-BFT(22.1) * θ)
    eb = exp(b / T)
    βₘ = (B1 / T) * eb / (eb - one(FT))^2 + B2 * f^2
    δβ = exp(-BFT(9.963) + BFT(0.0372) * (T - BFT(273.16)))
    β  = βₘ + δβ

    return Complex(BFT(3.1884) + BFT(9.1e-4) * (T - BFT(273)), α / f + β * f)
end

"""
    dielectric(mod::SoilMW, T, f)

Compute the complex relative permittivity `ε(T, f)` of moist soil following the
Dobson semi-empirical mixing model in Ulaby & Long (2014, Chapter 4, Eqs. 4.66–4.70).

# Arguments
- `mod::SoilMW`: soil texture, moisture, and bulk-density parameters.
- `T::Real`: temperature in `K` (`265 <= T <= 310`).
- `f::Real`: frequency in `GHz`.

# Examples
```jldoctest
julia> s = SoilMW(sand_frac = 0.2, clay_frac = 0.2, mᵥ = 0.3, ρ = 1.7);

julia> real(dielectric(s, 283.0, 10.0)) > 10
true
```
"""
function dielectric(mod::SoilMW{FTS}, T_in::FTT, f_in::FTF) where {FTS<:Real,FTT<:Real,FTF<:Real}
    FT = float(promote_type(FTS, FTT, FTF))
    BFT = _leaf_base_float(FT)
    sand_frac = FT(mod.sand_frac)
    clay_frac = FT(mod.clay_frac)
    mᵥ = FT(mod.mᵥ)
    ρ = FT(mod.ρ)
    T = FT(T_in)
    f = FT(f_in)
    f_p = prospect_primal(f)

    f_hz = f * BFT(1e9)
    α  = BFT(0.65)
    β₁ = BFT(1.27) - BFT(0.519) * sand_frac - BFT(0.152) * clay_frac
    β₂ = BFT(2.06) - BFT(0.928) * sand_frac - BFT(0.255) * clay_frac
    ϵ₀ = BFT(8.854e-12)

    σ = if f_p > BFT(1.3)
        BFT(-1.645) + BFT(1.939) * ρ - BFT(2.256) * sand_frac + BFT(1.594) * clay_frac
    elseif BFT(0.3) <= f_p <= BFT(1.3)
        BFT(0.0467) + BFT(0.22) * ρ - BFT(0.411) * sand_frac + BFT(0.661) * clay_frac
    else
        zero(FT)
    end

    if iszero(prospect_primal(mᵥ)) && _same_val_and_partials(mᵥ, zero(FT))
        return Complex((one(FT) + BFT(0.66) * ρ)^inv(α), zero(FT))
    end

    ϵw = dielectric(LiquidPureWater(), T, f)
    cond_term = ((BFT(2.65) - ρ) / BFT(2.65)) * σ / ((BFT(2) * FT(π)) * ϵ₀ * f_hz)
    epsr = (one(FT) + BFT(0.66) * ρ + mᵥ^β₁ * real(ϵw)^α - mᵥ)^inv(α)
    epsi = mᵥ^β₂ * imag(ϵw) + mᵥ^(β₂ - one(FT)) * cond_term
    return Complex(epsr, epsi)
end

"""
    dielectric(mod::LeafUlabyElRayes1987, T, f)

Compute the complex relative permittivity `ε(f)` of a fresh leaf with
gravimetric moisture `mod.M_g` at frequency `f` in `GHz` (`0.2 <= f <= 20`)
following Ulaby & El-Rayes (1987):

```math
\\varepsilon = \\varepsilon_r + v_{fw}\\,\\varepsilon_{fw}
             + v_b\\,\\varepsilon_b .
```

Temperature `T` in `K` is accepted for signature consistency with other
[`dielectric`](@ref) methods. Uses the `+i ε''` loss sign convention.

# Examples
```jldoctest
julia> leaf = LeafUlabyElRayes1987(M_g = 0.5);

julia> round(dielectric(leaf, 295.0, 5.0); digits = 4)
14.4008 + 4.6895im
```

# Reference
Ulaby & El-Rayes (1987), IEEE TGRS 25(5), 550–557 (Eqs. 11–14);
Ulaby & Long (2014), Chapter 11, Section 11-9.
"""
function dielectric(mod::LeafUlabyElRayes1987{FTM}, ::Real, f_in::FTF) where {FTM<:Real,FTF<:Real}
    FT = float(promote_type(FTM, FTF))
    BFT = _leaf_base_float(FT)
    M_g = FT(mod.M_g)
    σ = FT(mod.σ)
    f = FT(f_in)
    M_p = prospect_primal(M_g)
    f_p = prospect_primal(f)
    zero(BFT) <= M_p <= BFT(0.7) ||
        throw(ArgumentError("Gravimetric moisture must be ∈ [0, 0.7] (got $M_p)"))
    BFT(0.2) <= f_p <= BFT(20) ||
        throw(ArgumentError("Frequency must be ∈ [0.2, 20] GHz for Ulaby & El-Rayes 1987 (got $f_p)"))

    v_fw = M_g * (BFT(0.55) * M_g - BFT(0.076))
    v_b  = BFT(4.64) * M_g^2 / (one(FT) + BFT(7.36) * M_g^2)
    ε_r  = BFT(1.7) - BFT(0.74) * M_g + BFT(6.16) * M_g^2

    ε_fw = BFT(4.9) + BFT(75.0) / Complex(one(FT), -f / BFT(18.0)) +
           Complex(zero(FT), BFT(17.9751) * σ / f)

    # sqrt(-i * x) = sqrt(x / 2) * (1 - i) for real x > 0 (AD-safe on Complex{Dual}!)
    s_half = sqrt(f / BFT(0.36))
    ε_b  = BFT(2.9) + BFT(55.0) / Complex(one(FT) + s_half, -s_half)

    return ε_r + v_fw * ε_fw + v_b * ε_b
end
