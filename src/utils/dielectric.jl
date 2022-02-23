
"""
$(FUNCTIONNAME)(mod::LiquidSaltWater, T::FT,f::FT,S::FT)

Computes the complex dieletric of liquid salty walter (Ulaby & Long book)
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct, store polynomial terms
- `T`  Temperature in `[⁰K]`
- `f`  Frequency in `[GHz]`
- `S`  Salinity in `[PSU]`

# Examples
```julia-repl
julia> w = CanopyOptics.LiquidSaltWater()     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,10.0,10.0,10.0)
51.79254073931596 + 38.32304044382495im
```
"""
function dielectric(mod::LiquidSaltWater, T::FT,f::FT,S::FT) where FT<:Real
    #@unpack TempPoly, SalPoly1, SalPoly2, α₀PolNom, α₀PolDenom, α₁, aU  = mod
    @assert 0 ≤ S ≤ 45 "Salinity should be ∈ [0,45] PSU"
    @assert 265 ≤ T ≤ 310 "Temperature should be ∈ [265,310] K"
    # Equations were tuned for ⁰C but want K input!
    T = T - FT(273.15)        
    # Compute σ (Eqs 4.21 - 4.21 in Ulaby and Long)
    σ₃₅ = TempPoly(T)
    P   = S * SalPoly1(S) / SalPoly2(S)
    Q   = FT(1) + (α₀PolNom(S)/α₀PolDenom(S) * (T - FT(15)))/(T + α₁(T));
        
    σ = σ₃₅ * P * Q
 
    ϵS   = FT(87.85306) * exp(FT(-0.00456992)*T - aU[1]*S - aU[2]*S^2 - aU[3]*S*T);
    ϵOne =  aU[4] * exp(-aU[5]*T - aU[6]*S - aU[7]*S*T);
    τ1   = (aU[8]  + aU[9] *S) * exp(aU[10] / (T+aU[11]));
    τ2   = (aU[12] + aU[13]*S) * exp(aU[14] / (T+aU[15]));
    ϵInf =  aU[16] + aU[17]*T  + aU[18]*S;

    # Complex Permittivity Calculation
    ϵ = ((ϵS-ϵOne) / (FT(1) - 1im * 2π * f * τ1)) + ((ϵOne-ϵInf)./(1-1im * 2π * f * τ2)) + (ϵInf) + 1im*((FT(17.9751)*σ) /f);
end

"""
$(FUNCTIONNAME)(mod::LiquidPureWater, T::FT,f::FT,S::FT=0)

Computes the complex dieletric of liquid pure walter (Ulaby & Long book)
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct, store polynomial terms
- `T`  Temperature in `[⁰K]`
- `f`  Frequency in `[GHz]`

# Examples
```julia-repl
julia> w = CanopyOptics.LiquidPureWater()     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,10.0,10.0)
53.45306674215052 + 38.068642430090044im
```
"""
function dielectric(mod::LiquidPureWater, T::FT,f::FT,S::FT=FT(0)) where FT<:Real
    @assert 0 ≤ S ≤ 45 "Salinity should be ∈ [0,45] PSU"
    @assert 265 ≤ T ≤ 310 "Temperature should be ∈ [265,310] K"
    # Equations were tuned for ⁰C but want K input!
    T = T - FT(273.15)

    ϵS   = FT(87.85306) * exp(FT(-0.00456992)*T);
    ϵOne =  aU[4]  * exp(-aU[5]*T);
    τ1   =  aU[8]  * exp(aU[10] / (T+aU[11]));
    τ2   =  aU[12] * exp(aU[14] / (T+aU[15]));
    ϵInf =  aU[16] + aU[17]*T;
    # Complex Permittivity Calculation
    ϵ = ((ϵS-ϵOne) / (FT(1) - 1im * 2π * f * τ1)) + ((ϵOne-ϵInf)./(1-1im * 2π * f * τ2)) + (ϵInf);
end

"""
$(FUNCTIONNAME)(mod::PureIce, T::FT,f::FT)

Computes the complex dieletric of liquid pure walter (Ulaby & Long book)
# Arguments
- `mod` an [`PureIce`](@ref) type struct, store polynomial terms
- `T`  Temperature in `[⁰K]`
- `f`  Frequency in `[GHz]`

# Examples
```julia-repl
julia> w = CanopyOptics.PureIce()     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,10.0,10.0)
53.45306674215052 + 38.068642430090044im
```
"""
function dielectric(mod::PureIce, T::FT,f::FT,_=nothing) where FT
    @assert 233 ≤ T ≤ 273.15 "Temperature should be ∈ [233,273.15] K"
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

