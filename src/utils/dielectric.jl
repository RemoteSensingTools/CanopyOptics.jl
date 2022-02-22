
"""
$(FUNCTIONNAME)(mod::LiquidSaltWater, T::FT,f::FT,S::FT)

Computes the complex dieletric of liquid salty walter
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct, store polynomial terms
- `T`  Temperature in `[⁰C]`
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

Computes the complex dieletric of liquid pure walter
# Arguments
- `mod` an [`LiquidSaltWater`](@ref) type struct, store polynomial terms
- `T`  Temperature in `[⁰C]`
- `f`  Frequency in `[GHz]`
- `S`  Salinity in `[PSU]`

# Examples
```julia-repl
julia> w = CanopyOptics.LiquidPureWater()     # Create struct for salty seawater
julia> CanopyOptics.dielectric(w,10.0,10.0,10.0)
53.45306674215052 + 38.068642430090044im
```
"""
function dielectric(mod::LiquidPureWater, T::FT,f::FT,S::FT=FT(0)) where FT<:Real
    ϵS   = FT(87.85306) * exp(FT(-0.00456992)*T);
    ϵOne =  aU[4]  * exp(-aU[5]*T);
    τ1   =  aU[8]  * exp(aU[10] / (T+aU[11]));
    τ2   =  aU[12] * exp(aU[14] / (T+aU[15]));
    ϵInf =  aU[16] + aU[17]*T;
    # Complex Permittivity Calculation
    ϵ = ((ϵS-ϵOne) / (FT(1) - 1im * 2π * f * τ1)) + ((ϵOne-ϵInf)./(1-1im * 2π * f * τ2)) + (ϵInf);
end

