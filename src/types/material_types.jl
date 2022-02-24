"Abstract Material Type (can be all)"
abstract type AbstractMaterial end

"Abstract water type"
abstract type AbstractWater <: AbstractMaterial end
"Abstract soil type"
abstract type AbstractSoil  <: AbstractMaterial end

"Pure liquid water"
struct LiquidPureWater <: AbstractWater end

"Salty liquid water"
Base.@kwdef struct LiquidSaltWater{FT} <: AbstractWater
    "Salinity in `[PSU]`"
    S::FT = FT(10)
end

"Pure Ice"
struct PureIce <: AbstractWater end

"Soil MW properties"
Base.@kwdef struct SoilMW{FT} <: AbstractSoil 
    "Sand Fraction ∈ [0,1]"
    sand_frac::FT  = FT(0.2)
    "Clay Fraction ∈ [0,1]"
    clay_frac::FT  = FT(0.1)
    "Volumetric water content ∈ [0,1]"
    mᵥ::FT = FT(0.35)
    "Bulk density, `g/cm³` (typical value is 1.7 gcm³)"
    ρ::FT          = FT(1.7)
end



