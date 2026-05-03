"Abstract Material Type (can be all)"
abstract type AbstractMaterial end

"Abstract water type"
abstract type AbstractWater      <: AbstractMaterial end
"Abstract soil type"
abstract type AbstractSoil       <: AbstractMaterial end
"Abstract vegetation type (leaves, needles, branches, trunks)"
abstract type AbstractVegetation <: AbstractMaterial end

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

"""
    LeafUlabyElRayes1987(; M_g=0.5, σ=1.27)

Microwave dielectric model for leaves following Ulaby & El-Rayes (1987).
A dual-Debye mixing model combining a non-dispersive residual term, a
free-water Debye component, and a bound-water Cole–Cole-like component.
Valid for fresh leaves at 0.2–20 GHz with gravimetric moisture
`M_g ∈ [0, 0.7]`.

Pass to [`dielectric`](@ref) together with temperature `T` [K] (unused
by this model — kept for signature consistency) and frequency `f` [GHz].

# Fields
- `M_g::FT` — gravimetric moisture content ∈ [0, 0.7]
- `σ::FT`  — effective ionic conductivity of free water `[S/m]` (typ. 1.27)

# Reference
Ulaby, F. T. & El-Rayes, M. A. (1987). *Microwave Dielectric Spectrum
of Vegetation – Part II: Dual-Dispersion Model.* IEEE TGRS 25(5),
550–557. Reproduced in Ulaby & Long (2014), §11-9.
"""
struct LeafUlabyElRayes1987{FT<:Real} <: AbstractVegetation
    M_g::FT
    σ::FT
end

LeafUlabyElRayes1987(; M_g = 0.5, σ = 1.27) =
    LeafUlabyElRayes1987(promote(M_g, σ)...)
