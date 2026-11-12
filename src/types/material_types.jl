"""
    AbstractMaterial

Supertype for materials supporting microwave complex permittivity evaluation via
[`dielectric`](@ref).
"""
abstract type AbstractMaterial end

"""
    AbstractWater <: AbstractMaterial

Supertype for liquid and frozen water permittivity models.
"""
abstract type AbstractWater <: AbstractMaterial end

"""
    AbstractSoil <: AbstractMaterial

Supertype for moist-soil microwave permittivity models.
"""
abstract type AbstractSoil <: AbstractMaterial end

"""
    AbstractVegetation <: AbstractMaterial

Supertype for leaf, needle, branch, and trunk microwave permittivity models.
"""
abstract type AbstractVegetation <: AbstractMaterial end

"""
    LiquidPureWater()

Singleton material tag for pure liquid water in [`dielectric`](@ref).
"""
struct LiquidPureWater <: AbstractWater end

"""
    LiquidSaltWater{FT<:Real}(S)
    LiquidSaltWater(; S = 10.0)

Saline liquid water material parameterized by salinity `S` in practical salinity
units (`PSU`, valid on `[0, 45]`).

# Fields
- `S::FT`: salinity in `PSU` (`0 <= S <= 45`).

# Examples
```jldoctest
julia> LiquidSaltWater(S = 35.0)
LiquidSaltWater{Float64}(35.0)
```
"""
struct LiquidSaltWater{FT<:Real} <: AbstractWater
    S::FT
end

LiquidSaltWater(; S::Real = 10.0) = LiquidSaltWater{typeof(float(S))}(float(S))
LiquidSaltWater{FT}(; S::Real = FT(10)) where {FT<:Real} = LiquidSaltWater{FT}(FT(S))

"""
    PureIce()

Singleton material tag for pure water ice in [`dielectric`](@ref).
"""
struct PureIce <: AbstractWater end

"""
    SoilMW{FT<:Real}(sand_frac, clay_frac, mᵥ, ρ)
    SoilMW(; sand_frac = 0.2, clay_frac = 0.1, mᵥ = 0.35, ρ = 1.7)

Dobson–Ulaby microwave soil material parameterized by texture fractions,
volumetric moisture, and bulk density.

# Fields
- `sand_frac::FT`: sand mass fraction in `[0, 1]`.
- `clay_frac::FT`: clay mass fraction in `[0, 1]`.
- `mᵥ::FT`: volumetric soil water content in `[0, 1]`.
- `ρ::FT`: dry soil bulk density in `g cm⁻³` (typical value `1.7`).

# Examples
```jldoctest
julia> s = SoilMW(sand_frac = 0.3, clay_frac = 0.2, mᵥ = 0.25, ρ = 1.6);

julia> s.mᵥ
0.25
```
"""
struct SoilMW{FT<:Real} <: AbstractSoil
    sand_frac::FT
    clay_frac::FT
    mᵥ::FT
    ρ::FT
end

function SoilMW(sand_frac::Real, clay_frac::Real, mᵥ::Real, ρ::Real)
    FT = float(promote_type(typeof(sand_frac), typeof(clay_frac), typeof(mᵥ), typeof(ρ)))
    return SoilMW{FT}(FT(sand_frac), FT(clay_frac), FT(mᵥ), FT(ρ))
end

function SoilMW(; sand_frac::Real = 0.2, clay_frac::Real = 0.1, mᵥ::Real = 0.35, ρ::Real = 1.7)
    return SoilMW(sand_frac, clay_frac, mᵥ, ρ)
end

function SoilMW{FT}(;
    sand_frac::Real = FT(0.2),
    clay_frac::Real = FT(0.1),
    mᵥ::Real = FT(0.35),
    ρ::Real = FT(1.7),
) where {FT<:Real}
    return SoilMW{FT}(FT(sand_frac), FT(clay_frac), FT(mᵥ), FT(ρ))
end

"""
    LeafUlabyElRayes1987{FT<:Real}(M_g, σ)
    LeafUlabyElRayes1987(; M_g = 0.5, σ = 1.27)

Dual-dispersion microwave leaf dielectric model of Ulaby & El-Rayes (1987),
combining a non-dispersive residual term, a free-water Debye component, and a
bound-water Cole–Cole-like component across `0.2–20 GHz`.

Pass to [`dielectric`](@ref) together with temperature `T` in `K` (accepted for
signature consistency) and frequency `f` in `GHz`.

# Fields
- `M_g::FT`: gravimetric leaf moisture content in `[0, 0.7]`.
- `σ::FT`: effective ionic conductivity of free water in `S m⁻¹` (typical `1.27`).

# Reference
Ulaby, F. T. & El-Rayes, M. A. (1987). *Microwave Dielectric Spectrum of
Vegetation – Part II: Dual-Dispersion Model.* IEEE TGRS 25(5), 550–557;
Ulaby & Long (2014), Chapter 11, Section 11-9.
"""
struct LeafUlabyElRayes1987{FT<:Real} <: AbstractVegetation
    M_g::FT
    σ::FT
end

function LeafUlabyElRayes1987(M_g::Real, σ::Real)
    FT = float(promote_type(typeof(M_g), typeof(σ)))
    return LeafUlabyElRayes1987{FT}(FT(M_g), FT(σ))
end

LeafUlabyElRayes1987(; M_g::Real = 0.5, σ::Real = 1.27) = LeafUlabyElRayes1987(M_g, σ)

function LeafUlabyElRayes1987{FT}(; M_g::Real = FT(0.5), σ::Real = FT(1.27)) where {FT<:Real}
    return LeafUlabyElRayes1987{FT}(FT(M_g), FT(σ))
end
