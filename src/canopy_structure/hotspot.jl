"""
    AbstractHotSpot{FT}

Model for bidirectional hotspot corrections to canopy gap probability.

Hotspot terms are propagation corrections, not leaf scattering matrices: they
modify the joint probability that the solar and viewing paths see the same
canopy gaps.  They should therefore be applied where direct-beam source terms
are assembled, not inside `compute_Z_matrices`.
"""
abstract type AbstractHotSpot{FT<:Real} end

_hotspot_parameter_type(args...) = float(promote_type(map(typeof, args)...))
@inline _hotspot_value(x) = x
@inline _hotspot_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)

"""
    NoHotSpot()
    NoHotSpot{FT}()

Default model with independent sun and view gaps. The joint gap probability is
`exp(-(k_s + k_o)L)`.
"""
struct NoHotSpot{FT<:Real} <: AbstractHotSpot{FT} end

NoHotSpot() = NoHotSpot{Float64}()

"""
    KuuskHotSpot(; h = 0.1)
    KuuskHotSpot{FT}(; h = FT(0.1))

Kuusk-style hotspot correction with dimensionless size parameter `h`, roughly
leaf size divided by canopy height. Larger `h` broadens the hotspot. `h <= 0`
disables the correction and returns the independent-gap probability.
"""
struct KuuskHotSpot{FT<:Real} <: AbstractHotSpot{FT}
    h::FT
end

function KuuskHotSpot(; h = 0.1)
    FT = _hotspot_parameter_type(h)
    return KuuskHotSpot{FT}(FT(h))
end

KuuskHotSpot{FT}(; h = FT(0.1)) where {FT<:Real} =
    KuuskHotSpot{FT}(FT(h))

"""
    canopy_extinction(G_value, μ)

Convert projected area `G(μ)` into the directional extinction coefficient
`k = G(μ) / μ` used by Beer-Lambert canopy gap probabilities.
"""
function canopy_extinction(G_value::Real, μ::Real)
    FT = _hotspot_parameter_type(G_value, μ)
    return FT(G_value) / FT(μ)
end

"""
    hotspot_separation(μ_s, μ_o, dϕ, h)

Angular separation parameter used by the Kuusk hotspot correction:

```math
α = \\frac{1}{h}\\sqrt{\\tan^2 θ_s + \\tan^2 θ_o
    - 2\\tan θ_s\\tan θ_o\\cos Δϕ}.
```
"""
function hotspot_separation(μ_s::Real, μ_o::Real, dϕ::Real, h::Real)
    FT = _hotspot_parameter_type(μ_s, μ_o, dϕ, h)
    μs = FT(μ_s)
    μo = FT(μ_o)
    h′ = FT(h)
    _hotspot_value(h′) <= 0 && return FT(Inf)

    tan_s = sqrt(max(zero(FT), one(FT) - μs^2)) / μs
    tan_o = sqrt(max(zero(FT), one(FT) - μo^2)) / μo
    Δ² = tan_s^2 + tan_o^2 - 2 * tan_s * tan_o * cos(FT(dϕ))
    return sqrt(max(zero(FT), Δ²)) / h′
end

"""
    hotspot_correction(model, k_s, k_o, μ_s, μ_o, dϕ, L)

Return the multiplicative correction `C_hs` applied to the independent joint
gap probability. `k_s` and `k_o` are directional extinction coefficients,
`μ_s` and `μ_o` are positive cosines for sun and view paths, `dϕ` is relative
azimuth in radians, and `L` is cumulative LAI or area index.
"""
function hotspot_correction(::NoHotSpot, k_s, k_o, μ_s, μ_o, dϕ, L)
    FT = _hotspot_parameter_type(k_s, k_o, μ_s, μ_o, dϕ, L)
    return one(FT)
end

function hotspot_correction(model::KuuskHotSpot, k_s, k_o, μ_s, μ_o, dϕ, L)
    FT = _hotspot_parameter_type(model.h, k_s, k_o, μ_s, μ_o, dϕ, L)
    h = FT(model.h)
    _hotspot_value(h) <= 0 && return one(FT)

    ks = FT(k_s)
    ko = FT(k_o)
    L′ = FT(L)
    α = hotspot_separation(FT(μ_s), FT(μ_o), FT(dϕ), h)
    αL = α * L′
    root_ksko = sqrt(ks * ko)

    exponent = if _hotspot_value(abs(αL)) < 1e-8
        root_ksko * L′
    else
        root_ksko / α * (one(FT) - exp(-αL))
    end
    return exp(exponent)
end

"""
    joint_gap_probability(model, k_s, k_o, μ_s, μ_o, dϕ, L)

Return the full correlated bidirectional gap probability

```math
P_{so}(L) = \\exp[-(k_s + k_o)L] C_{hs}(L).
```

Use [`NoHotSpot`](@ref) for the independent-path Beer-Lambert result.
"""
function joint_gap_probability(model::AbstractHotSpot, k_s, k_o,
                               μ_s, μ_o, dϕ, L)
    FT = _hotspot_parameter_type(k_s, k_o, μ_s, μ_o, dϕ, L)
    ks = FT(k_s)
    ko = FT(k_o)
    L′ = FT(L)
    base = exp(-(ks + ko) * L′)
    return base * hotspot_correction(model, ks, ko, FT(μ_s), FT(μ_o), FT(dϕ), L′)
end
