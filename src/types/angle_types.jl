"""
    dirVector{FT}(θ, ϕ)

Represent a unit spherical direction by polar/zenith angle `θ` and azimuth
angle `ϕ` (both in radians).

# Fields
- `θ::FT`: polar/zenith angle in radians.
- `ϕ::FT`: azimuth angle in radians.

# Examples
```jldoctest
julia> dirVector(0.5, 1.2)
dirVector{Float64}(0.5, 1.2)
```
"""
struct dirVector{FT}
    θ::FT
    ϕ::FT
end

function dirVector(θ::Real, ϕ::Real)
    FT = promote_type(typeof(θ), typeof(ϕ))
    return dirVector{FT}(FT(θ), FT(ϕ))
end

"""
    dirVector_μ{FT}(μ, ϕ)

Represent a unit spherical direction by polar cosine `μ = cos(θ) ∈ [-1, 1]` and
azimuth angle `ϕ` in radians.

# Fields
- `μ::FT`: polar cosine `cos(θ)` in `[-1, 1]`.
- `ϕ::FT`: azimuth angle in radians.

# Examples
```jldoctest
julia> dirVector_μ(0.8, 0.0)
dirVector_μ{Float64}(0.8, 0.0)
```
"""
struct dirVector_μ{FT}
    μ::FT
    ϕ::FT
end

function dirVector_μ(μ::Real, ϕ::Real)
    FT = promote_type(typeof(μ), typeof(ϕ))
    return dirVector_μ{FT}(FT(μ), FT(ϕ))
end

@inline function _unit_sin_from_cos(μ::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    u = (one(FT) - μ) * (one(FT) + μ)
    return prospect_primal(u) <= zero(BFT) ? zero(FT) : sqrt(u)
end

@inline function _planar_tan_distance(tan_s::FT, tan_o::FT, dϕ::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    s_dϕ = sin(dϕ)
    s_half = sin(dϕ / (one(FT) + one(FT)))
    dx = (tan_s - tan_o) + (one(FT) + one(FT)) * tan_o * (s_half * s_half)
    dy = tan_o * s_dϕ
    if iszero(prospect_primal(dy)) && _same_val_and_partials(dy, zero(FT))
        return abs(dx)
    end
    d2 = dx * dx + dy * dy
    return prospect_primal(d2) <= zero(BFT) ? zero(FT) : sqrt(d2)
end

dot(Ω₁::dirVector, Ω₂::dirVector) =
    cos(Ω₁.θ) * cos(Ω₂.θ) + sin(Ω₁.θ) * sin(Ω₂.θ) * cos(Ω₂.ϕ - Ω₁.ϕ)

function dot(Ω₁::dirVector_μ, Ω₂::dirVector_μ)
    FT = float(promote_type(typeof(Ω₁.μ), typeof(Ω₁.ϕ), typeof(Ω₂.μ), typeof(Ω₂.ϕ)))
    μ₁ = FT(Ω₁.μ)
    μ₂ = FT(Ω₂.μ)
    s₁ = _unit_sin_from_cos(μ₁)
    s₂ = _unit_sin_from_cos(μ₂)
    return μ₁ * μ₂ + s₁ * s₂ * cos(FT(Ω₂.ϕ) - FT(Ω₁.ϕ))
end

"""
    dot_product(μ₁, μ₂, dϕ)

Compute the Cartesian inner product `μ₁ μ₂ + √(1 - μ₁²) √(1 - μ₂²) cos(dϕ)`
between two unit directions parameterized by polar cosines `μ₁, μ₂ ∈ [-1, 1]`
and relative azimuth `dϕ` in radians.
"""
function dot_product(μ₁_in::Real, μ₂_in::Real, dϕ_in::Real)
    FT = float(promote_type(typeof(μ₁_in), typeof(μ₂_in), typeof(dϕ_in)))
    μ₁ = FT(μ₁_in)
    μ₂ = FT(μ₂_in)
    s₁ = _unit_sin_from_cos(μ₁)
    s₂ = _unit_sin_from_cos(μ₂)
    return μ₁ * μ₂ + s₁ * s₂ * cos(FT(dϕ_in))
end
