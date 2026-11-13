@doc raw"""
    DirectionCosine{FT}(μ, ϕ)
    DirectionCosine(; μ = nothing, θ = nothing, ϕ = 0)

Represent a unit spherical direction by polar cosine ``\mu = \cos\theta \in [-1, 1]``
and azimuth angle `ϕ` in radians.

Passing keyword `θ` instead of `μ` computes `μ = cos(θ)` (see also
[`direction_from_angles`](@ref)).

# Fields
- `μ::FT`: polar cosine ``\cos\theta`` in ``[-1, 1]``.
- `ϕ::FT`: azimuth angle in radians.

# Examples
```jldoctest
julia> DirectionCosine(0.8, 0.0)
DirectionCosine{Float64}(0.8, 0.0)
```
"""
struct DirectionCosine{FT<:Real}
    μ::FT
    ϕ::FT

    DirectionCosine{FT}(μ, ϕ) where {FT<:Real} = new{FT}(FT(μ), FT(ϕ))
end

function DirectionCosine(μ::Real, ϕ::Real)
    FT = _canopy_parameter_type(μ, ϕ)
    return DirectionCosine{FT}(FT(μ), FT(ϕ))
end

function DirectionCosine(; μ = nothing, θ = nothing, ϕ::Real = 0)
    (μ === nothing) == (θ === nothing) &&
        throw(ArgumentError("Specify exactly one of keyword `μ` or `θ`"))
    return μ !== nothing ? DirectionCosine(μ, ϕ) : direction_from_angles(θ, ϕ)
end

"""
    direction_from_angles(θ, ϕ)

Construct a [`DirectionCosine`](@ref) from polar/zenith angle `θ` and azimuth
angle `ϕ` (both in radians), setting `μ = cos(θ)`.

# Examples
```jldoctest
julia> direction_from_angles(0.0, 1.2)
DirectionCosine{Float64}(1.0, 1.2)
```
"""
function direction_from_angles(θ::Real, ϕ::Real)
    FT = _canopy_parameter_type(θ, ϕ)
    return DirectionCosine{FT}(cos(FT(θ)), FT(ϕ))
end

Base.@deprecate_binding dirVector_μ DirectionCosine
Base.Docs.doc!(
    @__MODULE__,
    Base.Docs.Binding(@__MODULE__, :dirVector_μ),
    Base.Docs.docstr(
        """
            dirVector_μ(μ, ϕ)

        Deprecated alias for [`DirectionCosine`](@ref).
        """,
        Dict{Symbol,Any}(:path => @__FILE__, :linenumber => @__LINE__, :module => @__MODULE__),
    ),
)

"""
    dirVector(θ, ϕ)

Deprecated constructor from polar angle `θ` and azimuth `ϕ`; use
[`direction_from_angles`](@ref) instead.
"""
function dirVector(θ::Real, ϕ::Real)
    Base.depwarn("`dirVector(θ, ϕ)` is deprecated, use `direction_from_angles(θ, ϕ)` instead.", :dirVector)
    return direction_from_angles(θ, ϕ)
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

function dot(Ω₁::DirectionCosine, Ω₂::DirectionCosine)
    FT = float(promote_type(typeof(Ω₁.μ), typeof(Ω₁.ϕ), typeof(Ω₂.μ), typeof(Ω₂.ϕ)))
    μ₁ = FT(Ω₁.μ)
    μ₂ = FT(Ω₂.μ)
    s₁ = _unit_sin_from_cos(μ₁)
    s₂ = _unit_sin_from_cos(μ₂)
    return μ₁ * μ₂ + s₁ * s₂ * cos(FT(Ω₂.ϕ) - FT(Ω₁.ϕ))
end

@doc raw"""
    dot_product(μ₁, μ₂, dϕ)

Compute the Cartesian inner product
``\mu_1 \mu_2 + \sqrt{1 - \mu_1^2}\sqrt{1 - \mu_2^2}\cos(\Delta\phi)``
between two unit directions parameterized by polar cosines ``\mu_1, \mu_2 \in [-1, 1]``
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
