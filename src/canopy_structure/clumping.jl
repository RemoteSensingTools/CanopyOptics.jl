"""
    AbstractClumping{FT}

Angular clumping model for canopy interception.

Clumping modifies the effective projected area used in gap probability and
extinction calculations:

```math
G_{eff}(μ) = Ω(μ) G(μ).
```

It does not change the leaf or wood scattering matrix `Z`; it changes how often
canopy elements are encountered along a path.
"""
abstract type AbstractClumping{FT<:Real} end

_clumping_parameter_type(args...) = float(promote_type(map(typeof, args)...))

"""
    NoClumping()
    NoClumping{FT}()

Default clumping model with `Ω(μ) = 1`.
"""
struct NoClumping{FT<:Real} <: AbstractClumping{FT} end

NoClumping() = NoClumping{Float64}()

"""
    ConstantClumping(; Ω = 1)
    ConstantClumping(; Ω₀ = 1)

Direction-independent clumping index.

Values `Ω < 1` reduce effective extinction relative to a random leaf
distribution. `Ω = 1` preserves the unclumped Beer-Lambert limit. Values above
one represent more regular spacing and are allowed for completeness.
"""
struct ConstantClumping{FT<:Real} <: AbstractClumping{FT}
    Ω₀::FT
end

function ConstantClumping(; Ω = nothing, Ω₀ = nothing)
    value = Ω === nothing ? (Ω₀ === nothing ? 1.0 : Ω₀) : Ω
    FT = _clumping_parameter_type(value)
    return ConstantClumping{FT}(FT(value))
end

ConstantClumping{FT}(; Ω = nothing, Ω₀ = nothing) where {FT<:Real} = begin
    value = Ω === nothing ? (Ω₀ === nothing ? one(FT) : Ω₀) : Ω
    ConstantClumping{FT}(FT(value))
end

"""
    ChenLeblancClumping(; Ω₀ = 0.7, c = 2, e = 2)

Angular clumping model after Chen and Leblanc-style empirical forms:

```math
Ω(θ) = \\frac{Ω_0}{Ω_0 + (1 - Ω_0)\\exp(-c θ^e)}, \\qquad θ = \\arccos μ.
```

The default shape gives `Ω(0) = Ω₀` at nadir and approaches one toward the
horizon, where clumping is visually less apparent along the slant path.
"""
struct ChenLeblancClumping{FT<:Real} <: AbstractClumping{FT}
    Ω₀::FT
    c::FT
    e::FT
end

function ChenLeblancClumping(; Ω₀ = 0.7, c = 2.0, e = 2.0)
    FT = _clumping_parameter_type(Ω₀, c, e)
    return ChenLeblancClumping{FT}(FT(Ω₀), FT(c), FT(e))
end

ChenLeblancClumping{FT}(; Ω₀ = FT(0.7), c = FT(2), e = FT(2)) where {FT<:Real} =
    ChenLeblancClumping{FT}(FT(Ω₀), FT(c), FT(e))

"""
    clumping_index(model, μ)

Evaluate `Ω(μ)` for one direction cosine or an array of direction cosines.
"""
clumping_index(::NoClumping, μ::Real) = one(μ)

clumping_index(model::ConstantClumping, μ::Real) =
    model.Ω₀ + zero(μ)

function clumping_index(model::ChenLeblancClumping, μ::Real)
    μₚ = μ + zero(model.Ω₀) + zero(model.c) + zero(model.e)
    oneμ = one(μₚ)
    μ_clamped = min(max(μₚ, -oneμ), oneμ)
    θ = acos(μ_clamped)
    return model.Ω₀ / (model.Ω₀ + (one(model.Ω₀) - model.Ω₀) *
                       exp(-model.c * θ^model.e))
end

clumping_index(model::AbstractClumping, μ::AbstractArray) =
    clumping_index.(Ref(model), μ)

"""
    effective_G(clumping, G_value, μ)

Apply clumping to a Ross projection factor: `G_eff = Ω(μ) G`.

This helper is intended for propagation/extinction code. Z-matrix construction
should generally keep using the unclumped leaf-angle `G` normalization.
"""
effective_G(model::AbstractClumping, G_value::Real, μ::Real) =
    clumping_index(model, μ) * G_value

effective_G(model::AbstractClumping, G_values::AbstractArray, μ::AbstractArray) =
    effective_G.(Ref(model), G_values, μ)
