"""
    planophile_leaves(FT = Float64)

Construct a planophile (predominantly horizontal) [`LeafDistribution`](@ref)
from a `Beta` distribution with mean inclination `26.76°` and standard
deviation `18.51°` (Bonan, 2019, Table 2.1).

# Examples
```jldoctest
julia> planophile_leaves(Float64) isa LeafDistribution{Float64}
true
```
"""
function planophile_leaves(::Type{FT} = Float64) where {FT<:Real}
    α, β = βparameters(deg2rad(FT(26.76)), deg2rad(FT(18.51)))
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2) / FT(π))
end

@doc raw"""
    planophile_leaves2(FT = Float64)

Construct a `Beta(1, 2.768)` planophile [`LeafDistribution`](@ref)
approximating Bunnik's ``g_L(\theta_L) = \frac{2}{\pi}(1 + \cos 2\theta_L)``
distribution (Shultis & Myneni, 1988, Eq. 59).
"""
function planophile_leaves2(::Type{FT} = Float64) where {FT<:Real}
    return LeafDistribution(Beta(one(FT), FT(2.768)), FT(2) / FT(π))
end

@doc raw"""
    uniform_leaves(FT = Float64)

Construct a uniform [`LeafDistribution`](@ref) (``g_L(\theta_L) = 2/\pi`` over
``\theta_L \in [0, \pi/2]``).
"""
function uniform_leaves(::Type{FT} = Float64) where {FT<:Real}
    return LeafDistribution(Uniform(zero(FT), one(FT)), FT(2) / FT(π))
end

"""
    plagiophile_leaves(FT = Float64)

Construct a plagiophile (predominantly `45°`) [`LeafDistribution`](@ref) from a
`Beta` distribution with mean inclination `45.0°` and standard deviation
`16.27°` (Bonan, 2019, Table 2.1).
"""
function plagiophile_leaves(::Type{FT} = Float64) where {FT<:Real}
    α, β = βparameters(deg2rad(FT(45)), deg2rad(FT(16.27)))
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2) / FT(π))
end

"""
    erectophile_leaves(FT = Float64)

Construct an erectophile (predominantly vertical) [`LeafDistribution`](@ref)
from a `Beta` distribution with mean inclination `63.24°` and standard
deviation `18.5°` (Bonan, 2019, Table 2.1).
"""
function erectophile_leaves(::Type{FT} = Float64) where {FT<:Real}
    α, β = βparameters(deg2rad(FT(63.24)), deg2rad(FT(18.5)))
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2) / FT(π))
end

"""
    spherical_leaves(FT = Float64)

Construct a spherical [`LeafDistribution`](@ref) (`G(μ) ≈ 0.5`) from a `Beta`
distribution with mean inclination `57.3°` and standard deviation `21.55°`
(Bonan, 2019, Table 2.1).

# Examples
```jldoctest
julia> isapprox(only(G([0.5], spherical_leaves())), 0.5; atol = 2e-2)
true
```
"""
function spherical_leaves(::Type{FT} = Float64) where {FT<:Real}
    α, β = βparameters(deg2rad(FT(57.3)), deg2rad(FT(21.55)))
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2) / FT(π))
end

"""
    flat_leaves(FT = Float64)

Construct a nearly horizontal [`LeafDistribution`](@ref) (`Beta(1, 100)`)
peaked sharply near `θ_L ≈ 0`.
"""
function flat_leaves(::Type{FT} = Float64) where {FT<:Real}
    return LeafDistribution(Beta(one(FT), FT(100)), FT(2) / FT(π))
end

"""
    beta_leaves(α, β)

Construct a [`LeafDistribution`](@ref) over `[0, π/2]` from arbitrary `Beta(α, β)`
shape parameters (corresponding to `(ν, μ)` in the RAMI benchmark convention).

# Examples
```jldoctest
julia> beta_leaves(2.0, 2.0) isa LeafDistribution{Float64}
true
```
"""
function beta_leaves(α_in::Real, β_in::Real)
    FT = float(promote_type(typeof(α_in), typeof(β_in)))
    return LeafDistribution(Beta(FT(α_in), FT(β_in)), FT(2) / FT(π))
end
