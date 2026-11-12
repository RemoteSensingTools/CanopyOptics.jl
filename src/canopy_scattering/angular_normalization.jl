"""
    AbstractAngularNormalization

Supertype for numerical angular normalization policies shared across bi-Lambertian
canopy scattering and leaf emission geometry.
"""
abstract type AbstractAngularNormalization end

"""
    ContinuousAngularNormalization()

Use the continuously normalized Lambertian angular law (the default).
Its integral on a finite outgoing-stream quadrature need not be exactly one.
"""
struct ContinuousAngularNormalization <: AbstractAngularNormalization end

"""
    DiscreteAngularNormalization(μ, weights; atol)

Normalize each leaf face's angular redistribution against the supplied positive-
hemisphere cosine quadrature. Both arrays must share a floating-point type;
`0 < μ ≤ 1`, nonnegative weights, and `abs(sum(weights)-1) ≤ atol` are required.
Zero-weight viewing/solar nodes are allowed. No weights are rescaled.

The constructor copies the host arrays. Treat the owned arrays as read-only and
reconstruct the object when angles or weights change. Scattering calls require
the exact same angle values, order and precision. Include the normalization and
both angular arrays in any downstream optical-cache identity.

This is an opt-in discretization of a bi-Lambertian angular PDF, not a change to
leaf albedo, projected area G, absorbed PAR or photon yield. It preserves the
scalar area kernel's reciprocity while making its zeroth-mode column integral
conservative to rounding accuracy. Finite angular errors and solver depth errors
remain. Specular/composite normalization is not provided by this policy.
"""
struct DiscreteAngularNormalization{FT<:AbstractFloat} <: AbstractAngularNormalization
    μ::Vector{FT}
    weights::Vector{FT}
    function DiscreteAngularNormalization(μ::AbstractVector{FT},weights::AbstractVector{FT};
                                          atol::FT) where {FT<:AbstractFloat}
        n = length(μ)
        n == length(weights) > 0 || throw(DimensionMismatch("Angular nodes and weights must align"))
        isfinite(atol)&&atol>=zero(FT) || throw(ArgumentError("Angular tolerance must be finite and nonnegative"))
        all(x->isfinite(x)&&zero(FT)<x<=one(FT),μ) || throw(ArgumentError("Positive-hemisphere cosines must lie in (0,1]"))
        all(x->isfinite(x)&&x>=zero(FT),weights) || throw(ArgumentError("Angular weights must be finite and nonnegative"))
        abs(sum(weights)-one(FT))<=atol || throw(ArgumentError("Hemisphere weights must integrate unity without rescaling"))
        μ_own = copyto!(Vector{FT}(undef, n), μ)
        w_own = Vector{FT}(undef, n)
        for (k, w) in enumerate(weights)
            w_own[k] = iszero(w) ? zero(FT) : w
        end
        new{FT}(μ_own, w_own)
    end
end

Base.:(==)(a::DiscreteAngularNormalization, b::DiscreteAngularNormalization) =
    a.μ == b.μ && a.weights == b.weights

Base.hash(a::DiscreteAngularNormalization, h::UInt) =
    hash(a.weights, hash(a.μ, hash(:DiscreteAngularNormalization, h)))

_check_angular_normalization(::ContinuousAngularNormalization,μ) = nothing
function _check_angular_normalization(normalization::DiscreteAngularNormalization{FT},μ) where {FT}
    eltype(μ)===FT && length(μ)==length(normalization.μ) &&
        all(splat(==), zip(μ, normalization.μ)) ||
        throw(ArgumentError("Scattering angles must match the normalization grid, order and precision"))
    nothing
end

"""
    leaf_angular_scale(normalization, μ_leaf)

Host geometry factor for one leaf inclination. The continuous convention returns
one. The discrete convention returns `1/H`, where
`H = 2sum(wᵢ * (P₀(μᵢ, μ_leaf) + N₀(μᵢ, μ_leaf)))` is the quadrature integral
of a single Lambertian leaf face. `μ_leaf` lies in `[0,1]` and uses the measure's
floating-point precision. Multiply each leaf's scattering contribution or source
projection by this factor once; do not change G, leaf area or physiological yield.

Precompute this host geometry before copying optical/source arrays to a GPU.
"""
function leaf_angular_scale(::ContinuousAngularNormalization,μ_leaf::FT) where {FT<:Real}
    isfinite(μ_leaf)&&zero(FT)<=μ_leaf<=one(FT) || throw(ArgumentError("Leaf inclination cosine must lie in [0,1]"))
    one(FT)
end

function leaf_angular_scale(normalization::DiscreteAngularNormalization{FT},μ_leaf::FT) where {FT}
    isfinite(μ_leaf)&&zero(FT)<=μ_leaf<=one(FT) || throw(ArgumentError("Leaf inclination cosine must lie in [0,1]"))
    integral = zero(FT)
    P,N = zeros(FT,1),zeros(FT,1)
    for i in eachindex(normalization.μ)
        _one_sided_projection_moments!(P,N,normalization.μ[i],μ_leaf,0)
        integral += normalization.weights[i]*(P[1]+N[1])
    end
    H = FT(2)*integral
    isfinite(H)&&H>zero(FT) || throw(ArgumentError("Leaf angular law has no finite positive discrete measure"))
    inv(H)
end

_leaf_angular_scale_for_Z(normalization, μ_leaf) = leaf_angular_scale(normalization, μ_leaf)

function _leaf_angular_scale_for_Z(
    normalization::DiscreteAngularNormalization{FTN},
    μ_leaf::FTL,
) where {FTN<:AbstractFloat,FTL<:Real}
    isfinite(μ_leaf) && zero(FTL) <= μ_leaf <= one(FTL) ||
        throw(ArgumentError("Leaf inclination cosine must lie in [0,1]"))
    T = float(promote_type(FTN, FTL))
    integral = zero(T)
    P, N = zeros(T, 1), zeros(T, 1)
    μ_l = T(μ_leaf)
    for i in eachindex(normalization.μ)
        _one_sided_projection_moments!(P, N, T(normalization.μ[i]), μ_l, 0)
        integral += T(normalization.weights[i]) * (P[1] + N[1])
    end
    H = T(2) * integral
    isfinite(H) && H > zero(T) ||
        throw(ArgumentError("Leaf angular law has no finite positive discrete measure"))
    return inv(H)
end

"""
    leaf_projection_moments(μ, μ_leaf, m_max;
                            normalization=ContinuousAngularNormalization())

Return `(P, N)`, cosine Fourier moments of the positive/negative leaf-face
projections for orders `0:m_max`. Ray cosine `μ` is signed; the leaf normal points
upwards with inclination cosine `μ_leaf ∈ [0,1]`. Both use the same precision.
The continuous moments contain `1/(2π) ∫ projection*cos(mϕ) dϕ`.

For a two-face Lambertian source, `2P*front + 2N*back` gives its azimuth-integrated
source coefficients before leaf-normal azimuth factors. The optional discrete
normalization applies [`leaf_angular_scale`](@ref) once to both arrays, consistently
with bi-Lambertian scattering on that measure. The azimuth convention and μ path
factor of a particular RT solver belong in its adapter, not in these moments.
"""
function leaf_projection_moments(μ::FT,μ_leaf::FT,m_max::Integer;
        normalization::AbstractAngularNormalization=ContinuousAngularNormalization()) where {FT<:Real}
    isfinite(μ)&&abs(μ)<=one(FT) || throw(ArgumentError("Ray cosine must lie in [-1,1]"))
    m_max>=0 || throw(ArgumentError("Fourier order must be nonnegative"))
    scale = leaf_angular_scale(normalization,μ_leaf)
    P,N = _one_sided_projection_moments(μ,μ_leaf,Int(m_max))
    P .*= scale
    N .*= scale
    P,N
end
