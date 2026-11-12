"""
    prospect_primal(x)

Extract the primal scalar value of `x`, stripping nested `ForwardDiff.Dual`
wrappers while returning non-dual numbers unchanged.

# Examples
```jldoctest
julia> CanopyOptics.prospect_primal(1.42)
1.42
```
"""
@inline prospect_primal(x) = x
@inline prospect_primal(x::ForwardDiff.Dual) = prospect_primal(ForwardDiff.value(x))

@inline _leaf_base_float(::Type{T}) where {T<:AbstractFloat} = T
@inline _leaf_base_float(::Type{ForwardDiff.Dual{Tag,V,N}}) where {Tag,V,N} = _leaf_base_float(V)
@inline _leaf_base_float(::Type{T}) where {T<:Real} = float(T)
@inline _leaf_base_float(x::Real) = _leaf_base_float(typeof(x))

struct _UnsetLeafKw end
const _UNSET_LEAF_KW = _UnsetLeafKw()
@inline _kw_type(::_UnsetLeafKw) = Union{}
@inline _kw_type(x::Real) = typeof(x)

@inline function _resolve_pro_dry_matter(::Type{FT}, Cm, Cprot, Ccbc, ::AbstractString) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    if Cm === _UNSET_LEAF_KW
        cm_v = zero(FT)
        cprot_v = Cprot === _UNSET_LEAF_KW ? FT(BFT(0.001)) : FT(Cprot)
        ccbc_v = Ccbc === _UNSET_LEAF_KW ? FT(BFT(0.009)) : FT(Ccbc)
    else
        cm_v = FT(Cm)
        cprot_v = Cprot === _UNSET_LEAF_KW ? zero(FT) : FT(Cprot)
        ccbc_v = Ccbc === _UNSET_LEAF_KW ? zero(FT) : FT(Ccbc)
    end
    if prospect_primal(cm_v) > zero(BFT) &&
       (prospect_primal(cprot_v) > zero(BFT) || prospect_primal(ccbc_v) > zero(BFT))
        throw(ArgumentError(
            "Cannot set both bulk dry matter Cm > 0 and Cprot/Ccbc > 0 (double-counts LMA); omit Cm when using Cprot + Ccbc or omit Cprot/Ccbc when using Cm",
        ))
    end
    return cm_v, cprot_v, ccbc_v
end

"""
    AbstractCanopyScatteringType{FT<:Real}

Supertype for leaf, wood, composite, and multi-component canopy angular
scattering models projected onto discrete-ordinate Fourier `Z` matrices.
"""
abstract type AbstractCanopyScatteringType{FT<:Real} end

"""
    CanopyQuadrature(n_leaf, n_azimuth = 48)
    CanopyQuadrature(; n_leaf = 64, n_azimuth = 48, nQuad = nothing)

Quadrature node counts used when projecting canopy scattering models onto
Fourier `Z` matrices.

`n_leaf` sets Gauss–Legendre nodes over leaf inclination `[0, π/2]`, and
`n_azimuth` sets outgoing-azimuth quadrature nodes for specular, coupled-BSDF,
and reference azimuth integrations.

# Examples
```jldoctest
julia> q = CanopyQuadrature(n_leaf = 32, n_azimuth = 24)
CanopyQuadrature(32, 24)
```
"""
struct CanopyQuadrature
    n_leaf::Int
    n_azimuth::Int

    function CanopyQuadrature(n_leaf::Integer, n_azimuth::Integer)
        n_leaf > 0 || throw(ArgumentError("n_leaf must be positive"))
        n_azimuth > 0 || throw(ArgumentError("n_azimuth must be positive"))
        return new(Int(n_leaf), Int(n_azimuth))
    end
end

CanopyQuadrature(n_leaf::Integer) = CanopyQuadrature(n_leaf, 48)

function CanopyQuadrature(; n_leaf::Integer = 64,
                          n_azimuth::Integer = 48,
                          nQuad = nothing)
    if nQuad !== nothing
        n_leaf = nQuad
        n_azimuth = nQuad
    end
    return CanopyQuadrature(n_leaf, n_azimuth)
end

_canopy_parameter_type(args...) = float(promote_type(map(typeof, args)...))

function _bilambertian_canopy_scattering(R::Real, T::Real)
    R_prom, T_prom = promote(R, T)
    FT = _canopy_parameter_type(R_prom, T_prom)
    return BiLambertianCanopyScattering{FT}(FT(R_prom), FT(T_prom))
end

function _specular_canopy_scattering(nᵣ::Real, κ::Real)
    nᵣ_prom, κ_prom = promote(nᵣ, κ)
    FT = _canopy_parameter_type(nᵣ_prom, κ_prom)
    return SpecularCanopyScattering{FT}(FT(nᵣ_prom), FT(κ_prom))
end

"""
    BiLambertianCanopyScattering{FT<:Real}(R, T)
    BiLambertianCanopyScattering(; R = 0.3, T = 0.1)

Face-symmetric diffuse bi-Lambertian canopy leaf scattering model parameterized
by hemispherical reflectance `R` and transmittance `T`.

# Fields
- `R::FT`: leaf hemispherical reflectance.
- `T::FT`: leaf hemispherical transmittance.

# Examples
```jldoctest
julia> BiLambertianCanopyScattering(0.42, 0.38)
BiLambertianCanopyScattering{Float64}(0.42, 0.38)
```
"""
struct BiLambertianCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}
    R::FT
    T::FT
end

BiLambertianCanopyScattering(R::Real, T::Real) =
    _bilambertian_canopy_scattering(R, T)

BiLambertianCanopyScattering(R::Integer, T::Integer) =
    _bilambertian_canopy_scattering(R, T)

function BiLambertianCanopyScattering(; R = 0.3, T = 0.1, nQuad = nothing)
    return _bilambertian_canopy_scattering(R, T)
end

function BiLambertianCanopyScattering{FT}(; R = FT(0.3), T = FT(0.1),
                                          nQuad = nothing) where {FT<:Real}
    return BiLambertianCanopyScattering{FT}(FT(R), FT(T))
end

"""
    SpecularCanopyScattering{FT<:Real}(nᵣ, κ)
    SpecularCanopyScattering(; nᵣ = 1.5, κ = 0.5)

Specular leaf-surface canopy scattering model parameterized by cuticular
refractive index `nᵣ` and Nilson–Kuusk roughness attenuation `κ`.

# Fields
- `nᵣ::FT`: leaf surface refractive index.
- `κ::FT`: surface roughness attenuation parameter (`κ >= 0`).

# Examples
```jldoctest
julia> SpecularCanopyScattering(nᵣ = 1.45, κ = 0.2)
SpecularCanopyScattering{Float64}(1.45, 0.2)
```
"""
struct SpecularCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}
    nᵣ::FT
    κ::FT
end

SpecularCanopyScattering(nᵣ::Real, κ::Real) =
    _specular_canopy_scattering(nᵣ, κ)

SpecularCanopyScattering(nᵣ::Integer, κ::Integer) =
    _specular_canopy_scattering(nᵣ, κ)

function SpecularCanopyScattering(; nᵣ = 1.5, κ = 0.5, nQuad = nothing)
    return _specular_canopy_scattering(nᵣ, κ)
end

function SpecularCanopyScattering{FT}(; nᵣ = FT(1.5), κ = FT(0.5),
                                      nQuad = nothing) where {FT<:Real}
    return SpecularCanopyScattering{FT}(FT(nᵣ), FT(κ))
end

"""
    CompositeCanopyScattering(components...)

Combine multiple [`AbstractCanopyScatteringType`](@ref) models additively on a
shared leaf-angle distribution. Components are evaluated independently and their
`Z` matrices are summed; `a + b` constructs `CompositeCanopyScattering(a, b)`.
"""
struct CompositeCanopyScattering{FT<:Real,T<:Tuple} <: AbstractCanopyScatteringType{FT}
    components::T
end

_canopy_scattering_ft(::AbstractCanopyScatteringType{FT}) where {FT} = FT
_flatten_component(c::CompositeCanopyScattering) = c.components
_flatten_component(c::AbstractCanopyScatteringType) = (c,)

function CompositeCanopyScattering(components::AbstractCanopyScatteringType...)
    isempty(components) && throw(ArgumentError("CompositeCanopyScattering needs at least one component"))
    flat = ()
    for component in components
        flat = (flat..., _flatten_component(component)...)
    end
    FT = promote_type(map(_canopy_scattering_ft, flat)...)
    return CompositeCanopyScattering{FT,typeof(flat)}(flat)
end

Base.:+(a::AbstractCanopyScatteringType, b::AbstractCanopyScatteringType) =
    CompositeCanopyScattering(a, b)

"""
    AbstractLeafDistribution{FT<:Real}

Supertype for canopy element inclination angle distributions over `[0, π/2]`.
"""
abstract type AbstractLeafDistribution{FT<:Real} end

"""
    LeafDistribution{FT<:Real,D<:UnivariateDistribution}(LD, scaling)
    LeafDistribution(LD, scaling)

Represent a leaf inclination angle distribution `g_L(θ_L)` on `θ_L ∈ [0, π/2]`
via a `Distributions.UnivariateDistribution` `LD` on `u = 2 θ_L / π ∈ [0, 1]`
and Jacobian normalization factor `scaling` (typically `2/π`).

# Fields
- `LD::D`: univariate probability distribution on `[0, 1]`.
- `scaling::FT`: density normalization factor mapping `[0, 1]` to `[0, π/2]`.
"""
struct LeafDistribution{FT<:Real,D<:UnivariateDistribution} <: AbstractLeafDistribution{FT}
    LD::D
    scaling::FT
end

function LeafDistribution{FT}(LD::D, scaling::Real) where {FT<:Real,D<:UnivariateDistribution}
    return LeafDistribution{FT,D}(LD, FT(scaling))
end

Base.:(==)(a::LeafDistribution, b::LeafDistribution) =
    a.LD == b.LD && a.scaling == b.scaling

function Base.hash(a::LeafDistribution, h::UInt)
    h = hash(:LeafDistribution, h)
    if applicable(params, a.LD)
        h = hash(nameof(typeof(a.LD)), h)
        for p in params(a.LD)
            p_canon = (p isa Real && iszero(p)) ? zero(p) : p
            h = hash(p_canon, h)
        end
    else
        h = hash(a.LD, h)
    end
    s_canon = iszero(a.scaling) ? zero(a.scaling) : a.scaling
    return hash(s_canon, h)
end

"""
    AbstractLeafProperties

Supertype for leaf biochemical and structural state containers accepted by
leaf optical budget solvers.
"""
abstract type AbstractLeafProperties end

"""
    LeafProspectProProperties{FT}(; N=1.4, Ccab=40.0, Ccar=10.0, Canth=0.5,
                                  Cbrown=0.0, Cw=0.012, Cm=..., Cprot=..., Ccbc=...)

Store leaf biochemical and structural state for PROSPECT-D (bulk dry matter
`Cm` mode) and PROSPECT-PRO (`Cprot + Ccbc` dry-matter partition mode).

When constructed with keywords:
- Omitting `Cm` selects PROSPECT-PRO mode (`Cm = 0`, defaulting `Cprot = 0.001`
  and `Ccbc = 0.009`).
- Supplying `Cm` explicitly (including `Cm = 0` or a `ForwardDiff.Dual`) selects
  PROSPECT-D mode (`Cprot = 0`, `Ccbc = 0`) so dry matter is not double-counted.
- Setting both `Cm > 0` and `Cprot > 0` (or `Ccbc > 0`) throws an `ArgumentError`.

# Fields
- `N::FT`: leaf mesophyll structure parameter (`N >= 1`).
- `Ccab::FT`: chlorophyll *a*+*b* content in `µg cm⁻²`.
- `Ccar::FT`: carotenoid content in `µg cm⁻²`.
- `Canth::FT`: anthocyanin content in `µg cm⁻²` (Féret et al., 2017/2021).
- `Cbrown::FT`: senescent brown pigment content (arbitrary units).
- `Cw::FT`: equivalent water thickness in `cm`.
- `Cm::FT`: bulk dry matter content (LMA) in `g cm⁻²` (`0` in `Cprot + Ccbc` mode).
- `Cprot::FT`: protein content in `g cm⁻²`.
- `Ccbc::FT`: carbon-based constituents content in `g cm⁻²`.

# Examples
```jldoctest
julia> leaf = LeafProspectProProperties(N = 1.6, Ccab = 45.0, Cm = 0.008);

julia> (leaf.Cm, leaf.Cprot, leaf.Ccbc)
(0.008, 0.0, 0.0)
```
"""
struct LeafProspectProProperties{FT} <: AbstractLeafProperties
    N::FT
    Ccab::FT
    Ccar::FT
    Canth::FT
    Cbrown::FT
    Cw::FT
    Cm::FT
    Cprot::FT
    Ccbc::FT
end

function LeafProspectProProperties{FT}(;
    N = _UNSET_LEAF_KW,
    Ccab = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
) where {FT}
    BFT = _leaf_base_float(FT)
    n_v = N === _UNSET_LEAF_KW ? FT(BFT(1.4)) : FT(N)
    cab_v = Ccab === _UNSET_LEAF_KW ? FT(BFT(40.0)) : FT(Ccab)
    car_v = Ccar === _UNSET_LEAF_KW ? FT(BFT(10.0)) : FT(Ccar)
    ant_v = Canth === _UNSET_LEAF_KW ? FT(BFT(0.5)) : FT(Canth)
    brn_v = Cbrown === _UNSET_LEAF_KW ? zero(FT) : FT(Cbrown)
    cw_v = Cw === _UNSET_LEAF_KW ? FT(BFT(0.012)) : FT(Cw)
    cm_v, cprot_v, ccbc_v = _resolve_pro_dry_matter(FT, Cm, Cprot, Ccbc, "LeafProspectProProperties")
    return LeafProspectProProperties{FT}(
        n_v, cab_v, car_v, ant_v, brn_v, cw_v, cm_v, cprot_v, ccbc_v,
    )
end

function LeafProspectProProperties(;
    N = _UNSET_LEAF_KW,
    Ccab = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
)
    T_raw = promote_type(
        _kw_type(N), _kw_type(Ccab), _kw_type(Ccar), _kw_type(Canth),
        _kw_type(Cbrown), _kw_type(Cw), _kw_type(Cm), _kw_type(Cprot), _kw_type(Ccbc),
    )
    FT = T_raw === Union{} ? Float64 : float(T_raw)
    return LeafProspectProProperties{FT}(;
        N, Ccab, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
    )
end

"""
    AbstractProspectProperties

Supertype for wavelength-gridded leaf refractive-index and specific-absorption
cross-section tables ([`LeafOpticalProperties`](@ref) and
[`MultiPigmentOpticalProperties`](@ref)).
"""
abstract type AbstractProspectProperties end

"""
    LeafOpticalProperties{FT,FT2}(λ, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)

Store spectral grid `λ`, aqueous mesophyll refractive index `nᵣ`, and
PROSPECT-PRO specific absorption coefficient vectors. Construct from the bundled
artifact via [`createLeafOpticalStruct`](@ref).

# Fields
- `λ::FT2`: spectral bin centers (with or without `Unitful` units).
- `nᵣ::Vector{FT}`: refractive index of aqueous leaf mesophyll.
- `Kcab::Vector{FT}`: pooled chlorophyll *a*+*b* specific absorption (`cm² µg⁻¹`).
- `Kcar::Vector{FT}`: carotenoid specific absorption (`cm² µg⁻¹`).
- `Kant::Vector{FT}`: anthocyanin specific absorption (`cm² µg⁻¹`).
- `Kb::Vector{FT}`: senescent brown pigment specific absorption (arbitrary units).
- `Kw::Vector{FT}`: liquid water specific absorption (`cm⁻¹`).
- `Km::Vector{FT}`: bulk dry matter specific absorption (`cm² g⁻¹`).
- `Kp::Vector{FT}`: protein specific absorption (`cm² g⁻¹`).
- `Kcbc::Vector{FT}`: carbon-based constituent specific absorption (`cm² g⁻¹`).
"""
struct LeafOpticalProperties{FT,FT2} <: AbstractProspectProperties
    λ::FT2
    nᵣ::Vector{FT}
    Kcab::Vector{FT}
    Kcar::Vector{FT}
    Kant::Vector{FT}
    Kb::Vector{FT}
    Kw::Vector{FT}
    Km::Vector{FT}
    Kp::Vector{FT}
    Kcbc::Vector{FT}

    function LeafOpticalProperties{FT,FT2}(
        λ::FT2,
        nᵣ::AbstractVector,
        Kcab::AbstractVector,
        Kcar::AbstractVector,
        Kant::AbstractVector,
        Kb::AbstractVector,
        Kw::AbstractVector,
        Km::AbstractVector,
        Kp::AbstractVector,
        Kcbc::AbstractVector,
    ) where {FT,FT2}
        n = length(nᵣ)
        (length(λ) == n && length(Kcab) == n && length(Kcar) == n && length(Kant) == n &&
         length(Kb) == n && length(Kw) == n && length(Km) == n && length(Kp) == n && length(Kcbc) == n) ||
            throw(DimensionMismatch("All LeafOpticalProperties spectral vectors must match length(nᵣ) = $n"))
        conv(v) = v isa Vector{FT} ? v : copyto!(Vector{FT}(undef, length(v)), v)
        return new{FT,FT2}(
            λ, conv(nᵣ), conv(Kcab), conv(Kcar), conv(Kant),
            conv(Kb), conv(Kw), conv(Km), conv(Kp), conv(Kcbc),
        )
    end
end

function LeafOpticalProperties(
    λ::FT2,
    nᵣ::AbstractVector,
    Kcab::AbstractVector,
    Kcar::AbstractVector,
    Kant::AbstractVector,
    Kb::AbstractVector,
    Kw::AbstractVector,
    Km::AbstractVector,
    Kp::AbstractVector,
    Kcbc::AbstractVector,
) where {FT2}
    FT = float(promote_type(
        eltype(nᵣ), eltype(Kcab), eltype(Kcar), eltype(Kant),
        eltype(Kb), eltype(Kw), eltype(Km), eltype(Kp), eltype(Kcbc),
    ))
    return LeafOpticalProperties{FT,FT2}(λ, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)
end

function Base.:(==)(a::LeafOpticalProperties, b::LeafOpticalProperties)
    return a.λ == b.λ && a.nᵣ == b.nᵣ && a.Kcab == b.Kcab && a.Kcar == b.Kcar &&
           a.Kant == b.Kant && a.Kb == b.Kb && a.Kw == b.Kw && a.Km == b.Km &&
           a.Kp == b.Kp && a.Kcbc == b.Kcbc
end

function Base.hash(a::LeafOpticalProperties, h::UInt)
    h = hash(:LeafOpticalProperties, h)
    for v in (a.λ, a.nᵣ, a.Kcab, a.Kcar, a.Kant, a.Kb, a.Kw, a.Km, a.Kp, a.Kcbc)
        h = hash(v, h)
    end
    return h
end
