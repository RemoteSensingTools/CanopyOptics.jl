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

"Abstract Type for canopy scattering"
abstract type AbstractCanopyScatteringType{FT<:Real} end

"""
    CanopyQuadrature(; n_leaf = 64, n_azimuth = 48)

Integration controls used when projecting canopy scattering models onto
Fourier Z matrices.

`n_leaf` controls the leaf-inclination quadrature used by the bi-Lambertian
kernel. `n_azimuth` controls the outgoing-azimuth quadrature used by the
specular kernel and legacy brute-force reference paths.
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

"Model for bi-lambertian canopy leaf scattering"
struct BiLambertianCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}
    "Lambertian Reflectance"
    R::FT
    "Lambertian Transmission"
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

"Model for specular canopy leaf scattering"
struct SpecularCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}
    "Refractive index"
    nᵣ::FT
    "Roughness parameter"
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

Additive canopy scattering model.  Components are evaluated independently and
their Z matrices are summed, so a leaf can carry both diffuse bi-Lambertian and
specular surface terms.
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

"Abstract Type for leaf distributions"
abstract type AbstractLeafDistribution{FT<:Real} end

"""
    struct LeafDistribution{FT<:Real}
A struct that defines the leaf angular distribution in radians (from 0->π/2; here scaled to [0,1])
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafDistribution{FT<:Real} <: AbstractLeafDistribution{FT}
    "Julia Univariate Distribution from Distributions.js"
    LD::UnivariateDistribution
    "Scaling factor to normalize distribution (here mostly 2/π as Beta distribution is from [0,1])"
    scaling::FT
end


####################################
"Abstract Type for leaf composition"
abstract type AbstractLeafProperties end

"""
    struct LeafProspectProProperties{FT} <: AbstractLeafProperties

Leaf biochemical and structural state for PROSPECT-D (bulk dry matter `Cm` mode)
and PROSPECT-PRO (`Cprot + Ccbc` dry-matter partition mode).

When constructed with keywords:
- Omitting `Cm` selects **PROSPECT-PRO mode** (`Cm = 0`, defaulting `Cprot = 0.001`, `Ccbc = 0.009`).
- Supplying `Cm` explicitly (including `Cm = 0` or `ForwardDiff.Dual`) selects **PROSPECT-D mode**,
  defaulting `Cprot = 0, Ccbc = 0` so bulk dry matter is never silently double-counted.
- Setting both `Cm > 0` and `Cprot > 0` (or `Ccbc > 0`) throws an `ArgumentError`.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafProspectProProperties{FT} <: AbstractLeafProperties
    ### Prospect related parameters
    "Leaf structure parameter N >= 1 (fractional plate count)"
    N::FT
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Ccab::FT
    "Carotenoid content `[µg cm⁻²]`"
    Ccar::FT
    "Anthocyanin content `[µg cm⁻²]` (mass basis matching Féret et al. 2017/2021 `OPTI_2021`)"
    Canth::FT
    "Brown pigments content in arbitrary units"
    Cbrown::FT
    "Equivalent water thickness `[cm]`, typical about 0.002–0.015"
    Cw::FT
    "Bulk dry matter content (LMA) `[g cm⁻²]`, typical about 0.003–0.016 (0 when `Cprot + Ccbc` are used)"
    Cm::FT
    "Protein content `[g cm⁻²]`"
    Cprot::FT
    "Carbon-based constituents content `[g cm⁻²]` (cellulose, lignin, sugars, ...)"
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

####################################
"Abstract type for Prospect model versions"
abstract type AbstractProspectProperties end

"""
    struct LeafOpticalProperties{FT,FT2} <: AbstractProspectProperties
A struct which stores important absorption cross sections of pigments, water, etc
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafOpticalProperties{FT,FT2} <: AbstractProspectProperties
    "Wavelength `[length]`"
    λ::FT2 #typeof(([1.0])u"nm")
    ### Prospect-PRO related parameters
    "Refractive index of leaf material "
    nᵣ::Array{FT, 1}
    "specific absorption coefficient of chlorophyll (a+b) `[cm² μg⁻¹]`" 
    Kcab::Array{FT, 1}
    "specific absorption coefficient of carotenoids `[cm² μg⁻¹]`"
    Kcar::Array{FT, 1}       
    "specific absorption coefficient of anthocyanins `[cm² μg⁻¹]`"       
    Kant::Array{FT, 1} 
    "specific absorption coefficient of brown pigments (arbitrary units)" 
    Kb::Array{FT, 1} 
    "specific absorption coefficient of water `[cm⁻¹]`"    
    Kw::Array{FT, 1}                  
    "specific absorption coefficient of dry matter `[cm² g⁻¹]`" 
    Km::Array{FT, 1}
    "specific absorption coefficient of proteins `[cm² g⁻¹]`"  
    Kp::Array{FT, 1}
    "specific absorption coefficient of carbon based constituents `[cm² g⁻¹]`" 
    Kcbc::Array{FT, 1} 
end
