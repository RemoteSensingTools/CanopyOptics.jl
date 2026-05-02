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
abstract type AbstractLeafDistribution{FT<:AbstractFloat} end

"""
    struct LeafDistribution{FT<:AbstractFloat}
A struct that defines the leaf angular distribution in radians (from 0->π/2; here scaled to [0,1])
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafDistribution{FT<:AbstractFloat} <: AbstractLeafDistribution{FT}
    "Julia Univariate Distribution from Distributions.js"
    LD::UnivariateDistribution
    "Scaling factor to normalize distribution (here mostly 2/π as Beta distribution is from [0,1])"
    scaling::FT
end


####################################
"Abstract Type for leaf composition"
abstract type AbstractLeafProperties end

"""
    struct LeafProperties{FT}
A struct which stores important variables of leaf chemistry and structure
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct LeafProspectProProperties{FT} <: AbstractLeafProperties
    ### Prospect related parameters
    "Leaf structure parameter [0-3]"
    N::FT       = FT(1.4  )
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Ccab::FT    = FT(40.0)#u"µg/cm^2"
    "Carotenoid content `[µg cm⁻²]`"
    Ccar::FT    = FT(10.0)#u"µg/cm^2"
    "Anthocynanin content `[nmol cm⁻²]`"
    Canth::FT   = FT(0.5)#u"nmol/cm^2"
    "Brown pigments content in arbitrary units"
    Cbrown::FT = FT(0.0)
    "Equivalent water thickness `[cm]`, typical about 0.002-0.015"
    Cw::FT    = FT(0.012)#u"cm"
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`, typical about 0.003-0.016"
    Cm::FT   = FT(0.0)#u"g/cm^2"
    "protein content `[g/cm]`"
    Cprot::FT  = FT(0.001)#u"g/cm^2"
    "Carbone-based constituents content in `[g/cm⁻²]` (cellulose, lignin, sugars...)"
    Ccbc::FT   =  FT(0.009)#u"g/cm^2"
end

####################################
"Abstract type for Prospect model versions"
abstract type AbstractProspectProperties end

"""
    struct PigmentOpticalProperties{FT}
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
    "specific absorption coefficient of Anthocyanins `[cm² nmol⁻¹]`"       
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
