"Abstract Type for canopy scattering"
abstract type AbstractCanopyScatteringType{FT<:AbstractFloat} end

"Model for bi-lambertian canopy leaf scattering"
Base.@kwdef struct BiLambertianCanopyScattering{FT<:AbstractFloat} <: AbstractCanopyScatteringType{FT}
    "Lambertian Reflectance"
    R::FT = FT(0.3)
    "Lambertian Transmission"
    T::FT  = FT(0.1)
    "Number of quadrature points in inclination angle"
    nQuad::Int = 10
end

"Model for specular canopy leaf scattering"
Base.@kwdef struct SpecularCanopyScattering{FT<:AbstractFloat} <: AbstractCanopyScatteringType{FT}
    "Refractive index"
    nᵣ::FT = FT(1.5)
    "Roughness parameter"
    κ::FT  = FT(0.5)
    "Number of quadrature points in azimuth"
    nQuad::Int = 20
end

"Abstract Type for leaf distributions"
abstract type AbstractLeafDistribution end

"""
    struct LeafDistribution{FT<:AbstractFloat}
A struct that defines the leaf angular distribution in radians (from 0->π/2; here scaled to [0,1])
# Fields
$(DocStringExtensions.FIELDS)
"""
struct LeafDistribution{FT<:AbstractFloat} <: AbstractLeafDistribution
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