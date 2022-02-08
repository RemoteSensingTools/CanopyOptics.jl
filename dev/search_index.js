var documenterSearchIndex = {"docs":
[{"location":"#CanopyOptics.jl","page":"Home","title":"CanopyOptics.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package to compute canopy scattering properties","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Use leaf angle distributions to compute bi-lambertian area scattering matrices\nCompute specular reflection\nCompute leaf reflectance and transmittance based on Prospect-PRO","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The latest release of CanopyOptics can be installed from the Julia REPL prompt with","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]add https://github.com/RemoteSensingTools/CanopyOptics.jl","category":"page"},{"location":"#Function-Documentation","page":"Home","title":"Function Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Author = \"Christian Frankenberg\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CanopyOptics]","category":"page"},{"location":"#CanopyOptics.AbstractCanopyScatteringType","page":"Home","title":"CanopyOptics.AbstractCanopyScatteringType","text":"Abstract Type for canopy scattering\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.AbstractLeafDistribution","page":"Home","title":"CanopyOptics.AbstractLeafDistribution","text":"Abstract Type for leaf distributions\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.AbstractLeafProperties","page":"Home","title":"CanopyOptics.AbstractLeafProperties","text":"Abstract Type for leaf composition\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.AbstractProspectProperties","page":"Home","title":"CanopyOptics.AbstractProspectProperties","text":"Abstract type for Prospect model versions\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.BiLambertianCanopyScattering","page":"Home","title":"CanopyOptics.BiLambertianCanopyScattering","text":"Model for bi-lambertian canopy leaf scattering\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.LeafDistribution","page":"Home","title":"CanopyOptics.LeafDistribution","text":"struct LeafDistribution{FT<:AbstractFloat}\n\nA struct that defines the leaf angular distribution in radians (from 0->π/2; here scaled to [0,1])\n\nFields\n\nLD\nJulia Univariate Distribution from Distributions.js\nscaling\nScaling factor to normalize distribution (here mostly 2/π as Beta distribution is from [0,1])\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.LeafOpticalProperties","page":"Home","title":"CanopyOptics.LeafOpticalProperties","text":"struct PigmentOpticalProperties{FT}\n\nA struct which stores important absorption cross sections of pigments, water, etc\n\nFields\n\nλ\nWavelength [length]\nnᵣ\nRefractive index of leaf material\nKcab\nspecific absorption coefficient of chlorophyll (a+b) [cm² μg⁻¹]\nKcar\nspecific absorption coefficient of carotenoids [cm² μg⁻¹]\nKant\nspecific absorption coefficient of Anthocyanins [cm² nmol⁻¹]\nKb\nspecific absorption coefficient of brown pigments (arbitrary units)\nKw\nspecific absorption coefficient of water [cm⁻¹]\nKm\nspecific absorption coefficient of dry matter [cm² g⁻¹]\nKp\nspecific absorption coefficient of proteins [cm² g⁻¹]\nKcbc\nspecific absorption coefficient of carbon based constituents [cm² g⁻¹]\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.LeafProspectProProperties","page":"Home","title":"CanopyOptics.LeafProspectProProperties","text":"struct LeafProperties{FT}\n\nA struct which stores important variables of leaf chemistry and structure\n\nFields\n\nN\nLeaf structure parameter [0-3]\nCcab\nChlorophyll a+b content [µg cm⁻²]\nCcar\nCarotenoid content [µg cm⁻²]\nCanth\nAnthocynanin content [nmol cm⁻²]\nCbrown\nBrown pigments content in arbitrary units\nCw\nEquivalent water thickness [cm], typical about 0.002-0.015\nCm\nDry matter content (dry leaf mass per unit area) [g cm⁻²], typical about 0.003-0.016\nCprot\nprotein content [g/cm]\nCcbc\nCarbone-based constituents content in [g/cm⁻²] (cellulose, lignin, sugars...)\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.SpecularCanopyScattering","page":"Home","title":"CanopyOptics.SpecularCanopyScattering","text":"Model for specular canopy leaf scattering\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.dirVector","page":"Home","title":"CanopyOptics.dirVector","text":"Struct for spherical coordinate directions in θ (elevation angle) and ϕ (azimuth angle)\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.dirVector_μ","page":"Home","title":"CanopyOptics.dirVector_μ","text":"Struct for spherical coordinate directions in θ (elevation angle) and ϕ (azimuth angle)\n\n\n\n\n\n","category":"type"},{"location":"#CanopyOptics.A-Union{Tuple{FT}, Tuple{FT, FT}} where FT<:Real","page":"Home","title":"CanopyOptics.A","text":"Integrated projection of leaf area for a single leaf inclination of θₗ, assumes azimuthally uniform distribution\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.Fᵣ-Union{Tuple{FT}, Tuple{FT, FT}} where FT","page":"Home","title":"CanopyOptics.Fᵣ","text":"Fresnel reflection, needs pol_type soon too!\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.G-Union{Tuple{FT}, Tuple{Array{FT}, CanopyOptics.AbstractLeafDistribution}} where FT","page":"Home","title":"CanopyOptics.G","text":"G(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20)\n\nReturns the integrated projection of leaf area in the direction of μ, assumes azimuthally uniform distribution and a LD distribution for leaf polar angle θ.  This function is often referred to as the function O(B) (Goudriaan 1977) or G(Ζ) (Ross 1975,1981), see Bonan modeling book, eqs. 14.21-14.26. \n\nArguments\n\nμ an array of cos(θ) (directions [0,1]) \nLD an AbstractLeafDistribution type struct, includes a leaf distribution function\nnLeg an optional parameter for the number of legendre polynomials to integrate over the leaf distribution (default=20)\n\nExamples\n\njulia> μ,w = CanopyOptics.gauleg(10,0.0,1.0);       # Create 10 quadrature points in μ      \njulia> LD  = CanopyOptics.spherical_leaves()        # Create a default spherical leaf distribution\njulia> G   = CanopyOptics.G(μ, LD)                  # Compute G(μ)\n10-element Vector{Float64}:\n 0.5002522783000879\n 0.5002715115149204\n 0.5003537989277846\n 0.5004432798701134\n 0.5005134448870893\n 0.5003026448466977\n 0.4999186257540982\n 0.4994511190721635\n 0.49907252201082375\n 0.49936166823681594\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.H-Union{Tuple{FT}, Tuple{FT, FT}} where FT","page":"Home","title":"CanopyOptics.H","text":"Eq 46 in Shultis and Myneni\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.K-Union{Tuple{FT}, Tuple{FT, FT}} where FT","page":"Home","title":"CanopyOptics.K","text":"The reduction factor proposed by Nilson and Kuusk, κ ≈ 0.1-0.3\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.calctav-Union{Tuple{FT}, Tuple{FT, FT}} where FT<:AbstractFloat","page":"Home","title":"CanopyOptics.calctav","text":"calctav(α::FT, nr::FT) where {FT<:AbstractFloat}\n\nComputes transmission of isotropic radiation across a dielectric surface  (Stern F., 1964; Allen W.A.,Appl. Opt., 3(1):111-113 1973)).  From calctav.m in PROSPECT-D\n\nArguments\n\nα angle of incidence [degrees]\nnr Index of refraction\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.compute_Z_matrices-Union{Tuple{FT}, Tuple{CanopyOptics.BiLambertianCanopyScattering, Vector{FT}, CanopyOptics.AbstractLeafDistribution, Int64}} where FT","page":"Home","title":"CanopyOptics.compute_Z_matrices","text":"compute_Z_matrices(mod::BiLambertianCanopyScattering, μ::Array{FT,1}, LD::AbstractLeafDistribution, m::Int)\n\nComputes the single scattering Z matrices (𝐙⁺⁺ for same incoming and outgoing sign of μ, 𝐙⁻⁺ for a change in direction). Internally computes the azimuthally-averaged area scattering transfer function following Shultis and Myneni (https://doi.org/10.1016/0022-4073(88)90079-9), Eq 43::\n\nΓ(μ - μ) = int_0^1 dμ_L g_L(μ_L)t_L Ψ(μ μ μ_L) + r_L Ψ(μ μ μ_L)\n\nassuming an azimuthally uniform leaf angle distribution. Normalized Γ as 𝐙 = 4Γ/(ϖ⋅G(μ)). Returns 𝐙⁺⁺, 𝐙⁻⁺ \n\nArguments\n\nmod : A bilambertian canopy scattering model BiLambertianCanopyScattering, uses R,T,nQuad from that model.\nμꜛ::Array{FT,1}: Quadrature points ∈ [0,1]\nLD a AbstractLeafDistribution struct that describes the leaf angular distribution function.\nm: Fourier moment (for azimuthally uniform leave distributions such as here, only m=0 returns non-zero matrices)\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.compute_lambertian_Γ-Union{Tuple{FT}, Tuple{Vector{FT}, Vector{FT}, Any, Any, CanopyOptics.AbstractLeafDistribution}} where FT","page":"Home","title":"CanopyOptics.compute_lambertian_Γ","text":"compute_lambertian_Γ(μ::Array{FT,1},μꜛ::Array{FT,1}, r,t, LD::AbstractLeafDistribution; nLeg = 20)\n\nComputes the azimuthally-averaged area scattering transfer function following Shultis and Myneni (https://doi.org/10.1016/0022-4073(88)90079-9), Eq 43:\n\nΓ(μ - μ) = int_0^1 dμ_L g_L(μ_L)t_L Ψ(μ μ μ_L) + r_L Ψ(μ μ μ_L)\n\nassuming an azimuthally uniform leaf angle distribution.\n\nArguments\n\nμ::Array{FT,1} : Quadrature points incoming direction (cos(θ))\nμꜛ::Array{FT,1}: Quadrature points outgoing direction (cos(θ))\nr : Leaf lambertian reflectance\nt : Leaf lambertian transmittance\nLD a AbstractLeafDistribution struct that describes the leaf angular distribution function.\nnLeg = 20: number of quadrature points used for integration over all leaf angles (default is 20).\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.compute_Ψ-Union{Tuple{FT}, Tuple{Vector{FT}, FT}} where FT","page":"Home","title":"CanopyOptics.compute_Ψ","text":"Eq 45 in Shultis and Myneni, fixed grid in μ for both μ and μ' \n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.compute_Ψ-Union{Tuple{FT}, Tuple{Vector{FT}, Vector{FT}, FT}} where FT","page":"Home","title":"CanopyOptics.compute_Ψ","text":"Eq 45 in Shultis and Myneni, fixed grid in μ for both μ and μ' \n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.createLeafOpticalStruct-Tuple{Any}","page":"Home","title":"CanopyOptics.createLeafOpticalStruct","text":"createLeafOpticalStruct(λ_bnds)\n\nLoads in the PROSPECT-PRO database of pigments (and other) absorption cross section in leaves, returns a [`LeafOpticalProperties`](@ref) type struct with spectral units attached.\n\nArguments\n\n- `λ_bnds` an array (with or without a spectral grid unit) that defines the upper and lower limits over which to average the absorption cross sections\n\nExamples\n\njulia> using Unitful                                                               # Need to include Unitful package \njulia> opti = createLeafOpticalStruct((400.0:5:2400)*u\"nm\");                       # in nm\njulia> opti = createLeafOpticalStruct((0.4:0.1:2.4)*u\"μm\");                        # in μm\njulia> opti = CanopyOptics.createLeafOpticalStruct((10000.0:100:25000.0)u\"1/cm\");  # in wavenumber (cm⁻¹)\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.createLeafOpticalStruct-Tuple{}","page":"Home","title":"CanopyOptics.createLeafOpticalStruct","text":"createLeafOpticalStruct()\nAs in createLeafOpticalStruct(λ_bnds) but reads in the in Prospect-PRO database at original resolution (400-2500nm in 1nm steps)\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.erectophile_leaves","page":"Home","title":"CanopyOptics.erectophile_leaves","text":"erectophile_leaves(FT=Float64)\n\nCreates a LeafDistribution following a erectophile (mostly vertical) distribution using a Beta distribution\n\nStandard values for θₗ (63.24) and s (18.5) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1\n\n\n\n\n\n","category":"function"},{"location":"#CanopyOptics.plagiophile_leaves","page":"Home","title":"CanopyOptics.plagiophile_leaves","text":"plagiophile_leaves(FT=Float64)\n\nCreates a LeafDistribution following a plagiophile (mostly 45degrees) distribution using a Beta distribution\n\nStandard values for θₗ (45.0) and s (16.27) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1\n\n\n\n\n\n","category":"function"},{"location":"#CanopyOptics.planophile_leaves","page":"Home","title":"CanopyOptics.planophile_leaves","text":"planophile_leaves(FT=Float64)\n\nCreates a LeafDistribution following a planophile (mostly horizontal) distribution using a Beta distribution\n\nStandard values for θₗ (26.76) and s (18.51) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1\n\n\n\n\n\n","category":"function"},{"location":"#CanopyOptics.prospect-Union{Tuple{FT}, Tuple{LeafProspectProProperties{FT}, Any}} where FT<:AbstractFloat","page":"Home","title":"CanopyOptics.prospect","text":"prospect(leaf::LeafProspectProProperties{FT},\n                optis) where {FT<:AbstractFloat}\n\nComputes leaf optical properties (reflectance and transittance) based on PROSPECT-PRO\n\nArguments\n\nleaf  LeafProspectProProperties type struct which provides leaf composition\noptis LeafOpticalProperties type struct, which provides absorption cross sections and spectral grid\n\nExamples\n\njulia-repl julia> opti = createLeafOpticalStruct((400.0:5:2400)*u\"nm\"); julia> leaf = LeafProspectProProperties{Float64}(Ccab=30.0); julia> T,R = prospect(leaf,opti);`\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.scattering_model_lambertian-Tuple{dirVector, dirVector, dirVector, Any, Any}","page":"Home","title":"CanopyOptics.scattering_model_lambertian","text":"Eq 37 in Shultis and Myneni \n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.spherical_leaves","page":"Home","title":"CanopyOptics.spherical_leaves","text":"spherical_leaves(FT=Float64)\n\nCreates a LeafDistribution following a spherical (distributed as facets on a sphere) distribution using a Beta distribution\n\nStandard values for θₗ (57.3) and s (21.55) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1\n\n\n\n\n\n","category":"function"},{"location":"#CanopyOptics.uniform_leaves","page":"Home","title":"CanopyOptics.uniform_leaves","text":"uniform_leaves(FT=Float64)\n\nCreates a LeafDistribution following a uniform distribution (all angles equally likely)\n\n\n\n\n\n","category":"function"},{"location":"#CanopyOptics.unit2nm-Tuple{Any}","page":"Home","title":"CanopyOptics.unit2nm","text":"Convert an input array with Spectral units (or none) to a grid in nm\n\n\n\n\n\n","category":"method"},{"location":"#CanopyOptics.βparameters-Tuple{Any, Any}","page":"Home","title":"CanopyOptics.βparameters","text":"Get Beta Distribution parameters using standard tabulated values θₗ,s, Eq. 2.9 and 2.10 in Bonan\n\n\n\n\n\n","category":"method"}]
}