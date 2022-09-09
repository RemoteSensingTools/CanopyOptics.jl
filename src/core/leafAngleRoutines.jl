
"Integrated projection of leaf area for a single leaf inclination of θₗ, assumes azimuthally uniform distribution"
function A(θ::FT, θₗ::FT) where FT<:Real # Suniti: why not use the expressions in Eq 35-36 from Schultis & Myneni (1987)?
    a = cos(θ) * cos(θₗ)
    # Eq. 14.24 in Bonan et al.
    if θₗ ≤ FT(π/2) - θ
        return a
    else
        b = sin(θ) * sin(θₗ)
        c = sqrt(sin(θₗ)^2  - cos(θ)^2)
        return FT(2/π)*(c + a * asin(a/b))
    end
end

"Integrated projection of leaf area for a single leaf inclination of θₗ, assumes azimuthally uniform distribution"
function Asm(θ::FT, θₗ::FT) where FT<:Real # Suniti: Eq 35-36 from Schultis & Myneni (1987)
    a = cos(θ) * cos(θₗ)
    # Eq. 14.24 in Bonan et al.
    if θₗ ≤ FT(π/2) - θ
        return a
    else
        b = sin(θ) * sin(θₗ)
        c = sqrt(FT(1)-(a/b)^2)
        return FT(2/π)*(a * (acos(-a/b)-π/2) + b * c) # (1/π) * (cosθ.cosθₗ.cos⁻¹(cotθ.cotθₗ) + sinθ.sinθₗ.√(1-cot²θ.cot²θₗ))
    end
end



"""
    $(FUNCTIONNAME)(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20)

Returns the integrated projection of leaf area in the direction of μ, assumes azimuthally uniform distribution and a LD distribution for leaf polar angle θ. 
This function is often referred to as the function O(B) (Goudriaan 1977) or G(Ζ) (Ross 1975,1981), see Bonan modeling book, eqs. 14.21-14.26. 

# Arguments
- `μ` an array of cos(θ) (directions [0,1]) 
- `LD` an [`AbstractLeafDistribution`](@ref) type struct, includes a leaf distribution function
- `nLeg` an optional parameter for the number of legendre polynomials to integrate over the leaf distribution (default=20)

# Examples
```julia-repl
julia> μ,w = CanopyOptics.gauleg(10,0.0,1.0);       # Create 10 quadrature points in μ      
julia> LD  = CanopyOptics.spherical_leaves()        # Create a default spherical leaf distribution
julia> G   = CanopyOptics.G(μ, LD)                  # Compute G(μ)
10-element Vector{Float64}:
 0.5002522783000879
 0.5002715115149204
 0.5003537989277846
 0.5004432798701134
 0.5005134448870893
 0.5003026448466977
 0.4999186257540982
 0.4994511190721635
 0.49907252201082375
 0.49936166823681594
```
"""
function G(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20) where FT
    θₗ,w = gauleg(nLeg,FT(0),FT(π/2))
    Fᵢ = pdf.(LD.LD,2θₗ/π) * LD.scaling
    @show Fᵢ' * w
    res = similar(μ);
    θ = acos.(μ)
    for i in eachindex(μ)
        res[i] =  sum(w .* Fᵢ .* A.(θ[i],θₗ))
    end
    return res
end

"Brute Force G calculation (for testing"
function bfG(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20) where FT
    nQuad = 580
    ϕ, w_azi = gauleg(nQuad,FT(0),FT(2π));
    # Reference angles to integrate over in both ϕ and μ
    
    μ_l, w = gauleg(180,0.0,1.0);
    Ω_l  = [dirVector_μ(a,b) for a in μ_l, b in ϕ];
    θₗ = acos.(μ_l)
    # Have to divide by sin(θ) again to get ∂θ/∂μ for integration (weights won't work)
    Fᵢ = pdf.(LD.LD,2θₗ/π)  * LD.scaling ./ abs.(sin.(θₗ))
    #@show Fᵢ' * w 
    #Fᵢ = Fᵢ ./ (Fᵢ' * w)
    #@show Fᵢ' * w
    res = similar(μ);
    
    for i in eachindex(μ)
        Ω = dirVector_μ(μ[i],0.0);
        #res[i] =  sum(w .* Fᵢ .* A.(θ[i],θₗ))
        # Double integration here:
        res[i] =  ((Fᵢ .* abs.(dot.((Ω,),Ω_l)))' * w)' * w_azi /(2π)
    end
    return res
end



"""
    $(FUNCTIONNAME)(μ::Array{FT,1},μꜛ::Array{FT,1}, r,t, LD::AbstractLeafDistribution; nLeg = 20)

Computes the azimuthally-averaged area scattering transfer function following Shultis and Myneni (https://doi.org/10.1016/0022-4073(88)90079-9), Eq 43:

``Γ(μ' -> μ) = \\int_0^1 dμ_L g_L(μ_L)[t_L Ψ⁺(μ, μ', μ_L) + r_L Ψ⁻(μ, μ', μ_L)]``

assuming an azimuthally uniform leaf angle distribution.
# Arguments
- `μ::Array{FT,1}` : Quadrature points incoming direction (cos(θ))
- `μꜛ::Array{FT,1}`: Quadrature points outgoing direction (cos(θ))
- `r` : Leaf lambertian reflectance
- `t` : Leaf lambertian transmittance
- `LD` a [`AbstractLeafDistribution`](@ref) struct that describes the leaf angular distribution function.
- `nLeg = 20`: number of quadrature points used for integration over all leaf angles (default is 20).
"""
function compute_lambertian_Γ(μ::Array{FT,1},μꜛ::Array{FT,1}, r,t, LD::AbstractLeafDistribution; nLeg = 20) where FT
    Γ = zeros(length(μ), length(μ))
    θₗ,w = gauleg(nLeg,FT(0),FT(π/2))
    for i in eachindex(θₗ)
        Ψ⁺, Ψ⁻ = compute_Ψ(μ,μꜛ, cos(θₗ[i]));
        Γ += pdf.(LD.LD,2θₗ[i]/π) * LD.scaling * w[i] * (t * Ψ⁺ + r * Ψ⁻)
    end
    return Γ
end

"""
    $(FUNCTIONNAME)(mod::BiLambertianCanopyScattering, μ::Array{FT,1}, LD::AbstractLeafDistribution, m::Int)

Computes the single scattering Z matrices (𝐙⁺⁺ for same incoming and outgoing sign of μ, 𝐙⁻⁺ for a change in direction). Internally computes the azimuthally-averaged area scattering transfer function following Shultis and Myneni (https://doi.org/10.1016/0022-4073(88)90079-9), Eq 43::

``Γ(μ' -> μ) = \\int_0^1 dμ_L g_L(μ_L)[t_L Ψ⁺(μ, μ', μ_L) + r_L Ψ⁻(μ, μ', μ_L)]``

assuming an azimuthally uniform leaf angle distribution. Normalized Γ as 𝐙 = 4Γ/(ϖ⋅G(μ)).
Returns 𝐙⁺⁺, 𝐙⁻⁺ 

# Arguments
- `mod` : A bilambertian canopy scattering model [`BiLambertianCanopyScattering`](@ref), uses R,T,nQuad from that model.
- `μ::Array{FT,1}`: Quadrature points ∈ [0,1]
- `LD` a [`AbstractLeafDistribution`](@ref) struct that describes the leaf angular distribution function.
- `m`: Fourier moment (for azimuthally uniform leave distributions such as here, only m=0 returns non-zero matrices)
"""
function compute_Z_matrices(mod::BiLambertianCanopyScattering, μ::Array{FT,1}, LD::AbstractLeafDistribution, m::Int) where FT
    (;R,T,nQuad) = mod
    # Transmission (same direction)
    𝐙⁺⁺ = zeros(length(μ), length(μ))
    # Reflection (change direction)
    𝐙⁻⁺ = zeros(length(μ), length(μ))
    
    # skip everything beyond m=0
    if m>0  
        return 𝐙⁺⁺, 𝐙⁻⁺
    end
    # Ross kernel
    G = CanopyOptics.G(μ, LD)
    # Single Scattering Albedo (should make this a vector too)
    ϖ = R+T

    θₗ,w = gauleg(nQuad,FT(0),FT(π/2));
    for i in eachindex(θₗ)
        Ψ⁺, Ψ⁻ = compute_Ψ(μ,μ, cos(θₗ[i]));
        𝐙⁺⁺ += pdf.(LD.LD,2θₗ[i]/π) * LD.scaling * w[i] * (T * Ψ⁺ + R * Ψ⁻) 
        Ψ⁺, Ψ⁻ = compute_Ψ(μ,-μ, cos(θₗ[i]));
        𝐙⁻⁺ += pdf.(LD.LD,2θₗ[i]/π) * LD.scaling * w[i] * (T * Ψ⁺ + R * Ψ⁻) 
    end
    return 4𝐙⁺⁺ ./(G*ϖ), 4𝐙⁻⁺ ./(G*ϖ)
end

# Page 20, top of Knyazikhin and Marshak
# Example 
# ϕ = range(0.0, 2π,  length=200)
# θ = range(0.0, π/2, length=150)
# dirs = [dirVector(a,b) for a in θ, b in ϕ];
# R = CanopyOptics.compute_specular_reflection.([dirs[10,1]],dirs, [1.5], [0.3], [LD])
function compute_reflection(mod::SpecularCanopyScattering, Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}, LD) where FT
    (;nᵣ,κ) = mod
    Ωstar = getSpecularΩ(Ωⁱⁿ, Ωᵒᵘᵗ)
    #θstar = min(abs(Ωstar.θ), (π-abs(Ωstar.θ))) # min(abs(Ωstar.θ), abs(π+Ωstar.θ))
    θstar = Ωstar.θ;
    #if Ωⁱⁿ.θ ≈ Ωᵒᵘᵗ.θ && Ωⁱⁿ.ϕ ≈ Ωᵒᵘᵗ.ϕ
    #    θstar = Ωⁱⁿ.θ
    #end
    # Still needs to be implemented!
    # incident angle on leaf surface (half of in and out angle):
    sa = Ωⁱⁿ ⋅ Ωᵒᵘᵗ 
    sa > 1 ? sa = FT(1) : nothing
    αstar = acos(abs(sa))/2
    #@show Ωstar.ϕ, Ωstar.θ
    #a = (Ωⁱⁿ ⋅ Ωstar) * (Ωᵒᵘᵗ ⋅ Ωstar)
    return FT(1/8) * pdf(LD.LD,2θstar/π) * LD.scaling * K(κ, αstar) * Fᵣ(nᵣ,αstar)
    
end

function compute_reflection(mod::SpecularCanopyScattering,Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}, LD) where FT
    (;nᵣ,κ) = mod
    Ωstar, αstar = getSpecularΩ(Ωⁱⁿ, Ωᵒᵘᵗ)
    # Can change this later as well do have the pdf in μ, not theta!
    θstar = acos(abs(Ωstar.μ));
    # Eq. 2.39 in "Discrete Ordinates Method for Photon Transport in Leaf Canopies", page 59
    return FT(1/8) * pdf(LD.LD,2θstar/π) * LD.scaling * K(κ, αstar) * Fᵣ(nᵣ,αstar)
end

function compute_Z_matrices(mod::SpecularCanopyScattering, μ::Array{FT,1}, LD::AbstractLeafDistribution, m::Int) where FT
    (;nᵣ, κ, nQuad) = mod
    # Transmission (same direction)
    𝐙⁺⁺ = zeros(length(μ), length(μ))
    # Reflection (change direction)
    𝐙⁻⁺ = zeros(length(μ), length(μ))
    
    # Quadrature points in the azimuth:
    ϕ, w_azi = gauleg(nQuad,FT(0),FT(2π));
    # Fourier weights (cosine decomposition)
    f_weights = cos.(m*ϕ)
    for i in eachindex(μ)
        # Incoming beam at ϕ = 0
        Ωⁱⁿ = dirVector_μ(μ[i], FT(0));
        # Create outgoing vectors in θ and ϕ
        dirOutꜛ = [dirVector_μ(a,b) for a in μ, b in ϕ];
        dirOutꜜ = [dirVector_μ(a,b) for a in -μ, b in ϕ];
        # Compute over μ and μ_azi:
        Zup   = compute_reflection.((mod,),(Ωⁱⁿ,),dirOutꜛ, (LD,));
        Zdown = compute_reflection.((mod,),(Ωⁱⁿ,),dirOutꜜ, (LD,));
        # integrate over the azimuth:
        𝐙⁻⁺[i,:] = Zup   * (w_azi .* f_weights)
        𝐙⁺⁺[i,:] = Zdown * (w_azi .* f_weights)
    end
    return 𝐙⁺⁺, 𝐙⁻⁺
end

"The reduction factor proposed by Nilson and Kuusk, κ ≈ 0.1-0.3, returns exp(-κ * tan(abs(α))"
function K(κ::FT, α::FT) where FT 
    exp(-κ * tan(abs(α)));
end
