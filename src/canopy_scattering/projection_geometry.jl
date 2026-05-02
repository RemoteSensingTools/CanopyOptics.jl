
"""
    A(θ::FT, θₗ::FT) where FT<:Real

Returns the azimuthally-integrated projected leaf area function for beam direction `θ`
and leaf inclination `θₗ` (both in radians), assuming azimuthal uniformity.
Follows Bonan, *Modeling Earth's Climate*, Eq. 14.24:

``A(θ, θₗ) = \\begin{cases}
  \\cos θ\\,\\cos θₗ & θₗ \\le π/2 - θ \\\\
  \\frac{2}{π}\\!\\left[\\sqrt{\\sin^2 θₗ - \\cos^2 θ} + \\cos θ\\,\\cos θₗ\\,
    \\arcsin\\!\\left(\\frac{\\cos θ\\,\\cos θₗ}{\\sin θ\\,\\sin θₗ}\\right)\\right]
  & \\text{otherwise}
\\end{cases}``

Used by [`G`](@ref) to compute the Ross G-function.
See also [`Asm`](@ref) for the algebraically equivalent Shultis & Myneni form.
"""
function A(θ::FT, θₗ::FT) where FT<:Real # Suniti: why not use the expressions in Eq 35-36 from Schultis & Myneni (1987)?
    a = cos(θ) * cos(θₗ)
    #@show θ
    # Eq. 14.24 in Bonan et al.
    if θₗ ≤ (FT(π/2) - θ)
    #@show "<"
    #@show θₗ, a
        return a
    else
    #@show ">"
        b = sin(θ) * sin(θₗ)
        c = sqrt(max(zero(FT), sin(θₗ)^2 - cos(θ)^2))
     #   @show θₗ, FT(2/π)*(c + a * asin(a/b))
        return FT(2/π)*(c + a * asin(a/b))
    end
end

"""
    Asm(θ::FT, θₗ::FT) where FT<:Real

Algebraically equivalent form of [`A`](@ref) following Shultis & Myneni (1987)
Eqs. 35–36.  Used for cross-validation; production code uses [`A`](@ref).
"""
function Asm(θ::FT, θₗ::FT) where FT<:Real # Suniti: Eq 35-36 from Schultis & Myneni (1987)
    a = cos(θ) * cos(θₗ)
    # Eq. 14.24 in Bonan et al.
    if θₗ ≤ (FT(π/2) - θ)
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
- `LD` an `AbstractLeafDistribution` type struct, includes a leaf distribution function
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
function G(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg=40) where FT
    θₗ, w = gauleg(nLeg, FT(0), FT(π/2))
    Fᵢ = FT.(pdf.(LD.LD, FT(2) .* θₗ ./ FT(π))) .* FT(LD.scaling)
    Fᵢ = Fᵢ / (w' * Fᵢ)   # normalize leaf angle distribution
    θ  = acos.(abs.(μ))
    G  = (w .* Fᵢ)' * A.(θ', θₗ)
    return FT.(G')
end

"""
    G2(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg=40)

Reference implementation of [`G`](@ref) that integrates in leaf-normal
cosine `μ_L` instead of leaf inclination angle `θ_L`.

This path is retained for cross-checking the production `θ_L` quadrature in
[`G`](@ref). New code should normally call `G`.
"""
function G2(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg=40) where FT
    μl,w = gauleg(nLeg,FT(0),FT(1))
    θₗ = acos.(μl)
    Fᵢ = pdf.(LD.LD,2θₗ/π) * LD.scaling * π/2 
    Fᵢ = Fᵢ ./ (w'*Fᵢ)
    θ = acos.(abs.(μ))
    G = (w .* Fᵢ)' * A.(θ',θₗ)
    return G'
end

"""
    bfG(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20)

Brute-force two-angle reference for the Ross projection factor `G(μ)`.

This directly integrates `abs(Ω ⋅ Ω_L)` over leaf inclination and leaf azimuth.
It is useful for tests or debugging quadrature changes, but it is intentionally
slower than [`G`](@ref).
"""
function bfG(μ::Array{FT}, LD::AbstractLeafDistribution; nLeg=20) where FT
    nQuad = 100
    ϕ, w_azi = gauleg(nQuad,FT(0),FT(2π));
    # Reference angles to integrate over in both ϕ and μ
    
    μ_l, w = gauleg(180,0.0,1.0);
    Ω_l  = [dirVector_μ(a,b) for a in μ_l, b in ϕ];
    θₗ = acos.(μ_l)
    # Have to divide by sin(θ) again to get ∂θ/∂μ for integration (weights won't work)
    Fᵢ = pdf.(LD.LD,2θₗ/π)  * LD.scaling * π/2#./ abs.(sin.(θₗ))
    Fᵢ = Fᵢ ./ (Fᵢ' * w)
    #@show Fᵢ' * w
    res = similar(μ);
    
    for i in eachindex(μ)
        Ω = dirVector_μ(abs(μ[i]),0.0);
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
- `LD` an `AbstractLeafDistribution` struct that describes the leaf angular distribution function.
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
