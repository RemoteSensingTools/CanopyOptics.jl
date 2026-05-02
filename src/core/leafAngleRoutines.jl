
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
        c = sqrt(sin(θₗ)^2  - cos(θ)^2)
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
    θ  = acos.(μ)
    G  = (w .* Fᵢ)' * A.(θ', θₗ)
    return FT.(G')
end

function G2(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg=40) where FT
    μl,w = gauleg(nLeg,FT(0),FT(1))
    θₗ = acos.(μl)
    Fᵢ = pdf.(LD.LD,2θₗ/π) * LD.scaling * π/2 
    Fᵢ = Fᵢ ./ (w'*Fᵢ)
    θ = acos.(μ)
    G = (w .* Fᵢ)' * A.(θ',θₗ)
    return G'
end

"Brute Force G calculation (for testing"
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

"""
    $(FUNCTIONNAME)(mod::BiLambertianCanopyScattering, μ::Array{FT,1}, LD::AbstractLeafDistribution, m::Int)

Computes the single scattering Z matrices (𝐙⁺⁺ for same incoming and outgoing sign of μ, 𝐙⁻⁺ for a change in direction). Internally computes the azimuthally-averaged area scattering transfer function following Shultis and Myneni (https://doi.org/10.1016/0022-4073(88)90079-9), Eq 43::

``Γ(μ' -> μ) = \\int_0^1 dμ_L g_L(μ_L)[t_L Ψ⁺(μ, μ', μ_L) + r_L Ψ⁻(μ, μ', μ_L)]``

assuming an azimuthally uniform leaf angle distribution. Normalized Γ as 𝐙 = 4Γ/(ϖ⋅G(μ)).
Returns 𝐙⁺⁺, 𝐙⁻⁺ 

# Arguments
- `mod` : A bilambertian canopy scattering model [`BiLambertianCanopyScattering`](@ref), uses R,T,nQuad from that model.
- `μ::Array{FT,1}`: Quadrature points ∈ [0,1]
- `LD` an `AbstractLeafDistribution` struct that describes the leaf angular distribution function.
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
    return 4𝐙⁺⁺ ./(G'*ϖ), 4𝐙⁻⁺ ./(G'*ϖ)
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

"""
    compute_Γ(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT},
              LD::AbstractLeafDistribution) where FT

Computes the azimuthally-resolved area scattering transfer function
`Γ(Ωⁱⁿ → Ωᵒᵘᵗ)` for a specific direction pair by direct double integration over
leaf polar and azimuthal angles, following Shultis & Myneni (1988) Eqs. (38)–(39):

``Γ = R\\,Γ^- + T\\,Γ^+``

where `Γ⁻` integrates the negative (reflection) part and `Γ⁺` the positive
(transmission) part of `(Ωⁱⁿ ⋅ Ωᴸ)(Ωᵒᵘᵗ ⋅ Ωᴸ)` over all leaf orientations.

Used by [`compute_Z_matrices_aniso`](@ref) to build Fourier-decomposed Z matrices.
See [`compute_lambertian_Γ`](@ref) for the azimuthally-averaged (m=0 only) equivalent.
"""
function compute_Γ(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}, LD::AbstractLeafDistribution) where FT
    (;R,T, nQuad) = mod
    #nQuad = 80
    #μ_l, w = gauleg(nQuad,0.0,1.0);
    θₗ, w = gauleg(nQuad,FT(0),FT(π/2));
    μ_l = cos.(θₗ)
    # Quadrature points in the azimuth (this has to go over 2π):
    ϕ, w_azi = gauleg(nQuad+1,FT(0),FT(2π));
    
    Fᵢ = pdf.(LD.LD,2θₗ/π) * LD.scaling  
    
    # Create leaf angles in μ and ϕ
    Ω_l  = [dirVector_μ(a,b) for a in μ_l, b in ϕ];
    
    # Compute the angles (in, leaf angles, out, leaf)
    integrand = ((Ωⁱⁿ,) .⋅ Ω_l) .* ((Ωᵒᵘᵗ,) .⋅ Ω_l)

    # Integrate over positive and negatives separately
    iPos = (integrand+abs.(integrand))./2
    iNeg = (integrand-abs.(integrand))./2
    
    # Eq 39 in Shultis and Myneni, double integration here
    Γ⁻ = -1/2π * (Fᵢ .* w)' * (iNeg * w_azi)
    Γ⁺ =  1/2π * (Fᵢ .* w)' * (iPos * w_azi)
    # Eq 38 in Shultis and Myneni
    return R * Γ⁻ + T * Γ⁺ 
end

"""
    compute_reflection(mod::SpecularCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT},
                       LD::AbstractLeafDistribution) where FT

Returns the specular bidirectional scattering coefficient for direction pair
`(Ωⁱⁿ, Ωᵒᵘᵗ)`, following Knyazikhin & Marshak, *Discrete Ordinates Method for
Photon Transport in Leaf Canopies*, Eq. 2.39:

``f_s(Ω' \\to Ω) = \\frac{1}{8}\\,g_L(θ^*)\\,K(κ, α^*)\\,F_r(n_r, α^*)``

where:
- `θ*` = polar angle of the specular leaf normal (from `getSpecularΩ`)
- `α*` = incidence half-angle = `arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ) / 2`
- `K(κ, α*)` = Nilson–Kuusk roughness factor (from `K`)
- `Fᵣ(nᵣ, α*)` = unpolarized Fresnel reflectance (from `Fᵣ`)

For full Stokes-vector propagation (vSmartMOM.jl), use `fresnel_components`
to obtain `r_s, r_p` and construct the 4×4 Mueller reflection matrix directly.

Note: only reflection is currently modelled; specular transmission is not yet implemented.
"""
function compute_reflection(mod::SpecularCanopyScattering,Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}, LD::AbstractLeafDistribution) where FT
    (;nᵣ,κ) = mod
    Ωstar, αstar = getSpecularΩ(Ωⁱⁿ, Ωᵒᵘᵗ)
    # Can change this later as well do have the pdf in μ, not theta!
    θstar = acos(abs(Ωstar.μ));
    # Eq. 2.39 in "Discrete Ordinates Method for Photon Transport in Leaf Canopies", page 59
    return FT(1/8) * pdf(LD.LD,2θstar/π) * LD.scaling * K(κ, αstar) * Fᵣ(nᵣ,αstar)
end

"""
    compute_Z_matrices(mod::SpecularCanopyScattering, μ::Array{FT,1},
                       LD::AbstractLeafDistribution, m::Int) where FT

Computes the Fourier-`m` component of the single-scattering phase matrices
`(𝐙⁺⁺, 𝐙⁻⁺)` for a specular leaf surface by integrating
[`compute_reflection`](@ref) over the azimuthal quadrature grid.

- `𝐙⁺⁺[i,j]`: same-hemisphere scattering (μ>0 → μ>0, forward scatter)
- `𝐙⁻⁺[i,j]`: opposite-hemisphere scattering (μ>0 → μ<0, backscatter)

Azimuth integration uses `nQuad` Gauss-Legendre points over `[0, 2π]` from `mod`;
Fourier weights are `cos(m ϕ)`.
"""
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
        # dirOutꜛ (same hemisphere, μ>0) → forward scatter → 𝐙⁺⁺
        # dirOutꜜ (opposite hemisphere, μ<0) → back scatter  → 𝐙⁻⁺
        𝐙⁺⁺[i,:] = Zup   * (w_azi .* f_weights)
        𝐙⁻⁺[i,:] = Zdown * (w_azi .* f_weights)
    end
    return 𝐙⁺⁺, 𝐙⁻⁺
end

function compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering, μ::AbstractArray{FT,1}, LD::AbstractLeafDistribution, m::Int) where FT
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices_aniso_analytic(mod, μ, LD, m)
    return Z⁺⁺[:, :, m + 1], Z⁻⁺[:, :, m + 1]
end

@inline function _psi_same(Pᵢ::FT, Nᵢ::FT, Pₒ::FT, Nₒ::FT) where {FT}
    return FT(2) * (Pᵢ * Pₒ + Nᵢ * Nₒ)
end

@inline function _psi_opposite(Pᵢ::FT, Nᵢ::FT, Pₒ::FT, Nₒ::FT) where {FT}
    return FT(2) * (Pᵢ * Nₒ + Nᵢ * Pₒ)
end

function _normalise_Z!(Z⁺⁺::AbstractArray{FT,3},
                       Z⁻⁺::AbstractArray{FT,3},
                       G::AbstractVector,
                       ϖ::FT) where {FT}
    nμ = size(Z⁺⁺, 1)
    nm = size(Z⁺⁺, 3)
    for k in 1:nm
        ff = k == 1 ? FT(2) : FT(4)
        for j in 1:nμ
            scale = ff / (ϖ * FT(G[j]))
            @inbounds for i in 1:nμ
                Z⁺⁺[i, j, k] *= scale
                Z⁻⁺[i, j, k] *= scale
            end
        end
    end
    return Z⁺⁺, Z⁻⁺
end

function _leaf_inclination_quadrature(LD::AbstractLeafDistribution, nQuad::Int, ::Type{FT}) where {FT}
    θₗ, wθ = gauleg(nQuad, FT(0), FT(π / 2))
    Fₗ = FT.(pdf.(LD.LD, 2θₗ / FT(π))) .* FT(LD.scaling)
    return (θ = θₗ, w_measure = wθ .* Fₗ)
end

"""
    _clipped_projection_moments!(P, N, μ::FT, μ_L::FT, m_max::Int) -> (P, N)
    _clipped_projection_moments(μ::FT, μ_L::FT, m_max::Int) -> (P, N)

Closed-form cosine Fourier moments of the clipped leaf projection for
orders `m = 0:m_max`.

For a viewing direction with signed cosine `μ` and a leaf normal with
inclination cosine `μ_L`, write

```math
x(ψ) = a + b\\cos ψ, \\qquad
a = μ μ_L, \\qquad
b = \\sqrt{1 - μ^2}\\sqrt{1 - μ_L^2}.
```

This routine returns

```math
P_m = \\frac{1}{2π}\\int_0^{2π} \\max(x(ψ), 0)\\cos(mψ)\\,dψ,
```

and

```math
N_m = \\frac{1}{2π}\\int_0^{2π} \\max(-x(ψ), 0)\\cos(mψ)\\,dψ.
```

The `m = 0` value `P_0` is the Shultis and Myneni (1988) projected-area
function `H` (their Eq. 35/46 in this codebase).  For `|a| < b`, with
`ψ* = acos(-a / b)`, the antiderivative gives

```math
P_m =
\\frac{1}{π}\\left[
\\frac{a\\sin(mψ*)}{m} +
\\frac{b}{2}\\left(
\\frac{\\sin((m-1)ψ*)}{m-1} +
\\frac{\\sin((m+1)ψ*)}{m+1}
\\right)\\right],
```

with the removable limits

```math
P_0 = \\frac{aψ* + b\\sin ψ*}{π}, \\qquad
P_1 = \\frac{a\\sin ψ* + \\frac{b}{2}ψ* + \\frac{b}{4}\\sin(2ψ*)}{π}.
```

`N_m` is derived from the exact identity

```math
P_m - N_m = a\\,δ_{m0} + \\frac{b}{2}\\,δ_{m1},
```

which follows because `max(x,0) - max(-x,0) = x`.
"""
function _clipped_projection_moments!(P::AbstractVector{FT},
                                      N::AbstractVector{FT},
                                      μ::FT, μ_L::FT,
                                      m_max::Int) where {FT<:Real}
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))
    length(P) >= m_max + 1 || throw(DimensionMismatch("P must have length at least m_max + 1"))
    length(N) >= m_max + 1 || throw(DimensionMismatch("N must have length at least m_max + 1"))
    fill!(P, zero(FT))
    fill!(N, zero(FT))

    a = μ * μ_L
    b = sqrt(max(zero(FT), one(FT) - μ * μ)) *
        sqrt(max(zero(FT), one(FT) - μ_L * μ_L))

    if b == zero(FT)
        P[1] = max(a, zero(FT))
    elseif a >= b
        P[1] = a
        if m_max >= 1
            P[2] = b / FT(2)
        end
    elseif a <= -b
        # P remains zero; N is filled from P - N below.
    else
        πFT = FT(π)
        ψ = acos(clamp(-a / b, -one(FT), one(FT)))
        P[1] = (a * ψ + b * sin(ψ)) / πFT
        if m_max >= 1
            P[2] = (a * sin(ψ) + (b / FT(2)) * ψ +
                    (b / FT(4)) * sin(FT(2) * ψ)) / πFT
        end
        for m in 2:m_max
            mFT = FT(m)
            P[m + 1] = (
                a * sin(mFT * ψ) / mFT +
                (b / FT(2)) * (
                    sin(FT(m - 1) * ψ) / FT(m - 1) +
                    sin(FT(m + 1) * ψ) / FT(m + 1)
                )
            ) / πFT
        end
    end

    for m in 0:m_max
        diff = m == 0 ? a : (m == 1 ? b / FT(2) : zero(FT))
        N[m + 1] = P[m + 1] - diff
    end

    return P, N
end

function _clipped_projection_moments(μ::FT, μ_L::FT, m_max::Int) where {FT<:Real}
    P = zeros(FT, m_max + 1)
    N = zeros(FT, m_max + 1)
    return _clipped_projection_moments!(P, N, μ, μ_L, m_max)
end

"""
    compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                      μ::AbstractVector{FT},
                                      LD::AbstractLeafDistribution,
                                      m_max::Int) -> (Z⁺⁺, Z⁻⁺)

Compute all scalar Fourier moments `m = 0:m_max` of the bi-Lambertian
canopy phase matrices without azimuthal quadrature.

The returned arrays have shape `(length(μ), length(μ), m_max + 1)` and
use the vSmartMOM convention

```math
Z[i_{out}, j_{in}, m+1],
```

where `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is sign-change reflection.
The single-scattering albedo `ϖ = R + T` is divided out, so for
conservative leaves the `m = 0` column integral of `Z⁺⁺ + Z⁻⁺` targets
`≈ 2` under the same quadrature and `G(μ)` convention used by vSmartMOM.

# Origin trace

Shultis and Myneni (1988), Eq. 45 writes the azimuthally averaged
Lambertian canopy kernels as products of one-direction clipped projection
functions.  For each leaf inclination this implementation generalizes that
factorization from `m = 0` to arbitrary cosine Fourier order by using the
closed-form moments returned by `_clipped_projection_moments`.

For signed incoming/outgoing directions, let `(Pᵢ, Nᵢ)` and `(Pₒ, Nₒ)` be
the clipped projection moments.  The circular-convolution theorem gives the
per-leaf Fourier coefficients

```math
Ψ^+_m = P_{i,m}P_{o,m} + N_{i,m}N_{o,m},
\\qquad
Ψ^-_m = P_{i,m}N_{o,m} + N_{i,m}P_{o,m}.
```

The code multiplies these products by `2` before applying the existing
`f_0 = 2`, `f_{m>0} = 4` factors.  That leading `2` is not leaf physics; it
converts the `1/(2π)` Fourier coefficients above to the half-range cosine
moment used by the historical `compute_Z_matrices_aniso` implementation:
`(2/π)∫_0^π Γ(Δϕ)cos(mΔϕ)dΔϕ`.  This is the normalization expected by
vSmartMOM's elemental kernels, where `ϖ` is applied outside `Z`.
"""
function compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int) where {FT<:AbstractFloat}
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))

    (; R, T, nQuad) = mod
    R_leaf = FT(R)
    T_leaf = FT(T)
    ϖ = R_leaf + T_leaf

    nμ = length(μ)
    nm = m_max + 1
    Z⁺⁺ = zeros(FT, nμ, nμ, nm)
    Z⁻⁺ = zeros(FT, nμ, nμ, nm)
    ϖ <= zero(FT) && return Z⁺⁺, Z⁻⁺

    leaf_quad = _leaf_inclination_quadrature(LD, nQuad, FT)
    θₗ = leaf_quad.θ
    w_measure = leaf_quad.w_measure
    G = vec(CanopyOptics.G(Array(μ), LD))

    Pꜜ = zeros(FT, nμ, nm)
    Nꜜ = zeros(FT, nμ, nm)
    Pꜛ = zeros(FT, nμ, nm)
    Nꜛ = zeros(FT, nμ, nm)
    P = zeros(FT, nm)
    N = zeros(FT, nm)
    Pm = zeros(FT, nm)
    Nm = zeros(FT, nm)

    for l in eachindex(θₗ)
        μ_L = cos(θₗ[l])
        leaf_weight = w_measure[l]

        for i in eachindex(μ)
            _clipped_projection_moments!(P, N, μ[i], μ_L, m_max)
            _clipped_projection_moments!(Pm, Nm, -μ[i], μ_L, m_max)
            @inbounds for k in 1:nm
                Pꜜ[i, k] = P[k]
                Nꜜ[i, k] = N[k]
                Pꜛ[i, k] = Pm[k]
                Nꜛ[i, k] = Nm[k]
            end
        end

        for k in 1:nm, j in 1:nμ
            @inbounds begin
                Pj = Pꜜ[j, k]
                Nj = Nꜜ[j, k]
                for i in 1:nμ
                    Ψpp_same = _psi_same(Pj, Nj, Pꜜ[i, k], Nꜜ[i, k])
                    Ψpp_opp  = _psi_opposite(Pj, Nj, Pꜜ[i, k], Nꜜ[i, k])
                    Ψmp_same = _psi_same(Pj, Nj, Pꜛ[i, k], Nꜛ[i, k])
                    Ψmp_opp  = _psi_opposite(Pj, Nj, Pꜛ[i, k], Nꜛ[i, k])

                    Z⁺⁺[i, j, k] += leaf_weight * (T_leaf * Ψpp_same + R_leaf * Ψpp_opp)
                    Z⁻⁺[i, j, k] += leaf_weight * (T_leaf * Ψmp_same + R_leaf * Ψmp_opp)
                end
            end
        end
    end

    return _normalise_Z!(Z⁺⁺, Z⁻⁺, G, ϖ)
end

"""
    compute_Γ_isotropic(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT},
                        Ωᵒᵘᵗ::dirVector_μ{FT}) where FT

Analytic area scattering transfer function assuming an isotropic (spherical) leaf
angle distribution, following Shultis & Myneni (1988) Eq. (40):

``Γ_{\\mathrm{iso}}(β) = \\frac{ω}{3π}(\\sin β - β \\cos β) + \\frac{T}{3}\\cos β``

where `β = arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ)` is the scattering angle and `ω = R + T`.

Used for validation against the general anisotropic [`compute_Γ`](@ref).
"""
function compute_Γ_isotropic(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT
    (;R,T, nQuad) = mod
    β = acos( Ωᵒᵘᵗ ⋅ Ωⁱⁿ)
    ω = R + T
    
    # Eq 40 in Shultis and Myneni
    Γ = (ω/3π) * (sin(β) - β * cos(β)) + T/3*cos(β)
    return Γ
end

function compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering,μ::AbstractArray{FT,1},LD::AbstractLeafDistribution, Zup, Zdown, m::Int) where FT
    # Zup/Zdown are retained only for API compatibility with callers that
    # previously used the brute-force precomputed-azimuth path.
    return compute_Z_matrices_aniso(mod, μ, LD, m)
end


function precompute_Zazi(mod::BiLambertianCanopyScattering, μ::AbstractArray{FT,1}, LD::AbstractLeafDistribution) where FT
    (;R,T, nQuad) = mod
    nQuad = nQuad
    # Quadrature points in the azimuth:
    ϕ, w_azi = gauleg(nQuad,FT(0),FT(π));
    # Fourier weights (cosine decomposition)
    
    # Transmission (same direction)
    Zup = zeros(length(μ),length(μ), nQuad)
    # Reflection (change direction)
    Zdown = zeros(length(μ),length(μ), nQuad)
    
    # Create outgoing vectors in θ and ϕ
    dirOutꜛ = [dirVector_μ(a,b) for a in -μ, b in ϕ];
    dirOutꜜ = [dirVector_μ(a,b) for a in μ, b in ϕ];

    Threads.@threads for i in eachindex(μ)
        # Incoming beam at ϕ = 0
        Ωⁱⁿ = dirVector_μ(μ[i], FT(0));
        # Compute over μ and μ_azi:
        Zup[i,:,:]   = compute_Γ.((mod,),(Ωⁱⁿ,),dirOutꜛ, (LD,));
        Zdown[i,:,:] = compute_Γ.((mod,),(Ωⁱⁿ,),dirOutꜜ, (LD,));
        #Zup[:,:,i]   = compute_Γ_isotropic.((mod,),(Ωⁱⁿ,),dirOutꜛ);
        #Zdown[:,:,i] = compute_Γ_isotropic.((mod,),(Ωⁱⁿ,),dirOutꜜ);
    end
    return Zup, Zdown
end

function precompute_Zazi_(mod::BiLambertianCanopyScattering, μ::AbstractArray{FT,1}, LD::AbstractLeafDistribution) where FT
    (;R,T, nQuad) = mod
    # Quadrature points in μ
    n_μ  = length(μ);
    arr_type = typeof(μ);

    #μ,w         = CanopyOptics.gauleg(n_μ,   FT(0),  FT(1.0));
    dϕ,  w_azi  = CanopyOptics.gauleg(nQuad, FT(0),  FT(π));
    dϕᴸ, w_aziᴸ = CanopyOptics.gauleg(nQuad+1, FT(0),  FT(2π));
    θᴸ,wᴸ       = CanopyOptics.gauleg(nQuad, FT(0),FT(π/2));
    
    μᴸ = cos.(θᴸ)
    Fᵢ = pdf.(LD.LD,2θᴸ/π) * LD.scaling
    # Reshape stuff:
    μⁱⁿ   = reshape(arr_type(μ),  n_μ,  1,     1,     1,     1   );
    μᵒᵘᵗ  = reshape(arr_type(deepcopy(μ)), 1,   n_μ,    1,     1,     1   );
    _dϕ   = reshape(arr_type(dϕ), 1,    1,   nQuad,   1,     1   );
    _μᴸ   = reshape(arr_type(μᴸ), 1,    1,     1,   nQuad,   1   );
    _dϕᴸ  = reshape(arr_type(dϕᴸ),1,    1,     1,     1,   nQuad+1 );

    # Quadrature points
    wᴸ  = wᴸ .*  Fᵢ
    #_w_azi  = reshape(arr_type(w_azi),  1,  1,   nQuad,   1,     1   );
    _wᴸ     = reshape(arr_type(wᴸ),     1,  1,    1,    nQuad,   1   );
    _w_aziᴸ = reshape(arr_type(w_aziᴸ), 1,  1,    1,      1,    nQuad+1);

    integrand  =  CanopyOptics.leaf_dot_products.(μⁱⁿ, -μᵒᵘᵗ, _dϕ,_μᴸ, _dϕᴸ);
    iPos       = (integrand+abs.(integrand))./2;
    iNeg       = (integrand-abs.(integrand))./2;
    Γ⁻         = -1/2π * sum(sum(iNeg.*_w_aziᴸ, dims=5).*_wᴸ,dims=4); 
    Γ⁺         =  1/2π * sum(sum(iPos.*_w_aziᴸ, dims=5).*_wᴸ,dims=4);

    Γ = R .* Γ⁻ .+ T .* Γ⁺; 
    Γup = reshape(Γ, n_μ,n_μ, nQuad);
    
    integrand  =  CanopyOptics.leaf_dot_products.(μⁱⁿ, μᵒᵘᵗ, _dϕ,_μᴸ, _dϕᴸ);
    iPos       = (integrand+abs.(integrand))./2;
    iNeg       = (integrand-abs.(integrand))./2;
    Γ⁻         = -1/2π * sum(sum(iNeg.*_w_aziᴸ, dims=5).*_wᴸ,dims=4); 
    Γ⁺         =  1/2π * sum(sum(iPos.*_w_aziᴸ, dims=5).*_wᴸ,dims=4);
    
    Γ = R .* Γ⁻ .+ T .* Γ⁺; 
    Γdown = reshape(Γ, n_μ,n_μ, nQuad);
    return Γup,Γdown
end


"""
    K(κ::FT, α::FT) where FT

Returns the Nilson–Kuusk leaf-surface roughness reduction factor:

``K(κ, α) = e^{-κ \\tan |α|}``

- `κ ≈ 0.1–0.3` controls surface roughness (`κ = 0` → smooth Fresnel surface)
- `α` is the incidence half-angle in radians

Used in [`compute_reflection`](@ref) to attenuate specular reflectance for rough leaves.
"""
function K(κ::FT, α::FT) where FT
    exp(-κ * tan(abs(α)));
end

function leaf_dot_products(μⁱⁿ::FT, μᵒᵘᵗ::FT, dϕᵒᵘᵗ::FT, μᴸ::FT, dϕᴸ::FT) where FT
    (μⁱⁿ  * μᴸ + sqrt(1-μⁱⁿ^2)   * sqrt(1-μᴸ^2) * cos(dϕᴸ)) *
    (μᵒᵘᵗ * μᴸ + sqrt(1-μᵒᵘᵗ^2)  * sqrt(1-μᴸ^2) * cos(dϕᴸ - dϕᵒᵘᵗ))
end
