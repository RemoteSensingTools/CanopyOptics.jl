"""
    compute_reflection(mod::SpecularCanopyScattering,
                       Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}, LD) where FT

Angle-based variant of the specular leaf bidirectional scattering coefficient.

This uses the same Knyazikhin-Marshak specular-leaf model as the `dirVector_μ`
method below, but accepts polar-angle direction vectors. The `dirVector_μ`
method is the one used by Z-matrix assembly; this method is retained for
interactive angle-grid checks.
"""
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
                       LD::AbstractLeafDistribution, m::Int;
                       quadrature = CanopyQuadrature()) where FT

Computes the Fourier-`m` component of the single-scattering phase matrices
`(𝐙⁺⁺, 𝐙⁻⁺)` for a specular leaf surface by integrating
[`compute_reflection`](@ref) over the azimuthal quadrature grid.

- `𝐙⁺⁺[i,j]`: same-hemisphere scattering (μ>0 → μ>0, forward scatter)
- `𝐙⁻⁺[i,j]`: opposite-hemisphere scattering (μ>0 → μ<0, backscatter)

Azimuth integration uses `quadrature.n_azimuth` Gauss-Legendre points over `[0, 2π]`;
Fourier weights are `cos(m ϕ)`.
"""
function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::Array{FT,1},
                            LD::AbstractLeafDistribution,
                            m::Int;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    (;nᵣ, κ) = mod
    q = _resolve_quadrature(quadrature, nQuad)
    ZFT = promote_type(FT, typeof(nᵣ), typeof(κ))
    # Transmission (same direction)
    𝐙⁺⁺ = zeros(ZFT, length(μ), length(μ))
    # Reflection (change direction)
    𝐙⁻⁺ = zeros(ZFT, length(μ), length(μ))
    
    # Quadrature points in the azimuth:
    ϕ, w_azi = gauleg(q.n_azimuth,FT(0),FT(2π));
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

"""
    compute_Z_matrices(mod::SpecularCanopyScattering,
                       μ::AbstractVector, LD::AbstractLeafDistribution,
                       m::Integer; quadrature = CanopyQuadrature())

Convenience method for non-`Array` vector inputs. It materializes `μ` with
`collect` and delegates to the specular azimuth-integration method.
"""
function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, collect(μ), LD, Int(m); quadrature = q)
end
"""
    K(κ::FT, α::FT) where FT

Returns the Nilson–Kuusk leaf-surface roughness reduction factor:

``K(κ, α) = e^{-κ \\tan |α|}``

- `κ ≈ 0.1–0.3` controls surface roughness (`κ = 0` → smooth Fresnel surface)
- `α` is the incidence half-angle in radians

Used in [`compute_reflection`](@ref) to attenuate specular reflectance for rough leaves.
"""
function K(κ::FTκ, α::FTα) where {FTκ<:Real,FTα<:Real}
    FT = promote_type(FTκ, FTα)
    exp(-FT(κ) * tan(abs(FT(α))))
end
