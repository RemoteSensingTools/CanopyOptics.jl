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

For Stokes-vector propagation, use [`compute_reflection_mueller`](@ref), which
constructs the Fresnel Mueller matrix and rotates it into the scattering-plane
basis used by the canopy Fourier Z matrices.

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
                       quadrature = CanopyQuadrature(), npol = 1) where FT

Computes the Fourier-`m` component of the single-scattering phase matrices
`(𝐙⁺⁺, 𝐙⁻⁺)` for a specular leaf surface by integrating
[`compute_reflection`](@ref), or [`compute_reflection_mueller`](@ref) when
`npol > 1`, over the azimuthal quadrature grid.

- `𝐙⁺⁺[i,j]`: same-hemisphere scattering (μ>0 → μ>0, forward scatter)
- `𝐙⁻⁺[i,j]`: opposite-hemisphere scattering (μ>0 → μ<0, backscatter)

Azimuth integration uses `quadrature.n_azimuth` Gauss-Legendre points over `[0, 2π]`;
Fourier weights are `cos(m ϕ)` for scalar I and the vSmartMOM Stokes
cos/sin kernel for vector cases. `npol = 1` is the scalar default;
`npol = 3` and `npol = 4` return Stokes-block matrices.

The raw Knyazikhin-Marshak specular coefficient is divided by
`G(μ_in)` column-by-column so it uses the same projected-area convention as the
bi-Lambertian canopy kernels. The Fresnel/roughness strength remains folded
into the returned kernel because a direction-independent specular
single-scattering albedo is not yet part of the public API.
"""
function _normalise_specular_projected_area!(Z⁺⁺::AbstractMatrix{FT},
                                             Z⁻⁺::AbstractMatrix{FT},
                                             G::AbstractVector,
                                             m::Int,
                                             npol::Int) where {FT}
    ff = m == 0 ? FT(2) : FT(4)
    nμ = length(G)
    @inbounds for j in 1:nμ, sj in 1:npol
        col = _stokes_index(j, sj, npol)
        g = FT(G[j])
        if _real_value(g) <= 0
            Z⁺⁺[:, col] .= zero(FT)
            Z⁻⁺[:, col] .= zero(FT)
        else
            scale = ff / g
            for row in axes(Z⁺⁺, 1)
                Z⁺⁺[row, col] *= scale
                Z⁻⁺[row, col] *= scale
            end
        end
    end
    return Z⁺⁺, Z⁻⁺
end

function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::Array{FT,1},
                            LD::AbstractLeafDistribution,
                            m::Int;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    (;nᵣ, κ) = mod
    q = _resolve_quadrature(quadrature, nQuad)
    n = _validate_npol(npol)
    ZFT = promote_type(FT, typeof(nᵣ), typeof(κ))

    # Quadrature points in the azimuth:
    ϕ, w_azi = gauleg(q.n_azimuth, FT(0), FT(2π))
    G_in = ZFT.(vec(G(μ, LD)))

    if n == 1
        𝐙⁺⁺ = zeros(ZFT, length(μ), length(μ))
        𝐙⁻⁺ = zeros(ZFT, length(μ), length(μ))

        for j in eachindex(μ)
            Ωⁱⁿ = dirVector_μ(μ[j], FT(0))
            for i in eachindex(μ)
                acc_pp = zero(ZFT)
                acc_mp = zero(ZFT)
                for ia in eachindex(ϕ)
                    weight = ZFT(w_azi[ia]) * ZFT(cos(m * ϕ[ia]))
                    acc_pp += weight * compute_reflection(mod, Ωⁱⁿ,
                                                          dirVector_μ(μ[i], ϕ[ia]), LD)
                    acc_mp += weight * compute_reflection(mod, Ωⁱⁿ,
                                                          dirVector_μ(-μ[i], ϕ[ia]), LD)
                end
                𝐙⁺⁺[i, j] = acc_pp
                𝐙⁻⁺[i, j] = acc_mp
            end
        end
        return _normalise_specular_projected_area!(𝐙⁺⁺, 𝐙⁻⁺, G_in, m, n)
    end

    nμ = length(μ)
    𝐙⁺⁺ = zeros(ZFT, n * nμ, n * nμ)
    𝐙⁻⁺ = zeros(ZFT, n * nμ, n * nμ)

    for j in eachindex(μ)
        Ωⁱⁿ = dirVector_μ(μ[j], FT(0))
        for i in eachindex(μ)
            for ia in eachindex(ϕ)
                dϕ = ϕ[ia]
                wϕ = ZFT(w_azi[ia])
                Mpp = compute_reflection_mueller(mod, Ωⁱⁿ, dirVector_μ(μ[i], dϕ), LD, n)
                Mmp = compute_reflection_mueller(mod, Ωⁱⁿ, dirVector_μ(-μ[i], dϕ), LD, n)

                @inbounds for sj in 1:n, si in 1:n
                    az = ZFT(_azimuthal_kernel(si, sj, m, dϕ))
                    row = _stokes_index(i, si, n)
                    col = _stokes_index(j, sj, n)
                    𝐙⁺⁺[row, col] += wϕ * az * Mpp[si, sj]
                    𝐙⁻⁺[row, col] += wϕ * az * Mmp[si, sj]
                end
            end
        end
    end
    return _normalise_specular_projected_area!(𝐙⁺⁺, 𝐙⁻⁺, G_in, m, n)
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
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, collect(μ), LD, Int(m);
                              quadrature = q, npol = npol)
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
