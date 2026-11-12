"""
    compute_reflection(mod::SpecularCanopyScattering, Ωⁱⁿ, Ωᵒᵘᵗ, LD)

Compute the specular bidirectional scattering coefficient for direction pair
`(Ωⁱⁿ, Ωᵒᵘᵗ)`, following Knyazikhin & Marshak, *Discrete Ordinates Method for
Photon Transport in Leaf Canopies*, Eq. 2.39:

``f_s(Ω' \\to Ω) = \\frac{1}{8}\\,g_L(θ^*)\\,K(κ, α^*)\\,F_r(n_r, α^*)``

where `θ*` is the polar angle of the specular leaf normal (from [`getSpecularΩ`](@ref)),
`α* = arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ) / 2` is the incidence half-angle, `K(κ, α*)` is the
Nilson–Kuusk roughness factor, and `Fᵣ(nᵣ, α*)` is unpolarized Fresnel reflectance.
Accepts [`DirectionCosine`](@ref) directions (construct from `(θ, ϕ)` via
[`direction_from_angles`](@ref)).

For Stokes-vector propagation, use [`compute_reflection_mueller`](@ref). Note that
`Ωⁱⁿ` points toward the illumination source (`-k_in`) and `Ωᵒᵘᵗ` points along the
reflected ray (`+k_out`); only reflection is modelled here (for coupled surface,
volume, and directional transmission, see [`CoupledLeafBSDF`](@ref)).
"""
function compute_reflection(
    mod::SpecularCanopyScattering,
    Ωⁱⁿ::DirectionCosine{FT1},
    Ωᵒᵘᵗ::DirectionCosine{FT2},
    LD::AbstractLeafDistribution,
) where {FT1<:Real,FT2<:Real}
    (; nᵣ, κ) = mod
    Ωstar, αstar = getSpecularΩ(Ωⁱⁿ, Ωᵒᵘᵗ)
    FT = typeof(αstar)
    BFT = _leaf_base_float(FT)
    u_L = one(FT) - Ωstar.μ * Ωstar.μ
    s_L = prospect_primal(u_L) <= zero(BFT) ? zero(FT) : sqrt(u_L)
    θstar = atan(s_L, abs(Ωstar.μ))
    u_arg = (FT(2) * θstar) / FT(π)
    # Eq. 2.39 in "Discrete Ordinates Method for Photon Transport in Leaf Canopies", page 59
    return FT(1 / 8) * _lad_unit_pdf(LD.LD, u_arg) * LD.scaling * K(κ, αstar) * Fᵣ(nᵣ, αstar)
end

function _normalise_specular_projected_area!(Z⁺⁺::AbstractMatrix{FT},
                                             Z⁻⁺::AbstractMatrix{FT},
                                             G_vec::AbstractVector,
                                             m::Int,
                                             npol::Int) where {FT}
    ff = m == 0 ? FT(2) : FT(4)
    nμ = length(G_vec)
    @inbounds for j in 1:nμ, sj in 1:npol
        col = _stokes_index(j, sj, npol)
        g = FT(G_vec[j])
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

"""
    compute_Z_matrices(mod::SpecularCanopyScattering, μ::AbstractVector,
                       LD::AbstractLeafDistribution, m::Integer;
                       quadrature = CanopyQuadrature(), npol = 1) -> (Z⁺⁺, Z⁻⁺)

Compute the Fourier-`m` component of the single-scattering phase matrices
`(Z⁺⁺, Z⁻⁺)` for a specular leaf surface by integrating [`compute_reflection`](@ref)
(or [`compute_reflection_mueller`](@ref) when `npol > 1`) over the azimuthal
quadrature grid:
- `Z⁺⁺[i, j]`: same-hemisphere scattering (`μ > 0 → μ > 0`, forward scatter).
- `Z⁻⁺[i, j]`: opposite-hemisphere scattering (`μ > 0 → μ < 0`, backscatter).

Azimuth integration uses `quadrature.n_azimuth` Gauss-Legendre points over `[0, 2π]`
with `cos(m ϕ)` weights for scalar `I` (`npol = 1`) and the vSmartMOM Stokes
`cos`/`sin` parity kernel for `npol = 3` or `4`. Columns are normalized by `G(μ_in)`
while retaining the Fresnel/roughness strength in the returned kernel.
"""
function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where {FT<:Real}
    (; nᵣ, κ) = mod
    m_int = Int(m)
    q = _resolve_quadrature(quadrature, nQuad)
    n = _validate_npol(npol)
    ZFT = float(promote_type(FT, typeof(nᵣ), typeof(κ), _lad_ft(LD)))
    QFT = _leaf_base_float(ZFT)
    DFT = float(promote_type(FT, QFT))

    μ_vec = μ isa Vector{FT} ? μ : copyto!(Vector{FT}(undef, length(μ)), μ)
    nμ = length(μ_vec)
    ϕ, w_azi = gauleg(q.n_azimuth, zero(QFT), QFT(2π))
    G_in = ZFT.(vec(G(μ_vec, LD; quadrature = q)))

    cos_w = Vector{ZFT}(undef, length(ϕ))
    @inbounds for ia in eachindex(ϕ)
        cos_w[ia] = ZFT(w_azi[ia]) * ZFT(cos(m_int * ϕ[ia]))
    end

    if n == 1
        𝐙⁺⁺ = zeros(ZFT, nμ, nμ)
        𝐙⁻⁺ = zeros(ZFT, nμ, nμ)

        @inbounds for j in 1:nμ
            # Downward incoming stream j travels along k_in = DirectionCosine(+μ[j], 0),
            # so the outward vector toward the illumination source is Ω_src = -k_in = DirectionCosine(-μ[j], π).
            Ωⁱⁿ = DirectionCosine(-DFT(μ_vec[j]), DFT(π))
            for i in 1:nμ
                μ_i = DFT(μ_vec[i])
                acc_pp = zero(ZFT)
                acc_mp = zero(ZFT)
                for ia in eachindex(ϕ)
                    weight = cos_w[ia]
                    dϕ = DFT(ϕ[ia])
                    acc_pp += weight * compute_reflection(mod, Ωⁱⁿ, DirectionCosine(μ_i, dϕ), LD)
                    acc_mp += weight * compute_reflection(mod, Ωⁱⁿ, DirectionCosine(-μ_i, dϕ), LD)
                end
                𝐙⁺⁺[i, j] = acc_pp
                𝐙⁻⁺[i, j] = acc_mp
            end
        end
        return _normalise_specular_projected_area!(𝐙⁺⁺, 𝐙⁻⁺, G_in, m_int, n)
    end

    sin_w = Vector{ZFT}(undef, length(ϕ))
    @inbounds for ia in eachindex(ϕ)
        sin_w[ia] = ZFT(w_azi[ia]) * ZFT(sin(m_int * ϕ[ia]))
    end

    𝐙⁺⁺ = zeros(ZFT, n * nμ, n * nμ)
    𝐙⁻⁺ = zeros(ZFT, n * nμ, n * nμ)

    @inbounds for j in 1:nμ
        Ωⁱⁿ = DirectionCosine(-DFT(μ_vec[j]), DFT(π))
        col1 = _stokes_index(j, 1, n)
        col2 = _stokes_index(j, 2, n)
        col3 = _stokes_index(j, 3, n)
        col4 = _stokes_index(j, 4, n)
        for i in 1:nμ
            μ_i = DFT(μ_vec[i])
            pp11 = zero(ZFT); pp12 = zero(ZFT); pp13 = zero(ZFT)
            pp21 = zero(ZFT); pp22 = zero(ZFT); pp23 = zero(ZFT)
            pp31 = zero(ZFT); pp32 = zero(ZFT); pp33 = zero(ZFT); pp44 = zero(ZFT)

            mp11 = zero(ZFT); mp12 = zero(ZFT); mp13 = zero(ZFT)
            mp21 = zero(ZFT); mp22 = zero(ZFT); mp23 = zero(ZFT)
            mp31 = zero(ZFT); mp32 = zero(ZFT); mp33 = zero(ZFT); mp44 = zero(ZFT)

            for ia in eachindex(ϕ)
                dϕ = DFT(ϕ[ia])
                cw = cos_w[ia]
                sw = sin_w[ia]

                a11, a12, a13, a21, a22, a23, a31, a32, a33, a44 =
                    _specular_mueller_elements(mod, Ωⁱⁿ, DirectionCosine(μ_i, dϕ), LD)
                pp11 += cw * a11; pp12 += cw * a12; pp13 += sw * a13
                pp21 += cw * a21; pp22 += cw * a22; pp23 += sw * a23
                pp31 += sw * a31; pp32 += sw * a32; pp33 += cw * a33
                pp44 += cw * a44

                b11, b12, b13, b21, b22, b23, b31, b32, b33, b44 =
                    _specular_mueller_elements(mod, Ωⁱⁿ, DirectionCosine(-μ_i, dϕ), LD)
                mp11 += cw * b11; mp12 += cw * b12; mp13 += sw * b13
                mp21 += cw * b21; mp22 += cw * b22; mp23 += sw * b23
                mp31 += sw * b31; mp32 += sw * b32; mp33 += cw * b33
                mp44 += cw * b44
            end

            row1 = _stokes_index(i, 1, n)
            row2 = _stokes_index(i, 2, n)
            row3 = _stokes_index(i, 3, n)

            𝐙⁺⁺[row1, col1] = pp11; 𝐙⁺⁺[row1, col2] = pp12; 𝐙⁺⁺[row1, col3] = pp13
            𝐙⁺⁺[row2, col1] = pp21; 𝐙⁺⁺[row2, col2] = pp22; 𝐙⁺⁺[row2, col3] = pp23
            𝐙⁺⁺[row3, col1] = pp31; 𝐙⁺⁺[row3, col2] = pp32; 𝐙⁺⁺[row3, col3] = pp33

            𝐙⁻⁺[row1, col1] = mp11; 𝐙⁻⁺[row1, col2] = mp12; 𝐙⁻⁺[row1, col3] = mp13
            𝐙⁻⁺[row2, col1] = mp21; 𝐙⁻⁺[row2, col2] = mp22; 𝐙⁻⁺[row2, col3] = mp23
            𝐙⁻⁺[row3, col1] = mp31; 𝐙⁻⁺[row3, col2] = mp32; 𝐙⁻⁺[row3, col3] = mp33

            if n == 4
                row4 = _stokes_index(i, 4, n)
                𝐙⁺⁺[row4, col4] = pp44
                𝐙⁻⁺[row4, col4] = mp44
            end
        end
    end
    return _normalise_specular_projected_area!(𝐙⁺⁺, 𝐙⁻⁺, G_in, m_int, n)
end
"""
    K(κ, α)

Compute the Nilson–Kuusk leaf-surface roughness reduction factor:

``K(κ, α) = e^{-κ \\tan |α|}``

- `κ ≈ 0.1–0.3` controls surface roughness (`κ = 0` gives a smooth Fresnel surface).
- `α` is the incidence half-angle in radians, with `abs(α) ≤ π/2`.
- At grazing incidence (`abs(α) = π/2`), returns `0` for `κ > 0` and `1` for `κ = 0`
  without evaluating `tan(Float32(π/2))`.

Used in [`compute_reflection`](@ref) to attenuate specular reflectance for rough leaves.
"""
function K(κ::FTκ, α::FTα) where {FTκ<:Real,FTα<:Real}
    FT = float(promote_type(FTκ, FTα))
    BFT = _leaf_base_float(FT)
    roughness, angle = FT(κ), abs(FT(α))
    if prospect_primal(angle) >= BFT(π) / BFT(2)
        return (iszero(prospect_primal(roughness)) && _same_val_and_partials(roughness, zero(FT))) ?
               one(FT) : zero(FT)
    end
    return exp(-roughness * tan(angle))
end
