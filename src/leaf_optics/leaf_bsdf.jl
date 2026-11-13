@doc raw"""
    CoupledLeafBSDF(leaf, surface; η_t = 0, γ_t = 4)
    CoupledLeafBSDF(dorsiventral_leaf; η_t = 0, γ_t = 4)
    CoupledLeafBSDF(asymmetric_surface_leaf; η_t = 0, γ_t = 4)

Energy-conserving, Helmholtz-reciprocal, Stokes-polarized leaf bidirectional
scattering distribution function (BRDF + BTDF + Mueller matrix) coupling rough
dielectric microfacet surface optics consistently with internal plate transport:

1. **First-surface reflection** (``\mathbf{M}_{r,\mathrm{surf}}``): Beckmann normal
   distribution (``\tan^2\theta_h = (H_x^2 + H_y^2) / H_z^2``), `erfcx`-based Smith
   shadowing–masking ``G_2 = G_1(\mu_i)\,G_1(\mu_o)``, trigonometric-free Fresnel
   amplitudes, and algebraic meridian-plane Stokes rotations (`npol` in ``\{1, 3, 4\}``),
   obeying exact Hovenier reciprocity
   ``\mathbf{M}(\boldsymbol{\omega}_o \to \boldsymbol{\omega}_i) = \mathbf{P}\,\mathbf{M}(\boldsymbol{\omega}_i \to \boldsymbol{\omega}_o)^\top \mathbf{P}``
   with ``\mathbf{P} = \operatorname{diag}(1, 1, -1, 1)``.
2. **Internal volume BRDF** (``f_{r,\mathrm{vol}}``): Coupled through external entry
   ``T_{\mathrm{ext},s}(\mu_i) = 1 - R_{\mathrm{surf},s}(\mu_i)`` and directional
   escape ``E_s(\mu_o) = T_{\mathrm{ext},s}(\mu_o) / T_{\mathrm{hem,ext},s}`` using
   the *same* [`RoughDielectricSurface`](@ref) (never a mismatched `ConeSurface`).
3. **Directional transmission** (``f_t``): Symmetric Kelemen–Szirmay-Kalos /
   Kulla–Conty closure combining diffuse escape (``1 - \eta_t``) with a normalized
   forward lobe (``\eta_t \in [0, 1)``, concentration ``\gamma_t \ge 0``),
   guaranteeing strict nonnegativity, exact hemispherical conservation
   ``R(\mu_i) + T(\mu_i) + A(\mu_i) = 1``, and exact directional Helmholtz
   reciprocity ``f_t(\boldsymbol{\Omega}_{\mathrm{in}} \to \boldsymbol{\Omega}_{\mathrm{out}}) = f_t(-\boldsymbol{\Omega}_{\mathrm{out}} \to -\boldsymbol{\Omega}_{\mathrm{in}})``.
"""
struct CoupledLeafBSDF{
    L<:AbstractLeafProperties,
    SAd<:RoughDielectricSurface,
    SAb<:RoughDielectricSurface,
    FT<:Real,
}
    leaf::L
    surface_adaxial::SAd
    surface_abaxial::SAb
    η_t::FT
    γ_t::FT

    function CoupledLeafBSDF(
        leaf::L,
        surface_adaxial::SAd,
        surface_abaxial::SAb,
        η_t::FT,
        γ_t::FT,
    ) where {
        L<:Union{
            LeafProspectProProperties,
            LeafMultiPigmentProperties,
            DorsiventralLeafProperties,
            AsymmetricSurfaceLeafProperties,
        },
        SAd<:RoughDielectricSurface,
        SAb<:RoughDielectricSurface,
        FT<:Real,
    }
        BFT = _leaf_base_float(FT)
        eta_v = prospect_primal(η_t)
        gam_v = prospect_primal(γ_t)
        isfinite(eta_v) && zero(BFT) <= eta_v < one(BFT) ||
            throw(ArgumentError("CoupledLeafBSDF requires finite η_t in [0, 1)"))
        isfinite(gam_v) && gam_v >= zero(BFT) ||
            throw(ArgumentError("CoupledLeafBSDF requires finite nonnegative γ_t"))

        if leaf isa Union{DorsiventralLeafProperties,AsymmetricSurfaceLeafProperties}
            (leaf.surface_adaxial isa RoughDielectricSurface &&
             leaf.surface_abaxial isa RoughDielectricSurface &&
             _same_val_and_partials(surface_adaxial.σ, leaf.surface_adaxial.σ) &&
             _same_val_and_partials(surface_adaxial.n_scale, leaf.surface_adaxial.n_scale) &&
             _same_val_and_partials(surface_abaxial.σ, leaf.surface_abaxial.σ) &&
             _same_val_and_partials(surface_abaxial.n_scale, leaf.surface_abaxial.n_scale)) ||
                throw(ArgumentError(
                    "CoupledLeafBSDF surface_adaxial and surface_abaxial must match leaf.surface_adaxial and leaf.surface_abaxial",
                ))
        else
            (_same_val_and_partials(surface_adaxial.σ, surface_abaxial.σ) &&
             _same_val_and_partials(surface_adaxial.n_scale, surface_abaxial.n_scale)) ||
                throw(ArgumentError(
                    "CoupledLeafBSDF with single-zone LeafProspectProProperties or LeafMultiPigmentProperties requires identical adaxial and abaxial surfaces; wrap the leaf in AsymmetricSurfaceLeafProperties or DorsiventralLeafProperties for distinct surfaces",
                ))
        end

        return new{L,SAd,SAb,FT}(leaf, surface_adaxial, surface_abaxial, η_t, γ_t)
    end
end

function CoupledLeafBSDF(
    leaf::Union{LeafProspectProProperties{FTL},LeafMultiPigmentProperties{FTL}},
    surface::RoughDielectricSurface{FTS} = RoughDielectricSurface{float(FTL)}();
    η_t::Real = zero(_leaf_base_float(promote_type(FTL, FTS))),
    γ_t::Real = _leaf_base_float(promote_type(FTL, FTS))(4),
) where {FTL<:Real,FTS<:Real}
    FT = float(promote_type(FTL, FTS, typeof(η_t), typeof(γ_t)))
    return CoupledLeafBSDF(leaf, surface, surface, FT(η_t), FT(γ_t))
end

function CoupledLeafBSDF(
    dors::DorsiventralLeafProperties;
    η_t::Real = zero(_leaf_base_float(typeof(dors.palisade.N))),
    γ_t::Real = _leaf_base_float(typeof(dors.palisade.N))(4),
)
    (dors.surface_adaxial isa RoughDielectricSurface &&
     dors.surface_abaxial isa RoughDielectricSurface) ||
        throw(ArgumentError("CoupledLeafBSDF(dorsiventral) requires RoughDielectricSurface on both faces"))
    FTS1 = _surface_ft(dors.surface_adaxial)
    FTS2 = _surface_ft(dors.surface_abaxial)
    FT = float(promote_type(typeof(dors.palisade.N), typeof(dors.spongy.N), FTS1, FTS2, typeof(η_t), typeof(γ_t)))
    return CoupledLeafBSDF(dors, dors.surface_adaxial, dors.surface_abaxial, FT(η_t), FT(γ_t))
end

function CoupledLeafBSDF(
    asym::AsymmetricSurfaceLeafProperties;
    η_t::Real = zero(_leaf_base_float(typeof(asym.leaf.N))),
    γ_t::Real = _leaf_base_float(typeof(asym.leaf.N))(4),
)
    (asym.surface_adaxial isa RoughDielectricSurface &&
     asym.surface_abaxial isa RoughDielectricSurface) ||
        throw(ArgumentError("CoupledLeafBSDF(asym) requires RoughDielectricSurface on both faces"))
    FTS1 = _surface_ft(asym.surface_adaxial)
    FTS2 = _surface_ft(asym.surface_abaxial)
    FT = float(promote_type(typeof(asym.leaf.N), FTS1, FTS2, typeof(η_t), typeof(γ_t)))
    return CoupledLeafBSDF(asym, asym.surface_adaxial, asym.surface_abaxial, FT(η_t), FT(γ_t))
end

@inline function _write_surface_brdf_mueller!(
    M_dest::AbstractArray{FT},
    idx::Int,
    μi::FT,
    ci_ϕ::FT,
    si_ϕ::FT,
    sθi::FT,
    μo::FT,
    co_ϕ::FT,
    so_ϕ::FT,
    sθo::FT,
    n_eff::FT,
    σ::FT,
    n::Int,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    ωix = sθi * ci_ϕ; ωiy = sθi * si_ϕ; ωiz = μi
    ωox = sθo * co_ϕ; ωoy = sθo * so_ϕ; ωoz = μo

    # Unnormalized half-vector H = ω_i + ω_o (H_z = μ_i + μ_o > 0 strictly everywhere!)
    Hx = ωix + ωox
    Hy = ωiy + ωoy
    Hz = ωiz + ωoz
    H2 = Hx * Hx + Hy * Hy + Hz * Hz
    Hnorm = sqrt(H2)
    Hz2 = Hz * Hz

    two = one(FT) + one(FT)
    half = inv(two)
    four = two + two

    tan2_th = (Hx * Hx + Hy * Hy) / Hz2
    cos2_th = Hz2 / H2
    D = beckmann_D(tan2_th, cos2_th, σ)
    G2 = smith_beckmann_G1(μi, σ) * smith_beckmann_G1(μo, σ)

    # Local incidence cosine c_i = ω_i ⋅ h = |H| / 2 and sin^2 α_inc = |ω_i - ω_o|^2 / 4
    c_loc = Hnorm * half
    dx = ωix - ωox; dy = ωiy - ωoy; dz = ωiz - ωoz
    s_loc2 = (dx * dx + dy * dy + dz * dz) / four

    n2 = n_eff * n_eff
    g = sqrt(n2 - s_loc2)
    rs = (c_loc - g) / (c_loc + g)
    rp = (n2 * c_loc - g) / (n2 * c_loc + g)

    prefactor = (D * G2) / (four * μi * μo)
    Rs = rs * rs
    Rp = rp * rp
    f_plus = (Rs + Rp) * half

    @inbounds if n == 1
        M_dest[1, 1, idx] = prefactor * f_plus
        return prefactor * f_plus
    end

    f_minus = (Rs - Rp) * half
    f_sp = rs * rp

    es_in_x = -si_ϕ; es_in_y = ci_ϕ; es_in_z = zero(FT)
    ep_in_x = μi * ci_ϕ; ep_in_y = μi * si_ϕ; ep_in_z = -sθi

    es_out_x = -so_ϕ; es_out_y = co_ϕ; es_out_z = zero(FT)
    ep_out_x = -μo * co_ϕ; ep_out_y = -μo * so_ϕ; ep_out_z = sθo

    inv_H = inv(Hnorm)
    hx = Hx * inv_H; hy = Hy * inv_H; hz = Hz * inv_H

    kin_x = -ωix; kin_y = -ωiy; kin_z = -ωiz
    vs_x = kin_y * hz - kin_z * hy
    vs_y = kin_z * hx - kin_x * hz
    vs_z = kin_x * hy - kin_y * hx
    vs_norm2 = vs_x * vs_x + vs_y * vs_y + vs_z * vs_z

    us_x, us_y, us_z = if prospect_primal(vs_norm2) > eps(BFT)^2
        inv_vs = inv(sqrt(vs_norm2))
        (vs_x * inv_vs, vs_y * inv_vs, vs_z * inv_vs)
    else
        (es_in_x, es_in_y, es_in_z)
    end

    c1 = es_in_x * us_x + es_in_y * us_y + es_in_z * us_z
    s1 = ep_in_x * us_x + ep_in_y * us_y + ep_in_z * us_z
    c2 = es_out_x * us_x + es_out_y * us_y + es_out_z * us_z
    s2 = ep_out_x * us_x + ep_out_y * us_y + ep_out_z * us_z

    cos2α1 = c1 * c1 - s1 * s1
    sin2α1 = two * c1 * s1
    cos2α2 = c2 * c2 - s2 * s2
    sin2α2 = two * c2 * s2

    @inbounds begin
        M_dest[1, 1, idx] = prefactor * f_plus
        M_dest[1, 2, idx] = prefactor * f_minus * cos2α1
        M_dest[1, 3, idx] = prefactor * f_minus * sin2α1

        M_dest[2, 1, idx] = prefactor * cos2α2 * f_minus
        M_dest[2, 2, idx] = prefactor * (cos2α2 * f_plus * cos2α1 + sin2α2 * f_sp * sin2α1)
        M_dest[2, 3, idx] = prefactor * (cos2α2 * f_plus * sin2α1 - sin2α2 * f_sp * cos2α1)

        M_dest[3, 1, idx] = prefactor * sin2α2 * f_minus
        M_dest[3, 2, idx] = prefactor * (sin2α2 * f_plus * cos2α1 - cos2α2 * f_sp * sin2α1)
        M_dest[3, 3, idx] = prefactor * (sin2α2 * f_plus * sin2α1 + cos2α2 * f_sp * cos2α1)

        if n == 4
            M_dest[4, 4, idx] = prefactor * f_sp
        end
    end
    return prefactor * f_plus
end

@doc raw"""
    surface_brdf_mueller(μ_i, ϕ_i, μ_o, ϕ_o, n_eff, σ, npol = 1)

Evaluate the single-bounce polarized microfacet BRDF Mueller matrix
``\mathbf{M}_{r,\mathrm{surf}}`` (size `npol × npol`, `npol` in ``\{1, 3, 4\}``)
in the leaf local frame (``\hat{\mathbf{z}} = (0, 0, 1)`` outward face normal),
where ``(\mu_i, \phi_i)`` points toward the illumination source
(``\mathbf{k}_{\mathrm{in}} = -\boldsymbol{\omega}_i``) and ``(\mu_o, \phi_o)``
points toward the detector (``\mathbf{k}_{\mathrm{out}} = +\boldsymbol{\omega}_o``),
with ``\mu_i, \mu_o \in (0, 1]``.

In the meridian-plane Stokes convention of `CanopyOptics`, ``+Q`` denotes
perpendicular (``s`` / horizontal) linear polarization and ``-Q`` denotes parallel
(``p`` / in-meridian) linear polarization. All geometry and Stokes rotations are
computed algebraically without `acos`, `asin`, or `atan` singularities.
"""
function surface_brdf_mueller(
    μ_i_in::FT1,
    ϕ_i_in::FT2,
    μ_o_in::FT3,
    ϕ_o_in::FT4,
    n_eff_in::FT5,
    σ_in::FT6,
    npol::Integer = 1,
) where {FT1<:Real,FT2<:Real,FT3<:Real,FT4<:Real,FT5<:Real,FT6<:Real}
    n = _validate_npol(npol)
    FT = float(promote_type(FT1, FT2, FT3, FT4, FT5, FT6))
    BFT = _leaf_base_float(FT)

    μi = FT(μ_i_in); ϕi = FT(ϕ_i_in)
    μo = FT(μ_o_in); ϕo = FT(ϕ_o_in)
    n_eff = FT(n_eff_in); σ = FT(σ_in)

    μi_v = prospect_primal(μi); μo_v = prospect_primal(μo)
    ne_v = prospect_primal(n_eff); sig_v = prospect_primal(σ)
    (isfinite(μi_v) && zero(BFT) < μi_v <= one(BFT) &&
     isfinite(μo_v) && zero(BFT) < μo_v <= one(BFT)) ||
        throw(ArgumentError("surface_brdf_mueller requires finite μ_i, μ_o in (0, 1]"))
    isfinite(ne_v) && ne_v > one(BFT) ||
        throw(ArgumentError("surface_brdf_mueller requires finite n_eff > 1"))
    isfinite(sig_v) && sig_v > zero(BFT) ||
        throw(ArgumentError("surface_brdf_mueller requires finite σ > 0"))

    ui = one(FT) - μi * μi
    uo = one(FT) - μo * μo
    sθi = prospect_primal(ui) <= zero(BFT) ? zero(FT) : sqrt(ui)
    sθo = prospect_primal(uo) <= zero(BFT) ? zero(FT) : sqrt(uo)

    si_ϕ, ci_ϕ = sincos(ϕi)
    so_ϕ, co_ϕ = sincos(ϕo)

    M = zeros(FT, n, n)
    _write_surface_brdf_mueller!(
        reshape(M, n, n, 1), 1, μi, ci_ϕ, si_ϕ, sθi, μo, co_ϕ, so_ϕ, sθo, n_eff, σ, n,
    )
    return M
end

@inline function _forward_lobe_sphere_scale(γ::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    two = one(FT) + one(FT)
    if prospect_primal(γ) <= BFT(1e-2)
        # 2π * (1/2 - γ/6 + γ^2/24 - γ^3/120 + γ^4/720 - γ^5/5040) (< 1 ulp at γ = 1e-2)
        g2 = γ * γ
        g4 = g2 * g2
        return two * FT(π) * (
            inv(two) - γ / BFT(6) + g2 / BFT(24) - (g2 * γ) / BFT(120) +
            g4 / BFT(720) - (g4 * γ) / BFT(5040)
        )
    end
    return two * FT(π) * (γ - one(FT) + exp(-γ)) / (γ * γ)
end

@inline function _btdf_lobe_escape_integral_precomp(
    μ::FT,
    sθ::FT,
    E_exit_gl12::AbstractVector,
    γ::FT,
    c_norm::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    pi_ft = FT(π)
    acc = zero(FT)
    @inbounds for kv in 1:12
        wv = BFT(_SURF_GL12_WV[kv])
        μo = FT(BFT(_SURF_GL12_MU[kv]))
        sθo = FT(BFT(_SURF_GL12_STH[kv]))
        E_o = FT(E_exit_gl12[kv])
        base_arg = γ * (μ * μo - one(FT))
        s_scale = -γ * sθ * sθo
        az_acc = zero(FT)
        for kp in 1:12
            wp = BFT(_SURF_GL12_WV[kp])
            cos_dϕ = FT(BFT(_SURF_GL12_COS_DPHI[kp]))
            az_acc += wp * exp(muladd(s_scale, cos_dϕ, base_arg))
        end
        acc += wv * E_o * (pi_ft * az_acc)
    end
    return acc / c_norm
end

@inline function _btdf_lobe_escape_integral(
    μ::FT,
    sθ::FT,
    surf_exit::RoughDielectricSurface,
    n_surf_exit::FT,
    t_hem_exit::FT,
    γ::FT,
    c_norm::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    n_eff = FT(surf_exit.n_scale) * n_surf_exit
    σ = FT(surf_exit.σ)
    pi_ft = FT(π)
    acc = zero(FT)

    @inbounds for kv in 1:12
        wv = BFT(_SURF_GL12_WV[kv])
        μo = FT(BFT(_SURF_GL12_MU[kv]))
        sθo = FT(BFT(_SURF_GL12_STH[kv]))
        r_ext_o = _microfacet_directional_reflectance_core(μo, sθo, n_eff, σ)
        E_o = (one(FT) - r_ext_o) / t_hem_exit
        base_arg = γ * (μ * μo - one(FT))
        s_scale = -γ * sθ * sθo

        az_acc = zero(FT)
        for kp in 1:12
            wp = BFT(_SURF_GL12_WV[kp])
            cos_dϕ = FT(BFT(_SURF_GL12_COS_DPHI[kp]))
            az_acc += wp * exp(muladd(s_scale, cos_dϕ, base_arg))
        end
        # ∫_0^{2π} dΔϕ = 2π * az_acc; ∫_0^1 μo dμo = 1/2 ∫_0^1 dv => factor π
        acc += wv * E_o * (pi_ft * az_acc)
    end
    return acc / c_norm
end

function _compute_gl12_escape_matrix(
    ::Type{FT},
    σ::FT,
    n_scale::FT,
    n_surf_vec::Vector{FT},
    t_hem_vec::Vector{FT},
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    nλ = length(n_surf_vec)
    e_mat = Matrix{FT}(undef, 12, nλ)
    w_gl12 = ntuple(Val(12)) do kv
        _microfacet_dir_cheb32_weights(FT, FT(BFT(_SURF_GL12_MU[kv])), FT(BFT(_SURF_GL12_STH[kv])), σ)
    end
    @inbounds for idx in 1:nλ
        n_eff = n_scale * n_surf_vec[idx]
        th = t_hem_vec[idx]
        if prospect_primal(n_eff) < BFT(1.08)
            for kv in 1:12
                μ = FT(BFT(_SURF_GL12_MU[kv]))
                sθ = FT(BFT(_SURF_GL12_STH[kv]))
                r_v = _microfacet_directional_reflectance_core(μ, sθ, n_eff, σ)
                e_mat[kv, idx] = (one(FT) - r_v) / th
            end
        else
            f_vals = ntuple(Val(_CHEB32_M)) do j
                fresnel_unpolarized_cos(FT(BFT(_CHEB32_C[j])), FT(BFT(_CHEB32_S2[j])), n_eff)
            end
            for kv in 1:12
                wkv = w_gl12[kv]
                r_v = zero(FT)
                for j in 1:_CHEB32_M
                    r_v = muladd(wkv[j], f_vals[j], r_v)
                end
                e_mat[kv, idx] = (one(FT) - r_v) / th
            end
        end
    end
    return e_mat
end

"""
    PreparedLeafBSDFSpectrum{FT<:Real}

Wavelength-resolved precomputation returned by [`prepare_leaf_bsdf`](@ref) storing
adaxial and abaxial cuticular surface parameters, outer refractive-index spectra
(`n_surf_ad`, `n_surf_ab`), hemispherical volume/transmission/absorption budgets,
12-node directional escape profiles, and Kelemen–Szirmay-Kalos directional
transmission normalization constants. Enables `O(1)` multi-angle BRDF, BTDF, and
Stokes Mueller evaluations via [`leaf_bsdf_spectrum`](@ref) and [`leaf_directional_budget`](@ref).
"""
struct PreparedLeafBSDFSpectrum{FT<:Real}
    nᵣ::Vector{FT}
    n_surf_ad::Vector{FT}
    n_surf_ab::Vector{FT}
    σ_ad::FT
    n_scale_ad::FT
    σ_ab::FT
    n_scale_ab::FT
    η_t::FT
    γ_t::FT
    t_hem_ext_ad::Vector{FT}
    t_hem_ext_ab::Vector{FT}
    r_vol_hem_ad::Vector{FT}
    r_vol_hem_ab::Vector{FT}
    t_hem::Vector{FT}
    a_hem_ad::Vector{FT}
    a_hem_ab::Vector{FT}
    c_norm::Vector{FT}
    d_bar::Vector{FT}
    e_gl12_ad::Matrix{FT}
    e_gl12_ab::Matrix{FT}
end

function PreparedLeafBSDFSpectrum{FT}(
    nᵣ::Vector{FT},
    n_surf_ad::Vector{FT},
    n_surf_ab::Vector{FT},
    σ_ad::FT,
    n_scale_ad::FT,
    σ_ab::FT,
    n_scale_ab::FT,
    η_t::FT,
    γ_t::FT,
    t_hem_ext_ad::Vector{FT},
    t_hem_ext_ab::Vector{FT},
    r_vol_hem_ad::Vector{FT},
    r_vol_hem_ab::Vector{FT},
    t_hem::Vector{FT},
    a_hem_ad::Vector{FT},
    a_hem_ab::Vector{FT},
    c_norm::Vector{FT},
    d_bar::Vector{FT},
) where {FT<:Real}
    e_ad = _compute_gl12_escape_matrix(FT, σ_ad, n_scale_ad, n_surf_ad, t_hem_ext_ad)
    e_ab = _compute_gl12_escape_matrix(FT, σ_ab, n_scale_ab, n_surf_ab, t_hem_ext_ab)
    return PreparedLeafBSDFSpectrum{FT}(
        nᵣ, n_surf_ad, n_surf_ab, σ_ad, n_scale_ad, σ_ab, n_scale_ab, η_t, γ_t,
        t_hem_ext_ad, t_hem_ext_ab, r_vol_hem_ad, r_vol_hem_ab, t_hem,
        a_hem_ad, a_hem_ab, c_norm, d_bar, e_ad, e_ab,
    )
end

function Base.:(==)(a::PreparedLeafBSDFSpectrum, b::PreparedLeafBSDFSpectrum)
    return a.nᵣ == b.nᵣ && a.n_surf_ad == b.n_surf_ad && a.n_surf_ab == b.n_surf_ab &&
           a.σ_ad == b.σ_ad && a.n_scale_ad == b.n_scale_ad &&
           a.σ_ab == b.σ_ab && a.n_scale_ab == b.n_scale_ab &&
           a.η_t == b.η_t && a.γ_t == b.γ_t &&
           a.t_hem_ext_ad == b.t_hem_ext_ad && a.t_hem_ext_ab == b.t_hem_ext_ab &&
           a.r_vol_hem_ad == b.r_vol_hem_ad && a.r_vol_hem_ab == b.r_vol_hem_ab &&
           a.t_hem == b.t_hem && a.a_hem_ad == b.a_hem_ad && a.a_hem_ab == b.a_hem_ab &&
           a.c_norm == b.c_norm && a.d_bar == b.d_bar &&
           a.e_gl12_ad == b.e_gl12_ad && a.e_gl12_ab == b.e_gl12_ab
end

function Base.hash(a::PreparedLeafBSDFSpectrum, h::UInt)
    h = hash(:PreparedLeafBSDFSpectrum, h)
    for f in (
        a.nᵣ, a.n_surf_ad, a.n_surf_ab, a.σ_ad, a.n_scale_ad, a.σ_ab, a.n_scale_ab,
        a.η_t, a.γ_t, a.t_hem_ext_ad, a.t_hem_ext_ab, a.r_vol_hem_ad, a.r_vol_hem_ab,
        a.t_hem, a.a_hem_ad, a.a_hem_ab, a.c_norm, a.d_bar, a.e_gl12_ad, a.e_gl12_ab,
    )
        h = hash(f, h)
    end
    return h
end

"""
    prepare_leaf_bsdf(model::CoupledLeafBSDF, optics::AbstractProspectProperties;
                      n_surface = nothing,
                      n_surface_adaxial = n_surface,
                      n_surface_abaxial = n_surface)

Precompute wavelength-dependent hemispherical internal budgets and symmetric
Kelemen–Szirmay-Kalos transmission normalization constants once per leaf state,
so subsequent multi-angle BRDF / BTDF / Stokes Mueller evaluations are O(1) per
wavelength and direction pair. Pass `n_surface = load_wax_refractive_index(λ_bnds)`
(or face-specific `n_surface_adaxial`, `n_surface_abaxial` on dorsiventral /
asymmetric-surface leaves) to include epicuticular wax dispersion.
"""
function prepare_leaf_bsdf(
    model::CoupledLeafBSDF{L,SAd,SAb,FTM},
    optics::AbstractProspectProperties;
    n_surface::Union{Nothing,AbstractVector{<:Real}} = nothing,
    n_surface_adaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
    n_surface_abaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
) where {L,SAd,SAb,FTM<:Real}
    _validate_optics_arrays(optics)
    nλ = length(optics.nᵣ)
    _validate_n_surface_vec(n_surface_adaxial, nλ, "n_surface_adaxial")
    _validate_n_surface_vec(n_surface_abaxial, nλ, "n_surface_abaxial")

    FTS1 = _surface_ft(model.surface_adaxial)
    FTS2 = _surface_ft(model.surface_abaxial)
    FTW1 = n_surface_adaxial === nothing ? eltype(optics.nᵣ) : eltype(n_surface_adaxial)
    FTW2 = n_surface_abaxial === nothing ? eltype(optics.nᵣ) : eltype(n_surface_abaxial)

    # Evaluate hemispherical internal budgets using the exact RoughDielectricSurface boundaries
    hem_ad, hem_ab, FT = if model.leaf isa Union{DorsiventralLeafProperties,AsymmetricSurfaceLeafProperties}
        b = dorsiventral_budget(
            model.leaf, optics;
            hemispherical = true,
            n_surface_adaxial = n_surface_adaxial,
            n_surface_abaxial = n_surface_abaxial,
        )
        FT0 = float(promote_type(eltype(b.adaxial.transmission), FTS1, FTS2, FTM, FTW1, FTW2))
        (b.adaxial, b.abaxial, FT0)
    else
        _validate_single_zone_leaf(model.leaf)
        if !(n_surface_adaxial === n_surface_abaxial)
            v_ad_chk = n_surface_adaxial === nothing ? optics.nᵣ : n_surface_adaxial
            v_ab_chk = n_surface_abaxial === nothing ? optics.nᵣ : n_surface_abaxial
            off_ad_c = firstindex(v_ad_chk) - 1
            off_ab_c = firstindex(v_ab_chk) - 1
            all(idx -> _same_val_and_partials(v_ad_chk[off_ad_c + idx], v_ab_chk[off_ab_c + idx]), 1:nλ) ||
                throw(ArgumentError(
                    "prepare_leaf_bsdf with single-zone LeafProspectProProperties or LeafMultiPigmentProperties requires identical adaxial and abaxial n_surface; wrap the leaf in AsymmetricSurfaceLeafProperties or DorsiventralLeafProperties for distinct surfaces",
                ))
        end
        K_tot = _total_optical_depth_vector(model.leaf, optics)
        FT0 = float(promote_type(eltype(K_tot), typeof(model.leaf.N), eltype(optics.nᵣ), FTS1, FTS2, FTM, FTW1, FTW2))
        T_v = Vector{FT0}(undef, nλ); Rv_v = Vector{FT0}(undef, nλ); A_v = Vector{FT0}(undef, nλ)
        N_ft = FT0(model.leaf.N)
        prep_hem0 = _prepare_surface_cheb(FT0, model.surface_adaxial, one(FT0), true)
        off_ad0 = n_surface_adaxial === nothing ? 0 : firstindex(n_surface_adaxial) - 1
        @inbounds for idx in 1:nλ
            nr = FT0(optics.nᵣ[idx])
            ns_idx = n_surface_adaxial === nothing ? nr : FT0(n_surface_adaxial[off_ad0 + idx])
            s_hem = _eval_prepared_surface(model.surface_adaxial, prep_hem0, one(FT0), nr, ns_idx, true)
            bh = _asymmetric_top_prospect_budget(
                FT0(K_tot[idx]) / N_ft, N_ft, s_hem.t_hem_ext, s_hem.t_int_diffuse, nr,
            )
            T_v[idx] = bh.transmission
            Rv_v[idx] = bh.reflectance_volume
            A_v[idx] = bh.absorptance
        end
        nt = (; transmission = T_v, reflectance_volume = Rv_v, absorptance = A_v)
        (nt, nt, FT0)
    end

    BFT = _leaf_base_float(FT)
    nr_vec = Vector{FT}(optics.nᵣ)
    n_surf_ad_vec = n_surface_adaxial === nothing ? copy(nr_vec) : copyto!(Vector{FT}(undef, nλ), n_surface_adaxial)
    n_surf_ab_vec = n_surface_abaxial === nothing ? copy(nr_vec) : copyto!(Vector{FT}(undef, nλ), n_surface_abaxial)
    σ_ad = FT(model.surface_adaxial.σ)
    ns_ad = FT(model.surface_adaxial.n_scale)
    σ_ab = FT(model.surface_abaxial.σ)
    ns_ab = FT(model.surface_abaxial.n_scale)
    η_t = FT(model.η_t)
    γ_t = FT(model.γ_t)

    t_hem_ad = Vector{FT}(undef, nλ)
    t_hem_ab = Vector{FT}(undef, nλ)
    c_norm_vec = Vector{FT}(undef, nλ)
    d_bar_vec = Vector{FT}(undef, nλ)

    two = one(FT) + one(FT)
    c_sph = _forward_lobe_sphere_scale(γ_t)
    same_nsurf = (n_surf_ad_vec === n_surf_ab_vec) ||
                 all(idx -> _same_val_and_partials(n_surf_ad_vec[idx], n_surf_ab_vec[idx]), 1:nλ)
    same_surf = (((model.surface_adaxial === model.surface_abaxial) ||
                  (_same_val_and_partials(σ_ad, σ_ab) && _same_val_and_partials(ns_ad, ns_ab))) &&
                 same_nsurf)
    has_lobe = prospect_primal(η_t) > zero(BFT)

    prep_hem_ad = _prepare_surface_cheb(FT, model.surface_adaxial, one(FT), true)
    prep_hem_ab = same_surf ? prep_hem_ad : _prepare_surface_cheb(FT, model.surface_abaxial, one(FT), true)

    @inbounds for idx in 1:nλ
        nr = nr_vec[idx]
        nsurf_ad = n_surf_ad_vec[idx]
        nsurf_ab = n_surf_ab_vec[idx]
        th_ad = _eval_prepared_surface(model.surface_adaxial, prep_hem_ad, one(FT), nr, nsurf_ad, true).t_hem_ext
        th_ab = same_surf ? th_ad :
                _eval_prepared_surface(model.surface_abaxial, prep_hem_ab, one(FT), nr, nsurf_ab, true).t_hem_ext
        t_hem_ad[idx] = th_ad
        t_hem_ab[idx] = th_ab
    end

    e_gl12_ad = _compute_gl12_escape_matrix(FT, σ_ad, ns_ad, n_surf_ad_vec, t_hem_ad)
    e_gl12_ab = same_surf ? copy(e_gl12_ad) : _compute_gl12_escape_matrix(FT, σ_ab, ns_ab, n_surf_ab_vec, t_hem_ab)

    @inbounds for idx in 1:nλ
        if has_lobe
            th_ad = t_hem_ad[idx]
            th_ab = t_hem_ab[idx]
            cn = two * c_sph * sqrt(inv(th_ad * th_ad) + inv(th_ab * th_ab))
            c_norm_vec[idx] = cn
            e_ab_col = view(e_gl12_ab, :, idx)
            acc_d = zero(FT)
            for kv in 1:12
                wv = BFT(_SURF_GL12_WV[kv])
                μ = FT(BFT(_SURF_GL12_MU[kv]))
                sθ = FT(BFT(_SURF_GL12_STH[kv]))
                L_ab_v = _btdf_lobe_escape_integral_precomp(μ, sθ, e_ab_col, γ_t, cn)
                acc_d += wv * e_gl12_ad[kv, idx] * (one(FT) - L_ab_v)
            end
            d_bar_vec[idx] = acc_d
        else
            c_norm_vec[idx] = one(FT)
            d_bar_vec[idx] = one(FT)
        end
    end

    return PreparedLeafBSDFSpectrum{FT}(
        nr_vec, n_surf_ad_vec, n_surf_ab_vec,
        σ_ad, ns_ad, σ_ab, ns_ab, η_t, γ_t,
        t_hem_ad, t_hem_ab,
        FT.(hem_ad.reflectance_volume), FT.(hem_ab.reflectance_volume),
        FT.(hem_ad.transmission),
        FT.(hem_ad.absorptance), FT.(hem_ab.absorptance),
        c_norm_vec, d_bar_vec, e_gl12_ad, e_gl12_ab,
    )
end

@inline function _prepare_dir_only_cheb(::Type{FT}, μi::FT, σ::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    u_sin2 = one(FT) - μi * μi
    if prospect_primal(u_sin2) <= zero(BFT)
        if u_sin2 isa ForwardDiff.Dual && !_same_val_and_partials(u_sin2, zero(FT))
            return nothing
        end
        return _microfacet_dir_cheb32_weights(FT, one(FT), zero(FT), σ)
    end
    return _microfacet_dir_cheb32_weights(FT, μi, sqrt(u_sin2), σ)
end

@inline function _eval_dir_only_cheb(
    w::Union{Nothing,NTuple{_CHEB32_M,FT}},
    μi::FT,
    n_eff::FT,
    σ::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    if w === nothing || prospect_primal(n_eff) < BFT(1.08)
        return microfacet_directional_reflectance(μi, n_eff, σ)
    end
    return _eval_cheb32_fresnel(w, n_eff)
end

@doc raw"""
    leaf_directional_budget(prep::PreparedLeafBSDFSpectrum, μ_i; face = :adaxial)
    leaf_directional_budget(model::CoupledLeafBSDF, optics, μ_i;
                            face = :adaxial, n_surface = nothing,
                            n_surface_adaxial = n_surface, n_surface_abaxial = n_surface)

Compute the exact directional-hemispherical budget
`(; reflectance, reflectance_surface, reflectance_volume, transmission, absorptance)`
for collimated incidence at cosine ``\mu_i \in (0, 1]`` on `face` (`:adaxial` or `:abaxial`).
Satisfies ``R(\mu_i) + T(\mu_i) + A(\mu_i) = 1`` across all ``\mu_i \in (0, 1]``.
"""
function leaf_directional_budget(
    prep::PreparedLeafBSDFSpectrum{FTP},
    μ_i_in::FTM;
    face::Symbol = :adaxial,
) where {FTP<:Real,FTM<:Real}
    face in (:adaxial, :abaxial) || throw(ArgumentError("face must be :adaxial or :abaxial"))
    FT = float(promote_type(FTP, FTM))
    BFT = _leaf_base_float(FT)
    μi = FT(μ_i_in)
    zero(BFT) < prospect_primal(μi) <= one(BFT) ||
        throw(ArgumentError("μ_i must lie in (0, 1]"))

    is_ad = face === :adaxial
    σ = FT(is_ad ? prep.σ_ad : prep.σ_ab)
    ns = FT(is_ad ? prep.n_scale_ad : prep.n_scale_ab)
    n_surf = is_ad ? prep.n_surf_ad : prep.n_surf_ab
    t_hem_ext = is_ad ? prep.t_hem_ext_ad : prep.t_hem_ext_ab
    r_vol_hem = is_ad ? prep.r_vol_hem_ad : prep.r_vol_hem_ab
    a_hem = is_ad ? prep.a_hem_ad : prep.a_hem_ab

    w_dir = _prepare_dir_only_cheb(FT, μi, σ)

    nλ = length(prep.nᵣ)
    R_tot = Vector{FT}(undef, nλ)
    R_surf = Vector{FT}(undef, nλ)
    R_vol = Vector{FT}(undef, nλ)
    T_dir = Vector{FT}(undef, nλ)
    A_dir = Vector{FT}(undef, nλ)

    @inbounds for idx in 1:nλ
        nsurf_idx = FT(n_surf[idx])
        r_s = _eval_dir_only_cheb(w_dir, μi, ns * nsurf_idx, σ)
        t_ext = one(FT) - r_s
        E_i = t_ext / t_hem_ext[idx]
        rv = E_i * r_vol_hem[idx]
        td = E_i * prep.t_hem[idx]
        ad = E_i * a_hem[idx]
        R_surf[idx] = r_s
        R_vol[idx] = rv
        R_tot[idx] = r_s + rv
        T_dir[idx] = td
        A_dir[idx] = ad
    end

    return (;
        reflectance = R_tot,
        reflectance_surface = R_surf,
        reflectance_volume = R_vol,
        transmission = T_dir,
        absorptance = A_dir,
    )
end

function leaf_directional_budget(
    model::CoupledLeafBSDF,
    optics::AbstractProspectProperties,
    μ_i::Real;
    face::Symbol = :adaxial,
    n_surface::Union{Nothing,AbstractVector{<:Real}} = nothing,
    n_surface_adaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
    n_surface_abaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
)
    prep = prepare_leaf_bsdf(model, optics; n_surface_adaxial, n_surface_abaxial)
    return leaf_directional_budget(prep, μ_i; face)
end

@doc raw"""
    leaf_bsdf_spectrum(prep::PreparedLeafBSDFSpectrum, μ_i, ϕ_i, μ_o, ϕ_o;
                       face = :adaxial, npol = 1)
    leaf_bsdf_spectrum(model::CoupledLeafBSDF, optics, μ_i, ϕ_i, μ_o, ϕ_o;
                       face = :adaxial, npol = 1)

Evaluate the spectral leaf BSDF at incident source direction ``(\mu_i, \phi_i)`` on
`face` (`:adaxial` or `:abaxial`), reflected detector direction ``(\mu_o, \phi_o)``
(``\mu_o > 0`` on the illuminated face), and transmitted detector direction
``(\mu_o, \phi_o)`` (``\mu_o > 0`` measured from the opposite face normal
``-\hat{\mathbf{z}}``, so straight-through forward continuation of
``(\theta_i, \phi_i = 0)`` is ``(\theta_o = \theta_i, \phi_o = \pi)``):
- `brdf`: total scalar BRDF ``[\mathrm{sr}^{-1}]`` (`nλ`-vector, equal to
  ``f_{r,\mathrm{surf}} + f_{r,\mathrm{vol}}``).
- `btdf`: total scalar BTDF ``[\mathrm{sr}^{-1}]`` (`nλ`-vector).
- `brdf_mueller`: `(npol, npol, nλ)` polarized BRDF Mueller array (``M_{11}`` equals
  `brdf`; diffuse internal multiple scattering populates ``(1, 1)`` unpolarized).
- `dolp`: reflected degree of linear polarization
  ``\sqrt{M_{21}^2 + M_{31}^2} / M_{11}`` when `npol >= 3` (for unpolarized
  incident light), or `zeros(FT, nλ)` when `npol == 1`.
"""
function leaf_bsdf_spectrum(
    prep::PreparedLeafBSDFSpectrum{FTP},
    μ_i_in::FT1,
    ϕ_i_in::FT2,
    μ_o_in::FT3,
    ϕ_o_in::FT4;
    face::Symbol = :adaxial,
    npol::Integer = 1,
) where {FTP<:Real,FT1<:Real,FT2<:Real,FT3<:Real,FT4<:Real}
    face in (:adaxial, :abaxial) || throw(ArgumentError("face must be :adaxial or :abaxial"))
    n = _validate_npol(npol)
    FT = float(promote_type(FTP, FT1, FT2, FT3, FT4))
    BFT = _leaf_base_float(FT)

    μi = FT(μ_i_in); ϕi = FT(ϕ_i_in)
    μo = FT(μ_o_in); ϕo = FT(ϕ_o_in)
    (isfinite(prospect_primal(μi)) && zero(BFT) < prospect_primal(μi) <= one(BFT) &&
     isfinite(prospect_primal(μo)) && zero(BFT) < prospect_primal(μo) <= one(BFT)) ||
        throw(ArgumentError("μ_i and μ_o must lie in (0, 1]"))

    is_ad = face === :adaxial
    σ_in = FT(is_ad ? prep.σ_ad : prep.σ_ab)
    ns_in = FT(is_ad ? prep.n_scale_ad : prep.n_scale_ab)
    n_surf_in = is_ad ? prep.n_surf_ad : prep.n_surf_ab
    σ_tr = FT(is_ad ? prep.σ_ab : prep.σ_ad)
    ns_tr = FT(is_ad ? prep.n_scale_ab : prep.n_scale_ad)
    n_surf_tr = is_ad ? prep.n_surf_ab : prep.n_surf_ad
    th_in = is_ad ? prep.t_hem_ext_ad : prep.t_hem_ext_ab
    th_tr = is_ad ? prep.t_hem_ext_ab : prep.t_hem_ext_ad
    rv_hem = is_ad ? prep.r_vol_hem_ad : prep.r_vol_hem_ab
    e_gl12_in = is_ad ? prep.e_gl12_ad : prep.e_gl12_ab
    e_gl12_tr = is_ad ? prep.e_gl12_ab : prep.e_gl12_ad

    ui = one(FT) - μi * μi
    uo = one(FT) - μo * μo
    sθi = prospect_primal(ui) <= zero(BFT) ? zero(FT) : sqrt(ui)
    sθo = prospect_primal(uo) <= zero(BFT) ? zero(FT) : sqrt(uo)
    si_ϕ, ci_ϕ = sincos(ϕi)
    so_ϕ, co_ϕ = sincos(ϕo)
    c_trans_dot = μi * μo - sθi * sθo * (co_ϕ * ci_ϕ + so_ϕ * si_ϕ)

    η_t = FT(prep.η_t)
    γ_t = FT(prep.γ_t)
    has_lobe = prospect_primal(η_t) > zero(BFT)
    inv_pi = inv(FT(π))
    same_sigma = _same_val_and_partials(σ_in, σ_tr)

    w_dir_i = _prepare_dir_only_cheb(FT, μi, σ_in)
    w_dir_o_in = _prepare_dir_only_cheb(FT, μo, σ_in)
    w_dir_o_tr = same_sigma ? w_dir_o_in : _prepare_dir_only_cheb(FT, μo, σ_tr)

    nλ = length(prep.nᵣ)
    brdf_vec = Vector{FT}(undef, nλ)
    btdf_vec = Vector{FT}(undef, nλ)
    M_arr = zeros(FT, n, n, nλ)
    dolp_vec = zeros(FT, nλ)

    @inbounds for idx in 1:nλ
        nsurf_i = FT(n_surf_in[idx])
        nsurf_t = FT(n_surf_tr[idx])
        n_eff_in = ns_in * nsurf_i
        n_eff_tr = ns_tr * nsurf_t

        # 1. Polarized surface BRDF Mueller matrix (written directly into M_arr[:, :, idx])
        f_surf = _write_surface_brdf_mueller!(
            M_arr, idx, μi, ci_ϕ, si_ϕ, sθi, μo, co_ϕ, so_ϕ, sθo, n_eff_in, σ_in, n,
        )

        # 2. Directional entry/escape factors
        r_ext_i = _eval_dir_only_cheb(w_dir_i, μi, n_eff_in, σ_in)
        r_ext_o_refl = _eval_dir_only_cheb(w_dir_o_in, μo, n_eff_in, σ_in)
        r_ext_o_tran = (same_sigma && _same_val_and_partials(n_eff_in, n_eff_tr)) ? r_ext_o_refl :
                       _eval_dir_only_cheb(w_dir_o_tr, μo, n_eff_tr, σ_tr)

        E_i = (one(FT) - r_ext_i) / FT(th_in[idx])
        E_o_refl = (one(FT) - r_ext_o_refl) / FT(th_in[idx])
        E_o_tran = (one(FT) - r_ext_o_tran) / FT(th_tr[idx])

        # 3. Diffuse unpolarized volume BRDF (symmetric in μ_i, μ_o!)
        f_vol = FT(rv_hem[idx]) * inv_pi * E_i * E_o_refl
        f_r = f_surf + f_vol
        brdf_vec[idx] = f_r
        M_arr[1, 1, idx] = f_r

        if n >= 3 && prospect_primal(f_r) > zero(BFT)
            q_pol = M_arr[2, 1, idx]
            u_pol = M_arr[3, 1, idx]
            # Avoid hypot(0, 0) derivative NaN at normal incidence or U == 0 in principal plane!
            pol_amp = if iszero(prospect_primal(q_pol)) && iszero(prospect_primal(u_pol))
                zero(FT)
            elseif _same_val_and_partials(u_pol, zero(FT))
                abs(q_pol)
            else
                hypot(q_pol, u_pol)
            end
            dolp_vec[idx] = pol_amp / f_r
        end

        # 4. Symmetric Kelemen-Szirmay-Kalos / Kulla-Conty coupled BTDF
        t_hem_val = FT(prep.t_hem[idx])
        if !has_lobe
            btdf_vec[idx] = t_hem_val * inv_pi * E_i * E_o_tran
        else
            cn = FT(prep.c_norm[idx])
            p_tilde = exp(γ_t * (c_trans_dot - one(FT))) / cn
            L_tr_i = _btdf_lobe_escape_integral_precomp(μi, sθi, view(e_gl12_tr, :, idx), γ_t, cn)
            L_in_o = _btdf_lobe_escape_integral_precomp(μo, sθo, view(e_gl12_in, :, idx), γ_t, cn)
            closure = (one(FT) - L_tr_i) * (one(FT) - L_in_o) * inv_pi / FT(prep.d_bar[idx])
            btdf_vec[idx] = t_hem_val * E_i * E_o_tran *
                            ((one(FT) - η_t) * inv_pi + η_t * (p_tilde + closure))
        end
    end

    return (;
        brdf = brdf_vec,
        btdf = btdf_vec,
        brdf_mueller = M_arr,
        dolp = dolp_vec,
    )
end

function leaf_bsdf_spectrum(
    model::CoupledLeafBSDF,
    optics::AbstractProspectProperties,
    μ_i::Real,
    ϕ_i::Real,
    μ_o::Real,
    ϕ_o::Real;
    face::Symbol = :adaxial,
    npol::Integer = 1,
    n_surface::Union{Nothing,AbstractVector{<:Real}} = nothing,
    n_surface_adaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
    n_surface_abaxial::Union{Nothing,AbstractVector{<:Real}} = n_surface,
)
    prep = prepare_leaf_bsdf(model, optics; n_surface_adaxial, n_surface_abaxial)
    return leaf_bsdf_spectrum(prep, μ_i, ϕ_i, μ_o, ϕ_o; face, npol)
end
