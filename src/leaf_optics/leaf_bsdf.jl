"""
    CoupledLeafBSDF(leaf, surface; η_t = 0, γ_t = 4)
    CoupledLeafBSDF(dorsiventral_leaf; η_t = 0, γ_t = 4)
    CoupledLeafBSDF(asymmetric_surface_leaf; η_t = 0, γ_t = 4)

Energy-conserving, Helmholtz-reciprocal, Stokes-polarized leaf bidirectional
scattering distribution function (BRDF + BTDF + Mueller matrix) coupling rough
dielectric microfacet surface optics consistently with internal plate transport:

1. **First-surface reflection (`M_{r,surf}`)**: Beckmann normal distribution
   (`tan² θ_h = (H_x² + H_y²) / H_z²`), `erfcx`-based Smith `G₂ = G₁(μ_i) G₁(μ_o)`,
   trigonometric-free Fresnel amplitudes, and algebraic meridian-plane Stokes
   rotations (`npol ∈ (1, 3, 4)`), obeying exact Hovenier reciprocity
   `M(ω_o -> ω_i) = P M(ω_i -> ω_o)ᵀ P` (`P = diag(1, 1, -1, 1)`).
2. **Internal volume BRDF (`f_{r,vol}`)**: Coupled through `T_{ext,s}(μ_i) = 1 - R_{surf,s}(μ_i)`
   and directional escape `E_s(μ_o) = T_{ext,s}(μ_o) / T_{hem,ext,s}` using the
   *same* [`RoughDielectricSurface`](@ref) (never a mismatched `ConeSurface`).
3. **Directional transmission (`f_t`)**: Symmetric Kelemen–Szirmay-Kalos /
   Kulla–Conty closure combining diffuse escape (`1 - η_t`) with a normalized
   forward lobe (`η_t ∈ [0, 1)`, concentration `γ_t >= 0`), guaranteeing strict
   nonnegativity, exact hemispherical conservation `R(μ_i) + T(μ_i) + A(μ_i) == 1`,
   and exact directional Helmholtz reciprocity `f_t(Ω_in -> Ω_out) == f_t(-Ω_out -> -Ω_in)`.
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

"""
    surface_brdf_mueller(μ_i, ϕ_i, μ_o, ϕ_o, n_eff, σ, npol = 1)

Single-bounce polarized microfacet BRDF Mueller matrix `M_{r,surf}` (`npol × npol`,
`npol ∈ (1, 3, 4)`) in the leaf local frame (`ẑ = (0, 0, 1)` outward face normal),
where `(μ_i, ϕ_i)` points toward the illumination source (`k_in = -ω_i`) and
`(μ_o, ϕ_o)` points toward the detector (`k_out = +ω_o`), with `μ_i, μ_o ∈ (0, 1]`.

In the meridian-plane Stokes convention of `CanopyOptics`, `+Q` denotes
perpendicular (`s` / horizontal) linear polarization and `-Q` denotes parallel
(`p` / in-meridian) linear polarization. All geometry and Stokes rotations are
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

    # Guard sqrt(1 - μ^2) at exact nadir μ == 1 so ForwardDiff wrt non-angle params never hits sqrt(0)
    ui = one(FT) - μi * μi
    uo = one(FT) - μo * μo
    sθi = prospect_primal(ui) <= zero(BFT) ? zero(FT) : sqrt(ui)
    sθo = prospect_primal(uo) <= zero(BFT) ? zero(FT) : sqrt(uo)

    ci_ϕ, si_ϕ = cos(ϕi), sin(ϕi)
    co_ϕ, so_ϕ = cos(ϕo), sin(ϕo)

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

    M = zeros(FT, n, n)
    if n == 1
        M[1, 1] = prefactor * f_plus
        return M
    end

    f_minus = (Rs - Rp) * half
    f_sp = rs * rp

    # Meridian polarization bases (smooth C^∞ even at θ_i = 0 or θ_o = 0!):
    # Incoming ray propagates along k_in = -ω_i; outgoing along k_out = +ω_o.
    # Right-handed triads (e_s, e_p = k × e_s, k):
    es_in_x = -si_ϕ; es_in_y = ci_ϕ; es_in_z = zero(FT)
    ep_in_x = μi * ci_ϕ; ep_in_y = μi * si_ϕ; ep_in_z = -sθi

    es_out_x = -so_ϕ; es_out_y = co_ϕ; es_out_z = zero(FT)
    ep_out_x = -μo * co_ϕ; ep_out_y = -μo * so_ϕ; ep_out_z = sθo

    inv_H = inv(Hnorm)
    hx = Hx * inv_H; hy = Hy * inv_H; hz = Hz * inv_H

    # Local reflection-plane normal u_s = (k_in × h) / |k_in × h|, where k_in = -ω_i
    kin_x = -ωix; kin_y = -ωiy; kin_z = -ωiz
    vs_x = kin_y * hz - kin_z * hy
    vs_y = kin_z * hx - kin_x * hz
    vs_z = kin_x * hy - kin_y * hx
    vs_norm2 = vs_x * vs_x + vs_y * vs_y + vs_z * vs_z

    us_x, us_y, us_z = if prospect_primal(vs_norm2) > eps(BFT)^2
        inv_vs = inv(sqrt(vs_norm2))
        (vs_x * inv_vs, vs_y * inv_vs, vs_z * inv_vs)
    else
        # Exact normal-incidence retroreflection limit (Rs == Rp, f_minus == 0):
        # M_loc = diag(Rs, Rs, -Rs, -Rs) is invariant to rotation of u_s ⊥ k_in!
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

    # Explicit 3x3 / 4x4 product prefactor * L2 * M_loc * L1 (zero intermediate matrices!)
    M[1, 1] = prefactor * f_plus
    M[1, 2] = prefactor * f_minus * cos2α1
    M[1, 3] = prefactor * f_minus * sin2α1

    M[2, 1] = prefactor * cos2α2 * f_minus
    M[2, 2] = prefactor * (cos2α2 * f_plus * cos2α1 + sin2α2 * f_sp * sin2α1)
    M[2, 3] = prefactor * (cos2α2 * f_plus * sin2α1 - sin2α2 * f_sp * cos2α1)

    M[3, 1] = prefactor * sin2α2 * f_minus
    M[3, 2] = prefactor * (sin2α2 * f_plus * cos2α1 - cos2α2 * f_sp * sin2α1)
    M[3, 3] = prefactor * (sin2α2 * f_plus * sin2α1 + cos2α2 * f_sp * cos2α1)

    if n == 4
        M[4, 4] = prefactor * f_sp
    end
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
    E_exit_gl12::AbstractVector{FT},
    γ::FT,
    c_norm::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    pi_ft = FT(π)
    acc = zero(FT)
    @inbounds for kv in 1:12
        v = (BFT(_SURF_GL12_X[kv]) + one(BFT)) * BFT(0.5)
        wv = BFT(_SURF_GL12_W[kv]) * BFT(0.5)
        μo = FT(sqrt(v))
        sθo = FT(sqrt(one(BFT) - v))
        E_o = E_exit_gl12[kv]
        az_acc = zero(FT)
        for kp in 1:12
            Δϕ = (BFT(_SURF_GL12_X[kp]) + one(BFT)) * (BFT(π) * BFT(0.5))
            wp = BFT(_SURF_GL12_W[kp]) * BFT(0.5)
            c_dot = μ * μo - sθ * sθo * cos(FT(Δϕ))
            az_acc += wp * exp(γ * (c_dot - one(FT)))
        end
        acc += wv * E_o * (pi_ft * az_acc)
    end
    return acc / c_norm
end

@inline function _btdf_lobe_escape_integral(
    μ::FT,
    sθ::FT,
    surf_exit::RoughDielectricSurface,
    nr::FT,
    t_hem_exit::FT,
    γ::FT,
    c_norm::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    n_eff = FT(surf_exit.n_scale) * nr
    σ = FT(surf_exit.σ)
    pi_ft = FT(π)
    acc = zero(FT)

    @inbounds for kv in 1:12
        v = (BFT(_SURF_GL12_X[kv]) + one(BFT)) * BFT(0.5)
        wv = BFT(_SURF_GL12_W[kv]) * BFT(0.5)
        μo = FT(sqrt(v))
        sθo = FT(sqrt(one(BFT) - v))
        r_ext_o = _microfacet_directional_reflectance_core(μo, sθo, n_eff, σ)
        E_o = (one(FT) - r_ext_o) / t_hem_exit

        az_acc = zero(FT)
        for kp in 1:12
            Δϕ = (BFT(_SURF_GL12_X[kp]) + one(BFT)) * (BFT(π) * BFT(0.5))
            wp = BFT(_SURF_GL12_W[kp]) * BFT(0.5)
            c_dot = μ * μo - sθ * sθo * cos(FT(Δϕ))
            az_acc += wp * exp(γ * (c_dot - one(FT)))
        end
        # ∫_0^{2π} dΔϕ = 2π * az_acc; ∫_0^1 μo dμo = 1/2 ∫_0^1 dv => factor π
        acc += wv * E_o * (pi_ft * az_acc)
    end
    return acc / c_norm
end

struct PreparedLeafBSDFSpectrum{FT<:Real}
    nᵣ::Vector{FT}
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
end

"""
    prepare_leaf_bsdf(model::CoupledLeafBSDF, optics::AbstractProspectProperties)

Precompute wavelength-dependent hemispherical internal budgets and symmetric
Kelemen–Szirmay-Kalos transmission normalization constants once per leaf state,
so subsequent multi-angle BRDF / BTDF / Stokes Mueller evaluations are O(1) per
wavelength and direction pair!
"""
function prepare_leaf_bsdf(
    model::CoupledLeafBSDF{L,SAd,SAb,FTM},
    optics::AbstractProspectProperties,
) where {L,SAd,SAb,FTM<:Real}
    _validate_optics_arrays(optics)
    nλ = length(optics.nᵣ)
    FTS1 = _surface_ft(model.surface_adaxial)
    FTS2 = _surface_ft(model.surface_abaxial)

    # Evaluate hemispherical internal budgets using the exact RoughDielectricSurface boundaries
    hem_ad, hem_ab, FT = if model.leaf isa Union{DorsiventralLeafProperties,AsymmetricSurfaceLeafProperties}
        b = dorsiventral_budget(model.leaf, optics; hemispherical = true)
        FT0 = float(promote_type(eltype(b.adaxial.transmission), FTS1, FTS2, FTM))
        (b.adaxial, b.abaxial, FT0)
    else
        _validate_single_zone_leaf(model.leaf)
        K_tot = _total_optical_depth_vector(model.leaf, optics)
        FT0 = float(promote_type(eltype(K_tot), typeof(model.leaf.N), eltype(optics.nᵣ), FTS1, FTS2, FTM))
        T_v = Vector{FT0}(undef, nλ); Rv_v = Vector{FT0}(undef, nλ); A_v = Vector{FT0}(undef, nλ)
        N_ft = FT0(model.leaf.N)
        @inbounds for idx in 1:nλ
            nr = FT0(optics.nᵣ[idx])
            s_hem = surface_interface_optics(model.surface_adaxial, one(FT0), nr; hemispherical = true)
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
    nr_vec = FT.(optics.nᵣ)
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
    same_surf = (model.surface_adaxial === model.surface_abaxial) ||
                (_same_val_and_partials(σ_ad, σ_ab) && _same_val_and_partials(ns_ad, ns_ab))
    has_lobe = prospect_primal(η_t) > zero(BFT)

    E_ad_gl12 = has_lobe ? Vector{FT}(undef, 12) : FT[]
    E_ab_gl12 = (has_lobe && !same_surf) ? Vector{FT}(undef, 12) : E_ad_gl12

    @inbounds for idx in 1:nλ
        nr = nr_vec[idx]
        s_ad = surface_interface_optics(model.surface_adaxial, one(FT), nr; hemispherical = true)
        th_ad = s_ad.t_hem_ext
        th_ab = if same_surf
            th_ad
        else
            surface_interface_optics(model.surface_abaxial, one(FT), nr; hemispherical = true).t_hem_ext
        end
        t_hem_ad[idx] = th_ad
        t_hem_ab[idx] = th_ab

        if has_lobe
            cn = two * c_sph * sqrt(inv(th_ad * th_ad) + inv(th_ab * th_ab))
            c_norm_vec[idx] = cn
            n_eff_ad = ns_ad * nr
            n_eff_ab = ns_ab * nr
            for kv in 1:12
                v = (BFT(_SURF_GL12_X[kv]) + one(BFT)) * BFT(0.5)
                μ = FT(sqrt(v))
                sθ = FT(sqrt(one(BFT) - v))
                r_ad_v = _microfacet_directional_reflectance_core(μ, sθ, n_eff_ad, σ_ad)
                E_ad_gl12[kv] = (one(FT) - r_ad_v) / th_ad
                if !same_surf
                    r_ab_v = _microfacet_directional_reflectance_core(μ, sθ, n_eff_ab, σ_ab)
                    E_ab_gl12[kv] = (one(FT) - r_ab_v) / th_ab
                end
            end
            acc_d = zero(FT)
            for kv in 1:12
                v = (BFT(_SURF_GL12_X[kv]) + one(BFT)) * BFT(0.5)
                wv = BFT(_SURF_GL12_W[kv]) * BFT(0.5)
                μ = FT(sqrt(v))
                sθ = FT(sqrt(one(BFT) - v))
                L_ab_v = _btdf_lobe_escape_integral_precomp(μ, sθ, E_ab_gl12, γ_t, cn)
                acc_d += wv * E_ad_gl12[kv] * (one(FT) - L_ab_v)
            end
            d_bar_vec[idx] = acc_d
        else
            c_norm_vec[idx] = one(FT)
            d_bar_vec[idx] = one(FT)
        end
    end

    return PreparedLeafBSDFSpectrum{FT}(
        nr_vec, σ_ad, ns_ad, σ_ab, ns_ab, η_t, γ_t,
        t_hem_ad, t_hem_ab,
        FT.(hem_ad.reflectance_volume), FT.(hem_ab.reflectance_volume),
        FT.(hem_ad.transmission),
        FT.(hem_ad.absorptance), FT.(hem_ab.absorptance),
        c_norm_vec, d_bar_vec,
    )
end

"""
    leaf_directional_budget(prep::PreparedLeafBSDFSpectrum, μ_i; face = :adaxial)
    leaf_directional_budget(model::CoupledLeafBSDF, optics, μ_i; face = :adaxial)

Exact directional-hemispherical `(; reflectance, reflectance_surface, reflectance_volume, transmission, absorptance)`
for collimated incidence at cosine `μ_i ∈ (0, 1]` on `face ∈ (:adaxial, :abaxial)`.
Satisfies `reflectance + transmission + absorptance == 1` identically across all `μ_i`.
"""
function leaf_directional_budget(
    prep::PreparedLeafBSDFSpectrum{FTP},
    μ_i_in::FTM;
    face::Symbol = :adaxial,
) where {FTP<:Real,FTM<:Real}
    face in (:adaxial, :abaxial) || throw(ArgumentError("face must be :adaxial or :abaxial"))
    FT = promote_type(FTP, FTM)
    BFT = _leaf_base_float(FT)
    μi = FT(μ_i_in)
    zero(BFT) < prospect_primal(μi) <= one(BFT) ||
        throw(ArgumentError("μ_i must lie in (0, 1]"))

    is_ad = face === :adaxial
    σ = is_ad ? prep.σ_ad : prep.σ_ab
    ns = is_ad ? prep.n_scale_ad : prep.n_scale_ab
    t_hem_ext = is_ad ? prep.t_hem_ext_ad : prep.t_hem_ext_ab
    r_vol_hem = is_ad ? prep.r_vol_hem_ad : prep.r_vol_hem_ab
    a_hem = is_ad ? prep.a_hem_ad : prep.a_hem_ab

    nλ = length(prep.nᵣ)
    R_tot = Vector{FT}(undef, nλ)
    R_surf = Vector{FT}(undef, nλ)
    R_vol = Vector{FT}(undef, nλ)
    T_dir = Vector{FT}(undef, nλ)
    A_dir = Vector{FT}(undef, nλ)

    @inbounds for idx in 1:nλ
        nr = prep.nᵣ[idx]
        r_s = microfacet_directional_reflectance(μi, ns * nr, σ)
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
)
    prep = prepare_leaf_bsdf(model, optics)
    return leaf_directional_budget(prep, μ_i; face)
end

"""
    leaf_bsdf_spectrum(prep::PreparedLeafBSDFSpectrum, μ_i, ϕ_i, μ_o, ϕ_o;
                       face = :adaxial, npol = 1)
    leaf_bsdf_spectrum(model::CoupledLeafBSDF, optics, μ_i, ϕ_i, μ_o, ϕ_o;
                       face = :adaxial, npol = 1)

Evaluate the full spectral leaf BSDF at incident source direction `(μ_i, ϕ_i)` on
`face ∈ (:adaxial, :abaxial)`, reflected detector direction `(μ_o, ϕ_o)` (`μ_o > 0`
on the illuminated face), and transmitted detector direction `(μ_o, ϕ_o)` (`μ_o > 0`
measured from the opposite face normal `-ẑ`, so straight-through forward continuation
of `(θ_i, ϕ_i = 0)` is `(θ_o = θ_i, ϕ_o = π)`):
- `brdf`: total scalar BRDF `[sr⁻¹]` (`nλ`-vector, `= brdf_surface + brdf_volume`)
- `btdf`: total scalar BTDF `[sr⁻¹]` (`nλ`-vector)
- `brdf_mueller`: `(npol, npol, nλ)` polarized BRDF Mueller array (`I -> I` equals `brdf`;
  diffuse internal multiple scattering populates `(1, 1)` unpolarized)
- `dolp`: reflected degree of linear polarization `hypot(M[2,1], M[3,1]) / M[1,1]` when
  `npol >= 3` (for unpolarized incident light), else `zeros(FT, nλ)`.
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
    FT = promote_type(FTP, FT1, FT2, FT3, FT4)
    BFT = _leaf_base_float(FT)

    μi = FT(μ_i_in); ϕi = FT(ϕ_i_in)
    μo = FT(μ_o_in); ϕo = FT(ϕ_o_in)
    (zero(BFT) < prospect_primal(μi) <= one(BFT) &&
     zero(BFT) < prospect_primal(μo) <= one(BFT)) ||
        throw(ArgumentError("μ_i and μ_o must lie in (0, 1]"))

    is_ad = face === :adaxial
    σ_in = is_ad ? prep.σ_ad : prep.σ_ab
    ns_in = is_ad ? prep.n_scale_ad : prep.n_scale_ab
    σ_tr = is_ad ? prep.σ_ab : prep.σ_ad
    ns_tr = is_ad ? prep.n_scale_ab : prep.n_scale_ad
    th_in = is_ad ? prep.t_hem_ext_ad : prep.t_hem_ext_ab
    th_tr = is_ad ? prep.t_hem_ext_ab : prep.t_hem_ext_ad
    rv_hem = is_ad ? prep.r_vol_hem_ad : prep.r_vol_hem_ab

    surf_in_obj = RoughDielectricSurface{FT}(σ_in, ns_in)
    surf_tr_obj = RoughDielectricSurface{FT}(σ_tr, ns_tr)

    ui = one(FT) - μi * μi
    uo = one(FT) - μo * μo
    sθi = prospect_primal(ui) <= zero(BFT) ? zero(FT) : sqrt(ui)
    sθo = prospect_primal(uo) <= zero(BFT) ? zero(FT) : sqrt(uo)
    c_trans_dot = μi * μo - sθi * sθo * cos(ϕo - ϕi)

    η_t = FT(prep.η_t)
    γ_t = FT(prep.γ_t)
    has_lobe = prospect_primal(η_t) > zero(BFT)
    inv_pi = inv(FT(π))

    nλ = length(prep.nᵣ)
    brdf_vec = Vector{FT}(undef, nλ)
    btdf_vec = Vector{FT}(undef, nλ)
    M_arr = zeros(FT, n, n, nλ)
    dolp_vec = zeros(FT, nλ)

    @inbounds for idx in 1:nλ
        nr = FT(prep.nᵣ[idx])
        n_eff_in = ns_in * nr
        n_eff_tr = ns_tr * nr

        # 1. Polarized surface BRDF Mueller matrix
        M_s = surface_brdf_mueller(μi, ϕi, μo, ϕo, n_eff_in, σ_in, n)

        # 2. Directional entry/escape factors
        r_ext_i = microfacet_directional_reflectance(μi, n_eff_in, σ_in)
        r_ext_o_refl = microfacet_directional_reflectance(μo, n_eff_in, σ_in)
        r_ext_o_tran = (_same_val_and_partials(σ_in, σ_tr) && _same_val_and_partials(ns_in, ns_tr)) ? r_ext_o_refl :
                       microfacet_directional_reflectance(μo, n_eff_tr, σ_tr)

        E_i = (one(FT) - r_ext_i) / FT(th_in[idx])
        E_o_refl = (one(FT) - r_ext_o_refl) / FT(th_in[idx])
        E_o_tran = (one(FT) - r_ext_o_tran) / FT(th_tr[idx])

        # 3. Diffuse unpolarized volume BRDF (symmetric in μ_i, μ_o!)
        f_vol = FT(rv_hem[idx]) * inv_pi * E_i * E_o_refl
        f_r = M_s[1, 1] + f_vol
        brdf_vec[idx] = f_r

        for cj in 1:n, ri in 1:n
            M_arr[ri, cj, idx] = M_s[ri, cj]
        end
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
            L_tr_i = _btdf_lobe_escape_integral(μi, sθi, surf_tr_obj, nr, FT(th_tr[idx]), γ_t, cn)
            L_in_o = _btdf_lobe_escape_integral(μo, sθo, surf_in_obj, nr, FT(th_in[idx]), γ_t, cn)
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
)
    prep = prepare_leaf_bsdf(model, optics)
    return leaf_bsdf_spectrum(prep, μ_i, ϕ_i, μ_o, ϕ_o; face, npol)
end
