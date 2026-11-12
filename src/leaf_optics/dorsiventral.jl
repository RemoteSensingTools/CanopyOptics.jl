# Two-sided Redheffer star-product and asymmetric/dorsiventral plate operators.
# Every operator tracks direct zone absorptance (palisade & spongy) without
# 1 - R - T subtraction, preserving conservation into the k -> 0 limit.

struct _TwoSidedTwoZoneLayer{FT<:Real}
    t_dn::FT
    r_top::FT
    a1_dn::FT
    a2_dn::FT
    t_up::FT
    r_bot::FT
    a1_up::FT
    a2_up::FT
end

@inline function _star_two_zone(L1::_TwoSidedTwoZoneLayer{FT}, L2::_TwoSidedTwoZoneLayer{FT}) where {FT<:Real}
    c = one(FT) - L1.r_bot * L2.r_top
    inv_c = inv(c)

    t1_c = L1.t_dn * inv_c
    t_dn = t1_c * L2.t_dn
    r_top = L1.r_top + t1_c * L2.r_top * L1.t_up
    a1_dn = L1.a1_dn + t1_c * (L2.a1_dn + L2.r_top * L1.a1_up)
    a2_dn = L1.a2_dn + t1_c * (L2.a2_dn + L2.r_top * L1.a2_up)

    t2_c = L2.t_up * inv_c
    t_up = t2_c * L1.t_up
    r_bot = L2.r_bot + t2_c * L1.r_bot * L2.t_dn
    a1_up = L2.a1_up + t2_c * (L1.a1_up + L1.r_bot * L2.a1_dn)
    a2_up = L2.a2_up + t2_c * (L1.a2_up + L1.r_bot * L2.a2_dn)

    return _TwoSidedTwoZoneLayer{FT}(t_dn, r_top, a1_dn, a2_dn, t_up, r_bot, a1_up, a2_up)
end

"""
Single homogeneous slab `k >= 0` bounded on top by `(t_top_in, t_top_out)` and on
bottom by `(t_bot_in, t_bot_out)`, assigning all absorption to `zone` (`1` or `2`).
Also returns `(r_vol_top, r_vol_bot)` (multiple-scattered volume reflectance
excluding first-surface `1 - t_in` reflection).
"""
@inline function _single_bounded_plate(
    k::FT,
    t_top_in::FT,
    t_top_out::FT,
    t_bot_in::FT,
    t_bot_out::FT,
    zone::Int,
) where {FT<:Real}
    att = prospect_internal_attenuation(k)
    τ, a = att.transmission, att.absorptance
    r_top_ext = one(FT) - t_top_in
    r_bot_ext = one(FT) - t_bot_in
    r_top_int = one(FT) - t_top_out
    r_bot_int = one(FT) - t_bot_out
    two = one(FT) + one(FT)

    # Strictly positive cancellation-free form of 1 - r_top_int * r_bot_int * τ^2
    d = (t_top_out + t_bot_out * r_top_int) + r_top_int * r_bot_int * a * (two - a)
    inv_d = inv(d)

    t_dn = t_top_in * τ * t_bot_out * inv_d
    r_vol_top = t_top_in * (τ * τ * r_bot_int * t_top_out) * inv_d
    r_top = r_top_ext + r_vol_top
    a_dn = t_top_in * a * (one(FT) + r_bot_int * τ) * inv_d

    t_up = t_bot_in * τ * t_top_out * inv_d
    r_vol_bot = t_bot_in * (τ * τ * r_top_int * t_bot_out) * inv_d
    r_bot = r_bot_ext + r_vol_bot
    a_up = t_bot_in * a * (one(FT) + r_top_int * τ) * inv_d

    z = zero(FT)
    layer = if zone == 1
        _TwoSidedTwoZoneLayer{FT}(t_dn, r_top, a_dn, z, t_up, r_bot, a_up, z)
    else
        _TwoSidedTwoZoneLayer{FT}(t_dn, r_top, z, a_dn, t_up, r_bot, z, a_up)
    end
    return layer, r_vol_top, r_vol_bot
end

"""
Symmetric 1-zone generalized PROSPECT spectral budget used by
`prospect_budget(leaf, optics, surface)`. Uses `(t_ext, t_top_out)` on the
illuminated exterior face and symmetric `(t_hem_surf, t_top_out)` on internal
and opposite boundaries (`t_hem_surf = t_top_out * nr²`): when
`surface == ConeSurface(40, 1)`, `t_top_out == calctav(90°, nr) / nr²` so this
reduces algebraically to `prospect_spectral_budget` to `< 3 eps(FT)`, while
preserving smoothness across all `N >= 1` (`N = 1.0`, `1.4`, `3.0`).
"""
@inline function _asymmetric_top_prospect_budget(
    k_unit::FT,
    N::FT,
    t_ext::FT,
    t_top_out::FT,
    nr::FT,
) where {FT<:Real}
    t_hem_surf = t_top_out * (nr * nr)
    p_sym, _, _ = _single_bounded_plate(k_unit, t_hem_surf, t_top_out, t_hem_surf, t_top_out, 1)
    p_top, r_vol_top, _ = _single_bounded_plate(k_unit, t_ext, t_top_out, t_hem_surf, t_top_out, 1)

    rem = prospect_plate_stack(p_sym.r_top, p_sym.t_dn, p_sym.a1_dn, N - one(FT))
    coupling = one(FT) - rem.reflectance * p_top.r_bot
    inv_c = inv(coupling)

    transmission = p_top.t_dn * rem.transmission * inv_c
    r_vol_coupled = p_top.t_dn * rem.reflectance * p_top.t_up * inv_c
    reflectance_volume = r_vol_top + r_vol_coupled
    reflectance = (one(FT) - t_ext) + reflectance_volume
    absorptance = p_top.a1_dn + p_top.t_dn * (rem.absorptance + rem.reflectance * p_top.a1_up) * inv_c
    return (; transmission, reflectance, absorptance, reflectance_volume)
end

@inline function _couple_external_interfaces_to_medium(
    r_top_ext::FT, t_top_in::FT, r_top_int::FT, t_top_out::FT,
    r_med_dn::FT, t_med_dn::FT, a_med_dn::FT,
    r_med_up::FT, t_med_up::FT, a_med_up::FT,
    r_bot_int::FT, t_bot_out::FT,
) where {FT<:Real}
    z = zero(FT)
    I_top = _TwoSidedTwoZoneLayer{FT}(t_top_in, r_top_ext, z, z, t_top_out, r_top_int, z, z)
    M = _TwoSidedTwoZoneLayer{FT}(t_med_dn, r_med_dn, a_med_dn, z, t_med_up, r_med_up, a_med_up, z)
    # Bottom interface only needs upward internal reflection r_bot_int and exit t_bot_out for top incidence
    I_bot = _TwoSidedTwoZoneLayer{FT}(t_bot_out, r_bot_int, z, z, z, z, z, z)
    tot = _star_two_zone(_star_two_zone(I_top, M), I_bot)
    return (;
        transmission = tot.t_dn,
        reflectance = tot.r_top,
        absorptance = tot.a1_dn,
    )
end

@inline function _sublayer_chord_split(k1::FT, k2::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    k12 = k1 + k2
    att1 = prospect_internal_attenuation(k1)
    att12 = prospect_internal_attenuation(k12)
    a1 = att1.absorptance
    # Avoid 1 - 1 subtraction when k1 > 1: a12 - a1 == τ1 - τ12, while preserving
    # ForwardDiff tangents at k2 == 0 (where diff == 0 with positive Dual derivative)
    diff = prospect_primal(k1) > one(BFT) ? (att1.transmission - att12.transmission) : (att12.absorptance - a1)
    a2_rem = prospect_primal(diff) < zero(BFT) ? zero(FT) : diff
    return att12.transmission, a1, a2_rem
end

"""
    single_slab_dorsiventral_budget(k_pal, k_spo, nr, surf_ad, surf_ab;
                                    μ_ad = 1, μ_ab = 1, hemispherical = false)

Exact `N = 1` unbroken two-sublayer dorsiventral slab (`I_ad ⋆ S₂(k_pal, k_spo) ⋆ I_ab`)
without any artificial internal air gap or `(2 E₃(k/2))²` splitting error.
"""
@inline function single_slab_dorsiventral_budget(
    k_pal_in::FT1,
    k_spo_in::FT2,
    nr_in::FT3,
    surf_ad::AbstractLeafSurface{FTS1},
    surf_ab::AbstractLeafSurface{FTS2};
    μ_ad::Real = 1,
    μ_ab::Real = 1,
    hemispherical::Bool = false,
) where {FT1<:Real,FT2<:Real,FT3<:Real,FTS1<:Real,FTS2<:Real}
    FT = promote_type(FT1, FT2, FT3, FTS1, FTS2, typeof(μ_ad), typeof(μ_ab))
    k_pal = FT(k_pal_in)
    k_spo = FT(k_spo_in)
    nr = FT(nr_in)

    s_ad = surface_interface_optics(surf_ad, FT(μ_ad), nr; hemispherical)
    s_ab = surface_interface_optics(surf_ab, FT(μ_ab), nr; hemispherical)

    τ_12, a1_dn, a2_dn = _sublayer_chord_split(k_pal, k_spo)
    _, a2_up, a1_up = _sublayer_chord_split(k_spo, k_pal)
    a_tot = a1_dn + a2_dn
    two = one(FT) + one(FT)

    r_ad_int = s_ad.r_int_diffuse
    r_ab_int = s_ab.r_int_diffuse
    d = (s_ad.t_int_diffuse + s_ab.t_int_diffuse * r_ad_int) +
        r_ad_int * r_ab_int * a_tot * (two - a_tot)
    inv_d = inv(d)

    t_ad = s_ad.t_ext * τ_12 * s_ab.t_int_diffuse * inv_d
    r_vol_ad = s_ad.t_ext * (τ_12 * τ_12 * r_ab_int * s_ad.t_int_diffuse) * inv_d
    r_ad = s_ad.r_ext + r_vol_ad
    ap_ad = s_ad.t_ext * (a1_dn + τ_12 * r_ab_int * a1_up) * inv_d
    as_ad = s_ad.t_ext * (a2_dn + τ_12 * r_ab_int * a2_up) * inv_d

    t_ab = s_ab.t_ext * τ_12 * s_ad.t_int_diffuse * inv_d
    r_vol_ab = s_ab.t_ext * (τ_12 * τ_12 * r_ad_int * s_ab.t_int_diffuse) * inv_d
    r_ab = s_ab.r_ext + r_vol_ab
    ap_ab = s_ab.t_ext * (a1_up + τ_12 * r_ad_int * a1_dn) * inv_d
    as_ab = s_ab.t_ext * (a2_up + τ_12 * r_ad_int * a2_dn) * inv_d

    return (;
        adaxial = (;
            transmission = t_ad,
            reflectance = r_ad,
            absorptance = ap_ad + as_ad,
            absorptance_palisade = ap_ad,
            absorptance_spongy = as_ad,
            reflectance_volume = r_vol_ad,
        ),
        abaxial = (;
            transmission = t_ab,
            reflectance = r_ab,
            absorptance = ap_ab + as_ab,
            absorptance_palisade = ap_ab,
            absorptance_spongy = as_ab,
            reflectance_volume = r_vol_ab,
        ),
    )
end

"""
Scalar two-zone dorsiventral mesophyll stack for `N_pal >= 1, N_spo >= 1`
(`N = N_pal + N_spo >= 2`):
`P_top(k_pal, surf_ad, int) ⋆ P_int(k_pal)^{N_pal - 1} ⋆ P_int(k_spo)^{N_spo - 1} ⋆ P_bot(k_spo, int, surf_ab)`
where `k_pal` and `k_spo` are per-plate optical depths (`K_pal_zone / N_pal` and
`K_spo_zone / N_spo`).
"""
@inline function _two_zone_dorsiventral_scalar(
    k_pal_unit::FT,
    N_pal::FT,
    k_spo_unit::FT,
    N_spo::FT,
    nr::FT,
    surf_ad::AbstractLeafSurface,
    surf_ab::AbstractLeafSurface,
    μ_ad::FT,
    μ_ab::FT,
    hemispherical::Bool,
) where {FT<:Real}
    s_ad = surface_interface_optics(surf_ad, μ_ad, nr; hemispherical)
    s_ab = surface_interface_optics(surf_ab, μ_ab, nr; hemispherical)
    t12_int = calctav_smooth(FT(90), nr)
    t21_int = t12_int / (nr * nr)
    z = zero(FT)

    # Top boundary plate (zone 1 = palisade)
    P_top, r_vol_top_direct, _ = _single_bounded_plate(
        k_pal_unit, s_ad.t_ext, s_ad.t_int_diffuse, t12_int, t21_int, 1,
    )
    # Internal palisade plates (m_pal = N_pal - 1 >= 0)
    P_pal_int, _, _ = _single_bounded_plate(
        k_pal_unit, t12_int, t21_int, t12_int, t21_int, 1,
    )
    rem_pal = prospect_plate_stack(P_pal_int.r_top, P_pal_int.t_dn, P_pal_int.a1_dn, N_pal - one(FT))
    L_pal_rem = _TwoSidedTwoZoneLayer{FT}(
        rem_pal.transmission, rem_pal.reflectance, rem_pal.absorptance, z,
        rem_pal.transmission, rem_pal.reflectance, rem_pal.absorptance, z,
    )

    # Internal spongy plates (m_spo = N_spo - 1 >= 0)
    P_spo_int, _, _ = _single_bounded_plate(
        k_spo_unit, t12_int, t21_int, t12_int, t21_int, 2,
    )
    rem_spo = prospect_plate_stack(P_spo_int.r_top, P_spo_int.t_dn, P_spo_int.a2_dn, N_spo - one(FT))
    L_spo_rem = _TwoSidedTwoZoneLayer{FT}(
        rem_spo.transmission, rem_spo.reflectance, z, rem_spo.absorptance,
        rem_spo.transmission, rem_spo.reflectance, z, rem_spo.absorptance,
    )

    # Bottom boundary plate (zone 2 = spongy)
    P_bot, _, r_vol_bot_direct = _single_bounded_plate(
        k_spo_unit, t12_int, t21_int, s_ab.t_ext, s_ab.t_int_diffuse, 2,
    )

    # Assemble interior + opposite boundary first so volume reflectance is computed
    # directly without subtracting r_ext from total reflectance!
    L_below_top = _star_two_zone(L_pal_rem, _star_two_zone(L_spo_rem, P_bot))
    c_ad = inv(one(FT) - P_top.r_bot * L_below_top.r_top)
    r_vol_ad = r_vol_top_direct + P_top.t_dn * L_below_top.r_top * P_top.t_up * c_ad
    L_full = _star_two_zone(P_top, L_below_top)

    L_above_bot = _star_two_zone(P_top, _star_two_zone(L_pal_rem, L_spo_rem))
    c_ab = inv(one(FT) - L_above_bot.r_bot * P_bot.r_top)
    r_vol_ab = r_vol_bot_direct + P_bot.t_up * L_above_bot.r_bot * P_bot.t_dn * c_ab

    return (;
        adaxial = (;
            transmission = L_full.t_dn,
            reflectance = L_full.r_top,
            absorptance = L_full.a1_dn + L_full.a2_dn,
            absorptance_palisade = L_full.a1_dn,
            absorptance_spongy = L_full.a2_dn,
            reflectance_volume = r_vol_ad,
        ),
        abaxial = (;
            transmission = L_full.t_up,
            reflectance = L_full.r_bot,
            absorptance = L_full.a1_up + L_full.a2_up,
            absorptance_palisade = L_full.a1_up,
            absorptance_spongy = L_full.a2_up,
            reflectance_volume = r_vol_ab,
        ),
    )
end

"""
    DorsiventralLeafProperties{LP<:AbstractLeafProperties, LS<:AbstractLeafProperties, SAd<:AbstractLeafSurface, SAb<:AbstractLeafSurface} <: AbstractLeafProperties

Two-zone dorsiventral leaf representation separating **palisade** (adaxial) and
**spongy** (abaxial) mesophyll layers (`palisade`, `spongy`, each a
[`LeafProspectProProperties`](@ref) or [`LeafMultiPigmentProperties`](@ref) with
zone plate count `N >= 1`) bounded by distinct adaxial (`surface_adaxial`) and
abaxial (`surface_abaxial`) dielectric surfaces.

When `single_slab == true`, evaluates the unbroken `N = 1` two-sublayer operator
[`single_slab_dorsiventral_budget`](@ref) using each zone's total constituent
optical depth without internal mesophyll air gaps.
"""
struct DorsiventralLeafProperties{
    LP<:AbstractLeafProperties,
    LS<:AbstractLeafProperties,
    SAd<:AbstractLeafSurface,
    SAb<:AbstractLeafSurface,
} <: AbstractLeafProperties
    palisade::LP
    spongy::LS
    surface_adaxial::SAd
    surface_abaxial::SAb
    single_slab::Bool

    function DorsiventralLeafProperties(
        palisade::LP,
        spongy::LS;
        surface_adaxial::SAd = ConeSurface(one(_leaf_base_float(typeof(palisade.N))) * 40),
        surface_abaxial::SAb = ConeSurface(one(_leaf_base_float(typeof(spongy.N))) * 40),
        single_slab::Bool = false,
    ) where {
        LP<:Union{LeafProspectProProperties,LeafMultiPigmentProperties},
        LS<:Union{LeafProspectProProperties,LeafMultiPigmentProperties},
        SAd<:AbstractLeafSurface,
        SAb<:AbstractLeafSurface,
    }
        np = prospect_primal(palisade.N)
        ns = prospect_primal(spongy.N)
        if !single_slab
            (isfinite(np) && np >= one(np) && isfinite(ns) && ns >= one(ns)) ||
                throw(ArgumentError("Two-zone DorsiventralLeafProperties requires N_pal >= 1 and N_spo >= 1 (or single_slab = true)"))
        else
            (isfinite(np) && np > zero(np) && isfinite(ns) && ns > zero(ns)) ||
                throw(ArgumentError("single_slab DorsiventralLeafProperties requires finite positive N"))
        end
        return new{LP,LS,SAd,SAb}(palisade, spongy, surface_adaxial, surface_abaxial, single_slab)
    end
end

"""
    AsymmetricSurfaceLeafProperties(leaf, surface_adaxial, surface_abaxial)

One-zone homogeneous mesophyll (`leaf.N >= 2`, or `leaf.N == 1` unbroken slab;
also supports all `leaf.N >= 1` when `surface_adaxial` and `surface_abaxial`
share internal diffuse exit transmittance, such as `ConeSurface(α_ad, 1)` and
`ConeSurface(α_ab, 1)`) bounded by distinct adaxial and abaxial surfaces.
"""
struct AsymmetricSurfaceLeafProperties{
    L<:Union{LeafProspectProProperties,LeafMultiPigmentProperties},
    SAd<:AbstractLeafSurface,
    SAb<:AbstractLeafSurface,
} <: AbstractLeafProperties
    leaf::L
    surface_adaxial::SAd
    surface_abaxial::SAb

    function AsymmetricSurfaceLeafProperties(
        leaf::L,
        surface_adaxial::SAd,
        surface_abaxial::SAb,
    ) where {
        L<:Union{LeafProspectProProperties,LeafMultiPigmentProperties},
        SAd<:AbstractLeafSurface,
        SAb<:AbstractLeafSurface,
    }
        nv = prospect_primal(leaf.N)
        isfinite(nv) && nv >= one(nv) ||
            throw(ArgumentError("AsymmetricSurfaceLeafProperties requires finite N >= 1"))
        return new{L,SAd,SAb}(leaf, surface_adaxial, surface_abaxial)
    end
end

function AsymmetricSurfaceLeafProperties(
    leaf::Union{LeafProspectProProperties,LeafMultiPigmentProperties};
    surface_adaxial::AbstractLeafSurface = ConeSurface(one(_leaf_base_float(typeof(leaf.N))) * 40),
    surface_abaxial::AbstractLeafSurface = ConeSurface(one(_leaf_base_float(typeof(leaf.N))) * 40),
)
    return AsymmetricSurfaceLeafProperties(leaf, surface_adaxial, surface_abaxial)
end

function _validate_single_zone_leaf(leaf::LeafProspectProProperties)
    _validate_leaf_contents((
        leaf.Ccab, leaf.Ccar, leaf.Canth, leaf.Cbrown,
        leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
    ))
end

function _validate_single_zone_leaf(leaf::LeafMultiPigmentProperties)
    _validate_leaf_contents((
        leaf.Cchla, leaf.Cchlb, leaf.Ccar, leaf.Canth, leaf.Cbrown,
        leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
    ))
end

"""
    dorsiventral_budget(leaf::DorsiventralLeafProperties, optics;
                        μ_ad = 1, μ_ab = 1, hemispherical = false)

Compute paired adaxial-illuminated and abaxial-illuminated leaf optical budgets
`(; adaxial, abaxial)`, where each face budget is a `NamedTuple` containing
`(; transmission, reflectance, absorptance, absorptance_palisade, absorptance_spongy, reflectance_volume)`.

Guarantees:
1. Exact conservation `R + T + A_pal + A_spo ≈ 1` without clipping or rescaling.
2. Exact reduction to [`prospect_budget`](@ref) when palisade and spongy zones
   have equal per-plate optical depth and `ConeSurface(40, 1)` on both faces.
3. Exact hemispherical (`hemispherical = true`) Helmholtz reciprocity
   `adaxial.transmission ≈ abaxial.transmission` even for arbitrary pigment,
   structural (`N_pal != N_spo`), and cuticular (`surf_ad != surf_ab`) asymmetry.
"""
function dorsiventral_budget(
    leaf::DorsiventralLeafProperties,
    optics::AbstractProspectProperties;
    μ_ad::Real = one(_leaf_base_float(eltype(optics.nᵣ))),
    μ_ab::Real = one(_leaf_base_float(eltype(optics.nᵣ))),
    hemispherical::Bool = false,
)
    _validate_single_zone_leaf(leaf.palisade)
    _validate_single_zone_leaf(leaf.spongy)
    _validate_optics_arrays(optics)

    K_pal_vec = _total_optical_depth_vector(leaf.palisade, optics)
    K_spo_vec = _total_optical_depth_vector(leaf.spongy, optics)

    FTS1 = _surface_ft(leaf.surface_adaxial)
    FTS2 = _surface_ft(leaf.surface_abaxial)
    FT = float(promote_type(
        eltype(K_pal_vec), eltype(K_spo_vec), typeof(leaf.palisade.N), typeof(leaf.spongy.N),
        eltype(optics.nᵣ), FTS1, FTS2, typeof(μ_ad), typeof(μ_ab),
    ))
    BFT = _leaf_base_float(FT)
    if !hemispherical
        (zero(BFT) < prospect_primal(μ_ad) <= one(BFT) &&
         zero(BFT) < prospect_primal(μ_ab) <= one(BFT)) ||
            throw(ArgumentError("μ_ad and μ_ab must lie in (0, 1]"))
    end

    nλ = length(optics.nᵣ)
    T_ad = Vector{FT}(undef, nλ); R_ad = Vector{FT}(undef, nλ); A_ad = Vector{FT}(undef, nλ)
    Ap_ad = Vector{FT}(undef, nλ); As_ad = Vector{FT}(undef, nλ); Rv_ad = Vector{FT}(undef, nλ)
    T_ab = Vector{FT}(undef, nλ); R_ab = Vector{FT}(undef, nλ); A_ab = Vector{FT}(undef, nλ)
    Ap_ab = Vector{FT}(undef, nλ); As_ab = Vector{FT}(undef, nλ); Rv_ab = Vector{FT}(undef, nλ)

    N_pal = FT(leaf.palisade.N)
    N_spo = FT(leaf.spongy.N)
    μ_ad_ft = FT(μ_ad)
    μ_ab_ft = FT(μ_ab)

    @inbounds for idx in 1:nλ
        nr = FT(optics.nᵣ[idx])
        res = if leaf.single_slab
            single_slab_dorsiventral_budget(
                FT(K_pal_vec[idx]), FT(K_spo_vec[idx]), nr,
                leaf.surface_adaxial, leaf.surface_abaxial;
                μ_ad = μ_ad_ft, μ_ab = μ_ab_ft, hemispherical,
            )
        else
            _two_zone_dorsiventral_scalar(
                FT(K_pal_vec[idx]) / N_pal, N_pal,
                FT(K_spo_vec[idx]) / N_spo, N_spo,
                nr, leaf.surface_adaxial, leaf.surface_abaxial,
                μ_ad_ft, μ_ab_ft, hemispherical,
            )
        end
        ad, ab = res.adaxial, res.abaxial
        T_ad[idx] = ad.transmission; R_ad[idx] = ad.reflectance; A_ad[idx] = ad.absorptance
        Ap_ad[idx] = ad.absorptance_palisade; As_ad[idx] = ad.absorptance_spongy; Rv_ad[idx] = ad.reflectance_volume
        T_ab[idx] = ab.transmission; R_ab[idx] = ab.reflectance; A_ab[idx] = ab.absorptance
        Ap_ab[idx] = ab.absorptance_palisade; As_ab[idx] = ab.absorptance_spongy; Rv_ab[idx] = ab.reflectance_volume
    end

    return (;
        adaxial = (;
            transmission = T_ad, reflectance = R_ad, absorptance = A_ad,
            absorptance_palisade = Ap_ad, absorptance_spongy = As_ad, reflectance_volume = Rv_ad,
        ),
        abaxial = (;
            transmission = T_ab, reflectance = R_ab, absorptance = A_ab,
            absorptance_palisade = Ap_ab, absorptance_spongy = As_ab, reflectance_volume = Rv_ab,
        ),
    )
end

function dorsiventral_budget(
    asym::AsymmetricSurfaceLeafProperties,
    optics::AbstractProspectProperties;
    μ_ad::Real = one(_leaf_base_float(eltype(optics.nᵣ))),
    μ_ab::Real = one(_leaf_base_float(eltype(optics.nᵣ))),
    hemispherical::Bool = false,
)
    _validate_single_zone_leaf(asym.leaf)
    _validate_optics_arrays(optics)

    K_tot = _total_optical_depth_vector(asym.leaf, optics)
    FTS1 = _surface_ft(asym.surface_adaxial)
    FTS2 = _surface_ft(asym.surface_abaxial)
    FT = float(promote_type(
        eltype(K_tot), typeof(asym.leaf.N), eltype(optics.nᵣ),
        FTS1, FTS2, typeof(μ_ad), typeof(μ_ab),
    ))
    BFT = _leaf_base_float(FT)
    if !hemispherical
        (zero(BFT) < prospect_primal(μ_ad) <= one(BFT) &&
         zero(BFT) < prospect_primal(μ_ab) <= one(BFT)) ||
            throw(ArgumentError("μ_ad and μ_ab must lie in (0, 1]"))
    end

    N_ft = FT(asym.leaf.N)
    two = one(FT) + one(FT)
    half = inv(two)
    n_val = prospect_primal(N_ft)
    μ_ad_ft = FT(μ_ad)
    μ_ab_ft = FT(μ_ab)

    # When both faces are ConeSurface(..., 1) with exact unit n_scale (including zero Dual partials),
    # internal exit t21 == calctav(90°, nr)/nr^2 on both faces, so _asymmetric_top_prospect_budget is C^∞ smooth wrt N for ALL N >= 1!
    can_use_cone_1zone = (asym.surface_adaxial isa ConeSurface &&
                          asym.surface_abaxial isa ConeSurface &&
                          _same_val_and_partials(asym.surface_adaxial.n_scale, one(typeof(asym.surface_adaxial.n_scale))) &&
                          _same_val_and_partials(asym.surface_abaxial.n_scale, one(typeof(asym.surface_abaxial.n_scale))))
    if one(BFT) < n_val < BFT(2) && !can_use_cone_1zone
        throw(ArgumentError(
            "AsymmetricSurfaceLeafProperties with non-unit n_scale or RoughDielectricSurface requires N == 1 or N >= 2 (use DorsiventralLeafProperties with N_pal >= 1, N_spo >= 1, or symmetric prospect_budget for 1 < N < 2)",
        ))
    end

    nλ = length(optics.nᵣ)
    T_ad = Vector{FT}(undef, nλ); R_ad = Vector{FT}(undef, nλ); A_ad = Vector{FT}(undef, nλ); Rv_ad = Vector{FT}(undef, nλ)
    T_ab = Vector{FT}(undef, nλ); R_ab = Vector{FT}(undef, nλ); A_ab = Vector{FT}(undef, nλ); Rv_ab = Vector{FT}(undef, nλ)

    @inbounds for idx in 1:nλ
        nr = FT(optics.nᵣ[idx])
        K = FT(K_tot[idx])
        if can_use_cone_1zone
            s_ad = surface_interface_optics(asym.surface_adaxial, μ_ad_ft, nr; hemispherical)
            s_ab = surface_interface_optics(asym.surface_abaxial, μ_ab_ft, nr; hemispherical)
            k_unit = K / N_ft
            b_ad = _asymmetric_top_prospect_budget(k_unit, N_ft, s_ad.t_ext, s_ad.t_int_diffuse, nr)
            b_ab = _asymmetric_top_prospect_budget(k_unit, N_ft, s_ab.t_ext, s_ab.t_int_diffuse, nr)
            T_ad[idx] = b_ad.transmission; R_ad[idx] = b_ad.reflectance
            A_ad[idx] = b_ad.absorptance; Rv_ad[idx] = b_ad.reflectance_volume
            T_ab[idx] = b_ab.transmission; R_ab[idx] = b_ab.reflectance
            A_ab[idx] = b_ab.absorptance; Rv_ab[idx] = b_ab.reflectance_volume
        elseif n_val == one(BFT)
            res = single_slab_dorsiventral_budget(
                K * half, K * half, nr, asym.surface_adaxial, asym.surface_abaxial;
                μ_ad = μ_ad_ft, μ_ab = μ_ab_ft, hemispherical,
            )
            T_ad[idx] = res.adaxial.transmission; R_ad[idx] = res.adaxial.reflectance
            A_ad[idx] = res.adaxial.absorptance; Rv_ad[idx] = res.adaxial.reflectance_volume
            T_ab[idx] = res.abaxial.transmission; R_ab[idx] = res.abaxial.reflectance
            A_ab[idx] = res.abaxial.absorptance; Rv_ab[idx] = res.abaxial.reflectance_volume
        else
            # Exact semigroup split N_pal = N/2 >= 1, N_spo = N/2 >= 1, k_unit = K / N
            N_half = N_ft * half
            k_unit = K / N_ft
            res = _two_zone_dorsiventral_scalar(
                k_unit, N_half, k_unit, N_half, nr,
                asym.surface_adaxial, asym.surface_abaxial,
                μ_ad_ft, μ_ab_ft, hemispherical,
            )
            T_ad[idx] = res.adaxial.transmission; R_ad[idx] = res.adaxial.reflectance
            A_ad[idx] = res.adaxial.absorptance; Rv_ad[idx] = res.adaxial.reflectance_volume
            T_ab[idx] = res.abaxial.transmission; R_ab[idx] = res.abaxial.reflectance
            A_ab[idx] = res.abaxial.absorptance; Rv_ab[idx] = res.abaxial.reflectance_volume
        end
    end

    return (;
        adaxial = (; transmission = T_ad, reflectance = R_ad, absorptance = A_ad, reflectance_volume = Rv_ad),
        abaxial = (; transmission = T_ab, reflectance = R_ab, absorptance = A_ab, reflectance_volume = Rv_ab),
    )
end
