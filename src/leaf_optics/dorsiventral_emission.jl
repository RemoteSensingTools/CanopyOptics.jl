"""
    DorsiventralPlateEscape(adaxial_depth, abaxial_depth, adaxial_transmission, abaxial_transmission[, interface_transmission = 1])

Two-zone dorsiventral internal emission escape kernel with finite nonnegative zone
emission optical depths `adaxial_depth ≥ 0` (palisade / zone 1) and `abaxial_depth ≥ 0`
(spongy / zone 2), internal-to-outside diffuse exit transmission probabilities
`0 < adaxial_transmission ≤ 1` and `0 < abaxial_transmission ≤ 1`, and internal
palisade–spongy boundary transmission `0 < interface_transmission ≤ 1` (default `1`
for an index-matched mesophyll boundary, or `< 1` for an intercellular air-space
interface).

Within-zone single-pass survival is `exp(-b_i)` with isotropic half-space initial
emission (`1/2` toward either boundary of the birth zone). Reabsorption is terminal
with no secondary re-emission, wavelength redistribution, physiological yield, or
directional radiance assumption. Photon conservation `front + back + reabsorbed ≈ 1`
holds to floating-point roundoff without clipping or post-hoc rescaling across normal
and subnormal interface transmissions.

This `isbits` scalar object owns five floating-point values, allocates zero memory
on evaluation, and runs identically on CPU and CUDA (`Float32`/`Float64`).
"""
struct DorsiventralPlateEscape{FT<:AbstractFloat}
    adaxial_depth::FT
    abaxial_depth::FT
    adaxial_transmission::FT
    abaxial_transmission::FT
    interface_transmission::FT
    function DorsiventralPlateEscape{FT}(
        b_ad::FT,
        b_ab::FT,
        t_ad::FT,
        t_ab::FT,
        t_int::FT = one(FT),
    ) where {FT<:AbstractFloat}
        isfinite(b_ad) && b_ad >= zero(FT) && isfinite(b_ab) && b_ab >= zero(FT) ||
            throw(ArgumentError("zone emission depths must be finite and nonnegative"))
        isfinite(t_ad) &&
        zero(FT) < t_ad <= one(FT) &&
        isfinite(t_ab) &&
        zero(FT) < t_ab <= one(FT) ||
            throw(ArgumentError("outer interface transmissions must be finite and in (0, 1]"))
        isfinite(t_int) && zero(FT) < t_int <= one(FT) ||
            throw(ArgumentError("internal interface transmission must be finite and in (0, 1]"))
        new{FT}(b_ad, b_ab, t_ad, t_ab, t_int)
    end
end

DorsiventralPlateEscape(
    b_ad::FT,
    b_ab::FT,
    t_ad::FT,
    t_ab::FT,
    t_int::FT = one(FT),
) where {FT<:AbstractFloat} = DorsiventralPlateEscape{FT}(b_ad, b_ab, t_ad, t_ab, t_int)

@inline function _normalize_zone_birth_amplitudes(
    forward::FT,
    backward::FT,
    zone_weight::FT,
) where {FT<:AbstractFloat}
    isfinite(forward) && isfinite(backward) && min(forward, backward) >= zero(FT) ||
        throw(ArgumentError("zone birth amplitudes must be finite and nonnegative"))
    scale = max(forward, backward)
    if iszero(scale)
        iszero(zone_weight) && return (one(FT), zero(FT))
        throw(ArgumentError("active zone with positive weight must have a positive birth amplitude"))
    end
    f, b = forward / scale, backward / scale
    total = f + b
    return (f / total, b / total)
end

"""
    DorsiventralBirthProfile(
        adaxial_weight, abaxial_weight,
        adaxial_excitation_depth, adaxial_forward, adaxial_backward,
        abaxial_excitation_depth, abaxial_forward, abaxial_backward,
    )
    DorsiventralBirthProfile(adaxial_weight, adaxial_profile, abaxial_weight, abaxial_profile)
    DorsiventralBirthProfile(adaxial_weight, abaxial_weight)

Two-zone normalized internal photon birth-depth distribution for [`DorsiventralPlateEscape`](@ref).

Local dimensionless depth `ζ_i ∈ [0, 1]` runs **anatomically adaxial → abaxial in both zones**:
- Zone 1 (adaxial / palisade): `ζ_1 = 0` at the outer adaxial epidermis and `ζ_1 = 1` at the
  internal palisade–spongy boundary, with conditional density proportional to
  `adaxial_forward * exp(-a_1 * ζ_1) + adaxial_backward * exp(-a_1 * (1 - ζ_1))`.
- Zone 2 (abaxial / spongy): `ζ_2 = 0` at the internal palisade–spongy boundary and `ζ_2 = 1`
  at the outer abaxial epidermis, with conditional density proportional to
  `abaxial_forward * exp(-a_2 * ζ_2) + abaxial_backward * exp(-a_2 * (1 - ζ_2))`.

Reversing the anatomical leaf orientation (`adaxial ↔ abaxial`) therefore maps
`(w_1, a_1, f_1, g_1) ↔ (w_2, a_2, g_2, f_2)`, swapping `front` and `back` escape
while leaving `reabsorbed` invariant.

The two-argument `DorsiventralBirthProfile(adaxial_weight, abaxial_weight)` constructor
creates a piecewise-uniform birth profile (`a_1 = a_2 = 0`) with normalized zone weights.
Passing [`UniformBirthProfile`](@ref) directly to [`leaf_emission_fates`](@ref) on a
[`DorsiventralPlateEscape`](@ref) weights both zones equally (`1/2, 1/2`).
"""
struct DorsiventralBirthProfile{FT<:AbstractFloat} <: AbstractLeafBirthProfile
    adaxial_weight::FT
    abaxial_weight::FT
    adaxial_excitation_depth::FT
    adaxial_forward_fraction::FT
    adaxial_backward_fraction::FT
    abaxial_excitation_depth::FT
    abaxial_forward_fraction::FT
    abaxial_backward_fraction::FT
    function DorsiventralBirthProfile{FT}(
        w_ad::FT,
        w_ab::FT,
        a_ad::FT,
        f_ad::FT,
        b_ad::FT,
        a_ab::FT,
        f_ab::FT,
        b_ab::FT,
    ) where {FT<:AbstractFloat}
        isfinite(w_ad) && isfinite(w_ab) && min(w_ad, w_ab) >= zero(FT) ||
            throw(ArgumentError("zone weights must be finite and nonnegative"))
        w_scale = max(w_ad, w_ab)
        w_scale > zero(FT) ||
            throw(ArgumentError("at least one zone birth weight must be positive"))
        w1_s, w2_s = w_ad / w_scale, w_ab / w_scale
        w_tot = w1_s + w2_s
        w1, w2 = w1_s / w_tot, w2_s / w_tot

        isfinite(a_ad) && a_ad >= zero(FT) && isfinite(a_ab) && a_ab >= zero(FT) ||
            throw(ArgumentError("zone excitation depths must be finite and nonnegative"))

        f1, g1 = _normalize_zone_birth_amplitudes(f_ad, b_ad, w1)
        f2, g2 = _normalize_zone_birth_amplitudes(f_ab, b_ab, w2)
        new{FT}(w1, w2, a_ad, f1, g1, a_ab, f2, g2)
    end
end

DorsiventralBirthProfile(
    w_ad::FT,
    w_ab::FT,
    a_ad::FT,
    f_ad::FT,
    b_ad::FT,
    a_ab::FT,
    f_ab::FT,
    b_ab::FT,
) where {FT<:AbstractFloat} =
    DorsiventralBirthProfile{FT}(w_ad, w_ab, a_ad, f_ad, b_ad, a_ab, f_ab, b_ab)

DorsiventralBirthProfile{FT}(w_ad::FT, w_ab::FT) where {FT<:AbstractFloat} =
    DorsiventralBirthProfile{FT}(
        w_ad,
        w_ab,
        zero(FT),
        one(FT),
        zero(FT),
        zero(FT),
        one(FT),
        zero(FT),
    )

DorsiventralBirthProfile(w_ad::FT, w_ab::FT) where {FT<:AbstractFloat} =
    DorsiventralBirthProfile{FT}(w_ad, w_ab)

DorsiventralBirthProfile(
    w_ad::FT,
    adaxial::ExponentialBirthProfile{FT},
    w_ab::FT,
    abaxial::ExponentialBirthProfile{FT},
) where {FT<:AbstractFloat} = DorsiventralBirthProfile{FT}(
    w_ad,
    w_ab,
    adaxial.excitation_depth,
    adaxial.forward_fraction,
    adaxial.backward_fraction,
    abaxial.excitation_depth,
    abaxial.forward_fraction,
    abaxial.backward_fraction,
)

DorsiventralBirthProfile(
    w_ad::FT,
    ::UniformBirthProfile,
    w_ab::FT,
    ::UniformBirthProfile,
) where {FT<:AbstractFloat} = DorsiventralBirthProfile(w_ad, w_ab)

# Single-zone capped reflection/termination step for a photon entering zone i from the
# internal interface toward its outer boundary (transmission t, optical depth b).
# Forms conditional termination branches (e_branch, a_branch) before excursion scaling
# so subnormal outer transmissions t ∈ (0, floatmin(FT)) retain full precision.
@inline function _zone_interface_step(b::FT, t::FT) where {FT}
    r = one(FT) - t
    τ = exp(-b)
    ℓ = -expm1(-b)
    R = r * τ * τ
    A = ℓ * (one(FT) + r * τ)
    E = τ * t
    T = E + A
    e_branch, a_branch = if τ >= inv(FT(2))
        a_scaled = expm1(b) * (one(FT) + r * τ)
        denom = t + a_scaled
        (t / denom, a_scaled / denom)
    elseif t < floatmin(FT)
        up = inv(eps(FT))
        e_up = τ * (t * up)
        denom = e_up + A * up
        (e_up / denom, (A * up) / denom)
    else
        (E / T, A / T)
    end
    return (; τ, ℓ, r, R, T, e_branch, a_branch)
end

@inline _scaled_pair_product(x::FT, y::FT, s_med::FT, s_max::FT) where {FT} =
    (max(x, y) / s_max) * (min(x, y) / s_med)

@inline function dorsiventral_plate_fates(
    escape::DorsiventralPlateEscape{FT},
    w_1::FT,
    w_2::FT,
    a_1::FT,
    f_1::FT,
    g_1::FT,
    a_2::FT,
    f_2::FT,
    g_2::FT,
) where {FT}
    b_1, b_2 = escape.adaxial_depth, escape.abaxial_depth
    t_1, t_2 = escape.adaxial_transmission, escape.abaxial_transmission
    t_int = escape.interface_transmission
    r_int = one(FT) - t_int
    half = inv(FT(2))

    # Zone 1 (ζ_1 = 0 adaxial outer face, ζ_1 = 1 internal interface)
    same_1, opp_1 = emission_overlaps(a_1, b_1)
    F_10 = (f_1 * same_1 + g_1 * opp_1) * half
    I_10 = (g_1 * same_1 + f_1 * opp_1) * half
    surv_1 = (same_1 + opp_1) * half
    L_10 = emission_first_loss(a_1, b_1, surv_1)
    z1 = _zone_interface_step(b_1, t_1)

    # Zone 2 (ζ_2 = 0 internal interface, ζ_2 = 1 abaxial outer face)
    same_2, opp_2 = emission_overlaps(a_2, b_2)
    I_20 = (f_2 * same_2 + g_2 * opp_2) * half
    B_20 = (g_2 * same_2 + f_2 * opp_2) * half
    surv_2 = (same_2 + opp_2) * half
    L_20 = emission_first_loss(a_2, b_2, surv_2)
    z2 = _zone_interface_step(b_2, t_2)

    # Pre-interface first-excursion arrivals at the internal boundary from each zone
    J_1 = w_1 * (I_10 + z1.r * z1.τ * F_10)
    J_2 = w_2 * (I_20 + z2.r * z2.τ * B_20)

    # Sort (T_1, T_2, t_int) symmetrically via 3 minmax operations so s_max >= s_med >= s_min > 0
    s_a, s_b = minmax(z1.T, z2.T)
    s_b, s_max = minmax(s_b, t_int)
    _, s_med = minmax(s_a, s_b)

    q_12 = _scaled_pair_product(z1.T, z2.T, s_med, s_max)
    q_t1 = _scaled_pair_product(t_int, z1.T, s_med, s_max)
    q_t2 = _scaled_pair_product(t_int, z2.T, s_med, s_max)

    # Subtraction-free termination probabilities in zone 1 vs zone 2; denominators ∈ [1, 3]
    w11 = r_int * q_12 + q_t1 * z2.R
    w12 = q_t2
    denom_1 = w11 + w12

    w21 = q_t1
    w22 = r_int * q_12 + q_t2 * z1.R
    denom_2 = w21 + w22

    W_1 = J_1 * (w11 / denom_1) + J_2 * (w21 / denom_2)
    W_2 = J_1 * (w12 / denom_1) + J_2 * (w22 / denom_2)

    front = w_1 * t_1 * F_10 + W_1 * z1.e_branch
    back = w_2 * t_2 * B_20 + W_2 * z2.e_branch
    reabsorbed =
        (w_1 * (L_10 + z1.r * z1.ℓ * F_10) + w_2 * (L_20 + z2.r * z2.ℓ * B_20)) +
        (W_1 * z1.a_branch + W_2 * z2.a_branch)
    return (; front, back, reabsorbed)
end

@inline function leaf_emission_fates(
    escape::DorsiventralPlateEscape{FT},
    ::UniformBirthProfile,
) where {FT}
    half = inv(FT(2))
    return dorsiventral_plate_fates(
        escape,
        half,
        half,
        zero(FT),
        one(FT),
        zero(FT),
        zero(FT),
        one(FT),
        zero(FT),
    )
end

@inline leaf_emission_fates(
    escape::DorsiventralPlateEscape{FT},
    birth::DorsiventralBirthProfile{FT},
) where {FT} = dorsiventral_plate_fates(
    escape,
    birth.adaxial_weight,
    birth.abaxial_weight,
    birth.adaxial_excitation_depth,
    birth.adaxial_forward_fraction,
    birth.adaxial_backward_fraction,
    birth.abaxial_excitation_depth,
    birth.abaxial_forward_fraction,
    birth.abaxial_backward_fraction,
)

@inline function _single_zone_plate_attenuation(k::FT) where {FT<:AbstractFloat}
    att = prospect_internal_attenuation(k)
    tau, loss = att.transmission, att.absorptance
    isfinite(tau) && zero(FT) < tau <= one(FT) || throw(
        DomainError(
            k,
            "unresolved zone traversal transmission; use higher precision or an independently resolved effective depth",
        ),
    )
    isfinite(loss) && zero(FT) <= loss <= one(FT) ||
        throw(DomainError(k, "invalid direct zone traversal absorption"))
    b = loss <= inv(FT(2)) ? -log1p(-loss) : -log(tau)
    return tau, loss, b
end

"""
    dorsiventral_plate_optics(
        k_adaxial, k_abaxial, nr;
        surface_adaxial = ConeSurface(40, 1),
        surface_abaxial = ConeSurface(40, 1),
        n_surface_adaxial = nr,
        n_surface_abaxial = nr,
        hemispherical = false,
        internal_air_gap = false,
    )

Host preparation helper mapping finite nonnegative palisade (`k_adaxial`) and spongy
(`k_abaxial`) material optical depths and aqueous mesophyll refractive index `nr > 1`
into effective exponential emission depths, outer entry/exit transmissions, internal
palisade–spongy boundary transmission, and a prepared [`DorsiventralPlateEscape`](@ref).

Outer diffuse internal exit transmissions `t_exit = t_hem_ext(n_eff, σ) / nr²` obey
isothermal detailed-balance étendue conservation (`n_eff = n_scale * n_surface_face`),
supporting both [`ConeSurface`](@ref) and [`RoughDielectricSurface`](@ref) as well as
epicuticular wax refractive indices (`n_surface_adaxial`, `n_surface_abaxial`).
`hemispherical` selects hemispherical (`90°`) vs `40°`-cone external **entry**
transmission (`adaxial_entry_transmission`, `abaxial_entry_transmission`) for excitation
transport; diffuse internal **exit** transmission is always hemispherical.

For the internal boundary between zones:
- `internal_air_gap = false` (default): `interface_transmission = 1` (unbroken
  index-matched mesophyll boundary, exact two-sublayer slab limit).
- `internal_air_gap = true`: symmetric diffuse mesophyll (`nr`) → intercellular air space
  (`n = 1`) → mesophyll (`nr`) transmission `t_gap = t90 / (nr² * (2 - t90))` where
  `t90 = calctav(90, nr)`, accounting for multiple internal reflections within the air gap.
"""
function dorsiventral_plate_optics(
    k_adaxial::FT,
    k_abaxial::FT,
    nr::FT;
    surface_adaxial::AbstractLeafSurface = ConeSurface{FT}(FT(40), one(FT)),
    surface_abaxial::AbstractLeafSurface = ConeSurface{FT}(FT(40), one(FT)),
    n_surface_adaxial::FT = nr,
    n_surface_abaxial::FT = nr,
    hemispherical::Bool = false,
    internal_air_gap::Bool = false,
) where {FT<:AbstractFloat}
    isfinite(k_adaxial) && k_adaxial >= zero(FT) && isfinite(k_abaxial) && k_abaxial >= zero(FT) ||
        throw(ArgumentError("zone material optical depths must be finite and nonnegative"))
    isfinite(nr) && nr > one(FT) ||
        throw(ArgumentError("mesophyll refractive index must be finite and greater than one"))
    isfinite(n_surface_adaxial) &&
    n_surface_adaxial > one(FT) &&
    isfinite(n_surface_abaxial) &&
    n_surface_abaxial > one(FT) ||
        throw(ArgumentError("surface refractive indices must be finite and greater than one"))

    tau_ad, loss_ad, b_ad = _single_zone_plate_attenuation(k_adaxial)
    tau_ab, loss_ab, b_ab = _single_zone_plate_attenuation(k_abaxial)

    surf_ad = surface_interface_optics(
        surface_adaxial,
        one(FT),
        nr;
        hemispherical = hemispherical,
        n_surface = n_surface_adaxial,
    )
    surf_ab = surface_interface_optics(
        surface_abaxial,
        one(FT),
        nr;
        hemispherical = hemispherical,
        n_surface = n_surface_abaxial,
    )
    entry_ad = FT(surf_ad.t_ext)
    exit_ad = FT(surf_ad.t_int_diffuse)
    entry_ab = FT(surf_ab.t_ext)
    exit_ab = FT(surf_ab.t_int_diffuse)

    t_int = if internal_air_gap
        t90 = FT(calctav_smooth(FT(90), nr))
        t90 / (nr * nr * (FT(2) - t90))
    else
        one(FT)
    end

    isfinite(entry_ad) &&
    zero(FT) < entry_ad <= one(FT) &&
    isfinite(exit_ad) &&
    zero(FT) < exit_ad <= one(FT) &&
    isfinite(entry_ab) &&
    zero(FT) < entry_ab <= one(FT) &&
    isfinite(exit_ab) &&
    zero(FT) < exit_ab <= one(FT) &&
    isfinite(t_int) &&
    zero(FT) < t_int <= one(FT) || throw(
        DomainError(
            nr,
            "invalid computed interface transmission; nearly matched indices are unsupported",
        ),
    )

    escape = DorsiventralPlateEscape{FT}(b_ad, b_ab, exit_ad, exit_ab, t_int)
    return (;
        adaxial_traversal_transmission = tau_ad,
        adaxial_traversal_absorptance = loss_ad,
        adaxial_effective_depth = b_ad,
        abaxial_traversal_transmission = tau_ab,
        abaxial_traversal_absorptance = loss_ab,
        abaxial_effective_depth = b_ab,
        adaxial_entry_transmission = entry_ad,
        adaxial_exit_transmission = exit_ad,
        abaxial_entry_transmission = entry_ab,
        abaxial_exit_transmission = exit_ab,
        interface_transmission = t_int,
        escape,
    )
end
