# `_leaf_base_float` and `prospect_primal` are defined in `src/types/canopy_types.jl`.

# Precomputed Gauss-Legendre nodes/weights and Chebyshev-Gauss-Lobatto tables (zero heap allocation)
const _SURF_GL5_X, _SURF_GL5_W = gausslegendre(5)
const _SURF_GL12_X, _SURF_GL12_W = gausslegendre(12)
const _SURF_GL16_X, _SURF_GL16_W = gausslegendre(16)
const _SURF_GL20_X, _SURF_GL20_W = gausslegendre(20)
const _SURF_GL24_X, _SURF_GL24_W = gausslegendre(24)

const _SURF_GL12_V = Tuple((_SURF_GL12_X .+ 1.0) .* 0.5)
const _SURF_GL12_WV = Tuple(_SURF_GL12_W .* 0.5)
const _SURF_GL12_MU = Tuple(sqrt.(_SURF_GL12_V))
const _SURF_GL12_STH = Tuple(sqrt.(max.(0.0, 1.0 .- _SURF_GL12_V)))
const _SURF_GL12_COS_DPHI = Tuple(cos.((_SURF_GL12_X .+ 1.0) .* (π * 0.5)))

const _SURF_GL20_S = Tuple((_SURF_GL20_X .+ 1.0) .* 0.5)
const _SURF_GL20_W_TUP = Tuple(_SURF_GL20_W)
const _SURF_GL24_COS_PHIH = Tuple(cos.((_SURF_GL24_X .+ 1.0) .* (π * 0.5)))
const _SURF_GL24_WP = Tuple(_SURF_GL24_W .* 0.5)

const _CHEB32_M = 32
const _CHEB32_C = ntuple(j -> (1.0 + cos(π * (j - 1) / (_CHEB32_M - 1))) * 0.5, _CHEB32_M)
const _CHEB32_S2 = ntuple(j -> max(0.0, 1.0 - _CHEB32_C[j]^2), _CHEB32_M)
const _CHEB32_CT = ntuple(_CHEB32_M) do j
    pj = (j == 1 || j == _CHEB32_M) ? 2.0 : 1.0
    ntuple(_CHEB32_M) do kp1
        k = kp1 - 1
        pk = (k == 0 || k == _CHEB32_M - 1) ? 2.0 : 1.0
        (2.0 / ((_CHEB32_M - 1) * pk * pj)) * cos(π * k * (j - 1) / (_CHEB32_M - 1))
    end
end

"""
    AbstractLeafSurface{FT<:Real}

Interface contract for the external air–cuticle–mesophyll dielectric boundary of
a leaf face. Concrete subtypes implement [`surface_interface_optics`](@ref).
"""
abstract type AbstractLeafSurface{FT<:Real} end

@inline _surface_ft(::AbstractLeafSurface{FT}) where {FT} = FT
@inline _same_val_and_partials(a::Real, b::Real) = a == b
@inline _same_val_and_partials(::ForwardDiff.Dual, ::Real) = false
@inline _same_val_and_partials(::Real, ::ForwardDiff.Dual) = false
@inline _same_val_and_partials(::ForwardDiff.Dual, ::ForwardDiff.Dual) = false
@inline function _same_val_and_partials(
    a::ForwardDiff.Dual{Tag,VA,N},
    b::ForwardDiff.Dual{Tag,VB,N},
) where {Tag,VA,VB,N}
    _same_val_and_partials(ForwardDiff.value(a), ForwardDiff.value(b)) || return false
    pa = ForwardDiff.partials(a).values
    pb = ForwardDiff.partials(b).values
    @inbounds for i in 1:N
        _same_val_and_partials(pa[i], pb[i]) || return false
    end
    return true
end

"""
    ConeSurface(α_deg = 40, n_scale = 1)

Isotropic incidence-cone dielectric boundary generalizing PROSPECT's Stern–Allen
surface (`α_deg = 40` reproduces legacy PROSPECT; `α_deg = 0` models collimated
nadir integrating-sphere illumination; `α_deg = 8` models 8°-port integrating
spheres). `n_scale` scales the aqueous refractive index `nr` at the external
epicuticular film (`n_eff = n_scale * nr`), while internal diffuse exit into the
common `nr` mesophyll obeys exact `1 / nr²` étendue invariance.
"""
struct ConeSurface{FT<:Real} <: AbstractLeafSurface{FT}
    α_deg::FT
    n_scale::FT

    function ConeSurface{FT}(α_deg::FT, n_scale::FT) where {FT<:Real}
        BFT = _leaf_base_float(FT)
        a_val = prospect_primal(α_deg)
        s_val = prospect_primal(n_scale)
        isfinite(a_val) && zero(BFT) <= a_val <= BFT(90) ||
            throw(ArgumentError("ConeSurface requires finite α_deg in [0, 90]"))
        isfinite(s_val) && s_val > zero(BFT) ||
            throw(ArgumentError("ConeSurface requires finite positive n_scale"))
        return new{FT}(α_deg, n_scale)
    end
end

ConeSurface(α_deg::FT, n_scale::FT = one(FT)) where {FT<:Real} =
    ConeSurface{float(FT)}(float(FT)(α_deg), float(FT)(n_scale))

function ConeSurface(α_deg::Real, n_scale::Real)
    FT = float(promote_type(typeof(α_deg), typeof(n_scale)))
    return ConeSurface{FT}(FT(α_deg), FT(n_scale))
end

function ConeSurface(; α_deg::Real = 40.0, n_scale::Real = one(typeof(α_deg)))
    return ConeSurface(α_deg, n_scale)
end

function ConeSurface{FT}(; α_deg = FT(40), n_scale = one(FT)) where {FT<:Real}
    return ConeSurface{FT}(FT(α_deg), FT(n_scale))
end

"""
    RoughDielectricSurface(σ = 0.15, n_scale = 1)

Rough Beckmann microfacet dielectric leaf surface with RMS facet slope
`σ >= 0.02` and effective epicuticular refractive-index factor `n_scale`
(`n_eff = n_scale * nr`). Uses exact `erfcx`-based Smith separable
masking-shadowing (`0 <= R_surf < 1` strictly) and thermodynamic `1 / nr²`
étendue scaling into the `nr` aqueous mesophyll.
"""
struct RoughDielectricSurface{FT<:Real} <: AbstractLeafSurface{FT}
    σ::FT
    n_scale::FT

    function RoughDielectricSurface{FT}(σ::FT, n_scale::FT) where {FT<:Real}
        BFT = _leaf_base_float(FT)
        sig_val = prospect_primal(σ)
        ns_val = prospect_primal(n_scale)
        isfinite(sig_val) && BFT(0.02) <= sig_val <= BFT(1.0) ||
            throw(ArgumentError("RoughDielectricSurface requires finite σ in [0.02, 1.0]"))
        isfinite(ns_val) && ns_val > zero(BFT) ||
            throw(ArgumentError("RoughDielectricSurface requires finite positive n_scale"))
        return new{FT}(σ, n_scale)
    end
end

RoughDielectricSurface(σ::FT, n_scale::FT = one(FT)) where {FT<:Real} =
    RoughDielectricSurface{float(FT)}(float(FT)(σ), float(FT)(n_scale))

function RoughDielectricSurface(σ::Real, n_scale::Real)
    FT = float(promote_type(typeof(σ), typeof(n_scale)))
    return RoughDielectricSurface{FT}(FT(σ), FT(n_scale))
end

function RoughDielectricSurface(; σ::Real = 0.15, n_scale::Real = one(typeof(σ)))
    return RoughDielectricSurface(σ, n_scale)
end

function RoughDielectricSurface{FT}(; σ = FT(0.15), n_scale = one(FT)) where {FT<:Real}
    return RoughDielectricSurface{FT}(FT(σ), FT(n_scale))
end

"""
    calctav_smooth(α_deg, nr)

Cancellation-free Stern (1964) / Allen (1973) cone-averaged dielectric surface
transmittance for `α_deg ∈ [0, 90]` and `nr > 1`, supporting `Float32`,
`Float64`, and `ForwardDiff.Dual` without precision promotion or `0/0`
singularities at `α_deg = 0` or `α_deg = 90`:
- For `sind(α_deg)^2 <= 0.04` (`α_deg <= 11.54°`, including `0°` and `8°`),
  evaluates the defining strictly-positive Fresnel integral `∫_0^1 T_F(sa2*ξ, nr) dξ`
  by 5-point Gauss–Legendre quadrature (`< 1-2 ulp` in `Float32`/`Float64`, `C^∞` at `0°`).
- Above `0.04` (including `40°` and `90°`), evaluates the exact Stern–Allen
  closed form with `_b1 = cosd(α_deg) * sqrt(nr^2 - sa2)`, avoiding `sqrt(0)`
  derivative blowup at `α_deg = 90°`.
"""
@inline function calctav_smooth(α_deg::FTA, nr::FTN) where {FTA<:Real,FTN<:Real}
    FT = float(promote_type(FTA, FTN))
    BFT = _leaf_base_float(FT)
    a_val = prospect_primal(α_deg)
    n_val = prospect_primal(nr)
    isfinite(a_val) && zero(BFT) <= a_val <= BFT(90) ||
        throw(ArgumentError("calctav_smooth requires α_deg in [0, 90]"))
    isfinite(n_val) && n_val > one(BFT) ||
        throw(ArgumentError("calctav_smooth requires refractive index nr > 1"))

    α = FT(α_deg)
    n = FT(nr)
    n2 = n * n
    sa = sind(α)
    sa2 = sa * sa
    two = one(FT) + one(FT)
    half = inv(two)

    if prospect_primal(sa2) <= BFT(0.04)
        acc = zero(FT)
        @inbounds for k in 1:5
            ξ = (BFT(_SURF_GL5_X[k]) + one(BFT)) * BFT(0.5)
            w = BFT(_SURF_GL5_W[k]) * BFT(0.5)
            u = sa2 * ξ
            ci = sqrt(one(FT) - u)
            g = sqrt(n2 - u)
            # Cancellation-free Fresnel transmittance: (1 - rs^2 + 1 - rp^2)/2
            ts_half = two * ci * g / ((ci + g) * (ci + g))
            tp_half = two * n2 * ci * g / ((n2 * ci + g) * (n2 * ci + g))
            acc += w * (ts_half + tp_half)
        end
        return acc
    end

    np1 = n + one(FT)
    a = (np1 * np1) * half
    a3 = a * a * a
    n4 = n2 * n2
    n6 = n4 * n2
    np = n2 + one(FT)
    np2 = np * np
    np3 = np2 * np
    nm1 = n2 - one(FT)
    nm2 = nm1 * nm1
    k_c = -nm2 / BFT(4)
    k2 = k_c * k_c

    _b1 = cosd(α) * sqrt(n2 - sa2)
    _b2 = sa2 - np * half
    b = _b1 - _b2
    b3 = b * b * b

    ts = (k2 / (BFT(6) * b3) + k_c / b - b * half) -
         (k2 / (BFT(6) * a3) + k_c / a - a * half)
    tp1 = -two * n2 * (b - a) / np2
    tp2 = -two * n2 * np * log(b / a) / nm2
    tp3 = n2 * (inv(b) - inv(a)) * half
    tp4 = BFT(16) * n4 * (n4 + one(FT)) *
          log((two * np * b - nm2) / (two * np * a - nm2)) / (np3 * nm2)
    tp5 = BFT(16) * n6 * (inv(two * np * b - nm2) - inv(two * np * a - nm2)) / np3
    tp = tp1 + tp2 + tp3 + tp4 + tp5
    return (ts + tp) / (two * sa2)
end

"""
    beckmann_D(tan2_theta_h, cos2_theta_h, σ)

Beckmann isotropic microfacet normal distribution `D(h)` evaluated directly from
`tan² θ_h = (H_x² + H_y²) / H_z²` and `cos² θ_h = H_z² / |H|²` without inverse
trigonometric functions (`C^∞` smooth at `θ_h = 0`).
"""
@inline function beckmann_D(tan2_th::FT1, cos2_th::FT2, σ::FT3) where {FT1<:Real,FT2<:Real,FT3<:Real}
    FT = float(promote_type(FT1, FT2, FT3))
    BFT = _leaf_base_float(FT)
    t_val = prospect_primal(tan2_th)
    c_val = prospect_primal(cos2_th)
    s_val = prospect_primal(σ)
    (isfinite(t_val) && t_val >= zero(BFT) &&
     isfinite(c_val) && c_val > zero(BFT) &&
     isfinite(s_val) && s_val > zero(BFT)) ||
        throw(ArgumentError("beckmann_D requires finite tan2_th >= 0, cos2_th > 0, and σ > 0"))
    sig2 = FT(σ) * FT(σ)
    c4 = FT(cos2_th) * FT(cos2_th)
    exponent = -FT(tan2_th) / sig2
    # Guard extreme exponent underflow cleanly while preserving AD
    prospect_primal(exponent) < BFT(-80) && return zero(FT)
    return exp(exponent) / (FT(π) * sig2 * c4)
end

"""
    smith_beckmann_G1(μ, σ)

Exact Beckmann 1D Smith masking-shadowing factor `G₁(μ; σ) = 1 / (1 + Λ(a))`
evaluated via `SpecialFunctions.erfcx`. Guards the flat `C^∞` nadir neighborhood
`(1 - μ²) σ² <= μ² / a_max²` before forming `sqrt(1 - μ²)`, eliminating `0/0`
and `sqrt(0)` `ForwardDiff.Dual` singularities at `μ = 1`.
"""
@inline function smith_beckmann_G1(μ_in::FT1, σ_in::FT2) where {FT1<:Real,FT2<:Real}
    FT = float(promote_type(FT1, FT2))
    BFT = _leaf_base_float(FT)
    μ = FT(μ_in)
    σ = FT(σ_in)
    μ_val = prospect_primal(μ)
    s_val = prospect_primal(σ)
    isfinite(s_val) && s_val > zero(BFT) ||
        throw(ArgumentError("smith_beckmann_G1 requires finite σ > 0"))
    isfinite(μ_val) && μ_val <= one(BFT) + BFT(16) * eps(BFT) ||
        throw(ArgumentError("smith_beckmann_G1 requires finite μ <= 1"))
    μ_val <= zero(BFT) && return zero(FT)

    # For a >= a_cut (5.0 in Float32, 8.5 in Float64), Λ(a) < 1e-13 (F32) or 1e-34 (F64),
    # so 1 / (1 + Λ) == 1 to 0 ulps and all derivatives vanish exponentially.
    a_cut2 = BFT === Float32 ? BFT(25) : BFT(72.25)
    u_sin2 = one(FT) - μ * μ
    if prospect_primal(u_sin2) * (s_val * s_val) * a_cut2 <= μ_val * μ_val
        return one(FT)
    end

    a = μ / (σ * sqrt(u_sin2))
    two = one(FT) + one(FT)
    Λ = exp(-a * a) * (inv(a * sqrt(FT(π))) - oftype(a, erfcx(a))) / two
    return inv(one(FT) + Λ)
end

"""
    fresnel_unpolarized_cos(c_i, s_i2, n_eff)

Algebraic unpolarized Fresnel reflectance `(r_s² + r_p²) / 2` from local
incidence cosine `c_i = ω_i ⋅ h >= 0` and `s_i² = sin² α_inc` without `acos`/`asin`.
"""
@inline function fresnel_unpolarized_cos(c_i::FT, s_i2::FT, n_eff::FT) where {FT<:Real}
    n2 = n_eff * n_eff
    g = sqrt(n2 - s_i2)
    rs = (c_i - g) / (c_i + g)
    rp = (n2 * c_i - g) / (n2 * c_i + g)
    two = one(FT) + one(FT)
    return (rs * rs + rp * rp) / two
end

@inline function _microfacet_directional_reflectance_core(
    μi::FT,
    sθi::FT,
    n_eff::FT,
    σ::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    g1i = smith_beckmann_G1(μi, σ)
    two = one(FT) + one(FT)

    if iszero(prospect_primal(sθi)) && _same_val_and_partials(sθi, zero(FT))
        t_max = inv(μi)
        arg = -(t_max * t_max) / (σ * σ)
        u_min = prospect_primal(arg) < BFT(-80) ? zero(FT) : exp(arg)
        du_span = one(FT) - u_min
        acc_nadir = zero(FT)
        @inbounds for ku in 1:20
            s = BFT(_SURF_GL20_S[ku])
            u = max(FT(floatmin(BFT)), one(FT) - du_span * (s * s))
            wu = BFT(_SURF_GL20_W_TUP[ku]) * s * du_span
            t = σ * sqrt(-log(u))
            μh = inv(sqrt(one(FT) + t * t))
            ci = μi * μh
            μo = clamp(two * ci * μh - μi, FT(eps(BFT)), one(FT))
            si2 = max(zero(FT), one(FT) - ci * ci)
            F = fresnel_unpolarized_cos(ci, si2, n_eff)
            g1o = smith_beckmann_G1(μo, σ)
            acc_nadir += wu * g1o * F * (ci / (μi * μh))
        end
        return g1i * acc_nadir
    end

    acc = zero(FT)
    @inbounds for kp in 1:24
        wp = BFT(_SURF_GL24_WP[kp])
        b = sθi * BFT(_SURF_GL24_COS_PHIH[kp])
        disc = sqrt(b * b + μi * μi)
        # Exact upper horizon root t_max > 0 of μ_o = (-μ_i t^2 + 2 b t + μ_i)/(1 + t^2) = 0,
        # paired with u = 1 - (1 - u_min) s^2 so sqrt(-log u) ~ s is C^∞ smooth at u = 1!
        t_max = prospect_primal(b) >= zero(BFT) ? (b + disc) / μi : μi / (disc - b)
        arg = -(t_max * t_max) / (σ * σ)
        u_min = prospect_primal(arg) < BFT(-80) ? zero(FT) : exp(arg)
        du_span = one(FT) - u_min

        acc_u = zero(FT)
        for ku in 1:20
            s = BFT(_SURF_GL20_S[ku])
            u = max(FT(floatmin(BFT)), one(FT) - du_span * (s * s))
            wu = BFT(_SURF_GL20_W_TUP[ku]) * s * du_span
            t = σ * sqrt(-log(u))
            μh = inv(sqrt(one(FT) + t * t))
            sθh = t * μh
            ci = b * sθh + μi * μh
            μo = clamp(two * ci * μh - μi, FT(eps(BFT)), one(FT))
            si2 = max(zero(FT), one(FT) - ci * ci)
            F = fresnel_unpolarized_cos(ci, si2, n_eff)
            g1o = smith_beckmann_G1(μo, σ)
            acc_u += wu * g1o * F * (ci / (μi * μh))
        end
        acc += wp * acc_u
    end
    return g1i * acc
end

@inline function _accum_cheb32_moments!(
    moms::Vector{FT},
    μi::FT,
    sθi::FT,
    σ::FT,
    weight_outer::FT,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    g1i = smith_beckmann_G1(μi, σ) * weight_outer
    two = one(FT) + one(FT)

    if iszero(prospect_primal(sθi)) && _same_val_and_partials(sθi, zero(FT))
        t_max = inv(μi)
        arg = -(t_max * t_max) / (σ * σ)
        u_min = prospect_primal(arg) < BFT(-80) ? zero(FT) : exp(arg)
        du_span = one(FT) - u_min
        @inbounds for ku in 1:20
            s = BFT(_SURF_GL20_S[ku])
            u = max(FT(floatmin(BFT)), one(FT) - du_span * (s * s))
            wu = BFT(_SURF_GL20_W_TUP[ku]) * s * du_span
            t = σ * sqrt(-log(u))
            μh = inv(sqrt(one(FT) + t * t))
            ci = μi * μh
            μo = clamp(two * ci * μh - μi, FT(eps(BFT)), one(FT))
            g1o = smith_beckmann_G1(μo, σ)
            wq = g1i * wu * g1o * (ci / (μi * μh))
            x = two * ci - one(FT)
            t0 = one(FT)
            t1 = x
            moms[1] += wq * t0
            moms[2] += wq * t1
            for k in 3:_CHEB32_M
                tk = two * x * t1 - t0
                moms[k] += wq * tk
                t0 = t1
                t1 = tk
            end
        end
        return moms
    end

    @inbounds for kp in 1:24
        wp = BFT(_SURF_GL24_WP[kp])
        b = sθi * BFT(_SURF_GL24_COS_PHIH[kp])
        disc = sqrt(b * b + μi * μi)
        t_max = prospect_primal(b) >= zero(BFT) ? (b + disc) / μi : μi / (disc - b)
        arg = -(t_max * t_max) / (σ * σ)
        u_min = prospect_primal(arg) < BFT(-80) ? zero(FT) : exp(arg)
        du_span = one(FT) - u_min

        for ku in 1:20
            s = BFT(_SURF_GL20_S[ku])
            u = max(FT(floatmin(BFT)), one(FT) - du_span * (s * s))
            wu = BFT(_SURF_GL20_W_TUP[ku]) * s * du_span
            t = σ * sqrt(-log(u))
            μh = inv(sqrt(one(FT) + t * t))
            sθh = t * μh
            ci = b * sθh + μi * μh
            μo = clamp(two * ci * μh - μi, FT(eps(BFT)), one(FT))
            g1o = smith_beckmann_G1(μo, σ)
            wq = g1i * wp * wu * g1o * (ci / (μi * μh))
            x = two * ci - one(FT)
            t0 = one(FT)
            t1 = x
            moms[1] += wq * t0
            moms[2] += wq * t1
            for k in 3:_CHEB32_M
                tk = two * x * t1 - t0
                moms[k] += wq * tk
                t0 = t1
                t1 = tk
            end
        end
    end
    return moms
end

@inline function _moments_to_cheb32_weights(moms::Vector{FT}) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    return ntuple(Val(_CHEB32_M)) do j
        cj = _CHEB32_CT[j]
        acc = zero(FT)
        @inbounds for k in 1:_CHEB32_M
            acc = muladd(BFT(cj[k]), moms[k], acc)
        end
        acc
    end
end

@inline function _microfacet_dir_cheb32_weights(::Type{FT}, μi::FT, sθi::FT, σ::FT) where {FT<:Real}
    moms = zeros(FT, _CHEB32_M)
    _accum_cheb32_moments!(moms, μi, sθi, σ, one(FT))
    return _moments_to_cheb32_weights(moms)
end

@inline function _eval_cheb32_fresnel(w::NTuple{_CHEB32_M,FT}, n_eff::FT) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    acc = zero(FT)
    @inbounds for j in 1:_CHEB32_M
        Fj = fresnel_unpolarized_cos(FT(BFT(_CHEB32_C[j])), FT(BFT(_CHEB32_S2[j])), n_eff)
        acc = muladd(w[j], Fj, acc)
    end
    return acc
end

struct _PreparedRoughSurfaceCheb{FT<:Real}
    σ::FT
    n_scale::FT
    hemispherical::Bool
    fallback_direct::Bool
    w_hem::NTuple{_CHEB32_M,FT}
    w_dir::NTuple{_CHEB32_M,FT}
end

@inline _prepare_surface_cheb(::Type{FT}, ::ConeSurface, ::FT, ::Bool) where {FT<:Real} = nothing

function _prepare_surface_cheb(
    ::Type{FT},
    surface::RoughDielectricSurface,
    μ_in::FT,
    hemispherical::Bool,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    σ = FT(surface.σ)
    n_scale = FT(surface.n_scale)

    moms_hem = zeros(FT, _CHEB32_M)
    @inbounds for k in 1:12
        wv = FT(BFT(_SURF_GL12_WV[k]))
        μ = FT(BFT(_SURF_GL12_MU[k]))
        sθ = FT(BFT(_SURF_GL12_STH[k]))
        _accum_cheb32_moments!(moms_hem, μ, sθ, σ, wv)
    end
    w_hem = _moments_to_cheb32_weights(moms_hem)

    if hemispherical
        return _PreparedRoughSurfaceCheb{FT}(σ, n_scale, true, false, w_hem, w_hem)
    end

    u_sin2 = one(FT) - μ_in * μ_in
    if prospect_primal(u_sin2) <= zero(BFT)
        if u_sin2 isa ForwardDiff.Dual && !_same_val_and_partials(u_sin2, zero(FT))
            return _PreparedRoughSurfaceCheb{FT}(σ, n_scale, false, true, w_hem, w_hem)
        end
        w_dir = _microfacet_dir_cheb32_weights(FT, one(FT), zero(FT), σ)
        return _PreparedRoughSurfaceCheb{FT}(σ, n_scale, false, false, w_hem, w_dir)
    end

    w_dir = _microfacet_dir_cheb32_weights(FT, μ_in, sqrt(u_sin2), σ)
    return _PreparedRoughSurfaceCheb{FT}(σ, n_scale, false, false, w_hem, w_dir)
end

@inline function _eval_prepared_surface(
    surface::ConeSurface,
    ::Nothing,
    μ_in::FT,
    nr::FT,
    n_surf::FT,
    hemispherical::Bool,
) where {FT<:Real}
    return surface_interface_optics(surface, μ_in, nr; hemispherical, n_surface = n_surf)
end

@inline function _eval_prepared_surface(
    surface::RoughDielectricSurface,
    prep::_PreparedRoughSurfaceCheb{FT},
    μ_in::FT,
    nr::FT,
    n_surf::FT,
    hemispherical::Bool,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    n_eff = prep.n_scale * n_surf
    (isfinite(prospect_primal(nr)) && prospect_primal(nr) > one(BFT)) ||
        throw(ArgumentError("Mesophyll refractive index nr must be finite and exceed 1"))
    (isfinite(prospect_primal(n_eff)) && prospect_primal(n_eff) > one(BFT)) ||
        throw(ArgumentError("Effective surface refractive index n_eff = n_scale * n_surface must be finite and exceed 1"))

    if prep.fallback_direct || prospect_primal(n_eff) < BFT(1.08)
        return surface_interface_optics(surface, μ_in, nr; hemispherical, n_surface = n_surf)
    end

    r_hem_ext = zero(FT)
    r_dir = zero(FT)
    wh = prep.w_hem
    wd = prep.w_dir
    @inbounds for j in 1:_CHEB32_M
        Fj = fresnel_unpolarized_cos(FT(BFT(_CHEB32_C[j])), FT(BFT(_CHEB32_S2[j])), n_eff)
        r_hem_ext = muladd(wh[j], Fj, r_hem_ext)
        if !hemispherical
            r_dir = muladd(wd[j], Fj, r_dir)
        end
    end
    r_ext = hemispherical ? r_hem_ext : r_dir
    t_hem_ext = one(FT) - r_hem_ext
    t_ext = one(FT) - r_ext
    t_int_diffuse = t_hem_ext / (nr * nr)
    r_int_diffuse = one(FT) - t_int_diffuse
    return (; r_ext, t_ext, r_int_diffuse, t_int_diffuse, t_hem_ext, n_eff)
end

"""
    microfacet_directional_reflectance(μ_in, n_eff, σ)

Single-bounce directional-hemispherical reflectance `R_{surf,1}(μ_in, n_eff, σ)`
of a Beckmann–Smith dielectric surface (`0 < R_{surf,1} < 1` strictly), evaluated
via exact Beckmann slope-cdf `(u, ϕ_h)` change of variables (`D(h) μ_h dω_h = du dϕ_h / 2π`),
which eliminates peak-width aliasing across all `σ >= 0.02`. Preserves exact
`ForwardDiff` tangents even at nadir `μ_in = 1`.
"""
@inline function microfacet_directional_reflectance(
    μ_in::FT1,
    n_eff_in::FT2,
    σ_in::FT3,
) where {FT1<:Real,FT2<:Real,FT3<:Real}
    FT = float(promote_type(FT1, FT2, FT3))
    BFT = _leaf_base_float(FT)
    μi = FT(μ_in)
    n_eff = FT(n_eff_in)
    σ = FT(σ_in)

    μ_p = prospect_primal(μi)
    n_p = prospect_primal(n_eff)
    sig_p = prospect_primal(σ)
    (isfinite(μ_p) && zero(BFT) < μ_p <= one(BFT)) ||
        throw(ArgumentError("microfacet_directional_reflectance requires finite μ_in in (0, 1]"))
    (isfinite(n_p) && n_p > one(BFT)) ||
        throw(ArgumentError("microfacet_directional_reflectance requires finite n_eff > 1"))
    (isfinite(sig_p) && sig_p > zero(BFT)) ||
        throw(ArgumentError("microfacet_directional_reflectance requires finite σ > 0"))

    u_sin2 = one(FT) - μi * μi
    if prospect_primal(u_sin2) <= zero(BFT)
        r0 = _microfacet_directional_reflectance_core(one(FT), zero(FT), n_eff, σ)
        if u_sin2 isa ForwardDiff.Dual && !_same_val_and_partials(u_sin2, zero(FT))
            # Smooth dependence on u_sin2 = 1 - μ_in^2 at nadir by azimuthal symmetry;
            # evaluate r_δ on primal (n_p, sig_p) so Dual perturbations in (n_eff, σ) are not divided by δ^2!
            δ = BFT === Float32 ? min(BFT(5e-2), sig_p * BFT(0.25)) : min(BFT(1e-3), sig_p * BFT(0.05))
            μ_δ = sqrt(one(BFT) - δ * δ)
            r_δ = _microfacet_directional_reflectance_core(μ_δ, δ, n_p, sig_p)
            slope_u = (r_δ - prospect_primal(r0)) / (δ * δ)
            return r0 + FT(slope_u) * u_sin2
        end
        return r0
    end

    sθi = sqrt(u_sin2)
    return _microfacet_directional_reflectance_core(μi, sθi, n_eff, σ)
end

"""
    microfacet_hemispherical_reflectance(n_eff, σ)

Hemispherical-hemispherical single-bounce reflectance
`2 ∫_0^1 R_{surf,1}(μ, n_eff, σ) μ dμ = ∫_0^1 R_{surf,1}(√v, n_eff, σ) dv`.
"""
@inline function microfacet_hemispherical_reflectance(n_eff_in::FT1, σ_in::FT2) where {FT1<:Real,FT2<:Real}
    FT = float(promote_type(FT1, FT2))
    BFT = _leaf_base_float(FT)
    n_eff = FT(n_eff_in)
    σ = FT(σ_in)
    n_p = prospect_primal(n_eff)
    sig_p = prospect_primal(σ)
    (isfinite(n_p) && n_p > one(BFT)) ||
        throw(ArgumentError("microfacet_hemispherical_reflectance requires finite n_eff > 1"))
    (isfinite(sig_p) && sig_p > zero(BFT)) ||
        throw(ArgumentError("microfacet_hemispherical_reflectance requires finite σ > 0"))

    acc = zero(FT)
    @inbounds for k in 1:12
        w = BFT(_SURF_GL12_WV[k])
        μ = BFT(_SURF_GL12_MU[k])
        sθ = BFT(_SURF_GL12_STH[k])
        acc += w * _microfacet_directional_reflectance_core(FT(μ), FT(sθ), n_eff, σ)
    end
    return acc
end

"""
    surface_interface_optics(surface::AbstractLeafSurface, μ_in, nr;
                             hemispherical = false, n_surface = nr)

Evaluate external entry and internal diffuse exit probabilities across a leaf's
external dielectric boundary:
- `r_ext`, `t_ext = 1 - r_ext`: external reflectance and entry transmittance
  (directional at `μ_in` when `hemispherical == false`, or hemispherical `2π -> 2π`
  when `hemispherical == true`), evaluated at outer cuticular refractive index
  `n_eff = surface.n_scale * n_surface`.
- `t_hem_ext`: hemispherical (`2π -> 2π`) external transmittance across `surface`.
- `t_int_diffuse = t_hem_ext / nr²`, `r_int_diffuse = 1 - t_int_diffuse`:
  diffuse internal-to-air exit transmittance and internal reflectance from the
  aqueous mesophyll (`nr > 1`) into air (`n = 1`).

Physical thin-cuticle approximation when `n_surface != nr` (e.g., Stuckens et al.
2009 epicuticular wax `n_wax(λ)` via [`load_wax_refractive_index`](@ref)): across
air (`n = 1`) — lossless thin wax cuticle (`n_eff`) — aqueous mesophyll (`nr`),
smooth internal wax–mesophyll Fresnel reflection is `(n_eff - nr)² / (n_eff + nr)²`
(`≈ 0.03–0.05%` in VNIR and `≤ 0.5%` at the `2450 nm` SWIR `nr ≈ 1.27` minimum),
with zero wax→mesophyll total internal reflection for externally entering rays
since `sin θ_wax ≤ 1 / n_eff < nr / n_eff` for all `nr > 1`. Neglecting that `< 0.5%`
secondary internal wax–mesophyll bounce, external entry is governed by the high-contrast
air–wax interface (`n_eff`), while diffuse internal exit from the aqueous mesophyll
(`nr`) into air (`n = 1`) obeys exact isothermal detailed-balance étendue invariance
`t_int_diffuse = t_hem_ext(n_eff) / nr² ∈ (0, 1)` strictly.
"""
@inline function surface_interface_optics(
    surface::ConeSurface{FTS},
    μ_in::FTM,
    nr_in::FTN;
    hemispherical::Bool = false,
    n_surface::FTW = nr_in,
) where {FTS<:Real,FTM<:Real,FTN<:Real,FTW<:Real}
    FT = float(promote_type(FTS, FTM, FTN, FTW))
    BFT = _leaf_base_float(FT)
    nr = FT(nr_in)
    n_surf = FT(n_surface)
    n_eff = FT(surface.n_scale) * n_surf
    (isfinite(prospect_primal(nr)) && prospect_primal(nr) > one(BFT)) ||
        throw(ArgumentError("Mesophyll refractive index nr must be finite and exceed 1"))
    (isfinite(prospect_primal(n_eff)) && prospect_primal(n_eff) > one(BFT)) ||
        throw(ArgumentError("Effective surface refractive index n_eff = n_scale * n_surface must be finite and exceed 1"))

    t_hem_ext = calctav_smooth(FT(90), n_eff)
    t_ext = hemispherical ? t_hem_ext : calctav_smooth(FT(surface.α_deg), n_eff)
    r_ext = one(FT) - t_ext
    t_int_diffuse = t_hem_ext / (nr * nr)
    r_int_diffuse = one(FT) - t_int_diffuse
    return (; r_ext, t_ext, r_int_diffuse, t_int_diffuse, t_hem_ext, n_eff)
end

@inline function surface_interface_optics(
    surface::RoughDielectricSurface{FTS},
    μ_in::FTM,
    nr_in::FTN;
    hemispherical::Bool = false,
    n_surface::FTW = nr_in,
) where {FTS<:Real,FTM<:Real,FTN<:Real,FTW<:Real}
    FT = float(promote_type(FTS, FTM, FTN, FTW))
    BFT = _leaf_base_float(FT)
    nr = FT(nr_in)
    n_surf = FT(n_surface)
    σ = FT(surface.σ)
    n_eff = FT(surface.n_scale) * n_surf
    (isfinite(prospect_primal(nr)) && prospect_primal(nr) > one(BFT)) ||
        throw(ArgumentError("Mesophyll refractive index nr must be finite and exceed 1"))
    (isfinite(prospect_primal(n_eff)) && prospect_primal(n_eff) > one(BFT)) ||
        throw(ArgumentError("Effective surface refractive index n_eff = n_scale * n_surface must be finite and exceed 1"))

    r_hem_ext = microfacet_hemispherical_reflectance(n_eff, σ)
    t_hem_ext = one(FT) - r_hem_ext
    if hemispherical
        r_ext = r_hem_ext
    else
        μ_val = prospect_primal(μ_in)
        isfinite(μ_val) && zero(BFT) < μ_val <= one(BFT) ||
            throw(ArgumentError("RoughDielectricSurface requires μ_in in (0, 1]"))
        r_ext = microfacet_directional_reflectance(FT(μ_in), n_eff, σ)
    end
    t_ext = one(FT) - r_ext
    t_int_diffuse = t_hem_ext / (nr * nr)
    r_int_diffuse = one(FT) - t_int_diffuse
    return (; r_ext, t_ext, r_int_diffuse, t_int_diffuse, t_hem_ext, n_eff)
end
