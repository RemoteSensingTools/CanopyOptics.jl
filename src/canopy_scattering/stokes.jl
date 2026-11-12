"""
    _validate_npol(npol) -> Int

Validate the requested Stokes-vector length for canopy Z matrices.

`npol = 1` is the scalar Stokes-I default.  `npol = 3` and `npol = 4`
match vSmartMOM's `Stokes_IQU` and `Stokes_IQUV` kernels.  A two-component
`I,Q` basis is intentionally not accepted because rotations of the Stokes
reference plane mix `Q` and `U`.
"""
function _validate_npol(npol::Integer)
    n = Int(npol)
    n in (1, 3, 4) || throw(ArgumentError("npol must be 1, 3, or 4"))
    return n
end

@inline _stokes_index(iμ::Integer, ipol::Integer, npol::Integer) =
    (Int(iμ) - 1) * Int(npol) + Int(ipol)

"""
    _azimuthal_kernel(si, sj, m, dϕ)

Fourier kernel for a Stokes matrix element at azimuth separation `dϕ`.

The convention matches vSmartMOM's surface Fourier decomposition: I/Q to I/Q
and U/V to U/V elements use cosine moments; cross-family elements use sine
moments.
"""
@inline function _azimuthal_kernel(si::Int, sj::Int, m::Int, dϕ::FT) where {FT<:Real}
    row_is_iq = si <= 2
    col_is_iq = sj <= 2
    return row_is_iq == col_is_iq ? cos(m * dϕ) : sin(m * dϕ)
end

"""
    _expand_intensity_Z(Z, npol)

Expand a scalar Stokes-I canopy Z matrix into a block Stokes matrix.

Diffuse Lambertian canopy scattering is unpolarized in this model, so only the
I→I entry of each stream-to-stream block is populated.  `npol = 1` returns
the input unchanged.
"""
function _expand_intensity_Z(Z::AbstractMatrix, npol::Integer)
    n = _validate_npol(npol)
    n == 1 && return Z

    nμ_out, nμ_in = size(Z)
    ZN = zeros(eltype(Z), n * nμ_out, n * nμ_in)
    @inbounds for j in 1:nμ_in, i in 1:nμ_out
        ZN[_stokes_index(i, 1, n), _stokes_index(j, 1, n)] = Z[i, j]
    end
    return ZN
end

function _expand_intensity_Z(Z::AbstractArray{<:Any,3}, npol::Integer)
    n = _validate_npol(npol)
    n == 1 && return Z

    nμ_out, nμ_in, nm = size(Z)
    ZN = zeros(eltype(Z), n * nμ_out, n * nμ_in, nm)
    @inbounds for k in 1:nm, j in 1:nμ_in, i in 1:nμ_out
        ZN[_stokes_index(i, 1, n), _stokes_index(j, 1, n), k] = Z[i, j, k]
    end
    return ZN
end

_expand_intensity_Z_pair(Zpp, Zmp, npol::Integer) =
    (_expand_intensity_Z(Zpp, npol), _expand_intensity_Z(Zmp, npol))

function _fresnel_mueller_real(nᵣ, α, npol::Integer)
    n = _validate_npol(npol)
    r_s, r_p = fresnel_components(nᵣ, α)
    FT = promote_type(typeof(r_s), typeof(r_p))
    M = zeros(FT, n, n)

    Rs = r_s^2
    Rp = r_p^2
    half_sum = (Rs + Rp) / 2
    half_diff = (Rs - Rp) / 2
    rsp = r_s * r_p

    M[1, 1] = half_sum
    if n >= 3
        M[1, 2] = half_diff
        M[2, 1] = half_diff
        M[2, 2] = half_sum
        M[3, 3] = rsp
        if n == 4
            M[4, 4] = rsp
        end
    end
    return M
end

@inline function _cartesian(Ω::dirVector_μ{FT}) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    u = one(FT) - Ω.μ * Ω.μ
    sθ = prospect_primal(u) <= zero(BFT) ? zero(FT) : sqrt(u)
    return (sθ * cos(Ω.ϕ), sθ * sin(Ω.ϕ), Ω.μ)
end

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

@inline function _cross3(a, b)
    return (a[2] * b[3] - a[3] * b[2],
            a[3] * b[1] - a[1] * b[3],
            a[1] * b[2] - a[2] * b[1])
end

@inline _norm3(a) = sqrt(_dot3(a, a))
@inline _real_value(x) = x
@inline _real_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)

@inline function _scale3(a, c)
    return (a[1] * c, a[2] * c, a[3] * c)
end

"""
    compute_reflection_mueller(mod, Ωin, Ωout, LD, npol)

Specular leaf bidirectional scattering Mueller matrix for one direction pair,
where `Ωin` points toward the illumination source (`k_in = -_cartesian(Ωin)`)
and `Ωout` points along the reflected ray (`k_out = +_cartesian(Ωout)`).

The scalar `I → I` element reduces identically to [`compute_reflection`](@ref).
For `npol = 3` or `4`, the local Fresnel Mueller matrix (`+Q` perpendicular to the
local reflection plane, `-Q` parallel) is rotated algebraically from the incoming
canopy meridian-plane triad `(e_s,in, e_p,in, k_in)` into the local reflection
plane and out to the outgoing canopy meridian-plane triad `(e_s,out, e_p,out, k_out)`,
matching the vSmartMOM Fourier `cos(m dϕ)` / `sin(m dϕ)` parity convention and
satisfying exact 4×4 Hovenier reciprocity `M(Ω2, Ω1) = P M(Ω1, Ω2)ᵀ P` (`P = diag(1, 1, -1, 1)`).
"""
function compute_reflection_mueller(mod::SpecularCanopyScattering,
                                    Ωin::dirVector_μ{FTΩ},
                                    Ωout::dirVector_μ{FTΩ},
                                    LD::AbstractLeafDistribution,
                                    npol::Integer = 1) where {FTΩ<:Real}
    n = _validate_npol(npol)
    (; nᵣ, κ) = mod
    BFTΩ = _leaf_base_float(FTΩ)
    Ωleaf, αstar = getSpecularΩ(Ωin, Ωout)
    u_L = one(FTΩ) - Ωleaf.μ * Ωleaf.μ
    s_L = prospect_primal(u_L) <= zero(BFTΩ) ? zero(FTΩ) : sqrt(u_L)
    θstar = atan(s_L, abs(Ωleaf.μ))

    M_loc = _fresnel_mueller_real(nᵣ, αstar, n)
    ZFT = eltype(M_loc)
    BFT = _leaf_base_float(ZFT)
    prefactor = ZFT(1 / 8) * ZFT(pdf(LD.LD, 2θstar / π)) *
                ZFT(LD.scaling) * K(κ, αstar)

    if n == 1
        return prefactor .* M_loc
    end

    # Smooth C^∞ right-handed canopy meridian polarization triads (e_s, e_p = k × e_s, k):
    # Incoming ray propagates along k_in = -_cartesian(Ωin); outgoing along k_out = +_cartesian(Ωout).
    μi = clamp(ZFT(Ωin.μ), -one(ZFT), one(ZFT))
    μo = clamp(ZFT(Ωout.μ), -one(ZFT), one(ZFT))
    ϕi = ZFT(Ωin.ϕ)
    ϕo = ZFT(Ωout.ϕ)

    ui = one(ZFT) - μi * μi
    uo = one(ZFT) - μo * μo
    si = prospect_primal(ui) <= zero(BFT) ? zero(ZFT) : sqrt(ui)
    so = prospect_primal(uo) <= zero(BFT) ? zero(ZFT) : sqrt(uo)

    ci_ϕ, si_ϕ = cos(ϕi), sin(ϕi)
    co_ϕ, so_ϕ = cos(ϕo), sin(ϕo)

    es_in_x, es_in_y, es_in_z = -si_ϕ, ci_ϕ, zero(ZFT)
    ep_in_x, ep_in_y, ep_in_z = μi * ci_ϕ, μi * si_ϕ, -si

    es_out_x, es_out_y, es_out_z = -so_ϕ, co_ϕ, zero(ZFT)
    ep_out_x, ep_out_y, ep_out_z = -μo * co_ϕ, -μo * so_ϕ, so

    kin_x, kin_y, kin_z = -si * ci_ϕ, -si * si_ϕ, -μi
    nleaf = _cartesian(dirVector_μ(ZFT(Ωleaf.μ), ZFT(Ωleaf.ϕ)))

    vs_x = kin_y * nleaf[3] - kin_z * nleaf[2]
    vs_y = kin_z * nleaf[1] - kin_x * nleaf[3]
    vs_z = kin_x * nleaf[2] - kin_y * nleaf[1]
    vs_norm2 = vs_x * vs_x + vs_y * vs_y + vs_z * vs_z

    us_x, us_y, us_z = if prospect_primal(vs_norm2) > eps(BFT)^2
        inv_vs = inv(sqrt(vs_norm2))
        (vs_x * inv_vs, vs_y * inv_vs, vs_z * inv_vs)
    else
        # Collinear normal-retroreflection / grazing limit: R_s == R_p (or K(κ, π/2)) makes M_loc invariant
        (es_in_x, es_in_y, es_in_z)
    end

    c1 = es_in_x * us_x + es_in_y * us_y + es_in_z * us_z
    s1 = ep_in_x * us_x + ep_in_y * us_y + ep_in_z * us_z
    c2 = es_out_x * us_x + es_out_y * us_y + es_out_z * us_z
    s2 = ep_out_x * us_x + ep_out_y * us_y + ep_out_z * us_z

    two = one(ZFT) + one(ZFT)
    cos2α1 = c1 * c1 - s1 * s1
    sin2α1 = two * c1 * s1
    cos2α2 = c2 * c2 - s2 * s2
    sin2α2 = two * c2 * s2

    f_plus = M_loc[1, 1]
    f_minus = M_loc[1, 2]
    f_sp = M_loc[3, 3]

    M = zeros(ZFT, n, n)
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
