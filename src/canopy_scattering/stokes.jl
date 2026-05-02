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

function _stokes_rotation_matrix(φ::FTφ, npol::Integer, ::Type{FT}) where {FTφ<:Real,FT}
    n = _validate_npol(npol)
    M = zeros(FT, n, n)
    M[1, 1] = one(FT)
    n == 1 && return M

    c2 = FT(cos(2 * φ))
    s2 = FT(sin(2 * φ))
    M[2, 2] = c2
    M[2, 3] = -s2
    M[3, 2] = s2
    M[3, 3] = c2
    if n == 4
        M[4, 4] = one(FT)
    end
    return M
end

@inline function _cartesian(Ω::dirVector_μ{FT}) where {FT<:Real}
    sθ = sqrt(max(zero(FT), one(FT) - Ω.μ^2))
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

function _plane_normal(a, b, ::Type{FT}) where {FT}
    n = _cross3(a, b)
    mag = _norm3(n)
    _real_value(mag) <= 1e-14 && return nothing
    return _scale3(n, inv(mag))
end

function _signed_plane_angle(axis, from_plane, to_plane, ::Type{FT}) where {FT}
    from_plane === nothing && return zero(FT)
    to_plane === nothing && return zero(FT)

    c = clamp(FT(_dot3(from_plane, to_plane)), -one(FT), one(FT))
    s = FT(_dot3(axis, _cross3(from_plane, to_plane)))
    return atan(s, c)
end

"""
    _specular_stokes_rotation_angles(Ωin, Ωout, Ωleaf, FT)

Return rotations from the scattering plane to the local leaf incidence plane
and from the local reflection plane back to the scattering plane.
"""
function _specular_stokes_rotation_angles(Ωin::dirVector_μ, Ωout::dirVector_μ,
                                          Ωleaf::dirVector_μ, ::Type{FT}) where {FT}
    kin = _cartesian(Ωin)
    kout = _cartesian(Ωout)
    nleaf = _cartesian(Ωleaf)

    scattering_plane = _plane_normal(kin, kout, FT)
    incidence_plane = _plane_normal(kin, nleaf, FT)
    reflection_plane = _plane_normal(kout, nleaf, FT)

    α1 = _signed_plane_angle(kin, scattering_plane, incidence_plane, FT)
    α2 = _signed_plane_angle(kout, scattering_plane, reflection_plane, FT)
    return α1, α2
end

"""
    compute_reflection_mueller(mod, Ωin, Ωout, LD, npol)

Specular leaf bidirectional scattering Mueller matrix for one direction pair.

The scalar I→I element reduces to [`compute_reflection`](@ref).  For
`npol = 3` or `4`, the Fresnel Mueller matrix is rotated between the
scattering plane and the local leaf incidence/reflection planes before the
Nilson-Kuusk roughness and leaf-normal PDF factors are applied.
"""
function compute_reflection_mueller(mod::SpecularCanopyScattering,
                                    Ωin::dirVector_μ{FTΩ},
                                    Ωout::dirVector_μ{FTΩ},
                                    LD::AbstractLeafDistribution,
                                    npol::Integer = 1) where {FTΩ<:Real}
    n = _validate_npol(npol)
    (; nᵣ, κ) = mod
    Ωleaf, αstar = getSpecularΩ(Ωin, Ωout)
    θstar = acos(abs(Ωleaf.μ))

    M = _fresnel_mueller_real(nᵣ, αstar, n)
    ZFT = eltype(M)
    prefactor = ZFT(1 / 8) * ZFT(pdf(LD.LD, 2θstar / π)) *
                ZFT(LD.scaling) * K(κ, αstar)

    if n > 1
        α1, α2 = _specular_stokes_rotation_angles(Ωin, Ωout, Ωleaf, ZFT)
        L1 = _stokes_rotation_matrix(-α1, n, ZFT)
        L2 = _stokes_rotation_matrix(α2, n, ZFT)
        M = L2 * M * L1
    end

    return prefactor .* M
end
