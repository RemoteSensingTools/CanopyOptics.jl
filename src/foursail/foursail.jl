"""
    FourSAILGeometry(LAD; sza_deg, vza_deg, raa_deg, quadrature=CanopyQuadrature())

Geometry and leaf-angle coefficients for a 4SAIL canopy solve.

The expensive geometry/LAD integration is wavelength independent and is done
once here. [`foursail!`](@ref) then evaluates one wavelength per kernel work
item from leaf reflectance, leaf transmittance, soil reflectance, and LAI.
"""
struct FourSAILGeometry{FT<:Real}
    μs::FT
    μo::FT
    dϕ::FT
    dso::FT
    ks::FT
    ko::FT
    sob::FT
    sof::FT
    sdb::FT
    sdf::FT
    dob::FT
    dof::FT
    ddb::FT
    ddf::FT
end

"""
    FourSAILGeometrySet(LAD; sza_deg, vza_deg, raa_deg, quadrature=CanopyQuadrature())
    FourSAILGeometrySet(geometries)

Structure-of-arrays storage for many [`FourSAILGeometry`](@ref) values. This is
the preferred input for phase-space/LUT solves because the kernel can evaluate
all `(wavelength, angle)` pairs in one launch.
"""
struct FourSAILGeometrySet{A<:AbstractVector}
    μs::A
    μo::A
    dϕ::A
    dso::A
    ks::A
    ko::A
    sob::A
    sof::A
    sdb::A
    sdf::A
    dob::A
    dof::A
    ddb::A
    ddf::A
end

"""
    FourSAILResult(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)

Primary 4SAIL reflectance/transmittance factors:

- `rddt`: bi-hemispherical reflectance over canopy plus soil
- `rsdt`: directional-hemispherical reflectance for solar illumination
- `rdot`: hemispherical-directional reflectance for the view direction
- `rsot`: bidirectional reflectance factor
- `taus`: total solar-direction canopy transmission (`tsd + tss`)
- `rsos`: direct single-scattering canopy contribution
- `rsost`: direct/direct contribution (`rsos` plus directly transmitted soil)
- `rsodt`: diffuse/multiple-scattering contribution, so `rsot = rsost + rsodt`
"""
struct FourSAILResult{A}
    rddt::A
    rsdt::A
    rdot::A
    rsot::A
    taus::A
    rsos::A
    rsost::A
    rsodt::A
end

function FourSAILResult(rddt::Number, rsdt::Number, rdot::Number, rsot::Number,
                        taus::Number, rsos::Number, rsost::Number,
                        rsodt::Number)
    vals = promote(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)
    return FourSAILResult{typeof(vals[1])}(vals...)
end

_foursail_zero_like(x::Number) = zero(x)
function _foursail_zero_like(x::AbstractArray)
    out = similar(x)
    fill!(out, zero(eltype(x)))
    return out
end

_foursail_copy_like(x::Number) = x
_foursail_copy_like(x::AbstractArray) = copy(x)

FourSAILResult(rddt, rsdt, rdot, rsot, taus) =
    FourSAILResult(rddt, rsdt, rdot, rsot, taus,
                   _foursail_zero_like(rsot),
                   _foursail_copy_like(rsot),
                   _foursail_zero_like(rsot))

function FourSAILGeometry(LAD::AbstractLeafDistribution;
                          sza_deg::Real,
                          vza_deg::Real,
                          raa_deg::Real,
                          quadrature::CanopyQuadrature=CanopyQuadrature())
    FT = float(promote_type(typeof(sza_deg), typeof(vza_deg), typeof(raa_deg)))
    μs = cosd(FT(sza_deg))
    μo = cosd(FT(vza_deg))
    μs > zero(FT) || throw(ArgumentError("sza_deg must be less than 90 degrees"))
    μo > zero(FT) || throw(ArgumentError("vza_deg must be less than 90 degrees"))

    dϕ = _foursail_wrap_azimuth(deg2rad(FT(raa_deg)))
    tan_s = sqrt(max(zero(FT), one(FT) - μs^2)) / μs
    tan_o = sqrt(max(zero(FT), one(FT) - μo^2)) / μo
    dso = sqrt(max(zero(FT), tan_s^2 + tan_o^2 - 2 * tan_s * tan_o * cos(dϕ)))

    ks = zero(FT); ko = zero(FT); bf = zero(FT)
    sob = zero(FT); sof = zero(FT)

    θl, wθ = gauleg(quadrature.n_leaf, zero(FT), FT(π / 2))
    lidf = FT.(pdf.(LAD.LD, FT(2) .* θl ./ FT(π))) .* FT(LAD.scaling)
    lidf ./= sum(wθ .* lidf)

    for i in eachindex(θl)
        weight = wθ[i] * lidf[i]
        χs, χo, frho, ftau = _foursail_volscatt(μs, μo, dϕ, θl[i])
        ks += weight * χs / μs
        ko += weight * χo / μo
        bf += weight * cos(θl[i])^2
        sob += weight * frho * FT(π) / (μs * μo)
        sof += weight * ftau * FT(π) / (μs * μo)
    end

    sdb = FT(0.5) * (ks + bf)
    sdf = FT(0.5) * (ks - bf)
    dob = FT(0.5) * (ko + bf)
    dof = FT(0.5) * (ko - bf)
    ddb = FT(0.5) * (one(FT) + bf)
    ddf = FT(0.5) * (one(FT) - bf)

    return FourSAILGeometry{FT}(μs, μo, dϕ, dso, ks, ko, sob, sof,
                                sdb, sdf, dob, dof, ddb, ddf)
end

FourSAILGeometry(; LAD=spherical_leaves(), kwargs...) =
    FourSAILGeometry(LAD; kwargs...)

Base.length(geoms::FourSAILGeometrySet) = length(geoms.ks)

@inline function Base.getindex(geoms::FourSAILGeometrySet, i::Integer)
    return FourSAILGeometry(geoms.μs[i], geoms.μo[i], geoms.dϕ[i],
                            geoms.dso[i], geoms.ks[i], geoms.ko[i],
                            geoms.sob[i], geoms.sof[i], geoms.sdb[i],
                            geoms.sdf[i], geoms.dob[i], geoms.dof[i],
                            geoms.ddb[i], geoms.ddf[i])
end

@inline _foursail_len(x::AbstractVector) = length(x)
@inline _foursail_len(x) = 1

@inline function _foursail_param_at(x::AbstractVector, i::Int, n::Int)
    return length(x) == 1 ? x[1] : x[i]
end
@inline _foursail_param_at(x, i::Int, n::Int) = x

function _foursail_geometry_count(sza_deg, vza_deg, raa_deg)
    n = max(_foursail_len(sza_deg), _foursail_len(vza_deg), _foursail_len(raa_deg))
    for (name, x) in ((:sza_deg, sza_deg), (:vza_deg, vza_deg), (:raa_deg, raa_deg))
        lx = _foursail_len(x)
        (lx == 1 || lx == n) ||
            throw(DimensionMismatch("$name must be scalar, length 1, or length $n"))
    end
    return n
end

function FourSAILGeometrySet(geoms::AbstractVector{<:FourSAILGeometry})
    isempty(geoms) && throw(ArgumentError("FourSAILGeometrySet needs at least one geometry"))
    FT = promote_type(map(g -> typeof(g.ks), geoms)...)
    return FourSAILGeometrySet(FT[g.μs for g in geoms],
                               FT[g.μo for g in geoms],
                               FT[g.dϕ for g in geoms],
                               FT[g.dso for g in geoms],
                               FT[g.ks for g in geoms],
                               FT[g.ko for g in geoms],
                               FT[g.sob for g in geoms],
                               FT[g.sof for g in geoms],
                               FT[g.sdb for g in geoms],
                               FT[g.sdf for g in geoms],
                               FT[g.dob for g in geoms],
                               FT[g.dof for g in geoms],
                               FT[g.ddb for g in geoms],
                               FT[g.ddf for g in geoms])
end

function FourSAILGeometrySet(LAD::AbstractLeafDistribution;
                             sza_deg,
                             vza_deg,
                             raa_deg,
                             quadrature::CanopyQuadrature=CanopyQuadrature())
    n = _foursail_geometry_count(sza_deg, vza_deg, raa_deg)
    geoms = [FourSAILGeometry(LAD;
                              sza_deg = _foursail_param_at(sza_deg, i, n),
                              vza_deg = _foursail_param_at(vza_deg, i, n),
                              raa_deg = _foursail_param_at(raa_deg, i, n),
                              quadrature = quadrature)
             for i in 1:n]
    return FourSAILGeometrySet(geoms)
end

FourSAILGeometrySet(; LAD=spherical_leaves(), kwargs...) =
    FourSAILGeometrySet(LAD; kwargs...)

@inline _foursail_value(x) = x
@inline _foursail_value(x::ForwardDiff.Dual) = ForwardDiff.value(x)
@inline _foursail_tol(x) = sqrt(eps(float(typeof(_foursail_value(x)))))
@inline _foursail_nonnegative(x) = _foursail_value(x) < 0 ? zero(x) : x

@inline function _foursail_mode(mode::Symbol)
    mode === :full && return Val(:full)
    mode === :single && return Val(:single)
    throw(ArgumentError("Unknown FourSAIL mode: $mode. Use :full or :single."))
end

@inline function _foursail_wrap_azimuth(ϕ::FT) where {FT<:Real}
    twopi = FT(2π)
    ψ = mod(abs(ϕ), twopi)
    return ψ > FT(π) ? twopi - ψ : ψ
end

function _foursail_volscatt(μs::FT, μo::FT, dϕ::FT, θl::FT) where {FT<:Real}
    costs = μs
    costo = μo
    sints = sqrt(max(zero(FT), one(FT) - μs^2))
    sinto = sqrt(max(zero(FT), one(FT) - μo^2))
    cospsi = cos(dϕ)

    costl = cos(θl)
    sintl = sin(θl)
    cs = costl * costs
    co = costl * costo
    ss = sintl * sints
    so = sintl * sinto

    cosbts = abs(ss) > FT(1e-6) ? -cs / ss : FT(5)
    cosbto = abs(so) > FT(1e-6) ? -co / so : FT(5)

    if abs(cosbts) < one(FT)
        bts = acos(cosbts)
        ds = ss
    else
        bts = FT(π)
        ds = cs
    end
    χs = FT(2 / π) * ((bts - FT(π / 2)) * cs + sin(bts) * ss)

    if abs(cosbto) < one(FT)
        bto = acos(cosbto)
        doo = so
    else
        bto = FT(π)
        doo = co
    end
    χo = FT(2 / π) * ((bto - FT(π / 2)) * co + sin(bto) * so)

    btran1 = abs(bts - bto)
    btran2 = FT(π) - abs(bts + bto - FT(π))

    if dϕ <= btran1
        bt1 = dϕ
        bt2 = btran1
        bt3 = btran2
    else
        bt1 = btran1
        if dϕ <= btran2
            bt2 = dϕ
            bt3 = btran2
        else
            bt2 = btran2
            bt3 = dϕ
        end
    end

    t1 = 2 * cs * co + ss * so * cospsi
    t2 = bt2 > zero(FT) ?
         sin(bt2) * (2 * ds * doo + ss * so * cos(bt1) * cos(bt3)) :
         zero(FT)

    denom = FT(2) * FT(π)^2
    frho = ((FT(π) - bt2) * t1 + t2) / denom
    ftau = (-bt2 * t1 + t2) / denom
    return χs, χo, max(frho, zero(FT)), max(ftau, zero(FT))
end

@inline function _foursail_j1(k, l, lai)
    del = (k - l) * lai
    if _foursail_value(abs(del)) <= 1e-3
        return lai * FT_like(del, 0.5) * (exp(-k * lai) + exp(-l * lai)) *
               (one(del) - del^2 / 12)
    else
        return (exp(-l * lai) - exp(-k * lai)) / (k - l)
    end
end

@inline FT_like(x, value) = oftype(x + one(x), value)

@inline function _foursail_j2(k, l, lai)
    denom = k + l
    return -expm1(-denom * lai) / denom
end

function _foursail_hotspot(lai, hotspot, tss, ks, ko, dso)
    if _foursail_value(hotspot) <= 0
        return tss * exp(-ko * lai), _foursail_j2(ks, ko, lai) / lai
    end

    alf = dso / hotspot * 2 / (ks + ko)
    if _foursail_value(alf) > 200
        alf = oftype(alf, 200)
    end

    if _foursail_value(abs(alf)) <= _foursail_tol(alf)
        return tss, -expm1(-ks * lai) / (ks * lai)
    end

    fhot = lai * sqrt(ks * ko)
    x1 = zero(alf)
    y1 = zero(alf)
    f1 = one(alf)
    ca = exp(-alf)
    fint = (one(alf) - ca) / 20
    sumint = zero(alf)

    for i in 1:20
        x2 = i < 20 ? -log(one(alf) - i * fint) / alf : one(alf)
        y2 = -(ko + ks) * lai * x2 + fhot * (one(alf) - exp(-alf * x2)) / alf
        f2 = exp(y2)
        dy = y2 - y1
        sumint += _foursail_value(abs(dy)) <= _foursail_tol(dy) ?
                  f1 * (x2 - x1) :
                  (f2 - f1) * (x2 - x1) / dy
        x1 = x2
        y1 = y2
        f1 = f2
    end

    return f1, sumint
end

@inline _foursail_eval(ρ, τ, rsoil, geom::FourSAILGeometry, lai, hotspot) =
    _foursail_eval(ρ, τ, rsoil, geom, lai, hotspot, Val(:full))

function _foursail_eval(ρ, τ, rsoil, geom::FourSAILGeometry, lai, hotspot,
                        ::Val{:single})
    if _foursail_value(lai) <= _foursail_tol(lai)
        z = zero(ρ + τ + rsoil + lai)
        return FourSAILResult(rsoil, rsoil, rsoil, rsoil,
                              one(ρ + τ + rsoil + lai), z, rsoil, z)
    end

    w = geom.sob * ρ + geom.sof * τ
    tss = exp(-geom.ks * lai)
    tsstoo, sumint = _foursail_hotspot(lai, hotspot, tss, geom.ks, geom.ko, geom.dso)
    rsos = w * lai * sumint
    rsost = rsos + tsstoo * rsoil
    z = zero(rsost)
    return FourSAILResult(z, z, z, rsost, tss, rsos, rsost, z)
end

function _foursail_eval(ρ, τ, rsoil, geom::FourSAILGeometry, lai, hotspot,
                        ::Val{:full})
    if _foursail_value(lai) <= _foursail_tol(lai)
        z = zero(ρ + τ + rsoil + lai)
        return FourSAILResult(rsoil, rsoil, rsoil, rsoil,
                              one(ρ + τ + rsoil + lai), z, rsoil, z)
    end

    sigb = geom.ddb * ρ + geom.ddf * τ
    sigf = geom.ddf * ρ + geom.ddb * τ
    att = one(sigf) - sigf
    m2 = (att + sigb) * (att - sigb)
    # Conservative-leaf singularity: at ρ + τ = 1 we get m → 0, rinf → 1,
    # `denom = 1 − rinf²·e²` → 0, and rdd/tdd/tsd/… reduce to 0/0. Floor m
    # at ∛eps (well above the catastrophic-cancellation threshold for
    # `1 − e2` while still corresponding to an input perturbation
    # ~10⁻¹¹, so the answer matches the well-defined ρ + τ → 1⁻ limit to
    # several digits). Same workaround as canonical PROSAIL Fortran.
    m = max(sqrt(_foursail_nonnegative(m2)), cbrt(eps(float(typeof(_foursail_value(m2))))))

    sb = geom.sdb * ρ + geom.sdf * τ
    sf = geom.sdf * ρ + geom.sdb * τ
    vb = geom.dob * ρ + geom.dof * τ
    vf = geom.dof * ρ + geom.dob * τ
    w = geom.sob * ρ + geom.sof * τ

    rinf = _foursail_value(abs(sigb)) <= _foursail_tol(sigb) ?
           zero(sigb) : (att - m) / sigb
    e1 = exp(-m * lai)
    e2 = e1 * e1
    rinf2 = rinf * rinf
    re = rinf * e1
    denom = one(rinf2) - rinf2 * e2

    J1ks = _foursail_j1(geom.ks, m, lai)
    J2ks = _foursail_j2(geom.ks, m, lai)
    J1ko = _foursail_j1(geom.ko, m, lai)
    J2ko = _foursail_j2(geom.ko, m, lai)

    Ps = (sf + sb * rinf) * J1ks
    Qs = (sf * rinf + sb) * J2ks
    Pv = (vf + vb * rinf) * J1ko
    Qv = (vf * rinf + vb) * J2ko

    rdd = rinf * (one(e2) - e2) / denom
    tdd = (one(rinf2) - rinf2) * e1 / denom
    tsd = (Ps - re * Qs) / denom
    rsd = (Qs - re * Ps) / denom
    tdo = (Pv - re * Qv) / denom
    rdo = (Qv - re * Pv) / denom

    tss = exp(-geom.ks * lai)
    too = exp(-geom.ko * lai)

    z = _foursail_j2(geom.ks, geom.ko, lai)
    g1 = (z - J1ks * too) / (geom.ko + m)
    g2 = (z - J1ko * tss) / (geom.ks + m)
    Tv1 = (vf * rinf + vb) * g1
    Tv2 = (vf + vb * rinf) * g2
    T1 = Tv1 * (sf + sb * rinf)
    T2 = Tv2 * (sf * rinf + sb)
    T3 = (rdo * Qs + tdo * Ps) * rinf
    rsod = (T1 + T2 - T3) / (one(rinf2) - rinf2)

    tsstoo, sumint = _foursail_hotspot(lai, hotspot, tss, geom.ks, geom.ko, geom.dso)
    rsos = w * lai * sumint

    dn = one(rsoil) - rsoil * rdd
    rddt = rdd + tdd * rsoil * tdd / dn
    rsdt = rsd + (tsd + tss) * rsoil * tdd / dn
    rdot = rdo + tdd * rsoil * (tdo + too) / dn
    rsodt = rsod + ((tss + tsd) * tdo + (tsd + tss * rsoil * rdd) * too) *
                    rsoil / dn
    rsost = rsos + tsstoo * rsoil
    rsot = rsost + rsodt
    taus = tsd + tss

    return FourSAILResult(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)
end

@kernel function _foursail_kernel_vecsoil!(rddt, rsdt, rdot, rsot, taus,
                                           rsos, rsost, rsodt,
                                           leaf_ρ, leaf_τ, soil_ρ,
                                           geom, lai, hotspot, mode)
    i = @index(Global)
    @inbounds begin
        result = _foursail_eval(leaf_ρ[i], leaf_τ[i], soil_ρ[i],
                                geom, lai, hotspot, mode)
        rddt[i] = result.rddt
        rsdt[i] = result.rsdt
        rdot[i] = result.rdot
        rsot[i] = result.rsot
        taus[i] = result.taus
        rsos[i] = result.rsos
        rsost[i] = result.rsost
        rsodt[i] = result.rsodt
    end
end

@kernel function _foursail_kernel_scalarsoil!(rddt, rsdt, rdot, rsot, taus,
                                              rsos, rsost, rsodt,
                                              leaf_ρ, leaf_τ, soil_ρ,
                                              geom, lai, hotspot, mode)
    i = @index(Global)
    @inbounds begin
        result = _foursail_eval(leaf_ρ[i], leaf_τ[i], soil_ρ,
                                geom, lai, hotspot, mode)
        rddt[i] = result.rddt
        rsdt[i] = result.rsdt
        rdot[i] = result.rdot
        rsot[i] = result.rsot
        taus[i] = result.taus
        rsos[i] = result.rsos
        rsost[i] = result.rsost
        rsodt[i] = result.rsodt
    end
end

@kernel function _foursail_angle_kernel_scalarsoil!(rddt, rsdt, rdot, rsot, taus,
                                                    rsos, rsost, rsodt,
                                                    leaf_ρ, leaf_τ, soil_ρ,
                                                    μs, μo, dϕ, dso, ks, ko,
                                                    sob, sof, sdb, sdf, dob, dof,
                                                    ddb, ddf, lai, hotspot, mode)
    iλ, ia = @index(Global, NTuple)
    @inbounds begin
        geom = FourSAILGeometry(μs[ia], μo[ia], dϕ[ia], dso[ia], ks[ia], ko[ia],
                                sob[ia], sof[ia], sdb[ia], sdf[ia], dob[ia],
                                dof[ia], ddb[ia], ddf[ia])
        result = _foursail_eval(leaf_ρ[iλ], leaf_τ[iλ], soil_ρ,
                                geom, lai, hotspot, mode)
        rddt[iλ, ia] = result.rddt
        rsdt[iλ, ia] = result.rsdt
        rdot[iλ, ia] = result.rdot
        rsot[iλ, ia] = result.rsot
        taus[iλ, ia] = result.taus
        rsos[iλ, ia] = result.rsos
        rsost[iλ, ia] = result.rsost
        rsodt[iλ, ia] = result.rsodt
    end
end

@kernel function _foursail_angle_kernel_vecsoil!(rddt, rsdt, rdot, rsot, taus,
                                                 rsos, rsost, rsodt,
                                                 leaf_ρ, leaf_τ, soil_ρ,
                                                 μs, μo, dϕ, dso, ks, ko,
                                                 sob, sof, sdb, sdf, dob, dof,
                                                 ddb, ddf, lai, hotspot, mode)
    iλ, ia = @index(Global, NTuple)
    @inbounds begin
        geom = FourSAILGeometry(μs[ia], μo[ia], dϕ[ia], dso[ia], ks[ia], ko[ia],
                                sob[ia], sof[ia], sdb[ia], sdf[ia], dob[ia],
                                dof[ia], ddb[ia], ddf[ia])
        result = _foursail_eval(leaf_ρ[iλ], leaf_τ[iλ], soil_ρ[iλ],
                                geom, lai, hotspot, mode)
        rddt[iλ, ia] = result.rddt
        rsdt[iλ, ia] = result.rsdt
        rdot[iλ, ia] = result.rdot
        rsot[iλ, ia] = result.rsot
        taus[iλ, ia] = result.taus
        rsos[iλ, ia] = result.rsos
        rsost[iλ, ia] = result.rsost
        rsodt[iλ, ia] = result.rsodt
    end
end

@kernel function _foursail_angle_kernel_matsoil!(rddt, rsdt, rdot, rsot, taus,
                                                 rsos, rsost, rsodt,
                                                 leaf_ρ, leaf_τ, soil_ρ,
                                                 μs, μo, dϕ, dso, ks, ko,
                                                 sob, sof, sdb, sdf, dob, dof,
                                                 ddb, ddf, lai, hotspot, mode)
    iλ, ia = @index(Global, NTuple)
    @inbounds begin
        geom = FourSAILGeometry(μs[ia], μo[ia], dϕ[ia], dso[ia], ks[ia], ko[ia],
                                sob[ia], sof[ia], sdb[ia], sdf[ia], dob[ia],
                                dof[ia], ddb[ia], ddf[ia])
        result = _foursail_eval(leaf_ρ[iλ], leaf_τ[iλ], soil_ρ[iλ, ia],
                                geom, lai, hotspot, mode)
        rddt[iλ, ia] = result.rddt
        rsdt[iλ, ia] = result.rsdt
        rdot[iλ, ia] = result.rdot
        rsot[iλ, ia] = result.rsot
        taus[iλ, ia] = result.taus
        rsos[iλ, ia] = result.rsos
        rsost[iλ, ia] = result.rsost
        rsodt[iλ, ia] = result.rsodt
    end
end

_soil_eltype(x::AbstractArray) = eltype(x)
_soil_eltype(x::Number) = typeof(x)

function _check_foursail_soil_dims(soil_ρ, nλ::Int, nangle::Int)
    if soil_ρ isa AbstractMatrix
        size(soil_ρ) == (nλ, nangle) ||
            throw(DimensionMismatch("soil_reflectance matrix must have size ($nλ, $nangle)"))
    elseif soil_ρ isa AbstractVector
        length(soil_ρ) == nλ ||
            throw(DimensionMismatch("soil_reflectance vector must match leaf spectra length"))
    end
    return nothing
end

function _geometry_for_leaf(geoms::FourSAILGeometrySet, leaf_ρ::AbstractArray)
    if leaf_ρ isa CUDA.CuArray
        geoms.ks isa CUDA.CuArray && return geoms
        return FourSAILGeometrySet(CUDA.CuArray(geoms.μs),
                                   CUDA.CuArray(geoms.μo),
                                   CUDA.CuArray(geoms.dϕ),
                                   CUDA.CuArray(geoms.dso),
                                   CUDA.CuArray(geoms.ks),
                                   CUDA.CuArray(geoms.ko),
                                   CUDA.CuArray(geoms.sob),
                                   CUDA.CuArray(geoms.sof),
                                   CUDA.CuArray(geoms.sdb),
                                   CUDA.CuArray(geoms.sdf),
                                   CUDA.CuArray(geoms.dob),
                                   CUDA.CuArray(geoms.dof),
                                   CUDA.CuArray(geoms.ddb),
                                   CUDA.CuArray(geoms.ddf))
    else
        geoms.ks isa CUDA.CuArray || return geoms
        return FourSAILGeometrySet(Array(geoms.μs),
                                   Array(geoms.μo),
                                   Array(geoms.dϕ),
                                   Array(geoms.dso),
                                   Array(geoms.ks),
                                   Array(geoms.ko),
                                   Array(geoms.sob),
                                   Array(geoms.sof),
                                   Array(geoms.sdb),
                                   Array(geoms.sdf),
                                   Array(geoms.dob),
                                   Array(geoms.dof),
                                   Array(geoms.ddb),
                                   Array(geoms.ddf))
    end
end

_array_for_leaf(x::Number, leaf_ρ::AbstractArray) = x
function _array_for_leaf(x::AbstractArray, leaf_ρ::AbstractArray)
    if leaf_ρ isa CUDA.CuArray
        return x isa CUDA.CuArray ? x : CUDA.CuArray(x)
    else
        return x isa CUDA.CuArray ? Array(x) : x
    end
end

function _allocate_foursail_result(leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                                   soil_ρ, geom::FourSAILGeometry, lai, hotspot)
    n = length(leaf_ρ)
    length(leaf_τ) == n ||
        throw(DimensionMismatch("leaf_reflectance and leaf_transmittance must have the same length"))
    soil_ρ isa AbstractVector && length(soil_ρ) != n &&
        throw(DimensionMismatch("soil_reflectance must be scalar or match leaf spectra length"))
    OT = promote_type(eltype(leaf_ρ), eltype(leaf_τ), _soil_eltype(soil_ρ),
                      typeof(lai), typeof(hotspot), typeof(geom.ks))
    return FourSAILResult(similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n),
                          similar(leaf_ρ, OT, n))
end

function _allocate_foursail_result(leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                                   soil_ρ, geoms::FourSAILGeometrySet, lai, hotspot)
    nλ = length(leaf_ρ)
    nangle = length(geoms)
    length(leaf_τ) == nλ ||
        throw(DimensionMismatch("leaf_reflectance and leaf_transmittance must have the same length"))
    _check_foursail_soil_dims(soil_ρ, nλ, nangle)
    OT = promote_type(eltype(leaf_ρ), eltype(leaf_τ), _soil_eltype(soil_ρ),
                      typeof(lai), typeof(hotspot), eltype(geoms.ks))
    return FourSAILResult(similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle),
                          similar(leaf_ρ, OT, nλ, nangle))
end

"""
    foursail!(rddt, rsdt, rdot, rsot, taus,
              leaf_reflectance, leaf_transmittance, soil_reflectance,
              geometry, LAI; hotspot=0)
    foursail!(rddt, rsdt, rdot, rsot, taus,
              leaf_reflectance, leaf_transmittance, soil_reflectance,
              geometry_set, LAI; hotspot=0, mode=:full)

Evaluate 4SAIL over a wavelength vector using one KernelAbstractions work item
per wavelength. `soil_reflectance` may be a scalar or a vector matching the
leaf spectra.

For a [`FourSAILGeometrySet`](@ref), outputs are matrices of size
`(n_wavelength, n_angle)` and the kernel evaluates all wavelength-angle pairs.
Full-mode results also carry the bidirectional decomposition `rsost` and
`rsodt`, where `rsost` is the direct/direct term (leaf single scattering plus
directly transmitted soil) and `rsot = rsost + rsodt`. Use `mode=:single` only
when you want to skip the diffuse equations entirely; in that mode `rsot` equals
`rsost`, `rsodt` is zero, and diffuse-stream outputs are filled with zero.
"""
function _foursail!(rddt::AbstractVector, rsdt::AbstractVector,
                    rdot::AbstractVector, rsot::AbstractVector,
                    taus::AbstractVector, rsos::AbstractVector,
                    rsost::AbstractVector, rsodt::AbstractVector,
                    leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                    soil_ρ, geom::FourSAILGeometry, lai::Real;
                    hotspot::Real=0, mode::Symbol=:full)
    n = length(leaf_ρ)
    length(leaf_τ) == n ||
        throw(DimensionMismatch("leaf_reflectance and leaf_transmittance must have the same length"))
    all(length(out) == n for out in (rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)) ||
        throw(DimensionMismatch("all output arrays must match leaf spectra length"))

    backend = KernelAbstractions.get_backend(leaf_ρ)
    mode_val = _foursail_mode(mode)
    leaf_τ_dev = _array_for_leaf(leaf_τ, leaf_ρ)
    soil_dev = _array_for_leaf(soil_ρ, leaf_ρ)
    if soil_dev isa AbstractVector
        length(soil_dev) == n ||
            throw(DimensionMismatch("soil_reflectance vector must match leaf spectra length"))
        kernel! = _foursail_kernel_vecsoil!(backend)
        kernel!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                leaf_ρ, leaf_τ_dev, soil_dev,
                geom, lai, hotspot, mode_val; ndrange=n)
    else
        kernel! = _foursail_kernel_scalarsoil!(backend)
        kernel!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                leaf_ρ, leaf_τ_dev, soil_dev,
                geom, lai, hotspot, mode_val; ndrange=n)
    end
    KernelAbstractions.synchronize(backend)
    return FourSAILResult(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)
end

"""
    foursail!(rddt, rsdt, rdot, rsot, taus, leaf_ρ, leaf_τ, soil_ρ,
              geom, lai; hotspot=0, mode=:full)
    foursail!(result::FourSAILResult, leaf_ρ, leaf_τ, soil_ρ,
              geom, lai; hotspot=0, mode=:full)

In-place 4SAIL evaluation: writes canopy + soil reflectance and
transmittance factors into the supplied output arrays (or
[`FourSAILResult`](@ref)) for one geometry or a batched
[`FourSAILGeometrySet`](@ref). Outputs are dispatched on the shape of
the inputs:

- vector outputs + a [`FourSAILGeometry`](@ref) → one wavelength per
  vector element, one geometry.
- matrix outputs + a [`FourSAILGeometrySet`](@ref) → wavelengths along
  rows, geometries along columns; the GPU-friendly batched path.

The full 4SAIL bidirectional reflectance factor `rsot` is decomposed
as `rsot = rsost + rsodt`, where `rsost` is the direct/direct branch
(single scattering plus directly transmitted soil reflectance) and
`rsodt` is the diffuse / multiple-scattering remainder — see
[`FourSAILResult`](@ref).

Allocating wrappers: see [`foursail`](@ref).

# Reference
Verhoef, W. (1984). *Light scattering by leaf layers with application
to canopy reflectance modeling: the SAIL model.* Remote Sens. Env. 16,
125–141. Extended in Verhoef et al. (2007), 4SAIL.
"""
function foursail!(rddt::AbstractVector, rsdt::AbstractVector,
                   rdot::AbstractVector, rsot::AbstractVector,
                   taus::AbstractVector,
                   leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                   soil_ρ, geom::FourSAILGeometry, lai::Real;
                   hotspot::Real=0, mode::Symbol=:full)
    rsos = similar(rsot)
    rsost = similar(rsot)
    rsodt = similar(rsot)
    return _foursail!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                      leaf_ρ, leaf_τ, soil_ρ, geom, lai;
                      hotspot=hotspot, mode=mode)
end

function _foursail!(rddt::AbstractMatrix, rsdt::AbstractMatrix,
                    rdot::AbstractMatrix, rsot::AbstractMatrix,
                    taus::AbstractMatrix, rsos::AbstractMatrix,
                    rsost::AbstractMatrix, rsodt::AbstractMatrix,
                    leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                    soil_ρ, geoms::FourSAILGeometrySet, lai::Real;
                    hotspot::Real=0, mode::Symbol=:full)
    nλ = length(leaf_ρ)
    nangle = length(geoms)
    length(leaf_τ) == nλ ||
        throw(DimensionMismatch("leaf_reflectance and leaf_transmittance must have the same length"))
    all(size(out) == (nλ, nangle) for out in (rddt, rsdt, rdot, rsot, taus,
                                              rsos, rsost, rsodt)) ||
        throw(DimensionMismatch("all output arrays must have size ($nλ, $nangle)"))
    _check_foursail_soil_dims(soil_ρ, nλ, nangle)

    backend = KernelAbstractions.get_backend(leaf_ρ)
    geoms_dev = _geometry_for_leaf(geoms, leaf_ρ)
    leaf_τ_dev = _array_for_leaf(leaf_τ, leaf_ρ)
    soil_dev = _array_for_leaf(soil_ρ, leaf_ρ)
    mode_val = _foursail_mode(mode)
    if soil_dev isa AbstractMatrix
        kernel! = _foursail_angle_kernel_matsoil!(backend)
        kernel!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                leaf_ρ, leaf_τ_dev, soil_dev,
                geoms_dev.μs, geoms_dev.μo, geoms_dev.dϕ, geoms_dev.dso,
                geoms_dev.ks, geoms_dev.ko, geoms_dev.sob, geoms_dev.sof,
                geoms_dev.sdb, geoms_dev.sdf, geoms_dev.dob, geoms_dev.dof,
                geoms_dev.ddb, geoms_dev.ddf, lai, hotspot, mode_val;
                ndrange=(nλ, nangle))
    elseif soil_dev isa AbstractVector
        kernel! = _foursail_angle_kernel_vecsoil!(backend)
        kernel!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                leaf_ρ, leaf_τ_dev, soil_dev,
                geoms_dev.μs, geoms_dev.μo, geoms_dev.dϕ, geoms_dev.dso,
                geoms_dev.ks, geoms_dev.ko, geoms_dev.sob, geoms_dev.sof,
                geoms_dev.sdb, geoms_dev.sdf, geoms_dev.dob, geoms_dev.dof,
                geoms_dev.ddb, geoms_dev.ddf, lai, hotspot, mode_val;
                ndrange=(nλ, nangle))
    else
        kernel! = _foursail_angle_kernel_scalarsoil!(backend)
        kernel!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                leaf_ρ, leaf_τ_dev, soil_dev,
                geoms_dev.μs, geoms_dev.μo, geoms_dev.dϕ, geoms_dev.dso,
                geoms_dev.ks, geoms_dev.ko, geoms_dev.sob, geoms_dev.sof,
                geoms_dev.sdb, geoms_dev.sdf, geoms_dev.dob, geoms_dev.dof,
                geoms_dev.ddb, geoms_dev.ddf, lai, hotspot, mode_val;
                ndrange=(nλ, nangle))
    end
    KernelAbstractions.synchronize(backend)
    return FourSAILResult(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt)
end

function foursail!(rddt::AbstractMatrix, rsdt::AbstractMatrix,
                   rdot::AbstractMatrix, rsot::AbstractMatrix,
                   taus::AbstractMatrix,
                   leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                   soil_ρ, geoms::FourSAILGeometrySet, lai::Real;
                   hotspot::Real=0, mode::Symbol=:full)
    rsos = similar(rsot)
    rsost = similar(rsot)
    rsodt = similar(rsot)
    return _foursail!(rddt, rsdt, rdot, rsot, taus, rsos, rsost, rsodt,
                      leaf_ρ, leaf_τ, soil_ρ, geoms, lai;
                      hotspot=hotspot, mode=mode)
end

function foursail!(result::FourSAILResult,
                   leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                   soil_ρ, geom::FourSAILGeometry, lai::Real;
                   hotspot::Real=0, mode::Symbol=:full)
    return _foursail!(result.rddt, result.rsdt, result.rdot, result.rsot, result.taus,
                      result.rsos, result.rsost, result.rsodt,
                      leaf_ρ, leaf_τ, soil_ρ, geom, lai;
                      hotspot=hotspot, mode=mode)
end

function foursail!(result::FourSAILResult,
                   leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                   soil_ρ, geoms::FourSAILGeometrySet, lai::Real;
                   hotspot::Real=0, mode::Symbol=:full)
    return _foursail!(result.rddt, result.rsdt, result.rdot, result.rsot, result.taus,
                      result.rsos, result.rsost, result.rsodt,
                      leaf_ρ, leaf_τ, soil_ρ, geoms, lai;
                      hotspot=hotspot, mode=mode)
end

"""
    foursail(leaf_reflectance, leaf_transmittance, soil_reflectance,
             geometry, LAI; hotspot=0)
    foursail(leaf_reflectance, leaf_transmittance, soil_reflectance,
             geometry_set, LAI; hotspot=0, mode=:full)

Allocate outputs and evaluate 4SAIL over all wavelengths. For array inputs the
batched kernel path is used; for scalar inputs this returns a scalar
[`FourSAILResult`](@ref).
"""
function foursail(leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                  soil_ρ, geom::FourSAILGeometry, lai::Real;
                  hotspot::Real=0, mode::Symbol=:full)
    result = _allocate_foursail_result(leaf_ρ, leaf_τ, soil_ρ, geom, lai, hotspot)
    return foursail!(result, leaf_ρ, leaf_τ, soil_ρ, geom, lai;
                     hotspot=hotspot, mode=mode)
end

function foursail(leaf_ρ::AbstractVector, leaf_τ::AbstractVector,
                  soil_ρ, geoms::FourSAILGeometrySet, lai::Real;
                  hotspot::Real=0, mode::Symbol=:full)
    result = _allocate_foursail_result(leaf_ρ, leaf_τ, soil_ρ, geoms, lai, hotspot)
    return foursail!(result, leaf_ρ, leaf_τ, soil_ρ, geoms, lai;
                     hotspot=hotspot, mode=mode)
end

function foursail(leaf_ρ::Real, leaf_τ::Real, soil_ρ::Real,
                  geom::FourSAILGeometry, lai::Real; hotspot::Real=0,
                  mode::Symbol=:full)
    return _foursail_eval(leaf_ρ, leaf_τ, soil_ρ, geom, lai, hotspot,
                          _foursail_mode(mode))
end
