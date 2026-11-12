
"""
    A(θ::FT, θₗ::FT) where FT<:Real

Returns the azimuthally-integrated projected leaf area function for beam direction `θ`
and leaf inclination `θₗ` (both in radians), assuming azimuthal uniformity.
Follows Bonan, *Modeling Earth's Climate*, Eq. 14.24:

``A(θ, θₗ) = \\begin{cases}
  \\cos θ\\,\\cos θₗ & θₗ \\le π/2 - θ \\\\
  \\frac{2}{π}\\!\\left[\\sqrt{\\sin^2 θₗ - \\cos^2 θ} + \\cos θ\\,\\cos θₗ\\,
    \\arcsin\\!\\left(\\frac{\\cos θ\\,\\cos θₗ}{\\sin θ\\,\\sin θₗ}\\right)\\right]
  & \\text{otherwise}
\\end{cases}``

Used by [`G`](@ref) to compute the Ross G-function.
See also [`Asm`](@ref) for the algebraically equivalent Shultis & Myneni form.
"""
function A(θ::FT, θₗ::FT) where FT<:Real
    a = cos(θ) * cos(θₗ)
    # Bonan (2019), Eq. 14.24
    if θₗ ≤ (FT(π/2) - θ)
        return a
    else
        b = sin(θ) * sin(θₗ)
        c = sqrt(max(zero(FT), sin(θₗ)^2 - cos(θ)^2))
        return FT(2/π)*(c + a * asin(a/b))
    end
end

"""
    Asm(θ::FT, θₗ::FT) where FT<:Real

Algebraically equivalent form of [`A`](@ref) following Shultis & Myneni (1987)
Eqs. 35–36.  Used for cross-validation; production code uses [`A`](@ref).
"""
function Asm(θ::FT, θₗ::FT) where FT<:Real
    a = cos(θ) * cos(θₗ)
    if θₗ ≤ (FT(π/2) - θ)
        return a
    else
        b = sin(θ) * sin(θₗ)
        c = sqrt(FT(1)-(a/b)^2)
        return FT(2/π)*(a * (acos(-a/b)-π/2) + b * c)
    end
end



"""
    _leaf_inclination_quadrature(LD, n_leaf, ::Type{FT})

Compute leaf-inclination Gauss-Legendre nodes `θ` on `[0, π/2]` and normalized
leaf-area quadrature weights `w_measure` in precision `FT`.
"""
function _leaf_inclination_quadrature(LD::AbstractLeafDistribution, n_leaf::Int, ::Type{FT}) where {FT<:Real}
    TFT = float(FT)
    BFT = _leaf_base_float(TFT)
    θ, w = gauleg(n_leaf, zero(TFT), TFT(π) / TFT(2))
    w_measure = Vector{TFT}(undef, length(θ))
    inv_half_pi = TFT(2) / TFT(π)
    scale = TFT(LD.scaling)
    dist = LD.LD
    mass = zero(TFT)
    @inbounds for i in eachindex(θ, w, w_measure)
        wm = w[i] * TFT(pdf(dist, inv_half_pi * θ[i])) * scale
        w_measure[i] = wm
        mass += wm
    end
    isfinite(prospect_primal(mass)) && prospect_primal(mass) > zero(BFT) ||
        throw(ArgumentError("leaf-inclination quadrature must have finite positive mass"))
    inv_mass = inv(mass)
    @inbounds for i in eachindex(w_measure)
        w_measure[i] *= inv_mass
    end
    return (; θ, w_measure)
end

"""
    G(μ, LD; quadrature = CanopyQuadrature(), nLeg = nothing)

Compute the Ross projected leaf area factor in directions `μ` for an azimuthally
uniform leaf inclination distribution `LD`, preserving the container shape and
axes of `μ`.

Uses `quadrature.n_leaf` (or `nLeg` when explicitly passed) and the normalized
leaf-area measure shared with the canopy scattering kernels.

# Examples
```jldoctest
julia> isapprox(G([0.5, 1.0], spherical_leaves()), [0.5, 0.5]; atol = 1e-3)
true
```
"""
function G(
    μ::AbstractArray{FT},
    LD::AbstractLeafDistribution;
    quadrature::CanopyQuadrature = CanopyQuadrature(),
    nLeg = nothing,
) where {FT<:Real}
    LFT = _lad_ft(LD)
    OutFT = LFT <: ForwardDiff.Dual ? float(promote_type(FT, LFT)) : float(FT)
    WFT = LFT <: ForwardDiff.Dual ? OutFT : _leaf_base_float(OutFT)
    n_leaf = nLeg === nothing ? quadrature.n_leaf : Int(nLeg)
    leaf_quad = _leaf_inclination_quadrature(LD, n_leaf, WFT)
    θ_l = leaf_quad.θ
    w_l = leaf_quad.w_measure
    out = similar(μ, OutFT)
    fill!(out, zero(OutFT))
    @inbounds for l in eachindex(θ_l, w_l)
        μ_l = OutFT(cos(θ_l[l]))
        wl = OutFT(w_l[l])
        for idx in eachindex(out, μ)
            μ_abs = clamp(abs(OutFT(μ[idx])), zero(OutFT), one(OutFT))
            out[idx] = muladd(wl, H(μ_abs, μ_l) + H(-μ_abs, μ_l), out[idx])
        end
    end
    return out
end

"""
    G2(μ::AbstractArray, LD::AbstractLeafDistribution; nLeg = 40)

Compute the Ross projection factor `G(μ)` by integrating in leaf-normal cosine
`μ_L` instead of leaf inclination angle `θ_L` (retained for cross-checking [`G`](@ref)).
"""
function G2(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg = 40) where {FT<:Real}
    OutFT = float(FT)
    μl, w = gauleg(Int(nLeg), zero(OutFT), one(OutFT))
    θₗ = acos.(μl)
    Fᵢ = OutFT[OutFT(pdf(LD.LD, OutFT(2) * θ / OutFT(π))) * OutFT(LD.scaling) * (OutFT(π) / OutFT(2)) for θ in θₗ]
    w_norm = (w .* Fᵢ) ./ dot(w, Fᵢ)
    out = similar(μ, OutFT)
    @inbounds for idx in eachindex(out, μ)
        θ_i = acos(clamp(abs(OutFT(μ[idx])), zero(OutFT), one(OutFT)))
        acc = zero(OutFT)
        for l in eachindex(θₗ, w_norm)
            acc = muladd(w_norm[l], A(θ_i, θₗ[l]), acc)
        end
        out[idx] = acc
    end
    return out
end

"""
    bfG(μ::AbstractArray, LD::AbstractLeafDistribution; nLeg = 20)

Compute the Ross projection factor `G(μ)` by direct 2D quadrature of `|Ω ⋅ Ω_L|`
over leaf inclination and azimuth (reference implementation for [`G`](@ref)).
"""
function bfG(μ::AbstractArray{FT}, LD::AbstractLeafDistribution; nLeg = 20) where {FT<:Real}
    OutFT = float(FT)
    nQuad = 100
    ϕ, w_azi = gauleg(nQuad, zero(OutFT), OutFT(2π))
    μ_l, w = gauleg(180, zero(OutFT), one(OutFT))
    Ω_l = [DirectionCosine(a, b) for a in μ_l, b in ϕ]
    θₗ = acos.(μ_l)
    Fᵢ = OutFT[OutFT(pdf(LD.LD, OutFT(2) * θ / OutFT(π))) * OutFT(LD.scaling) * (OutFT(π) / OutFT(2)) for θ in θₗ]
    Fᵢ ./= dot(Fᵢ, w)
    res = similar(μ, OutFT)
    for i in eachindex(res, μ)
        Ω = DirectionCosine(abs(OutFT(μ[i])), zero(OutFT))
        res[i] = ((Fᵢ .* abs.(dot.((Ω,), Ω_l)))' * w)' * w_azi / OutFT(2π)
    end
    return res
end

"""
    compute_lambertian_Γ(μ, μꜛ, r, t, LD::AbstractLeafDistribution; nLeg = 20)

Compute the azimuthally averaged area scattering transfer matrix `Γ` of size
`(length(μ), length(μꜛ))` following Shultis & Myneni (1988), Eq. 43:

``Γ(μ_i, μ'_j) = \\int_0^{π/2} g_L(θ_L)\\,[t\\,Ψ^+(μ_i, μ'_j, \\cos θ_L) + r\\,Ψ^-(μ_i, μ'_j, \\cos θ_L)]\\,dθ_L``

for leaf reflectance `r`, leaf transmittance `t`, and azimuthally uniform leaf
angle distribution `LD`.

# Examples
```jldoctest
julia> Γ = CanopyOptics.compute_lambertian_Γ([0.4, 0.8], [-0.5], 0.4, 0.1, spherical_leaves());

julia> size(Γ)
(2, 1)
```
"""
function compute_lambertian_Γ(
    μ::AbstractVector{FT1},
    μꜛ::AbstractVector{FT2},
    r::Real,
    t::Real,
    LD::AbstractLeafDistribution{FTL};
    nLeg::Integer = 20,
) where {FT1<:Real,FT2<:Real,FTL<:Real}
    OutFT = float(promote_type(FT1, FT2, typeof(r), typeof(t), FTL))
    Γ = similar(μ, OutFT, axes(μ, 1), axes(μꜛ, 1))
    fill!(Γ, zero(OutFT))
    θₗ, w = gauleg(Int(nLeg), zero(OutFT), OutFT(π) / OutFT(2))
    r_ft = OutFT(r)
    t_ft = OutFT(t)
    scale = OutFT(LD.scaling)
    inv_half_pi = OutFT(2) / OutFT(π)
    dist = LD.LD

    @inbounds for k in eachindex(θₗ, w)
        θk = θₗ[k]
        μ_l = cos(θk)
        wk = OutFT(pdf(dist, inv_half_pi * θk)) * scale * w[k]
        for j in axes(μꜛ, 1)
            μo = OutFT(μꜛ[j])
            hp_o = H(μo, μ_l)
            hm_o = H(-μo, μ_l)
            for i in axes(μ, 1)
                μi = OutFT(μ[i])
                hp_i = H(μi, μ_l)
                hm_i = H(-μi, μ_l)
                ψ_plus = hp_i * hp_o + hm_i * hm_o
                ψ_minus = hp_i * hm_o + hm_i * hp_o
                Γ[i, j] = muladd(wk, t_ft * ψ_plus + r_ft * ψ_minus, Γ[i, j])
            end
        end
    end
    return Γ
end
