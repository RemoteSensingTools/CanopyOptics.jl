"""
    compute_Γ(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT},
              LD::AbstractLeafDistribution) where FT

Computes the azimuthally-resolved area scattering transfer function
`Γ(Ωⁱⁿ → Ωᵒᵘᵗ)` for a specific direction pair by direct double integration over
leaf polar and azimuthal angles, following Shultis & Myneni (1988) Eqs. (38)–(39):

``Γ = R\\,Γ^- + T\\,Γ^+``

where `Γ⁻` integrates the negative (reflection) part and `Γ⁺` the positive
(transmission) part of `(Ωⁱⁿ ⋅ Ωᴸ)(Ωᵒᵘᵗ ⋅ Ωᴸ)` over all leaf orientations.

Retained as a direct-quadrature reference path for [`precompute_Zazi`](@ref).
Production Z-matrix assembly uses [`compute_Z_matrices_aniso_analytic`](@ref).
See [`compute_lambertian_Γ`](@ref) for the azimuthally averaged (m=0 only)
equivalent.
"""
function compute_Γ(mod::BiLambertianCanopyScattering,
                   Ωⁱⁿ::dirVector_μ{FT},
                   Ωᵒᵘᵗ::dirVector_μ{FT},
                   LD::AbstractLeafDistribution;
                   quadrature::CanopyQuadrature = CanopyQuadrature()) where FT
    (;R,T) = mod
    n_leaf = quadrature.n_leaf
    n_azimuth = quadrature.n_azimuth
    #μ_l, w = gauleg(nQuad,0.0,1.0);
    θₗ, w = gauleg(n_leaf,FT(0),FT(π/2));
    μ_l = cos.(θₗ)
    # Quadrature points in the azimuth (this has to go over 2π):
    ϕ, w_azi = gauleg(n_azimuth,FT(0),FT(2π));
    
    Fᵢ = pdf.(LD.LD,2θₗ/π) * LD.scaling  
    
    # Create leaf angles in μ and ϕ
    Ω_l  = [dirVector_μ(a,b) for a in μ_l, b in ϕ];
    
    # Compute the angles (in, leaf angles, out, leaf)
    integrand = ((Ωⁱⁿ,) .⋅ Ω_l) .* ((Ωᵒᵘᵗ,) .⋅ Ω_l)

    # Integrate over positive and negatives separately
    iPos = (integrand+abs.(integrand))./2
    iNeg = (integrand-abs.(integrand))./2
    
    # Eq 39 in Shultis and Myneni, double integration here
    Γ⁻ = -1/2π * (Fᵢ .* w)' * (iNeg * w_azi)
    Γ⁺ =  1/2π * (Fᵢ .* w)' * (iPos * w_azi)
    # Eq 38 in Shultis and Myneni
    return R * Γ⁻ + T * Γ⁺ 
end

"""
    _psi_same(Pᵢ, Nᵢ, Pₒ, Nₒ)

Half-range Fourier product for same-sign bi-Lambertian scattering.

`P` and `N` are the one-sided leaf-projection moments for incoming (`ᵢ`) and
outgoing (`ₒ`) directions: `P` keeps the side of the leaf whose normal faces the
ray, while `N` keeps the opposite side. The leading factor of `2` matches the
normalization documented in
[`compute_Z_matrices_aniso_analytic`](@ref).
"""
@inline function _psi_same(Pᵢ::FT, Nᵢ::FT, Pₒ::FT, Nₒ::FT) where {FT}
    return FT(2) * (Pᵢ * Pₒ + Nᵢ * Nₒ)
end

"""
    _psi_opposite(Pᵢ, Nᵢ, Pₒ, Nₒ)

Half-range Fourier product for sign-changing bi-Lambertian scattering.
"""
@inline function _psi_opposite(Pᵢ::FT, Nᵢ::FT, Pₒ::FT, Nₒ::FT) where {FT}
    return FT(2) * (Pᵢ * Nₒ + Nᵢ * Pₒ)
end

"""
    _normalise_Z!(Z⁺⁺, Z⁻⁺, G, ϖ)

Apply the canopy Z-matrix convention in place.

For each incoming stream column, this divides by `ϖ * G(μ_in)` and applies the
Fourier normalization factor `2` for `m = 0` and `4` for `m > 0`. After this
normalization, conservative leaves satisfy an `m = 0` hemispheric column
integral of approximately `2`.
"""
function _normalise_Z!(Z⁺⁺::AbstractArray{FT,3},
                       Z⁻⁺::AbstractArray{FT,3},
                       G::AbstractVector,
                       ϖ::FT) where {FT}
    nμ = size(Z⁺⁺, 1)
    nm = size(Z⁺⁺, 3)
    for k in 1:nm
        ff = k == 1 ? FT(2) : FT(4)
        for j in 1:nμ
            scale = ff / (ϖ * FT(G[j]))
            @inbounds for i in 1:nμ
                Z⁺⁺[i, j, k] *= scale
                Z⁻⁺[i, j, k] *= scale
            end
        end
    end
    return Z⁺⁺, Z⁻⁺
end

"""
    _leaf_inclination_quadrature(LD, nQuad, FT)

Return leaf-inclination quadrature nodes and measure weights.

The returned `w_measure` is already the Gauss weight multiplied by the leaf
angle distribution density and scaling factor. Do not multiply by `pdf(LD, ...)`
again at the call site.
"""
function _leaf_inclination_quadrature(LD::AbstractLeafDistribution, nQuad::Int, ::Type{FT}) where {FT}
    θₗ, wθ = gauleg(nQuad, FT(0), FT(π / 2))
    Fₗ = FT.(pdf.(LD.LD, 2θₗ / FT(π))) .* FT(LD.scaling)
    return (θ = θₗ, w_measure = wθ .* Fₗ)
end

"""
    _one_sided_projection_moments!(P, N, μ::FT, μ_L::FT, m_max::Int) -> (P, N)
    _one_sided_projection_moments(μ::FT, μ_L::FT, m_max::Int) -> (P, N)

Closed-form cosine Fourier moments of the positive and negative sides of a
leaf projection for orders `m = 0:m_max`.

For a viewing direction with signed cosine `μ` and a leaf normal with
inclination cosine `μ_L`, write

```math
x(ψ) = a + b\\cos ψ, \\qquad
a = μ μ_L, \\qquad
b = \\sqrt{1 - μ^2}\\sqrt{1 - μ_L^2}.
```

Here `x(ψ)` is the signed cosine between the ray and the leaf normal as the
leaf azimuth `ψ` rotates around the canopy vertical. Positive `x` and negative
`x` are the two faces of the same flat leaf. This routine returns

```math
P_m = \\frac{1}{2π}\\int_0^{2π} \\max(x(ψ), 0)\\cos(mψ)\\,dψ,
```

and

```math
N_m = \\frac{1}{2π}\\int_0^{2π} \\max(-x(ψ), 0)\\cos(mψ)\\,dψ.
```

The older internal name used "clipped projection" because `max(x, 0)` clips
away the invisible/back-facing half of the signed projection. It is not a
clipped probability distribution.

The `m = 0` value `P_0` is the Shultis and Myneni (1988) projected-area
function `H` (their Eq. 35/46 in this codebase).  For `|a| < b`, with
`ψ* = acos(-a / b)`, the antiderivative gives

```math
P_m =
\\frac{1}{π}\\left[
\\frac{a\\sin(mψ*)}{m} +
\\frac{b}{2}\\left(
\\frac{\\sin((m-1)ψ*)}{m-1} +
\\frac{\\sin((m+1)ψ*)}{m+1}
\\right)\\right],
```

with the removable limits

```math
P_0 = \\frac{aψ* + b\\sin ψ*}{π}, \\qquad
P_1 = \\frac{a\\sin ψ* + \\frac{b}{2}ψ* + \\frac{b}{4}\\sin(2ψ*)}{π}.
```

`N_m` is derived from the exact identity

```math
P_m - N_m = a\\,δ_{m0} + \\frac{b}{2}\\,δ_{m1},
```

which follows because `max(x, 0) - max(-x, 0) = x`.
"""
function _one_sided_projection_moments!(P::AbstractVector{FT},
                                        N::AbstractVector{FT},
                                        μ::FT, μ_L::FT,
                                        m_max::Int) where {FT<:Real}
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))
    length(P) >= m_max + 1 || throw(DimensionMismatch("P must have length at least m_max + 1"))
    length(N) >= m_max + 1 || throw(DimensionMismatch("N must have length at least m_max + 1"))
    fill!(P, zero(FT))
    fill!(N, zero(FT))

    a = μ * μ_L
    b = sqrt(max(zero(FT), one(FT) - μ * μ)) *
        sqrt(max(zero(FT), one(FT) - μ_L * μ_L))

    if b == zero(FT)
        P[1] = max(a, zero(FT))
    elseif a >= b
        P[1] = a
        if m_max >= 1
            P[2] = b / FT(2)
        end
    elseif a <= -b
        # P remains zero; N is filled from P - N below.
    else
        πFT = FT(π)
        ψ = acos(clamp(-a / b, -one(FT), one(FT)))
        P[1] = (a * ψ + b * sin(ψ)) / πFT
        if m_max >= 1
            P[2] = (a * sin(ψ) + (b / FT(2)) * ψ +
                    (b / FT(4)) * sin(FT(2) * ψ)) / πFT
        end
        for m in 2:m_max
            mFT = FT(m)
            P[m + 1] = (
                a * sin(mFT * ψ) / mFT +
                (b / FT(2)) * (
                    sin(FT(m - 1) * ψ) / FT(m - 1) +
                    sin(FT(m + 1) * ψ) / FT(m + 1)
                )
            ) / πFT
        end
    end

    for m in 0:m_max
        diff = m == 0 ? a : (m == 1 ? b / FT(2) : zero(FT))
        N[m + 1] = P[m + 1] - diff
    end

    return P, N
end

function _one_sided_projection_moments(μ::FT, μ_L::FT, m_max::Int) where {FT<:Real}
    P = zeros(FT, m_max + 1)
    N = zeros(FT, m_max + 1)
    return _one_sided_projection_moments!(P, N, μ, μ_L, m_max)
end

_clipped_projection_moments!(P::AbstractVector{FT}, N::AbstractVector{FT},
                             μ::FT, μ_L::FT, m_max::Int) where {FT<:Real} =
    _one_sided_projection_moments!(P, N, μ, μ_L, m_max)

_clipped_projection_moments(μ::FT, μ_L::FT, m_max::Int) where {FT<:Real} =
    _one_sided_projection_moments(μ, μ_L, m_max)

"""
    compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                      μ::AbstractVector{FT},
                                      LD::AbstractLeafDistribution,
                                      m_max::Int;
                                      quadrature = CanopyQuadrature()) -> (Z⁺⁺, Z⁻⁺)

Compute all scalar Fourier moments `m = 0:m_max` of the bi-Lambertian
canopy phase matrices without azimuthal quadrature.

The returned arrays have shape `(length(μ), length(μ), m_max + 1)` and
use the vSmartMOM convention

```math
Z[i_{out}, j_{in}, m+1],
```

where `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is sign-change reflection.
The single-scattering albedo `ϖ = R + T` is divided out, so for
conservative leaves the `m = 0` column integral of `Z⁺⁺ + Z⁻⁺` targets
`≈ 2` under the same quadrature and `G(μ)` convention used by vSmartMOM.

# Origin trace

Shultis and Myneni (1988), Eq. 45 writes the azimuthally averaged
Lambertian canopy kernels as products of one-direction projected-area
functions. For each leaf inclination this implementation generalizes that
factorization from `m = 0` to arbitrary cosine Fourier order by using the
closed-form moments returned by `_one_sided_projection_moments`.

For signed incoming/outgoing directions, let `(Pᵢ, Nᵢ)` and `(Pₒ, Nₒ)` be
the positive-face and negative-face projection moments. The circular-convolution
theorem gives the per-leaf Fourier coefficients

```math
Ψ^+_m = P_{i,m}P_{o,m} + N_{i,m}N_{o,m},
\\qquad
Ψ^-_m = P_{i,m}N_{o,m} + N_{i,m}P_{o,m}.
```

The code multiplies these products by `2` before applying the existing
`f_0 = 2`, `f_{m>0} = 4` factors.  That leading `2` is not leaf physics; it
converts the `1/(2π)` Fourier coefficients above to the half-range cosine
moment used by the historical `compute_Z_matrices_aniso` implementation:
`(2/π)∫_0^π Γ(Δϕ)cos(mΔϕ)dΔϕ`.  This is the normalization expected by
vSmartMOM's elemental kernels, where `ϖ` is applied outside `Z`.
"""
function compute_Z_matrices_aniso_analytic(mod::BiLambertianCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int;
                                           quadrature::CanopyQuadrature = CanopyQuadrature(),
                                           nQuad = nothing) where {FT<:Real}
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))

    (; R, T) = mod
    q = _resolve_quadrature(quadrature, nQuad)
    ZFT = promote_type(FT, typeof(R), typeof(T))
    R_leaf = ZFT(R)
    T_leaf = ZFT(T)
    ϖ = R_leaf + T_leaf

    nμ = length(μ)
    nm = m_max + 1
    Z⁺⁺ = zeros(ZFT, nμ, nμ, nm)
    Z⁻⁺ = zeros(ZFT, nμ, nμ, nm)
    ϖ <= zero(ZFT) && return Z⁺⁺, Z⁻⁺

    μ_vec = collect(μ)
    μ_work = ZFT.(μ_vec)
    leaf_quad = _leaf_inclination_quadrature(LD, q.n_leaf, FT)
    θₗ = leaf_quad.θ
    w_measure = leaf_quad.w_measure
    G = ZFT.(vec(CanopyOptics.G(μ_vec, LD)))

    Pꜜ = zeros(ZFT, nμ, nm)
    Nꜜ = zeros(ZFT, nμ, nm)
    Pꜛ = zeros(ZFT, nμ, nm)
    Nꜛ = zeros(ZFT, nμ, nm)
    P = zeros(ZFT, nm)
    N = zeros(ZFT, nm)
    Pm = zeros(ZFT, nm)
    Nm = zeros(ZFT, nm)

    for l in eachindex(θₗ)
        μ_L = ZFT(cos(θₗ[l]))
        leaf_weight = ZFT(w_measure[l])

        for i in eachindex(μ_work)
            μ_i = μ_work[i]
            _one_sided_projection_moments!(P, N, μ_i, μ_L, m_max)
            _one_sided_projection_moments!(Pm, Nm, -μ_i, μ_L, m_max)
            @inbounds for k in 1:nm
                Pꜜ[i, k] = P[k]
                Nꜜ[i, k] = N[k]
                Pꜛ[i, k] = Pm[k]
                Nꜛ[i, k] = Nm[k]
            end
        end

        for k in 1:nm, j in 1:nμ
            @inbounds begin
                Pj = Pꜜ[j, k]
                Nj = Nꜜ[j, k]
                for i in 1:nμ
                    Ψpp_same = _psi_same(Pj, Nj, Pꜜ[i, k], Nꜜ[i, k])
                    Ψpp_opp  = _psi_opposite(Pj, Nj, Pꜜ[i, k], Nꜜ[i, k])
                    Ψmp_same = _psi_same(Pj, Nj, Pꜛ[i, k], Nꜛ[i, k])
                    Ψmp_opp  = _psi_opposite(Pj, Nj, Pꜛ[i, k], Nꜛ[i, k])

                    Z⁺⁺[i, j, k] += leaf_weight * (T_leaf * Ψpp_same + R_leaf * Ψpp_opp)
                    Z⁻⁺[i, j, k] += leaf_weight * (T_leaf * Ψmp_same + R_leaf * Ψmp_opp)
                end
            end
        end
    end

    return _normalise_Z!(Z⁺⁺, Z⁻⁺, G, ϖ)
end
