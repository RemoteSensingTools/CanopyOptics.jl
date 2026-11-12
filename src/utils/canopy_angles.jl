"""
    βparameters(θₗ, s)

Compute `Beta(α, β)` shape parameters `(α, β)` on `[0, 1]` from mean leaf
inclination angle `θₗ` and standard deviation `s` (both in radians), following
Bonan (2019, Table 2.1, Eqs. 2.9–2.10):

```math
\\alpha = \\frac{1 - (s^2 + \\theta_l^2)/(\\theta_l \\pi / 2)}{(s^2 + \\theta_l^2)/\\theta_l^2 - 1},
\\qquad
\\beta = \\left(\\frac{\\pi / 2}{\\theta_l} - 1\\right) \\alpha .
```

# Examples
```jldoctest
julia> α, β = βparameters(deg2rad(45.0), deg2rad(16.27));

julia> round.((α, β); digits = 3)
(3.325, 3.325)
```
"""
function βparameters(θₗ::Real, s::Real)
    FT = float(promote_type(typeof(θₗ), typeof(s)))
    θ = FT(θₗ)
    σ = FT(s)
    half_pi = FT(π) / FT(2)
    m2 = σ * σ + θ * θ
    α = (one(FT) - m2 / (θ * half_pi)) / (m2 / (θ * θ) - one(FT))
    β = (half_pi / θ - one(FT)) * α
    return α, β
end

"""
    H(μ, μₗ)

Evaluate the one-sided azimuthally averaged leaf projection factor `H(μ, μₗ)`
for stream polar cosine `μ = cos θ ∈ [-1, 1]` and leaf-normal cosine
`μₗ = cos θₗ ∈ [-1, 1]` (Shultis & Myneni, 1988, Eqs. 46–47):

```math
H(\\mu, \\mu_L) = \\frac{1}{2\\pi} \\int_0^{2\\pi} \\max\\!\\big(0,\\; \\mu\\mu_L + \\sqrt{1-\\mu^2}\\sqrt{1-\\mu_L^2}\\cos\\phi\\big)\\,d\\phi .
```

Evaluated algebraically without `cot(acos(⋅))` singularities so boundary inputs
(`|μ| == 1` or `μₗ == 0`) remain finite and `ForwardDiff`-compatible. Assembles
the azimuthally averaged `Ψ` kernels in [`compute_Ψ`](@ref).

# Examples
```jldoctest
julia> CanopyOptics.H(1.0, 0.0)
0.0
```
"""
@inline function H(μ_in::FT1, μₗ_in::FT2) where {FT1<:Real,FT2<:Real}
    FT = float(promote_type(FT1, FT2))
    BFT = _leaf_base_float(FT)
    μ = clamp(FT(μ_in), -one(FT), one(FT))
    μₗ = clamp(FT(μₗ_in), -one(FT), one(FT))
    s = _unit_sin_from_cos(μ)
    sₗ = _unit_sin_from_cos(μₗ)
    a = μ * μₗ
    b = s * sₗ
    a_p = prospect_primal(a)
    b_p = prospect_primal(b)
    if a_p >= b_p
        return a
    elseif a_p <= -b_p
        return zero(FT)
    else
        r = clamp(-a / b, -one(FT), one(FT))
        ϕ = acos(r)
        rad = (b - a) * (b + a)
        term_sin = prospect_primal(rad) <= zero(BFT) ? zero(FT) : sqrt(rad)
        return (a * ϕ + term_sin) / FT(π)
    end
end

"""
    compute_Ψ(μ::AbstractVector, μₗ::Real)
    compute_Ψ(μ::AbstractVector, μꜛ::AbstractVector, μₗ::Real)

Compute the azimuthally averaged `(Ψ⁺, Ψ⁻)` coupling matrices at leaf-normal
cosine `μₗ` following Shultis & Myneni (1988, Eq. 45):

```math
\\Psi^+(\\mu, \\mu') = H(\\mu, \\mu_L) H(\\mu', \\mu_L) + H(-\\mu, \\mu_L) H(-\\mu', \\mu_L),
```
```math
\\Psi^-(\\mu, \\mu') = H(\\mu, \\mu_L) H(-\\mu', \\mu_L) + H(-\\mu, \\mu_L) H(\\mu', \\mu_L).
```

Returns `(Ψ⁺, Ψ⁻)` of size `(length(μ), length(μ'))`, preserving the index axes
of `μ` and `μꜛ`. The two-argument method sets `μꜛ = μ`.
"""
function compute_Ψ(μ::AbstractVector, μₗ::Real)
    FT = float(promote_type(eltype(μ), typeof(μₗ)))
    μₗ_ft = FT(μₗ)
    H⁺ = H.(FT.(μ), μₗ_ft)
    H⁻ = H.(-FT.(μ), μₗ_ft)
    Ψ⁺ = H⁺ * H⁺' + H⁻ * H⁻'
    Ψ⁻ = H⁺ * H⁻' + H⁻ * H⁺'
    return Ψ⁺, Ψ⁻
end

function compute_Ψ(μ::AbstractVector, μꜛ::AbstractVector, μₗ::Real)
    FT = float(promote_type(eltype(μ), eltype(μꜛ), typeof(μₗ)))
    μₗ_ft = FT(μₗ)
    H₁⁺ = H.(FT.(μ), μₗ_ft)
    H₁⁻ = H.(-FT.(μ), μₗ_ft)
    H₂⁺ = H.(FT.(μꜛ), μₗ_ft)
    H₂⁻ = H.(-FT.(μꜛ), μₗ_ft)
    Ψ⁺ = H₁⁺ * H₂⁺' + H₁⁻ * H₂⁻'
    Ψ⁻ = H₁⁺ * H₂⁻' + H₁⁻ * H₂⁺'
    return Ψ⁺, Ψ⁻
end

"""
    scattering_model_lambertian(Ωⁱⁿ::DirectionCosine, Ωᵒᵘᵗ::DirectionCosine, Ωᴸ::DirectionCosine, r, t)

Evaluate the single-leaf bi-Lambertian phase function `γ(Ωⁱⁿ → Ωᵒᵘᵗ; Ωᴸ)` for
leaf reflectance `r` and transmittance `t` (Shultis & Myneni, 1988, Eq. 37):

```math
\\gamma(\\Omega' \\to \\Omega;\\,\\Omega_L) = \\begin{cases}
  r\\,|\\Omega \\cdot \\Omega_L| / \\pi & (\\Omega' \\cdot \\Omega_L)(\\Omega \\cdot \\Omega_L) < 0, \\\\
  t\\,|\\Omega \\cdot \\Omega_L| / \\pi & (\\Omega' \\cdot \\Omega_L)(\\Omega \\cdot \\Omega_L) \\ge 0 .
\\end{cases}
```
"""
function scattering_model_lambertian(
    Ωⁱⁿ::DirectionCosine,
    Ωᵒᵘᵗ::DirectionCosine,
    Ωᴸ::DirectionCosine,
    r::Real,
    t::Real,
)
    d_in = Ωⁱⁿ ⋅ Ωᴸ
    d_out = Ωᵒᵘᵗ ⋅ Ωᴸ
    FT = float(promote_type(typeof(d_in), typeof(d_out), typeof(r), typeof(t)))
    coeff = prospect_primal(d_in * d_out) < 0 ? FT(r) : FT(t)
    return coeff * abs(FT(d_out)) / FT(π)
end

@inline function _azimuth_0_2pi(
    y::FT,
    x::FT,
    ::Type{FT},
    ::Type{BFT},
) where {FT<:Real,BFT<:AbstractFloat}
    if iszero(prospect_primal(x)) && iszero(prospect_primal(y))
        return zero(FT)
    end
    θ = atan(y, x)
    return prospect_primal(θ) < zero(BFT) ? θ + FT(2) * FT(π) : θ
end

"""
    getSpecularΩ(Ωⁱⁿ::DirectionCosine, Ωᵒᵘᵗ::DirectionCosine) -> (Ωᴸ, α)

Compute the unit specular leaf normal `Ωᴸ` bisecting incoming source direction
`Ωⁱⁿ` (`-k_in`) and outgoing ray direction `Ωᵒᵘᵗ` (`+k_out`), together with the
local incidence half-angle `α = arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ) / 2` (Knyazikhin & Marshak,
1991, Eq. 18).

At the antipodal forward-grazing limit (`Ωᵒᵘᵗ ≈ -Ωⁱⁿ`, `α = π/2`), returns a
deterministic unit normal orthogonal to `Ωⁱⁿ` without coordinate singularity or
`ForwardDiff.Dual` `NaN` tangents.
"""
function getSpecularΩ(Ωⁱⁿ::DirectionCosine{FT1}, Ωᵒᵘᵗ::DirectionCosine{FT2}) where {FT1<:Real,FT2<:Real}
    FT = _canopy_parameter_type(Ωⁱⁿ.μ, Ωᵒᵘᵗ.μ)
    BFT = _leaf_base_float(FT)
    μ1 = clamp(FT(Ωⁱⁿ.μ), -one(FT), one(FT))
    μ2 = clamp(FT(Ωᵒᵘᵗ.μ), -one(FT), one(FT))
    ϕ1 = FT(Ωⁱⁿ.ϕ)
    ϕ2 = FT(Ωᵒᵘᵗ.ϕ)

    u1 = one(FT) - μ1 * μ1
    u2 = one(FT) - μ2 * μ2
    s1 = prospect_primal(u1) <= zero(BFT) ? zero(FT) : sqrt(u1)
    s2 = prospect_primal(u2) <= zero(BFT) ? zero(FT) : sqrt(u2)

    c1, sp1 = cos(ϕ1), sin(ϕ1)
    c2, sp2 = cos(ϕ2), sin(ϕ2)

    ω1x, ω1y, ω1z = s1 * c1, s1 * sp1, μ1
    ω2x, ω2y, ω2z = s2 * c2, s2 * sp2, μ2

    Hx, Hy, Hz = ω1x + ω2x, ω1y + ω2y, ω1z + ω2z
    Dx, Dy, Dz = ω1x - ω2x, ω1y - ω2y, ω1z - ω2z

    H2 = Hx * Hx + Hy * Hy + Hz * Hz
    D2 = Dx * Dx + Dy * Dy + Dz * Dz
    Hnorm = prospect_primal(H2) <= zero(BFT) ? zero(FT) : sqrt(H2)
    Dnorm = prospect_primal(D2) <= zero(BFT) ? zero(FT) : sqrt(D2)
    α = atan(Dnorm, Hnorm)

    if prospect_primal(H2) > (BFT(16) * eps(BFT))^2
        μ_L = clamp(Hz / Hnorm, -one(FT), one(FT))
        ϕ_L = _azimuth_0_2pi(Hy, Hx, FT, BFT)
        return DirectionCosine{FT}(μ_L, ϕ_L), α
    end

    # Antipodal grazing limit (Ωᵒᵘᵗ ≈ -Ωⁱⁿ, α ≈ π/2): deterministic unit normal ⊥ ω1
    nx = μ1 * c1
    ny = μ1 * sp1
    nz = -s1
    μ_L = clamp(nz, -one(FT), one(FT))
    ϕ_L = _azimuth_0_2pi(ny, nx, FT, BFT)
    return DirectionCosine{FT}(μ_L, ϕ_L), α
end
