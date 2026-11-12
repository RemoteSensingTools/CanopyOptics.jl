"""
    $(FUNCTIONNAME)(θₗ, s)

Returns Beta distribution parameters `(α, β)` from the mean leaf inclination angle `θₗ`
and its standard deviation `s` (both in radians), following Bonan Eqs. 2.9–2.10:

``α = \\frac{1 - (s^2 + θₗ^2)/(θₗ π/2)}{(s^2 + θₗ^2)/θₗ^2 - 1}``

``β = \\left(\\frac{π/2}{θₗ} - 1\\right) α``

The resulting `Beta(α, β)` distribution is defined on `[0, 1]`, which maps to leaf
zenith angles `[0, π/2]` via `θ_L = (π/2) u` for `u ~ Beta(α, β)`.
"""
function βparameters(θₗ,s)
    α = (1 - (s^2+θₗ^2) / (θₗ*π/2) ) / ( ((s^2+θₗ^2) / θₗ^2) - 1 )
    β = ( (π/2)/θₗ -1 ) * α
    return α,β
end

"""
    $(FUNCTIONNAME)(μ::FT, μₗ::FT) where FT

Computes the conditional scattering function `H(μ, μₗ)` for direction cosine `μ = cos θ`
and leaf inclination cosine `μₗ = cos θₗ`, following Shultis & Myneni (1988) Eq. (46).

Evaluated piecewise via Eq. (47): let `r = cot θ · cot θₗ`, then

``H(μ, μₗ) = \\begin{cases}
  μ μₗ & r > 1 \\\\
  0    & r < -1 \\\\
  \\frac{1}{π}\\!\\left[μ μₗ \\,ϕ_0 + \\sin θ \\sin θₗ \\sin ϕ_0\\right] & \\text{otherwise}
\\end{cases}``

where `ϕ_0 = arccos(−r)`.  Used to build the `Ψ` matrices via [`compute_Ψ`](@ref).
"""
function H(μ::FT,μₗ::FT) where FT
    # cot(θ) * cot(θₗ)
    r = cot(acos(μ)) * cot(acos(μₗ));
    if r > 1
        return μ * μₗ
    elseif r < -1
        return FT(0)
    else
        # Eq 47
        ϕ = acos(-r);
        return FT(1/π) * (μ*μₗ*ϕ + sqrt(FT(1)-μ^2) * sqrt(1-μₗ^2) * sin(ϕ)) 
    end
end

"""
    $(FUNCTIONNAME)(μ::AbstractArray{FT,1}, μₗ::FT) where FT
    $(FUNCTIONNAME)(μ::AbstractArray{FT,1}, μꜛ::AbstractArray{FT,1}, μₗ::FT) where FT

Computes the azimuthally-averaged Ψ matrices for a single leaf inclination cosine `μₗ`,
following Shultis & Myneni (1988) Eq. (45):

``Ψ^+(μ, μ') = H(μ, μₗ)\\,H(μ', μₗ) + H(-μ, μₗ)\\,H(-μ', μₗ)``

``Ψ^-(μ, μ') = H(μ, μₗ)\\,H(-μ', μₗ) + H(-μ, μₗ)\\,H(μ', μₗ)``

Returns `(Ψ⁺, Ψ⁻)` as matrices of size `(length(μ), length(μ'))`.
`Ψ⁺` couples same-side directions (transmission-like); `Ψ⁻` couples opposite-side
directions (reflection-like).  The one-array form uses `μ' = μ`.
"""
function compute_Ψ(μ::AbstractArray{FT,1}, μₗ::FT) where FT
    H⁺ = H.(μ,μₗ)
    H⁻ = H.(-μ,μₗ)
    Ψ⁺ = H⁺ * H⁺' + H⁻ * H⁻'
    Ψ⁻ = H⁺ * H⁻' + H⁻ * H⁺'
    return Ψ⁺, Ψ⁻
end

"""
    $(FUNCTIONNAME)(μ::AbstractArray{FT,1}, μꜛ::AbstractArray{FT,1}, μₗ::FT) where FT

Two-array form of [`compute_Ψ`](@ref): allows distinct incoming (`μ`) and outgoing
(`μꜛ`) quadrature grids.  See the one-array docstring for the full equations.
"""
function compute_Ψ(μ::AbstractArray{FT,1},μꜛ::AbstractArray{FT,1}, μₗ::FT) where FT
    Ψ⁺ = H.(μ,μₗ) * H.(μꜛ,μₗ)'  + H.(-μ,μₗ) * H.(-μꜛ,μₗ)'
    Ψ⁻ = H.(μ,μₗ) * H.(-μꜛ,μₗ)' + H.(-μ,μₗ) * H.(+μꜛ,μₗ)'
    return Ψ⁺, Ψ⁻
end

"""
    $(FUNCTIONNAME)(Ωⁱⁿ::dirVector, Ωᵒᵘᵗ::dirVector, Ωᴸ::dirVector, r, t)

Returns the single-leaf bi-Lambertian phase function for incoming direction `Ωⁱⁿ`,
outgoing direction `Ωᵒᵘᵗ`, leaf normal `Ωᴸ`, leaf reflectance `r`, and transmittance
`t`, following Shultis & Myneni (1988) Eq. (37):

``γ(Ω' \\to Ω;\\,Ωᴸ) = \\begin{cases}
  r\\,|Ω \\cdot Ωᴸ| / π & (Ω' \\cdot Ωᴸ)(Ω \\cdot Ωᴸ) < 0 \\quad (\\text{reflection}) \\\\
  t\\,|Ω \\cdot Ωᴸ| / π & (Ω' \\cdot Ωᴸ)(Ω \\cdot Ωᴸ) > 0 \\quad (\\text{transmission})
\\end{cases}``
"""
function scattering_model_lambertian(Ωⁱⁿ::dirVector, Ωᵒᵘᵗ::dirVector, Ωᴸ::dirVector, r, t)
    α = (Ωⁱⁿ ⋅ Ωᴸ) * (Ωᵒᵘᵗ ⋅ Ωᴸ)
    if α<0
        return r * abs(Ωᵒᵘᵗ ⋅ Ωᴸ) / π
    else
        return t * abs(Ωᵒᵘᵗ ⋅ Ωᴸ) / π
    end
end

"""
    $(FUNCTIONNAME)(Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}) where FT
    $(FUNCTIONNAME)(Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT

Returns the specular leaf normal direction `Ωᴸ` that satisfies the specular reflection
condition for incoming direction `Ωⁱⁿ` and outgoing direction `Ωᵒᵘᵗ`.

The leaf normal bisects `Ωⁱⁿ` and `Ωᵒᵘᵗ`, i.e. it lies along `Ωⁱⁿ + Ωᵒᵘᵗ` (normalised).
The incidence half-angle is `α = arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ) / 2`.
See Eq. 18 in Knyazikhin & Marshak, *Discrete Ordinates Method for Photon Transport
in Leaf Canopies*.

Here `Ωⁱⁿ` points toward the illumination source (`-k_in`) and `Ωᵒᵘᵗ` points toward
the outgoing ray (`+k_out`), so the specular leaf normal bisects `(Ωⁱⁿ, Ωᵒᵘᵗ)`
and satisfies `Ωⁱⁿ ⋅ Ωᴸ ≈ Ωᵒᵘᵗ ⋅ Ωᴸ`. At the antipodal grazing limit (`Ωᵒᵘᵗ ≈ -Ωⁱⁿ`,
`α = π/2`), returns a deterministic unit leaf normal orthogonal to `Ωⁱⁿ` without
singularity or `ForwardDiff.Dual` `NaN` tangents.

- The `dirVector` variant returns only `Ωᴸ`.
- The `dirVector_μ` variant returns `(Ωᴸ, α)` for direct use in [`compute_reflection`](@ref).
"""
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

function getSpecularΩ(Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}) where FT
    BFT = _leaf_base_float(FT)
    in_μ = dirVector_μ(cos(Ωⁱⁿ.θ), Ωⁱⁿ.ϕ)
    out_μ = dirVector_μ(cos(Ωᵒᵘᵗ.θ), Ωᵒᵘᵗ.ϕ)
    Ωᴸ_μ, _ = getSpecularΩ(in_μ, out_μ)
    u_L = one(FT) - Ωᴸ_μ.μ * Ωᴸ_μ.μ
    s_L = prospect_primal(u_L) <= zero(BFT) ? zero(FT) : sqrt(u_L)
    θ_L = atan(s_L, Ωᴸ_μ.μ)
    return dirVector(θ_L, Ωᴸ_μ.ϕ)
end

function getSpecularΩ(Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT
    BFT = _leaf_base_float(FT)
    μ1 = clamp(Ωⁱⁿ.μ, -one(FT), one(FT))
    μ2 = clamp(Ωᵒᵘᵗ.μ, -one(FT), one(FT))

    u1 = one(FT) - μ1 * μ1
    u2 = one(FT) - μ2 * μ2
    s1 = prospect_primal(u1) <= zero(BFT) ? zero(FT) : sqrt(u1)
    s2 = prospect_primal(u2) <= zero(BFT) ? zero(FT) : sqrt(u2)

    c1, sp1 = cos(Ωⁱⁿ.ϕ), sin(Ωⁱⁿ.ϕ)
    c2, sp2 = cos(Ωᵒᵘᵗ.ϕ), sin(Ωᵒᵘᵗ.ϕ)

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
        return dirVector_μ(μ_L, ϕ_L), α
    end

    # Antipodal grazing limit (Ωᵒᵘᵗ ≈ -Ωⁱⁿ, α ≈ π/2): deterministic unit normal ⊥ ω1
    nx = μ1 * c1
    ny = μ1 * sp1
    nz = -s1
    μ_L = clamp(nz, -one(FT), one(FT))
    ϕ_L = _azimuth_0_2pi(ny, nx, FT, BFT)
    return dirVector_μ(μ_L, ϕ_L), α
end
