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

An assertion verifies `Ωⁱⁿ ⋅ Ωᴸ ≈ Ωᵒᵘᵗ ⋅ Ωᴸ` (equal angle of incidence and exitance).

- The `dirVector` variant returns only `Ωᴸ`.
- The `dirVector_μ` variant returns `(Ωᴸ, α)` for direct use in [`compute_reflection`](@ref).
"""
function getSpecularΩ(Ωⁱⁿ::dirVector{FT}, Ωᵒᵘᵗ::dirVector{FT}) where FT
    sa = Ωⁱⁿ ⋅ Ωᵒᵘᵗ 
    sa > 1 ? sa = FT(1) : nothing
    
    α = acos(sa)/2
    # Relative azimuth angle:
    ϕ = (Ωᵒᵘᵗ.ϕ - Ωⁱⁿ.ϕ) 
    # Leaf polar angle:
    θₗ = acos((cos(Ωⁱⁿ.θ) + cos(Ωᵒᵘᵗ.θ))/2cos(α));
    
    t = (sin(Ωⁱⁿ.θ) + sin(Ωᵒᵘᵗ.θ) * cos(ϕ))/(2cos(α)*sin(θₗ))
    if isnan(t)
        t = FT(1.0)
    end
    # Leaf azimuth angle:
    ϕₗ = acos(max(FT(-1),min(t,FT(1))))
    ϕ > π   ? ϕₗ = 2π-ϕₗ-Ωⁱⁿ.ϕ : ϕₗ = ϕₗ+Ωⁱⁿ.ϕ
    #θₗ > π/2 ? θₗ = θₗ - π/2  : nothing 
    out = dirVector(θₗ,ϕₗ)
    #ϕ > π   ? out = dirVector(θₗ,2π-ϕₗ-Ωⁱⁿ.ϕ) : out = dirVector(θₗ,ϕₗ+Ωⁱⁿ.ϕ)
    #@show (Ωⁱⁿ ⋅ out) , (Ωᵒᵘᵗ ⋅ out)
    #@show abs(((Ωⁱⁿ ⋅ out) - (Ωᵒᵘᵗ ⋅ out)) / (Ωᵒᵘᵗ ⋅ out)) * 100
    #@show θₗ , ϕₗ 
    @assert abs(((Ωⁱⁿ ⋅ out) - (Ωᵒᵘᵗ ⋅ out)) / (Ωᵒᵘᵗ ⋅ out)) * 100 < 1 "Leaf angle wrong, $Ωⁱⁿ, $Ωᵒᵘᵗ, $out"
    return out 

end

function getSpecularΩ(Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT
    sa = Ωⁱⁿ ⋅ Ωᵒᵘᵗ 
    #@show sa
    sa > 1 ? sa = FT(1) : nothing
    # Scattering angle
    α = acos(sa)/FT(2)
    
    # Relative azimuth angle:
    ϕ = (Ωᵒᵘᵗ.ϕ - Ωⁱⁿ.ϕ) 
    # Leaf polar angle:
    μₗ = ((Ωⁱⁿ.μ + Ωᵒᵘᵗ.μ)/2cos(α));
    
    t = (sqrt(1-Ωⁱⁿ.μ^2) + sqrt(1-Ωᵒᵘᵗ.μ^2) * cos(ϕ))/(2cos(α)*sqrt(1-μₗ^2))
    #@show t, sqrt(1-Ωⁱⁿ.μ^2) + sqrt(1-Ωᵒᵘᵗ.μ)
    if isnan(t)
        t = FT(1.0)
    end
    # Leaf azimuth angle:
    ϕₗ = acos(max(FT(-1),min(t,FT(1))))
    ϕ > π   ? ϕₗ = 2π-ϕₗ-Ωⁱⁿ.ϕ : ϕₗ = ϕₗ+Ωⁱⁿ.ϕ
    out = dirVector_μ(μₗ,ϕₗ)
    @assert Ωⁱⁿ ⋅ out ≈  Ωᵒᵘᵗ ⋅ out "Leaf angle wrong in specular reflection, $Ωⁱⁿ, $Ωᵒᵘᵗ, $out"
    return out, α
end