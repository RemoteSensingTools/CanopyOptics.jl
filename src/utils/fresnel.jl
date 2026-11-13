"""
    fresnel_components(n, θᵢ)

Compute real-dielectric Fresnel amplitude reflectances `(r_s, r_p)` for an
air-to-medium interface (`n_air = 1`, medium refractive index `n >= 1`) at
incidence angle `θᵢ` in radians (Born & Wolf, *Principles of Optics*,
Eqs. 1.5.45–1.5.46):

```math
r_s = \\frac{\\cos\\theta_i - n\\cos\\theta_t}{\\cos\\theta_i + n\\cos\\theta_t},
\\qquad
r_p = \\frac{n\\cos\\theta_i - \\cos\\theta_t}{n\\cos\\theta_i + \\cos\\theta_t},
```

where ``\\sin\\theta_t = \\sin\\theta_i / n``.

# Examples
```jldoctest
julia> rs, rp = CanopyOptics.fresnel_components(1.5, 0.0);

julia> (round(rs; digits = 4), round(rp; digits = 4))
(-0.2, 0.2)
```
"""
@inline function fresnel_components(n::FTN, θᵢ::FTθ) where {FTN<:Real,FTθ<:Real}
    FT = float(promote_type(FTN, FTθ))
    BFT = _leaf_base_float(FT)
    n_medium = FT(n)
    θ = FT(θᵢ)
    sθ, cθ = sincos(θ)
    sin_t = sθ / n_medium
    u_t = (one(FT) - sin_t) * (one(FT) + sin_t)
    cθ_t = prospect_primal(u_t) <= zero(BFT) ? zero(FT) : sqrt(u_t)
    nc_t = n_medium * cθ_t
    nc_i = n_medium * cθ
    r_s = (cθ - nc_t) / (cθ + nc_t)
    r_p = (nc_i - cθ_t) / (nc_i + cθ_t)
    return r_s, r_p
end

"""
    Fᵣ(n, θᵢ)

Compute unpolarized Fresnel power reflectance ``(r_s^2 + r_p^2) / 2`` for an
air-to-medium interface with real refractive index ``n`` at incidence angle `θᵢ`
in radians. See [`fresnel_components`](@ref).

# Examples
```jldoctest
julia> round(CanopyOptics.Fᵣ(1.5, 0.0); digits = 4)
0.04
```
"""
@inline function Fᵣ(n::FTN, θᵢ::FTθ) where {FTN<:Real,FTθ<:Real}
    r_s, r_p = fresnel_components(n, θᵢ)
    FT = typeof(r_s)
    return (r_s * r_s + r_p * r_p) / FT(2)
end
