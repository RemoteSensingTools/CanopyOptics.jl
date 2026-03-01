"""
    fresnel_components(n::FT, θᵢ::FT) where FT

Returns the Fresnel amplitude reflectances `(r_s, r_p)` for an air–medium interface
(refractive index of air = 1) at incidence angle `θᵢ` (radians).

- `r_s`: s-polarization (TE, perpendicular) amplitude reflectance
- `r_p`: p-polarization (TM, parallel) amplitude reflectance

The refracted angle is computed via Snell's law: `sin θₜ = sin θᵢ / n`.

Fresnel equations (see Born & Wolf, *Principles of Optics*, Eqs. 1.5.45–1.5.46):

``r_s = \\frac{\\cos θᵢ - n \\cos θₜ}{\\cos θᵢ + n \\cos θₜ}``

``r_p = \\frac{n \\cos θᵢ - \\cos θₜ}{n \\cos θᵢ + \\cos θₜ}``

Note: `n` must be real. For full Stokes propagation in vSmartMOM.jl, use the
individual components to form the Mueller reflection matrix:

``M = \\begin{pmatrix} (r_s^2+r_p^2)/2 & (r_s^2-r_p^2)/2 & 0 & 0 \\\\
                       (r_s^2-r_p^2)/2 & (r_s^2+r_p^2)/2 & 0 & 0 \\\\
                       0 & 0 & r_s r_p & 0 \\\\
                       0 & 0 & 0 & r_s r_p \\end{pmatrix}``
"""
function fresnel_components(n::FT, θᵢ::FT) where FT
    nᵢ = FT(1)
    θₜ = asin(nᵢ * sin(θᵢ) / n)
    r_s = (nᵢ * cos(θᵢ) - n * cos(θₜ)) / (nᵢ * cos(θᵢ) + n * cos(θₜ))
    r_p = (n  * cos(θᵢ) - nᵢ * cos(θₜ)) / (n  * cos(θᵢ) + nᵢ * cos(θₜ))
    return r_s, r_p
end

"""
    Fᵣ(n::FT, θᵢ::FT) where FT

Returns the unpolarized Fresnel reflectance `(r_s² + r_p²) / 2` for an air–medium
interface at incidence angle `θᵢ` (radians), where `n` is the real refractive index
of the medium (refractive index of air = 1).

Calls [`fresnel_components`](@ref) internally. For polarization-resolved output
(needed for full Stokes vector propagation in vSmartMOM.jl), use `fresnel_components`
directly to obtain `r_s` and `r_p` separately.
"""
function Fᵣ(n::FT, θᵢ::FT) where FT
    r_s, r_p = fresnel_components(n, θᵢ)
    return (r_s^2 + r_p^2) / 2
end
