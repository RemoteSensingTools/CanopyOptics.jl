"Fresnel reflection"
function Fᵣ(n::FT,θᵢ::FT) where FT
    #@show rad2deg(θᵢ)
    # Refractive index of air
    nᵢ = FT(1)
    # Angle of transmitted light in medium:
    θₜ = asin(nᵢ * sin(θᵢ)/n)
    #(sin(α - θₛ)^2 / sin(α + θₛ)^2 +  tan(α - θₛ)^2 / tan(α + θₛ)^2)/2
    r_perpendicular = (nᵢ * cos(θᵢ) - n * cos(θₜ)) / (nᵢ * cos(θᵢ) + n * cos(θₜ))
    r_parallel      = (nᵢ * cos(θₜ) - n * cos(θᵢ)) / (nᵢ * cos(θₜ) + n * cos(θᵢ))
    return (r_perpendicular^2 +  r_parallel^2)/2
end