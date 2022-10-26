"""
    dirVector{FT}
Struct for spherical coordinate directions in θ (elevation angle) and ϕ (azimuth angle)
# Fields
$(DocStringExtensions.FIELDS)
"""
struct dirVector{FT}
    θ::FT
    ϕ::FT
end

"""
    dirVector_μ{FT}
Struct for spherical coordinate directions in θ (elevation angle) and ϕ (azimuth angle)"
# Fields
$(DocStringExtensions.FIELDS)
"""
struct dirVector_μ{FT}
    μ::FT
    ϕ::FT
end 

# Define dot product for the directional vector in spherical coordinates:
LinearAlgebra.dot(Ω₁::dirVector, Ω₂::dirVector) = cos(Ω₁.θ) * cos(Ω₂.θ) + sin(Ω₁.θ) * sin(Ω₂.θ) * cos(Ω₂.ϕ - Ω₁.ϕ)

# Define dot product for the directional vector in spherical coordinates:
LinearAlgebra.dot(Ω₁::dirVector_μ, Ω₂::dirVector_μ) = Ω₁.μ * Ω₂.μ + sqrt(1-Ω₁.μ^2)  * sqrt(1-Ω₂.μ^2) * cos(Ω₂.ϕ - Ω₁.ϕ)

function dot_product(μ₁::FT, μ₂::FT, dϕ::FT) where FT
    μ₁ * μ₂ + sqrt(1-μ₁^2)  * sqrt(1-μ₂^2) * cos(dϕ)
end
