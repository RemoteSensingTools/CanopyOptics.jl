# Copy of RadiativeTransfer.Scattering, had problems including it otherwise

"""
    type AbstractPolarizationType
Abstract Polarization type 
"""
abstract type AbstractPolarizationType  end

"""
    struct Stokes_IQUV{FT<:AbstractFloat}
A struct which defines full Stokes Vector ([I,Q,U,V]) RT code
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_IQUV{FT<:AbstractFloat} <: AbstractPolarizationType
    "Number of Stokes components (int)"
    n::Int = 4
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT}  = [1.0, 1.0, -1.0, -1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0, 0.0, 0.0, 0.0] #assuming completely unpolarized incident stellar radiation
end

"""
    struct Stokes_IQU{FT<:AbstractFloat}
A struct which defines Stokes Vector ([I,Q,U]) RT code
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_IQU{FT<:AbstractFloat} <: AbstractPolarizationType
    "Number of Stokes components (int)" 
    n::Int = 3
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT}  = [1.0, 1.0, -1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0, 0.0, 0.0] #assuming linearly unpolarized incident stellar radiation
end

"""
    struct Stokes_I{FT<:AbstractFloat}
A struct which define scalar I only RT code
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Stokes_I{FT<:AbstractFloat} <: AbstractPolarizationType 
    "Number of Stokes components (int)"
    n::Int = 1
    "Vector of length `n` for ... (see eq in Sanghavi )"
    D::Array{FT} = [1.0]
    "Incoming Stokes vector for scalar only"
    I₀::Array{FT} = [1.0]
end
