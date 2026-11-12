# Standalone Stokes configuration types (not included by default module load).

"""
    AbstractPolarizationType

Supertype for Stokes-vector configuration containers.
"""
abstract type AbstractPolarizationType end

"""
    Stokes_IQUV{FT<:AbstractFloat}(; n = 4, D = FT[1, 1, -1, -1], I₀ = FT[1, 0, 0, 0])
    Stokes_IQUV()

Configure 4-component `[I, Q, U, V]` polarized radiative transfer.
"""
struct Stokes_IQUV{FT<:AbstractFloat} <: AbstractPolarizationType
    n::Int
    D::Vector{FT}
    I₀::Vector{FT}
end

Stokes_IQUV{FT}(; n::Int = 4, D::Vector{FT} = FT[1, 1, -1, -1], I₀::Vector{FT} = FT[1, 0, 0, 0]) where {FT<:AbstractFloat} =
    Stokes_IQUV{FT}(n, D, I₀)
Stokes_IQUV(; kwargs...) = Stokes_IQUV{Float64}(; kwargs...)

"""
    Stokes_IQU{FT<:AbstractFloat}(; n = 3, D = FT[1, 1, -1], I₀ = FT[1, 0, 0])
    Stokes_IQU()

Configure 3-component `[I, Q, U]` linearly polarized radiative transfer.
"""
struct Stokes_IQU{FT<:AbstractFloat} <: AbstractPolarizationType
    n::Int
    D::Vector{FT}
    I₀::Vector{FT}
end

Stokes_IQU{FT}(; n::Int = 3, D::Vector{FT} = FT[1, 1, -1], I₀::Vector{FT} = FT[1, 0, 0]) where {FT<:AbstractFloat} =
    Stokes_IQU{FT}(n, D, I₀)
Stokes_IQU(; kwargs...) = Stokes_IQU{Float64}(; kwargs...)

"""
    Stokes_I{FT<:AbstractFloat}(; n = 1, D = FT[1], I₀ = FT[1])
    Stokes_I()

Configure scalar intensity-only `[I]` radiative transfer.
"""
struct Stokes_I{FT<:AbstractFloat} <: AbstractPolarizationType
    n::Int
    D::Vector{FT}
    I₀::Vector{FT}
end

Stokes_I{FT}(; n::Int = 1, D::Vector{FT} = FT[1], I₀::Vector{FT} = FT[1]) where {FT<:AbstractFloat} =
    Stokes_I{FT}(n, D, I₀)
Stokes_I(; kwargs...) = Stokes_I{Float64}(; kwargs...)
