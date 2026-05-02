"""
    AbstractWoodReflectance{FT}

Spectral reflectance model for bark or woody canopy elements.

Wood reflectance is kept separate from angular scattering. A reflectance model
answers "what is the bark albedo at this wavelength or wavenumber?", while a
canopy scattering model answers "how does that woody surface redirect light?".
"""
abstract type AbstractWoodReflectance{FT<:Real} end

"""
    AbstractLUTWoodReflectance{FT}

Base type for wood reflectance models backed by tabulated spectral data.
"""
abstract type AbstractLUTWoodReflectance{FT<:Real} <: AbstractWoodReflectance{FT} end

_check_spectral_grid_unit(grid_unit::Symbol) =
    grid_unit in (:nm, :cm_inv) ? grid_unit :
    throw(ArgumentError("grid_unit must be :nm or :cm_inv"))

_wood_parameter_type(args...) = float(promote_type(map(typeof, args)...))

function _spectral_coordinate(x::Real, from_unit::Symbol, to_unit::Symbol)
    _check_spectral_grid_unit(from_unit)
    _check_spectral_grid_unit(to_unit)
    x_model = x
    if from_unit != to_unit
        x > zero(x) || throw(ArgumentError("spectral coordinates must be positive for nm <-> cm_inv conversion"))
        x_model = oftype(x / one(x), 1e7 / x)
    end
    return x_model
end

"""
    ConstantWoodReflectance(R)
    ConstantWoodReflectance(; R = 0.25)

Wavelength-independent bark/wood reflectance.
"""
struct ConstantWoodReflectance{FT<:Real} <: AbstractWoodReflectance{FT}
    R::FT
end

ConstantWoodReflectance(; R = 0.25) = ConstantWoodReflectance(R)

ConstantWoodReflectance{FT}(; R = FT(0.25)) where {FT<:Real} =
    ConstantWoodReflectance{FT}(FT(R))

"""
    LUTWoodReflectance(grid, R; grid_unit = :nm, extrapolation = :clamp)

Piecewise-linear bark/wood reflectance lookup table.

`grid` and `R` must have the same length, with at least two points. The
constructor sorts the table by spectral coordinate. `grid_unit` is `:nm` for
wavelength in nanometers or `:cm_inv` for wavenumber. `extrapolation = :clamp`
uses the nearest endpoint outside the tabulated range; `:error` throws.
"""
struct LUTWoodReflectance{FT<:Real} <: AbstractLUTWoodReflectance{FT}
    grid::Vector{FT}
    R::Vector{FT}
    grid_unit::Symbol
    extrapolation::Symbol
end

function LUTWoodReflectance(grid::AbstractVector{<:Real},
                            R::AbstractVector{<:Real};
                            grid_unit::Symbol = :nm,
                            extrapolation::Symbol = :clamp)
    length(grid) == length(R) ||
        throw(DimensionMismatch("grid and R must have the same length"))
    length(grid) >= 2 ||
        throw(ArgumentError("LUTWoodReflectance needs at least two spectral points"))
    _check_spectral_grid_unit(grid_unit)
    extrapolation in (:clamp, :error) ||
        throw(ArgumentError("extrapolation must be :clamp or :error"))

    FT = _wood_parameter_type(grid..., R...)
    g = FT.(grid)
    r = FT.(R)
    sp = sortperm(g)
    g = g[sp]
    r = r[sp]
    all(diff(g) .> zero(FT)) ||
        throw(ArgumentError("grid values must be distinct"))
    return LUTWoodReflectance{FT}(g, r, grid_unit, extrapolation)
end

LUTWoodReflectance(; grid, R, grid_unit::Symbol = :nm,
                   extrapolation::Symbol = :clamp) =
    LUTWoodReflectance(grid, R; grid_unit, extrapolation)

"""
    PolynomialWoodReflectance(coeffs; grid_unit = :nm, x_offset = 0, x_scale = 1)

Polynomial bark/wood reflectance model.

Coefficients use ascending order:

```math
R(x) = c_0 + c_1 \\hat{x} + c_2 \\hat{x}^2 + \\cdots,
\\qquad \\hat{x} = (x - x_{offset}) / x_{scale}.
```

Use `x_offset` and `x_scale` to fit on a normalized wavelength or wavenumber
coordinate instead of raw nanometers.
"""
struct PolynomialWoodReflectance{FT<:Real} <: AbstractWoodReflectance{FT}
    coeffs::Vector{FT}
    grid_unit::Symbol
    x_offset::FT
    x_scale::FT
end

function PolynomialWoodReflectance(coeffs::AbstractVector{<:Real};
                                   grid_unit::Symbol = :nm,
                                   x_offset = nothing,
                                   x_scale = nothing)
    isempty(coeffs) && throw(ArgumentError("coeffs must be non-empty"))
    _check_spectral_grid_unit(grid_unit)
    FT = _wood_parameter_type(coeffs...,
                              x_offset === nothing ? 0 : x_offset,
                              x_scale === nothing ? 1 : x_scale)
    offset = x_offset === nothing ? zero(FT) : FT(x_offset)
    scale = x_scale === nothing ? one(FT) : FT(x_scale)
    scale != zero(FT) || throw(ArgumentError("x_scale must be non-zero"))
    return PolynomialWoodReflectance{FT}(FT.(coeffs), grid_unit, offset, scale)
end

PolynomialWoodReflectance(; coeffs, grid_unit::Symbol = :nm,
                          x_offset = nothing, x_scale = nothing) =
    PolynomialWoodReflectance(coeffs; grid_unit, x_offset, x_scale)

"""
    wood_reflectance(model, x; grid_unit = :nm)
    wood_reflectance(model, grid; grid_unit = :nm)

Evaluate a wood reflectance model at one spectral coordinate or a vector of
coordinates. `grid_unit` describes the input coordinate unit.
"""
wood_reflectance(model::ConstantWoodReflectance) = model.R

wood_reflectance(model::ConstantWoodReflectance, x::Real; grid_unit::Symbol = :nm) =
    model.R + zero(_spectral_coordinate(x, grid_unit, grid_unit))

function wood_reflectance(model::AbstractWoodReflectance,
                          grid::AbstractVector{<:Real};
                          grid_unit::Symbol = :nm)
    return [wood_reflectance(model, x; grid_unit) for x in grid]
end

function wood_reflectance(model::LUTWoodReflectance{FT},
                          x::Real;
                          grid_unit::Symbol = :nm) where {FT<:Real}
    x_model = _spectral_coordinate(x, grid_unit, model.grid_unit)
    g = model.grid
    r = model.R

    if x_model <= g[1]
        model.extrapolation == :clamp && return r[1] + zero(x_model)
        throw(ArgumentError("spectral coordinate is below the LUT range"))
    elseif x_model >= g[end]
        model.extrapolation == :clamp && return r[end] + zero(x_model)
        throw(ArgumentError("spectral coordinate is above the LUT range"))
    end

    idx = searchsortedlast(g, x_model)
    idx = min(max(idx, 1), length(g) - 1)
    x0 = g[idx]
    x1 = g[idx + 1]
    r0 = r[idx]
    r1 = r[idx + 1]
    t = (x_model - x0) / (x1 - x0)
    return r0 + t * (r1 - r0)
end

function wood_reflectance(model::PolynomialWoodReflectance{FT},
                          x::Real;
                          grid_unit::Symbol = :nm) where {FT<:Real}
    x_model = _spectral_coordinate(x, grid_unit, model.grid_unit)
    xhat = (x_model - model.x_offset) / model.x_scale
    value = zero(model.coeffs[end] + xhat)
    for c in Iterators.reverse(model.coeffs)
        value = value * xhat + c
    end
    return value
end

"""
    LambertianWoodCanopyScattering(reflectance)
    LambertianWoodCanopyScattering(R)
    LambertianWoodCanopyScattering(; R = 0.25, reflectance = nothing)

Opaque Lambertian woody element: reflection only, no transmission.

This is the angular-scattering counterpart to [`AbstractWoodReflectance`](@ref).
For scalar Z-matrix assembly it delegates to
`BiLambertianCanopyScattering(R, 0)`, preserving the existing canopy
normalization convention.
"""
struct LambertianWoodCanopyScattering{FT<:Real,W<:AbstractWoodReflectance{FT}} <:
       AbstractCanopyScatteringType{FT}
    reflectance::W
end

LambertianWoodCanopyScattering(R::Real) =
    LambertianWoodCanopyScattering(ConstantWoodReflectance(R))

function LambertianWoodCanopyScattering(; R = 0.25, reflectance = nothing)
    reflectance === nothing && return LambertianWoodCanopyScattering(R)
    return LambertianWoodCanopyScattering(reflectance)
end

const LambertianWood = LambertianWoodCanopyScattering

wood_reflectance(model::LambertianWoodCanopyScattering, args...; kwargs...) =
    wood_reflectance(model.reflectance, args...; kwargs...)

function _wood_scattering_reflectance(model::LambertianWoodCanopyScattering;
                                      spectral_coordinate = nothing,
                                      grid_unit::Symbol = :nm)
    if spectral_coordinate === nothing
        model.reflectance isa ConstantWoodReflectance &&
            return wood_reflectance(model.reflectance)
        throw(ArgumentError("spectral_coordinate is required for non-constant wood reflectance models"))
    end
    return wood_reflectance(model.reflectance, spectral_coordinate; grid_unit)
end
