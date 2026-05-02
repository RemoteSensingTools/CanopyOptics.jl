"""
    compute_Z_matrices(mod::BiLambertianCanopyScattering,
                       μ::AbstractVector, LD::AbstractLeafDistribution,
                       m::Integer; quadrature = CanopyQuadrature()) -> (Z⁺⁺, Z⁻⁺)

Compute one cosine Fourier moment of the bi-Lambertian canopy scattering
matrices.  Rows are outgoing streams and columns are incoming streams:
`Z[i_out, j_in]`.

This method uses the closed-form Fourier moments from
[`compute_Z_matrices_aniso_analytic`](@ref).  For all moments in one call, pass
a range such as `0:m_max`.
"""
function compute_Z_matrices(mod::BiLambertianCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices_aniso_analytic(mod, μ, LD, Int(m);
                                                 quadrature = q)
    return Z⁺⁺[:, :, m + 1], Z⁻⁺[:, :, m + 1]
end

"""
    compute_Z_matrices(mod::LambertianWoodCanopyScattering,
                       μ, LD, m; spectral_coordinate = nothing, grid_unit = :nm)

Compute one Fourier moment for opaque Lambertian wood.

The wood reflectance is evaluated, then the angular Z construction delegates to
`BiLambertianCanopyScattering(R, 0)`. For non-constant wood reflectance models,
pass `spectral_coordinate` in the unit specified by `grid_unit`.
"""
function compute_Z_matrices(mod::LambertianWoodCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            spectral_coordinate = nothing,
                            grid_unit::Symbol = :nm) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices(diffuse, μ, LD, m; quadrature = q)
end

"""
    _check_m_range(m_range::AbstractUnitRange{<:Integer})

Validate a Fourier-moment range for stacked Z-matrix assembly.

Moment ranges must be non-empty and non-negative because canopy Fourier
moments are stored with array index `m + 1`.
"""
function _check_m_range(m_range::AbstractUnitRange{<:Integer})
    isempty(m_range) && throw(ArgumentError("Fourier moment range must be non-empty"))
    first(m_range) < 0 && throw(ArgumentError("Fourier moments must be non-negative"))
    return m_range
end

"""
    _stack_Z_moments(compute_component_Z, mod, μ, LD, m_range; quadrature)

Build a three-dimensional `(i_out, j_in, m)` Z stack from a single-moment
component method.

This is used for components, such as specular scattering, that do not yet have
a native all-moments analytic implementation.
"""
function _stack_Z_moments(compute_component_Z, mod, μ, LD,
                          m_range::AbstractUnitRange{<:Integer};
                          quadrature::CanopyQuadrature = CanopyQuadrature())
    _check_m_range(m_range)
    Z⁺⁺₀, Z⁻⁺₀ = compute_component_Z(mod, μ, LD, first(m_range);
                                      quadrature = quadrature)
    nμ_out, nμ_in = size(Z⁺⁺₀)
    nm = length(m_range)
    Z⁺⁺ = similar(Z⁺⁺₀, nμ_out, nμ_in, nm)
    Z⁻⁺ = similar(Z⁻⁺₀, nμ_out, nμ_in, nm)
    Z⁺⁺[:, :, 1] .= Z⁺⁺₀
    Z⁻⁺[:, :, 1] .= Z⁻⁺₀

    k = 2
    for m in Iterators.drop(m_range, 1)
        Z⁺⁺ₘ, Z⁻⁺ₘ = compute_component_Z(mod, μ, LD, m;
                                          quadrature = quadrature)
        Z⁺⁺[:, :, k] .= Z⁺⁺ₘ
        Z⁻⁺[:, :, k] .= Z⁻⁺ₘ
        k += 1
    end
    return Z⁺⁺, Z⁻⁺
end

"""
    compute_Z_matrices(mod::AbstractCanopyScatteringType,
                       μ, LD, m_range::AbstractUnitRange; quadrature)

Compute a stack of canopy Z matrices for Fourier moments in `m_range`.

The returned arrays have shape `(length(μ), length(μ), length(m_range))` and
follow the package convention `Z[i_out, j_in, k]`, where
`k = m - first(m_range) + 1`.
"""
function compute_Z_matrices(mod::BiLambertianCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices_aniso_analytic(mod, μ, LD, last(m_range);
                                                 quadrature = q)
    return Z⁺⁺[:, :, m_range .+ 1], Z⁻⁺[:, :, m_range .+ 1]
end

function compute_Z_matrices(mod::LambertianWoodCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            spectral_coordinate = nothing,
                            grid_unit::Symbol = :nm) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices(diffuse, μ, LD, m_range; quadrature = q)
end

function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _stack_Z_moments(compute_Z_matrices, mod, μ, LD, m_range;
                            quadrature = q)
end

"""
    _sum_component_Z(compute_component_Z, components, μ, LD, m; quadrature)

Evaluate each component in a [`CompositeCanopyScattering`](@ref) model for one
Fourier moment and add the resulting Z matrices elementwise.
"""
function _sum_component_Z(compute_component_Z, components::Tuple, μ, LD, m::Int;
                          quadrature::CanopyQuadrature = CanopyQuadrature())
    Z⁺⁺, Z⁻⁺ = compute_component_Z(components[1], μ, LD, m;
                                    quadrature = quadrature)
    Z⁺⁺_sum = copy(Z⁺⁺)
    Z⁻⁺_sum = copy(Z⁻⁺)

    for i in 2:length(components)
        Z⁺⁺ᵢ, Z⁻⁺ᵢ = compute_component_Z(components[i], μ, LD, m;
                                          quadrature = quadrature)
        Z⁺⁺_sum = Z⁺⁺_sum .+ Z⁺⁺ᵢ
        Z⁻⁺_sum = Z⁻⁺_sum .+ Z⁻⁺ᵢ
    end
    return Z⁺⁺_sum, Z⁻⁺_sum
end

"""
    compute_Z_matrices(mod::CompositeCanopyScattering, μ, LD, m; quadrature)

Compute one Fourier moment for an additive canopy-scattering model.

Each component is evaluated independently and the resulting `(Z⁺⁺, Z⁻⁺)`
matrices are summed, so diffuse and specular leaf terms combine at the optical
property level.
"""
function compute_Z_matrices(mod::CompositeCanopyScattering,
                            μ::Array{FT,1},
                            LD::AbstractLeafDistribution,
                            m::Int;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices, mod.components, μ, LD, m;
                            quadrature = q)
end

function compute_Z_matrices(mod::CompositeCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices, mod.components, μ, LD, Int(m);
                            quadrature = q)
end

"""
    _sum_component_Z_stack(components, μ, LD, m_range; quadrature)

Stacked-moment analogue of [`_sum_component_Z`](@ref) for composite canopy
scattering models.
"""
function _sum_component_Z_stack(components::Tuple, μ, LD,
                                m_range::AbstractUnitRange{<:Integer};
                                quadrature::CanopyQuadrature = CanopyQuadrature())
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices(components[1], μ, LD, m_range;
                                  quadrature = quadrature)
    Z⁺⁺_sum = copy(Z⁺⁺)
    Z⁻⁺_sum = copy(Z⁻⁺)

    for i in 2:length(components)
        Z⁺⁺ᵢ, Z⁻⁺ᵢ = compute_Z_matrices(components[i], μ, LD, m_range;
                                          quadrature = quadrature)
        Z⁺⁺_sum = Z⁺⁺_sum .+ Z⁺⁺ᵢ
        Z⁻⁺_sum = Z⁻⁺_sum .+ Z⁻⁺ᵢ
    end
    return Z⁺⁺_sum, Z⁻⁺_sum
end

function compute_Z_matrices(mod::CompositeCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z_stack(mod.components, μ, LD, m_range;
                                  quadrature = q)
end

"""
    compute_Z_matrices_aniso(mod, μ, LD, m; quadrature)

Compatibility name for single-moment canopy Z assembly.

Historically this selected an azimuth-quadrature implementation for
bi-Lambertian scattering. It now delegates to the canonical
[`compute_Z_matrices`](@ref) API, which uses the analytic Fourier path where
available.
"""
function compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, m; quadrature = q)
end

function compute_Z_matrices_aniso(mod::LambertianWoodCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  spectral_coordinate = nothing,
                                  grid_unit::Symbol = :nm) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, m; quadrature = q,
                              spectral_coordinate = spectral_coordinate,
                              grid_unit = grid_unit)
end

function compute_Z_matrices_aniso(mod::SpecularCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, Array(μ), LD, m; quadrature = q)
end

function compute_Z_matrices_aniso(mod::CompositeCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices_aniso, mod.components, μ, LD, m;
                            quadrature = q)
end

"""
    compute_Z_matrices_aniso_analytic(mod::CompositeCanopyScattering,
                                      μ, LD, m_max; quadrature)

Compute all moments `0:m_max` for a composite canopy-scattering model by
delegating to the stacked [`compute_Z_matrices`](@ref) API.
"""
function compute_Z_matrices_aniso_analytic(mod::CompositeCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int;
                                           quadrature::CanopyQuadrature = CanopyQuadrature(),
                                           nQuad = nothing) where {FT<:Real}
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, 0:m_max; quadrature = q)
end

function compute_Z_matrices_aniso_analytic(mod::LambertianWoodCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int;
                                           quadrature::CanopyQuadrature = CanopyQuadrature(),
                                           nQuad = nothing,
                                           spectral_coordinate = nothing,
                                           grid_unit::Symbol = :nm) where {FT<:Real}
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices_aniso_analytic(diffuse, μ, LD, m_max;
                                             quadrature = q)
end
