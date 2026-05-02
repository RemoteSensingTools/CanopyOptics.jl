"""
    compute_Z_matrices(mod::BiLambertianCanopyScattering,
                       μ::AbstractVector, LD::AbstractLeafDistribution,
                       m::Integer; quadrature = CanopyQuadrature(), npol = 1) -> (Z⁺⁺, Z⁻⁺)

Compute one cosine Fourier moment of the bi-Lambertian canopy scattering
matrices.  Rows are outgoing streams and columns are incoming streams:
`Z[i_out, j_in]`.

`npol = 1` returns scalar Stokes-I matrices. `npol = 3` and `npol = 4`
expand the scalar diffuse result into vSmartMOM-style Stokes blocks with
only I→I populated.

This method uses the closed-form Fourier moments from
[`compute_Z_matrices_aniso_analytic`](@ref).  For all moments in one call, pass
a range such as `0:m_max`.
"""
function compute_Z_matrices(mod::BiLambertianCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices_aniso_analytic(mod, μ, LD, Int(m);
                                                 quadrature = q, npol = npol)
    return Z⁺⁺[:, :, m + 1], Z⁻⁺[:, :, m + 1]
end

"""
    compute_Z_matrices(mod::LambertianWoodCanopyScattering,
                       μ, LD, m; spectral_coordinate = nothing, grid_unit = :nm,
                       npol = 1)

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
                            grid_unit::Symbol = :nm,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices(diffuse, μ, LD, m; quadrature = q, npol = npol)
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
                          quadrature::CanopyQuadrature = CanopyQuadrature(),
                          npol::Integer = 1)
    _check_m_range(m_range)
    Z⁺⁺₀, Z⁻⁺₀ = compute_component_Z(mod, μ, LD, first(m_range);
                                      quadrature = quadrature, npol = npol)
    nμ_out, nμ_in = size(Z⁺⁺₀)
    nm = length(m_range)
    Z⁺⁺ = similar(Z⁺⁺₀, nμ_out, nμ_in, nm)
    Z⁻⁺ = similar(Z⁻⁺₀, nμ_out, nμ_in, nm)
    Z⁺⁺[:, :, 1] .= Z⁺⁺₀
    Z⁻⁺[:, :, 1] .= Z⁻⁺₀

    k = 2
    for m in Iterators.drop(m_range, 1)
        Z⁺⁺ₘ, Z⁻⁺ₘ = compute_component_Z(mod, μ, LD, m;
                                          quadrature = quadrature, npol = npol)
        Z⁺⁺[:, :, k] .= Z⁺⁺ₘ
        Z⁻⁺[:, :, k] .= Z⁻⁺ₘ
        k += 1
    end
    return Z⁺⁺, Z⁻⁺
end

"""
    compute_Z_matrices(mod::AbstractCanopyScatteringType,
                       μ, LD, m_range::AbstractUnitRange; quadrature, npol = 1)

Compute a stack of canopy Z matrices for Fourier moments in `m_range`.

The scalar default returns arrays with shape
`(length(μ), length(μ), length(m_range))`.  Passing `npol = 3` or `4`
returns Stokes-expanded arrays with the first two dimensions multiplied by
`npol`.  In both cases the package convention is `Z[i_out, j_in, k]`, where
`k = m - first(m_range) + 1`.
"""
function compute_Z_matrices(mod::BiLambertianCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices_aniso_analytic(mod, μ, LD, last(m_range);
                                                 quadrature = q, npol = npol)
    return Z⁺⁺[:, :, m_range .+ 1], Z⁻⁺[:, :, m_range .+ 1]
end

function compute_Z_matrices(mod::LambertianWoodCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            spectral_coordinate = nothing,
                            grid_unit::Symbol = :nm,
                            npol::Integer = 1) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices(diffuse, μ, LD, m_range; quadrature = q,
                              npol = npol)
end

function compute_Z_matrices(mod::SpecularCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m_range::AbstractUnitRange{<:Integer};
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _stack_Z_moments(compute_Z_matrices, mod, μ, LD, m_range;
                            quadrature = q, npol = npol)
end

"""
    _sum_component_Z(compute_component_Z, components, μ, LD, m; quadrature)

Evaluate each component in a [`CompositeCanopyScattering`](@ref) model for one
Fourier moment and add the resulting Z matrices elementwise.
"""
function _sum_component_Z(compute_component_Z, components::Tuple, μ, LD, m::Int;
                          quadrature::CanopyQuadrature = CanopyQuadrature(),
                          npol::Integer = 1)
    Z⁺⁺, Z⁻⁺ = compute_component_Z(components[1], μ, LD, m;
                                    quadrature = quadrature, npol = npol)
    Z⁺⁺_sum = copy(Z⁺⁺)
    Z⁻⁺_sum = copy(Z⁻⁺)

    for i in 2:length(components)
        Z⁺⁺ᵢ, Z⁻⁺ᵢ = compute_component_Z(components[i], μ, LD, m;
                                          quadrature = quadrature, npol = npol)
        Z⁺⁺_sum = Z⁺⁺_sum .+ Z⁺⁺ᵢ
        Z⁻⁺_sum = Z⁻⁺_sum .+ Z⁻⁺ᵢ
    end
    return Z⁺⁺_sum, Z⁻⁺_sum
end

"""
    compute_Z_matrices(mod::CompositeCanopyScattering, μ, LD, m; quadrature, npol = 1)

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
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices, mod.components, μ, LD, m;
                            quadrature = q, npol = npol)
end

function compute_Z_matrices(mod::CompositeCanopyScattering,
                            μ::AbstractVector{FT},
                            LD::AbstractLeafDistribution,
                            m::Integer;
                            quadrature::CanopyQuadrature = CanopyQuadrature(),
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices, mod.components, μ, LD, Int(m);
                            quadrature = q, npol = npol)
end

"""
    _sum_component_Z_stack(components, μ, LD, m_range; quadrature)

Stacked-moment analogue of [`_sum_component_Z`](@ref) for composite canopy
scattering models.
"""
function _sum_component_Z_stack(components::Tuple, μ, LD,
                                m_range::AbstractUnitRange{<:Integer};
                                quadrature::CanopyQuadrature = CanopyQuadrature(),
                                npol::Integer = 1)
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices(components[1], μ, LD, m_range;
                                  quadrature = quadrature, npol = npol)
    Z⁺⁺_sum = copy(Z⁺⁺)
    Z⁻⁺_sum = copy(Z⁻⁺)

    for i in 2:length(components)
        Z⁺⁺ᵢ, Z⁻⁺ᵢ = compute_Z_matrices(components[i], μ, LD, m_range;
                                          quadrature = quadrature, npol = npol)
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
                            nQuad = nothing,
                            npol::Integer = 1) where FT
    _check_m_range(m_range)
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z_stack(mod.components, μ, LD, m_range;
                                  quadrature = q, npol = npol)
end

function _component_G_weights(components::Tuple, μ::AbstractVector)
    Gs = map(component -> vec(G(μ, component.LAD)), components)
    FTG = eltype(Gs[1])
    for i in eachindex(components)
        FTG = promote_type(FTG, eltype(Gs[i]), typeof(components[i].area_index))
    end

    G_total = zeros(FTG, length(μ))
    for i in eachindex(components)
        area = FTG(components[i].area_index)
        Gi = Gs[i]
        @inbounds for j in eachindex(G_total)
            G_total[j] += area * FTG(Gi[j])
        end
    end
    return Gs, G_total
end

function _add_weighted_component_Z!(Zpp_sum::AbstractMatrix, Zmp_sum::AbstractMatrix,
                                    Zpp, Zmp, weights, npol::Int)
    nμ = length(weights)
    @inbounds for j in 1:nμ, sj in 1:npol
        col = _stokes_index(j, sj, npol)
        weight = weights[j]
        for row in axes(Zpp_sum, 1)
            Zpp_sum[row, col] += weight * Zpp[row, col]
            Zmp_sum[row, col] += weight * Zmp[row, col]
        end
    end
    return nothing
end

function _add_weighted_component_Z!(Zpp_sum::AbstractArray{<:Any,3},
                                    Zmp_sum::AbstractArray{<:Any,3},
                                    Zpp, Zmp, weights, npol::Int)
    nμ = length(weights)
    @inbounds for k in axes(Zpp_sum, 3), j in 1:nμ, sj in 1:npol
        col = _stokes_index(j, sj, npol)
        weight = weights[j]
        for row in axes(Zpp_sum, 1)
            Zpp_sum[row, col, k] += weight * Zpp[row, col, k]
            Zmp_sum[row, col, k] += weight * Zmp[row, col, k]
        end
    end
    return nothing
end

function _mixed_component_weights(component::CanopyComponent, G_component, G_total)
    weights = zeros(promote_type(eltype(G_component), eltype(G_total),
                                 typeof(component.area_index)), length(G_total))
    area = eltype(weights)(component.area_index)
    @inbounds for j in eachindex(weights)
        weights[j] = G_total[j] == zero(G_total[j]) ?
                     zero(eltype(weights)) :
                     area * G_component[j] / G_total[j]
    end
    return weights
end

function _compute_mixed_Z(canopy::MixedCanopy, μ::AbstractVector, m;
                          quadrature::CanopyQuadrature = CanopyQuadrature(),
                          nQuad = nothing,
                          npol::Integer = 1,
                          kwargs...)
    n = _validate_npol(npol)
    q = _resolve_quadrature(quadrature, nQuad)
    components = canopy.components
    Gs, G_total = _component_G_weights(components, μ)

    first_component = components[1]
    Zpp₁, Zmp₁ = compute_Z_matrices(first_component.scatterer, μ,
                                     first_component.LAD, m;
                                     quadrature = q, npol = n, kwargs...)
    Zpp_sum = zero.(Zpp₁)
    Zmp_sum = zero.(Zmp₁)
    weights = _mixed_component_weights(first_component, Gs[1], G_total)
    _add_weighted_component_Z!(Zpp_sum, Zmp_sum, Zpp₁, Zmp₁, weights, n)

    for i in 2:length(components)
        component = components[i]
        Zppᵢ, Zmpᵢ = compute_Z_matrices(component.scatterer, μ, component.LAD, m;
                                        quadrature = q, npol = n, kwargs...)
        weights = _mixed_component_weights(component, Gs[i], G_total)
        _add_weighted_component_Z!(Zpp_sum, Zmp_sum, Zppᵢ, Zmpᵢ, weights, n)
    end
    return Zpp_sum, Zmp_sum
end

"""
    compute_Z_matrices(canopy::MixedCanopy, μ, m; quadrature, npol = 1)
    compute_Z_matrices(canopy::MixedCanopy, μ, ignored_LD, m; quadrature, npol = 1)

Compute canopy Z matrices for a mixture of components with their own angle
distributions and area indices.  The `ignored_LD` compatibility form exists so
callers that already carry `(scatterer, μ, LAD, m)` can pass a `MixedCanopy` as
the scatterer without changing their call shape.
"""
compute_Z_matrices(canopy::MixedCanopy,
                   μ::AbstractVector,
                   m::Integer;
                   quadrature::CanopyQuadrature = CanopyQuadrature(),
                   nQuad = nothing,
                   npol::Integer = 1,
                   kwargs...) =
    _compute_mixed_Z(canopy, μ, Int(m);
                     quadrature = quadrature, nQuad = nQuad,
                     npol = npol, kwargs...)

compute_Z_matrices(canopy::MixedCanopy,
                   μ::AbstractVector,
                   m_range::AbstractUnitRange{<:Integer};
                   quadrature::CanopyQuadrature = CanopyQuadrature(),
                   nQuad = nothing,
                   npol::Integer = 1,
                   kwargs...) =
    _compute_mixed_Z(canopy, μ, m_range;
                     quadrature = quadrature, nQuad = nQuad,
                     npol = npol, kwargs...)

compute_Z_matrices(canopy::MixedCanopy,
                   μ::AbstractVector,
                   LD::AbstractLeafDistribution,
                   m;
                   quadrature::CanopyQuadrature = CanopyQuadrature(),
                   nQuad = nothing,
                   npol::Integer = 1,
                   kwargs...) =
    _compute_mixed_Z(canopy, μ, m;
                     quadrature = quadrature, nQuad = nQuad,
                     npol = npol, kwargs...)

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
                                  nQuad = nothing,
                                  npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, m; quadrature = q, npol = npol)
end

function compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  Zup, Zdown, m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices_aniso(mod, μ, LD, m; quadrature = q, npol = npol)
end

function compute_Z_matrices_aniso(mod::LambertianWoodCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  spectral_coordinate = nothing,
                                  grid_unit::Symbol = :nm,
                                  npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, m; quadrature = q,
                              spectral_coordinate = spectral_coordinate,
                              grid_unit = grid_unit, npol = npol)
end

function compute_Z_matrices_aniso(mod::SpecularCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, collect(μ), LD, m; quadrature = q, npol = npol)
end

function compute_Z_matrices_aniso(mod::CompositeCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  npol::Integer = 1) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return _sum_component_Z(compute_Z_matrices_aniso, mod.components, μ, LD, m;
                            quadrature = q, npol = npol)
end

function compute_Z_matrices_aniso(canopy::MixedCanopy,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing,
                                  npol::Integer = 1,
                                  kwargs...) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(canopy, μ, m;
                              quadrature = q, npol = npol, kwargs...)
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
                                           nQuad = nothing,
                                           npol::Integer = 1) where {FT<:Real}
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(mod, μ, LD, 0:m_max; quadrature = q, npol = npol)
end

function compute_Z_matrices_aniso_analytic(canopy::MixedCanopy,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int;
                                           quadrature::CanopyQuadrature = CanopyQuadrature(),
                                           nQuad = nothing,
                                           npol::Integer = 1,
                                           kwargs...) where {FT<:Real}
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices(canopy, μ, 0:m_max;
                              quadrature = q, npol = npol, kwargs...)
end

function compute_Z_matrices_aniso_analytic(mod::LambertianWoodCanopyScattering,
                                           μ::AbstractVector{FT},
                                           LD::AbstractLeafDistribution,
                                           m_max::Int;
                                           quadrature::CanopyQuadrature = CanopyQuadrature(),
                                           nQuad = nothing,
                                           spectral_coordinate = nothing,
                                           grid_unit::Symbol = :nm,
                                           npol::Integer = 1) where {FT<:Real}
    q = _resolve_quadrature(quadrature, nQuad)
    R = _wood_scattering_reflectance(mod;
                                     spectral_coordinate = spectral_coordinate,
                                     grid_unit = grid_unit)
    diffuse = BiLambertianCanopyScattering(R = R, T = zero(R))
    return compute_Z_matrices_aniso_analytic(diffuse, μ, LD, m_max;
                                             quadrature = q, npol = npol)
end
