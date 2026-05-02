"""
    CanopyComponent(scatterer, LAD, area_index; clumping = nothing)
    CanopyComponent(; scatterer, LAD, area_index, clumping = nothing)

One canopy scattering population.

`scatterer` is a canopy scattering model such as
[`BiLambertianCanopyScattering`](@ref) or
[`LambertianWoodCanopyScattering`](@ref), `LAD` is that population's element
angle distribution, and `area_index` is its one-sided area index: LAI for
leaves, WAI/BAI for woody material. Optional per-component clumping is used by
[`bulk_G`](@ref), but mixed Z-matrix normalization keeps using the unclumped
Ross `G` convention.
"""
struct CanopyComponent{FT<:Real,S<:AbstractCanopyScatteringType,
                       LD<:AbstractLeafDistribution,C<:AbstractClumping}
    scatterer::S
    LAD::LD
    area_index::FT
    clumping::C
end

_clumping_ft(::AbstractClumping{FT}) where {FT} = FT

function CanopyComponent(scatterer::AbstractCanopyScatteringType,
                         LAD::AbstractLeafDistribution,
                         area_index::Real;
                         clumping = nothing)
    FT = float(promote_type(_canopy_scattering_ft(scatterer), typeof(area_index),
                            clumping === nothing ? Float64 : _clumping_ft(clumping)))
    c = clumping === nothing ? NoClumping{FT}() : clumping
    return CanopyComponent{FT,typeof(scatterer),typeof(LAD),typeof(c)}(
        scatterer, LAD, FT(area_index), c)
end

function CanopyComponent(; scatterer, LAD, area_index = nothing,
                         AI = nothing, clumping = nothing)
    area = area_index === nothing ? AI : area_index
    area === nothing && throw(ArgumentError("CanopyComponent needs area_index or AI"))
    return CanopyComponent(scatterer, LAD, area; clumping = clumping)
end

"""
    MixedCanopy(components...)
    MixedCanopy((components...,))

Area-weighted canopy made from multiple scattering populations.

Each component may have its own scatterer, angle distribution, area index, and
clumping model. Z matrices are mixed column-by-column using the incoming-stream
projected area weights

```math
w_c(μ_j) =
\\frac{AI_c G_c(μ_j)}
     {\\sum_d AI_d G_d(μ_j)}.
```

The optional component clumping models are intentionally not included in this
Z normalization; use [`bulk_G`](@ref) with `clumped = true` for propagation.
"""
struct MixedCanopy{FT<:Real,T<:Tuple} <: AbstractCanopyScatteringType{FT}
    components::T
end

_component_ft(c::CanopyComponent{FT}) where {FT} = FT

function MixedCanopy(components::CanopyComponent...)
    isempty(components) && throw(ArgumentError("MixedCanopy needs at least one component"))
    FT = promote_type(map(_component_ft, components)...)
    return MixedCanopy{FT,typeof(components)}(components)
end

MixedCanopy(components::Tuple) = MixedCanopy(components...)

"""
    component_G(component, μ; clumped = false)

Projected area for one [`CanopyComponent`](@ref). Set `clumped = true` to apply
the component's clumping model.
"""
function component_G(component::CanopyComponent, μ::AbstractVector; clumped::Bool = false)
    G_raw = vec(G(μ, component.LAD))
    return clumped ? effective_G(component.clumping, G_raw, μ) : G_raw
end

"""
    bulk_G(canopy, μ; clumped = true)

Area-index-weighted canopy projected area:

```math
G_{bulk}(μ) = \\sum_c AI_c G_c(μ),
```

or `\\sum_c AI_c Ω_c(μ)G_c(μ)` when `clumped = true`.
"""
function bulk_G(canopy::MixedCanopy, μ::AbstractVector; clumped::Bool = true)
    first_G = component_G(canopy.components[1], μ; clumped = clumped)
    FT = promote_type(eltype(first_G), typeof(canopy.components[1].area_index))
    for component in Iterators.drop(canopy.components, 1)
        Gc = component_G(component, μ; clumped = clumped)
        FT = promote_type(FT, eltype(Gc), typeof(component.area_index))
    end

    total = zeros(FT, length(μ))
    for component in canopy.components
        Gc = component_G(component, μ; clumped = clumped)
        @inbounds for i in eachindex(total)
            total[i] += FT(component.area_index) * FT(Gc[i])
        end
    end
    return total
end

G(μ::AbstractArray, canopy::MixedCanopy; clumped::Bool = false) =
    bulk_G(canopy, μ; clumped = clumped)
