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

# Float type for a LeafDistribution — used in the CanopyComponent
# promotion below so a Float32 LAD + Float32 scatterer + Float32
# area_index does not silently get cast through Float64.
_lad_ft(::AbstractLeafDistribution{FT}) where {FT} = FT

function CanopyComponent(scatterer::AbstractCanopyScatteringType,
                         LAD::AbstractLeafDistribution,
                         area_index::Real;
                         clumping = nothing)
    # Promote without an artificial Float64 anchor: the resulting FT
    # is the natural promotion of the user-supplied types (scatterer,
    # LAD, area_index, and clumping if provided). The previous
    # `clumping === nothing ? Float64 : ...` term silently forced
    # everything to Float64 whenever the caller didn't pass an
    # explicit clumping model — defeating type-stable Float32 paths.
    FT = if clumping === nothing
        float(promote_type(_canopy_scattering_ft(scatterer),
                           _lad_ft(LAD),
                           typeof(area_index)))
    else
        float(promote_type(_canopy_scattering_ft(scatterer),
                           _lad_ft(LAD),
                           typeof(area_index),
                           _clumping_ft(clumping)))
    end
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

#####
##### Convenience constructors for the canonical tree-canopy populations
#####
##### Wraps the existing (scatterer, LAD, area_index) tuple into
##### named helpers — `LeafComponent`, `StemComponent`,
##### `BranchComponent` — that pick sensible defaults for each
##### population and accept overrides on every axis.  Then
##### `TreeCanopy(...)` composes them into a single `MixedCanopy`.
#####

"""
    LeafComponent(; LAI, scatterer = nothing, LAD = nothing,
                  clumping = nothing) -> CanopyComponent

Convenience constructor for the **leaf** scattering population of a
tree canopy.  Defaults to a `BiLambertianCanopyScattering` (reflection
+ transmission, the canonical leaf scatterer) and `spherical_leaves`
(Bonan/Norman default — projection function `G(μ) ≈ 0.5`).  Override
either via kwarg: e.g. `scatterer = some_prospect_driven_scatterer`
or `LAD = planophile_leaves(FT)`.

When `scatterer` and `LAD` are left at their defaults, they are
constructed in the float type of `LAI` — so `LeafComponent(LAI = 3.5f0)`
returns a `CanopyComponent{Float32, ...}` end-to-end, with no silent
Float64 promotion.

`LAI` is the one-sided leaf area index `[m² leaf / m² ground]`.  Must
be non-negative (throws `DomainError` otherwise).

The resulting `CanopyComponent` plugs directly into
[`MixedCanopy`](@ref) and the existing `compute_Z_matrices`,
`bulk_G`, `effective_G` paths — no other call sites change.
"""
function LeafComponent(; LAI::Real,
                      scatterer::Union{AbstractCanopyScatteringType, Nothing} = nothing,
                      LAD::Union{AbstractLeafDistribution, Nothing} = nothing,
                      clumping = nothing)
    LAI ≥ 0 || throw(DomainError(LAI, "LeafComponent: LAI must be non-negative"))
    FT = float(typeof(LAI))
    if scatterer === nothing
        scatterer = BiLambertianCanopyScattering{FT}()
    end
    if LAD === nothing
        LAD = spherical_leaves(_canopy_scattering_ft(scatterer))
    end
    return CanopyComponent(scatterer, LAD, LAI; clumping = clumping)
end

"""
    StemComponent(; SAI, R = 0.25, scatterer = nothing, LAD = nothing,
                  clumping = nothing) -> CanopyComponent

Convenience constructor for the **stem** scattering population.

Stems are modeled as opaque Lambertian (reflection only, no
transmission — the "no T, just R" point) with a default
**erectophile** (mostly vertical) angle distribution.  This captures
the fact that trunks intercept little overhead sun (small G at high
μ) but a lot of low-angle sun (large G at small μ).

`R` may be a scalar (broadband bark albedo) or any
[`AbstractWoodReflectance`](@ref) (e.g. `LUTWoodReflectance` for
band-resolved bark spectra).  The default `R = 0.25` is a broadband
typical mature-bark value (Bonan 2019, *Climate Change and Terrestrial
Ecosystem Modeling*, Table 14.1).  Pass `scatterer = ...` to override
the entire scattering model (rarely needed).

When `scatterer` and `LAD` are left at their defaults, they are
constructed in the float type of `SAI` — so `StemComponent(SAI = 0.9f0)`
returns a `CanopyComponent{Float32, ...}` end-to-end.

`SAI` is the one-sided stem area index `[m² stem / m² ground]`.  Must
be non-negative (throws `DomainError` otherwise).
"""
function StemComponent(; SAI::Real,
                      R = 0.25,
                      scatterer::Union{AbstractCanopyScatteringType, Nothing} = nothing,
                      LAD::Union{AbstractLeafDistribution, Nothing} = nothing,
                      clumping = nothing)
    SAI ≥ 0 || throw(DomainError(SAI, "StemComponent: SAI must be non-negative"))
    FT = float(typeof(SAI))
    if scatterer === nothing
        reflectance = if R isa AbstractWoodReflectance
            R
        else
            ConstantWoodReflectance{FT}(FT(R))
        end
        scatterer = LambertianWoodCanopyScattering(reflectance)
    end
    if LAD === nothing
        LAD = erectophile_leaves(_canopy_scattering_ft(scatterer))
    end
    return CanopyComponent(scatterer, LAD, SAI; clumping = clumping)
end

"""
    BranchComponent(; BAI, R = 0.30, scatterer = nothing, LAD = nothing,
                    clumping = nothing) -> CanopyComponent

Convenience constructor for the **branch** scattering population.

Similar to [`StemComponent`](@ref) but with a default
**plagiophile** (oblique, on-average ~45°-tilted) angle distribution
that represents branches hanging off the main stem at intermediate
inclinations.  The default `R = 0.30` is a broadband typical
younger-bark value, slightly higher than the mature-stem default to
reflect smoother branch bark; override with any
[`AbstractWoodReflectance`](@ref) for band-resolved spectra.

When `scatterer` and `LAD` are left at their defaults, they are
constructed in the float type of `BAI` — so `BranchComponent(BAI = 0.1f0)`
returns a `CanopyComponent{Float32, ...}` end-to-end.

`BAI` is the one-sided branch area index `[m² branch / m² ground]`.
Must be non-negative (throws `DomainError` otherwise).
"""
function BranchComponent(; BAI::Real,
                        R = 0.30,
                        scatterer::Union{AbstractCanopyScatteringType, Nothing} = nothing,
                        LAD::Union{AbstractLeafDistribution, Nothing} = nothing,
                        clumping = nothing)
    BAI ≥ 0 || throw(DomainError(BAI, "BranchComponent: BAI must be non-negative"))
    FT = float(typeof(BAI))
    if scatterer === nothing
        reflectance = if R isa AbstractWoodReflectance
            R
        else
            ConstantWoodReflectance{FT}(FT(R))
        end
        scatterer = LambertianWoodCanopyScattering(reflectance)
    end
    if LAD === nothing
        LAD = plagiophile_leaves(_canopy_scattering_ft(scatterer))
    end
    return CanopyComponent(scatterer, LAD, BAI; clumping = clumping)
end

"""
    TreeCanopy(; LAI, SAI = 0, BAI = 0,
               leaf_kwargs = (;), stem_kwargs = (;), branch_kwargs = (;))
        -> MixedCanopy

One-line constructor for "leaves + stems + branches" tree canopies.
Returns a [`MixedCanopy`](@ref) with up to three
[`CanopyComponent`](@ref)s.

A `SAI` or `BAI` of zero **omits the corresponding component
entirely** — `MixedCanopy` doesn't carry spurious zero-area
populations in its projected-area weighting.  When SAI = BAI = 0 the
result is a single-component `MixedCanopy` whose `bulk_G` matches
the leaf-only `G(μ, ::LeafDistribution)` exactly.

Each population's defaults can be overridden via the `*_kwargs`
NamedTuples, which are spliced into the corresponding
`LeafComponent` / `StemComponent` / `BranchComponent` call.  For
example, to use band-resolved bark reflectance:

```julia
using CanopyOptics
canopy = TreeCanopy(;
    LAI = 4.0,
    SAI = 0.9,
    BAI = 0.1,
    stem_kwargs = (; R = LUTWoodReflectance(
        grid = [400, 800, 2500], R = [0.10, 0.20, 0.45], grid_unit = :nm)),
)
```

The resulting `MixedCanopy` plugs into all existing paths
(`bulk_G`, `compute_Z_matrices`, etc.) and is the recommended way
to construct tree canopies in upstream consumers (RRTMGP,
ClimaLand multi-layer canopies, CanopyColumn).
"""
function TreeCanopy(; LAI::Real,
                   SAI::Real = zero(LAI),
                   BAI::Real = zero(LAI),
                   leaf_kwargs::NamedTuple = (;),
                   stem_kwargs::NamedTuple = (;),
                   branch_kwargs::NamedTuple = (;))
    # Per-component constructors throw `DomainError` on negative
    # inputs.  Reject up front so a negative SAI / BAI never silently
    # collapses to an LAI-only canopy via the `> zero(...)` gating
    # below.
    LAI ≥ 0 || throw(DomainError(LAI, "TreeCanopy: LAI must be non-negative"))
    SAI ≥ 0 || throw(DomainError(SAI, "TreeCanopy: SAI must be non-negative"))
    BAI ≥ 0 || throw(DomainError(BAI, "TreeCanopy: BAI must be non-negative"))
    components = CanopyComponent[]
    push!(components, LeafComponent(; LAI = LAI, leaf_kwargs...))
    if SAI > zero(SAI)
        push!(components, StemComponent(; SAI = SAI, stem_kwargs...))
    end
    if BAI > zero(BAI)
        push!(components, BranchComponent(; BAI = BAI, branch_kwargs...))
    end
    return MixedCanopy(Tuple(components))
end
