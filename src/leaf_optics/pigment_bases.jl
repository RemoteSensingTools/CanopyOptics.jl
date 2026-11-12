"""
    LeafMultiPigmentProperties{FT<:Real} <: AbstractLeafProperties

Leaf biochemical and structural state with separate chlorophyll *a* (`Cchla`)
and chlorophyll *b* (`Cchlb`) pools alongside carotenoids (`Ccar`), anthocyanins
(`Canth`), brown pigments (`Cbrown`), water (`Cw`), bulk dry matter (`Cm`),
protein (`Cprot`), and carbon-based constituents (`Ccbc`).
"""
struct LeafMultiPigmentProperties{FT<:Real} <: AbstractLeafProperties
    "Leaf structure parameter N >= 1 (fractional plate count)"
    N::FT
    "Chlorophyll a content `[µg cm⁻²]`"
    Cchla::FT
    "Chlorophyll b content `[µg cm⁻²]`"
    Cchlb::FT
    "Carotenoid content `[µg cm⁻²]`"
    Ccar::FT
    "Anthocyanin content (`[nmol cm⁻²]` with Féret `OPTI_2021` `Kant`, or `[µg cm⁻²]` with mass-based `Kant`)"
    Canth::FT
    "Brown pigments content in arbitrary units"
    Cbrown::FT
    "Equivalent water thickness `[cm]`"
    Cw::FT
    "Bulk dry matter content `[g cm⁻²]` (set to 0 when Cprot + Ccbc are used)"
    Cm::FT
    "Protein content `[g cm⁻²]`"
    Cprot::FT
    "Carbon-based constituents content `[g cm⁻²]`"
    Ccbc::FT
end

function LeafMultiPigmentProperties{FT}(;
    N = FT(1.4),
    Cchla = FT(30.0),
    Cchlb = FT(10.0),
    Ccar = FT(10.0),
    Canth = FT(0.5),
    Cbrown = zero(FT),
    Cw = FT(0.012),
    Cm = zero(FT),
    Cprot = FT(0.001),
    Ccbc = FT(0.009),
) where {FT<:Real}
    return LeafMultiPigmentProperties{FT}(
        FT(N), FT(Cchla), FT(Cchlb), FT(Ccar), FT(Canth),
        FT(Cbrown), FT(Cw), FT(Cm), FT(Cprot), FT(Ccbc),
    )
end

function LeafMultiPigmentProperties(;
    N::Real = 1.4,
    Cchla::Real = 30.0,
    Cchlb::Real = 10.0,
    Ccar::Real = 10.0,
    Canth::Real = 0.5,
    Cbrown::Real = 0.0,
    Cw::Real = 0.012,
    Cm::Real = 0.0,
    Cprot::Real = 0.001,
    Ccbc::Real = 0.009,
)
    FT = float(promote_type(
        typeof(N), typeof(Cchla), typeof(Cchlb), typeof(Ccar), typeof(Canth),
        typeof(Cbrown), typeof(Cw), typeof(Cm), typeof(Cprot), typeof(Ccbc),
    ))
    return LeafMultiPigmentProperties{FT}(
        FT(N), FT(Cchla), FT(Cchlb), FT(Ccar), FT(Canth),
        FT(Cbrown), FT(Cw), FT(Cm), FT(Cprot), FT(Ccbc),
    )
end

"""
    MultiPigmentOpticalProperties{FT<:Real, FT2} <: AbstractProspectProperties

Specific absorption cross-section spectra and refractive index `nᵣ` for leaf
models that resolve chlorophyll *a* (`Kchla`) and chlorophyll *b* (`Kchlb`)
separately across a user-supplied wavelength grid `λ`.
"""
struct MultiPigmentOpticalProperties{FT<:Real,FT2} <: AbstractProspectProperties
    "Wavelength grid"
    λ::FT2
    "Refractive index of aqueous leaf interior (> 1)"
    nᵣ::Vector{FT}
    "Specific absorption coefficient of chlorophyll a `[cm² μg⁻¹]`"
    Kchla::Vector{FT}
    "Specific absorption coefficient of chlorophyll b `[cm² μg⁻¹]`"
    Kchlb::Vector{FT}
    "Specific absorption coefficient of carotenoids `[cm² μg⁻¹]`"
    Kcar::Vector{FT}
    "Specific absorption coefficient of anthocyanins (`[cm² nmol⁻¹]` or `[cm² μg⁻¹]` matching `Canth`)"
    Kant::Vector{FT}
    "Specific absorption coefficient of brown pigments"
    Kb::Vector{FT}
    "Specific absorption coefficient of water `[cm⁻¹]`"
    Kw::Vector{FT}
    "Specific absorption coefficient of bulk dry matter `[cm² g⁻¹]`"
    Km::Vector{FT}
    "Specific absorption coefficient of proteins `[cm² g⁻¹]`"
    Kp::Vector{FT}
    "Specific absorption coefficient of carbon-based constituents `[cm² g⁻¹]`"
    Kcbc::Vector{FT}

    function MultiPigmentOpticalProperties{FT,FT2}(
        λ::FT2,
        nᵣ::Vector{FT},
        Kchla::Vector{FT},
        Kchlb::Vector{FT},
        Kcar::Vector{FT},
        Kant::Vector{FT},
        Kb::Vector{FT},
        Kw::Vector{FT},
        Km::Vector{FT},
        Kp::Vector{FT},
        Kcbc::Vector{FT},
    ) where {FT<:Real,FT2}
        nλ = length(nᵣ)
        nλ > 0 || throw(ArgumentError("MultiPigmentOpticalProperties requires non-empty spectra"))
        all(v -> length(v) == nλ, (Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)) ||
            throw(DimensionMismatch("All optical coefficient vectors must have length $(nλ)"))
        length(λ) == nλ ||
            throw(DimensionMismatch("Wavelength grid length $(length(λ)) must match optical vectors $(nλ)"))
        return new{FT,FT2}(λ, nᵣ, Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)
    end
end

function MultiPigmentOpticalProperties(
    λ,
    nᵣ::AbstractVector,
    Kchla::AbstractVector,
    Kchlb::AbstractVector,
    Kcar::AbstractVector,
    Kant::AbstractVector,
    Kb::AbstractVector,
    Kw::AbstractVector,
    Km::AbstractVector,
    Kp::AbstractVector,
    Kcbc::AbstractVector,
)
    FT = float(promote_type(
        eltype(nᵣ), eltype(Kchla), eltype(Kchlb), eltype(Kcar), eltype(Kant),
        eltype(Kb), eltype(Kw), eltype(Km), eltype(Kp), eltype(Kcbc),
    ))
    return MultiPigmentOpticalProperties{FT,typeof(λ)}(
        λ,
         Vector{FT}(nᵣ), Vector{FT}(Kchla), Vector{FT}(Kchlb), Vector{FT}(Kcar),
        Vector{FT}(Kant), Vector{FT}(Kb), Vector{FT}(Kw), Vector{FT}(Km),
        Vector{FT}(Kp), Vector{FT}(Kcbc),
    )
end

function MultiPigmentOpticalProperties(
    base::LeafOpticalProperties{FT,FT2};
    Kchla::AbstractVector = base.Kcab,
    Kchlb::AbstractVector = base.Kcab,
    Kcar::AbstractVector = base.Kcar,
    Kant::AbstractVector = base.Kant,
    Kb::AbstractVector = base.Kb,
    Kw::AbstractVector = base.Kw,
    Km::AbstractVector = base.Km,
    Kp::AbstractVector = base.Kp,
    Kcbc::AbstractVector = base.Kcbc,
    nᵣ::AbstractVector = base.nᵣ,
) where {FT<:Real,FT2}
    OFT = float(promote_type(
        FT, eltype(nᵣ), eltype(Kchla), eltype(Kchlb), eltype(Kcar), eltype(Kant),
        eltype(Kb), eltype(Kw), eltype(Km), eltype(Kp), eltype(Kcbc),
    ))
    return MultiPigmentOpticalProperties{OFT,FT2}(
        base.λ,
        Vector{OFT}(nᵣ), Vector{OFT}(Kchla), Vector{OFT}(Kchlb), Vector{OFT}(Kcar),
        Vector{OFT}(Kant), Vector{OFT}(Kb), Vector{OFT}(Kw), Vector{OFT}(Km),
        Vector{OFT}(Kp), Vector{OFT}(Kcbc),
    )
end

function prospect_constituent_coefficients(
    leaf::LeafMultiPigmentProperties,
    optics::MultiPigmentOpticalProperties,
)
    return (;
        chlorophyll_a = leaf.Cchla .* optics.Kchla,
        chlorophyll_b = leaf.Cchlb .* optics.Kchlb,
        carotenoids = leaf.Ccar .* optics.Kcar,
        anthocyanins = leaf.Canth .* optics.Kant,
        brown_material = leaf.Cbrown .* optics.Kb,
        water = leaf.Cw .* optics.Kw,
        dry_matter = leaf.Cm .* optics.Km,
        protein = leaf.Cprot .* optics.Kp,
        carbon_constituents = leaf.Ccbc .* optics.Kcbc,
    )
end

function prospect_constituent_coefficients(
    leaf::LeafMultiPigmentProperties,
    optics::LeafOpticalProperties,
)
    return (;
        chlorophyll_a = leaf.Cchla .* optics.Kcab,
        chlorophyll_b = leaf.Cchlb .* optics.Kcab,
        carotenoids = leaf.Ccar .* optics.Kcar,
        anthocyanins = leaf.Canth .* optics.Kant,
        brown_material = leaf.Cbrown .* optics.Kb,
        water = leaf.Cw .* optics.Kw,
        dry_matter = leaf.Cm .* optics.Km,
        protein = leaf.Cprot .* optics.Kp,
        carbon_constituents = leaf.Ccbc .* optics.Kcbc,
    )
end

function prospect_constituent_coefficients(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

"""
    constituent_optical_depths(leaf, optics, idx::Integer)

Scalar concentration × specific-absorption optical depths at wavelength index
`idx`, sharing a common `FT` for allocation-free use with
`partition_absorption(HomogeneousAbsorption(), absorbed[idx], coeffs)`.
"""
@inline function constituent_optical_depths(
    leaf::LeafProspectProProperties{FTL},
    optics::LeafOpticalProperties{FTO},
    idx::Integer,
) where {FTL<:Real,FTO<:Real}
    FT = promote_type(FTL, FTO)
    @boundscheck checkbounds(optics.nᵣ, idx)
    @inbounds return (;
        chlorophyll = FT(leaf.Ccab) * FT(optics.Kcab[idx]),
        carotenoids = FT(leaf.Ccar) * FT(optics.Kcar[idx]),
        anthocyanins = FT(leaf.Canth) * FT(optics.Kant[idx]),
        brown_material = FT(leaf.Cbrown) * FT(optics.Kb[idx]),
        water = FT(leaf.Cw) * FT(optics.Kw[idx]),
        dry_matter = FT(leaf.Cm) * FT(optics.Km[idx]),
        protein = FT(leaf.Cprot) * FT(optics.Kp[idx]),
        carbon_constituents = FT(leaf.Ccbc) * FT(optics.Kcbc[idx]),
    )
end

@inline function constituent_optical_depths(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties,
    ::Integer,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

@inline function constituent_optical_depths(
    leaf::LeafMultiPigmentProperties{FTL},
    optics::MultiPigmentOpticalProperties{FTO},
    idx::Integer,
) where {FTL<:Real,FTO<:Real}
    FT = promote_type(FTL, FTO)
    @boundscheck checkbounds(optics.nᵣ, idx)
    @inbounds return (;
        chlorophyll_a = FT(leaf.Cchla) * FT(optics.Kchla[idx]),
        chlorophyll_b = FT(leaf.Cchlb) * FT(optics.Kchlb[idx]),
        carotenoids = FT(leaf.Ccar) * FT(optics.Kcar[idx]),
        anthocyanins = FT(leaf.Canth) * FT(optics.Kant[idx]),
        brown_material = FT(leaf.Cbrown) * FT(optics.Kb[idx]),
        water = FT(leaf.Cw) * FT(optics.Kw[idx]),
        dry_matter = FT(leaf.Cm) * FT(optics.Km[idx]),
        protein = FT(leaf.Cprot) * FT(optics.Kp[idx]),
        carbon_constituents = FT(leaf.Ccbc) * FT(optics.Kcbc[idx]),
    )
end

@inline function constituent_optical_depths(
    leaf::LeafMultiPigmentProperties{FTL},
    optics::LeafOpticalProperties{FTO},
    idx::Integer,
) where {FTL<:Real,FTO<:Real}
    FT = promote_type(FTL, FTO)
    @boundscheck checkbounds(optics.nᵣ, idx)
    @inbounds return (;
        chlorophyll_a = FT(leaf.Cchla) * FT(optics.Kcab[idx]),
        chlorophyll_b = FT(leaf.Cchlb) * FT(optics.Kcab[idx]),
        carotenoids = FT(leaf.Ccar) * FT(optics.Kcar[idx]),
        anthocyanins = FT(leaf.Canth) * FT(optics.Kant[idx]),
        brown_material = FT(leaf.Cbrown) * FT(optics.Kb[idx]),
        water = FT(leaf.Cw) * FT(optics.Kw[idx]),
        dry_matter = FT(leaf.Cm) * FT(optics.Km[idx]),
        protein = FT(leaf.Cprot) * FT(optics.Kp[idx]),
        carbon_constituents = FT(leaf.Ccbc) * FT(optics.Kcbc[idx]),
    )
end

"""
    MultiChlorophyllAntenna(; chlb_transfer = 1, carotenoid_transfer = 0)

Antenna photon-selection rule for [`LeafMultiPigmentProperties`](@ref) partitions
containing `:chlorophyll_a` and `:chlorophyll_b` (and `:carotenoids` when
`carotenoid_transfer > 0`). Uses [`complementary_excitation`](@ref) to preserve
exact total-selection and zero-pigment limits.
"""
struct MultiChlorophyllAntenna{FT<:Real} <: AbstractAntennaSelection
    chlb_transfer::FT
    carotenoid_transfer::FT

    function MultiChlorophyllAntenna{FT}(chlb::FT, car::FT) where {FT<:Real}
        b_val = prospect_primal(chlb)
        c_val = prospect_primal(car)
        BFT = typeof(b_val)
        isfinite(b_val) && zero(BFT) <= b_val <= one(BFT) ||
            throw(ArgumentError("chlb_transfer must be finite and in [0, 1]"))
        isfinite(c_val) && zero(BFT) <= c_val <= one(BFT) ||
            throw(ArgumentError("carotenoid_transfer must be finite and in [0, 1]"))
        return new{FT}(chlb, car)
    end
end

function MultiChlorophyllAntenna(
    chlb_transfer::Real,
    carotenoid_transfer::Real = zero(typeof(chlb_transfer)),
)
    FT = float(promote_type(typeof(chlb_transfer), typeof(carotenoid_transfer)))
    return MultiChlorophyllAntenna{FT}(FT(chlb_transfer), FT(carotenoid_transfer))
end

function MultiChlorophyllAntenna{FT}(;
    chlb_transfer = one(FT),
    carotenoid_transfer = zero(FT),
) where {FT<:Real}
    return MultiChlorophyllAntenna{FT}(FT(chlb_transfer), FT(carotenoid_transfer))
end

function MultiChlorophyllAntenna(;
    chlb_transfer::Real = 1.0,
    carotenoid_transfer::Real = zero(typeof(chlb_transfer)),
)
    return MultiChlorophyllAntenna(chlb_transfer, carotenoid_transfer)
end

@inline function select_absorbed_photons(
    selection::MultiChlorophyllAntenna,
    partition::AbsorptionPartitionResult{FT},
) where {FT}
    fractions = partition.fractions
    haskey(fractions, :chlorophyll_a) && haskey(fractions, :chlorophyll_b) || return FT(NaN)
    f_b = FT(selection.chlb_transfer)
    f_car = FT(selection.carotenoid_transfer)
    has_car = haskey(fractions, :carotenoids)
    prospect_primal(f_car) > zero(prospect_primal(f_car)) && !has_car && return FT(NaN)
    partition.status == AbsorptionNoMaterial && return zero(FT)

    selected = fractions.chlorophyll_a + f_b * fractions.chlorophyll_b +
               (has_car ? f_car * fractions.carotenoids : zero(FT))
    unselected = sum(
        map(keys(fractions), values(fractions)) do name, value
            name === :chlorophyll_a ? zero(value) :
            name === :chlorophyll_b ? (one(f_b) - f_b) * value :
            name === :carotenoids ? (one(f_car) - f_car) * value : value
        end,
    )
    return complementary_excitation(partition, selected, unselected)
end

function _validate_leaf_contents(contents::Tuple)
    for x in contents
        v = prospect_primal(x)
        isfinite(v) && v >= zero(v) ||
            throw(ArgumentError("Leaf constituent contents must be finite and nonnegative"))
    end
    return nothing
end

function _validate_optics_arrays(optics::LeafOpticalProperties)
    (; Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc, nᵣ) = optics
    nonneg(x) = (v = prospect_primal(x); isfinite(v) && v >= zero(v))
    all(arr -> all(nonneg, arr), (Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)) ||
        throw(ArgumentError("Optical absorption coefficients must be finite and nonnegative"))
    all(x -> (v = prospect_primal(x); isfinite(v) && v > one(v)), nᵣ) ||
        throw(ArgumentError("Leaf refractive indices must be finite and > 1"))
    return nothing
end

function _validate_optics_arrays(optics::MultiPigmentOpticalProperties)
    (; Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc, nᵣ) = optics
    nonneg(x) = (v = prospect_primal(x); isfinite(v) && v >= zero(v))
    all(arr -> all(nonneg, arr), (Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)) ||
        throw(ArgumentError("Optical absorption coefficients must be finite and nonnegative"))
    all(x -> (v = prospect_primal(x); isfinite(v) && v > one(v)), nᵣ) ||
        throw(ArgumentError("Leaf refractive indices must be finite and > 1"))
    return nothing
end

@inline function _total_optical_depth_vector(
    leaf::LeafProspectProProperties,
    optics::LeafOpticalProperties,
)
    p = prospect_constituent_coefficients(leaf, optics)
    return p.chlorophyll + p.carotenoids + p.anthocyanins + p.brown_material +
           p.water + p.dry_matter + p.protein + p.carbon_constituents
end

function _total_optical_depth_vector(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

@inline function _total_optical_depth_vector(
    leaf::LeafMultiPigmentProperties,
    optics::AbstractProspectProperties,
)
    p = prospect_constituent_coefficients(leaf, optics)
    return p.chlorophyll_a + p.chlorophyll_b + p.carotenoids + p.anthocyanins +
           p.brown_material + p.water + p.dry_matter + p.protein + p.carbon_constituents
end

"""
    prospect_budget(leaf, optics, surface::AbstractLeafSurface; μ_in = one(FT))

Generalized PROSPECT plate budget accepting any [`AbstractLeafSurface`](@ref)
(`ConeSurface` or `RoughDielectricSurface`) and optional directional incidence
cosine `μ_in ∈ (0, 1]`, for either [`LeafProspectProProperties`](@ref) or
[`LeafMultiPigmentProperties`](@ref).

Note: First-surface entry uses `surface_interface_optics(surface, μ_in, nr)` and
symmetric diffuse plate interfaces use `t_int_diffuse = t_hem_ext / nr²` across all
`N >= 1` plates (smoothly supporting `1 <= N < 2`), reducing identically to legacy
`prospect_budget` when `surface == ConeSurface(40, 1)`. For distinct outer-epidermis
vs internal aqueous cell-wall interfaces, see [`DorsiventralLeafProperties`](@ref).
Returns `(; transmission, reflectance, absorptance, reflectance_volume)`.
"""
function prospect_budget(
    leaf::Union{LeafProspectProProperties{FTL},LeafMultiPigmentProperties{FTL}},
    optics::AbstractProspectProperties,
    surface::AbstractLeafSurface{FTS};
    μ_in::Real = one(_leaf_base_float(promote_type(FTL, eltype(optics.nᵣ), FTS))),
) where {FTL<:Real,FTS<:Real}
    n_val = prospect_primal(leaf.N)
    isfinite(n_val) && n_val >= one(n_val) ||
        throw(ArgumentError("PROSPECT requires finite N >= 1"))
    if leaf isa LeafProspectProProperties
        _validate_leaf_contents((
            leaf.Ccab, leaf.Ccar, leaf.Canth, leaf.Cbrown,
            leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
        ))
    else
        _validate_leaf_contents((
            leaf.Cchla, leaf.Cchlb, leaf.Ccar, leaf.Canth, leaf.Cbrown,
            leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
        ))
    end
    _validate_optics_arrays(optics)

    K_total = _total_optical_depth_vector(leaf, optics)
    FT = float(promote_type(eltype(K_total), typeof(leaf.N), eltype(optics.nᵣ), FTS, typeof(μ_in)))
    BFT = _leaf_base_float(FT)
    μ_val = prospect_primal(μ_in)
    isfinite(μ_val) && zero(BFT) < μ_val <= one(BFT) ||
        throw(ArgumentError("μ_in must be finite and in (0, 1]"))

    nλ = length(optics.nᵣ)
    T_vec = Vector{FT}(undef, nλ)
    R_vec = Vector{FT}(undef, nλ)
    A_vec = Vector{FT}(undef, nλ)
    Rvol_vec = Vector{FT}(undef, nλ)

    N_ft = FT(leaf.N)
    μ_ft = FT(μ_in)

    @inbounds for idx in 1:nλ
        nr = FT(optics.nᵣ[idx])
        k_unit = FT(K_total[idx]) / N_ft
        surf_opt = surface_interface_optics(surface, μ_ft, nr)
        b = _asymmetric_top_prospect_budget(
            k_unit, N_ft, surf_opt.t_ext, surf_opt.t_int_diffuse, nr,
        )
        T_vec[idx] = b.transmission
        R_vec[idx] = b.reflectance
        A_vec[idx] = b.absorptance
        Rvol_vec[idx] = b.reflectance_volume
    end

    return (;
        transmission = T_vec,
        reflectance = R_vec,
        absorptance = A_vec,
        reflectance_volume = Rvol_vec,
    )
end

function prospect_budget(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

function prospect(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

function prospect_budget(
    leaf::LeafMultiPigmentProperties{FTL},
    optics::AbstractProspectProperties,
) where {FTL<:Real}
    FT = float(promote_type(FTL, eltype(optics.nᵣ)))
    b = prospect_budget(leaf, optics, ConeSurface{FT}(FT(40), one(FT)))
    return (;
        transmission = b.transmission,
        reflectance = b.reflectance,
        absorptance = b.absorptance,
    )
end

function prospect(
    leaf::LeafMultiPigmentProperties,
    optics::AbstractProspectProperties,
)
    b = prospect_budget(leaf, optics)
    return b.transmission, b.reflectance
end

function prospect(
    leaf::Union{LeafProspectProProperties,LeafMultiPigmentProperties},
    optics::AbstractProspectProperties,
    surface::AbstractLeafSurface;
    kwargs...,
)
    b = prospect_budget(leaf, optics, surface; kwargs...)
    return b.transmission, b.reflectance
end

"""
    LeafSIPProperties{FT<:Real} <: AbstractLeafProperties

Spectral-invariants leaf model state following Wu et al. (2020, *Remote Sens.
Environ.* 250, 112131, `SIP_Model.m`), with explicit empirical structural slopes
(`alpha_lma`, `f_lma_slope`, `g_lma_slope`) so callers can supply TOML-configured
or calibrated invariant parameters without hardcoded magic numbers.
"""
Base.@kwdef struct LeafSIPProperties{FT<:Real} <: AbstractLeafProperties
    "Chlorophyll a+b content `[µg cm⁻²]` (set 0 when Cchla/Cchlb are used)"
    Ccab::FT        = FT(40.0)
    "Optional separate Chlorophyll a content `[µg cm⁻²]` (requires Ccab == 0)"
    Cchla::FT       = zero(FT)
    "Optional separate Chlorophyll b content `[µg cm⁻²]` (requires Ccab == 0)"
    Cchlb::FT       = zero(FT)
    "Carotenoid content `[µg cm⁻²]`"
    Ccar::FT        = FT(10.0)
    "Anthocyanin content (`[nmol cm⁻²]` or `[µg cm⁻²]` matching `Kant`)"
    Canth::FT       = FT(0.5)
    "Brown pigments content"
    Cbrown::FT      = zero(FT)
    "Equivalent water thickness `[cm]`"
    Cw::FT          = FT(0.012)
    "Bulk dry matter content `[g cm⁻²]` (set 0 when Cprot + Ccbc are used)"
    Cm::FT          = FT(0.008)
    "Protein content `[g cm⁻²]` (requires Cm == 0)"
    Cprot::FT       = zero(FT)
    "Carbon-based constituents content `[g cm⁻²]` (requires Cm == 0)"
    Ccbc::FT        = zero(FT)
    "Scattering element scaling factor `[cm² g⁻¹]` (Wu et al. 2020: 600)"
    alpha_lma::FT   = FT(600.0)
    "Recollision probability LMA slope `[cm² g⁻¹]` (Wu et al. 2020: 2765.0)"
    f_lma_slope::FT = FT(2765.0)
    "Asymmetry parameter LMA slope `[cm² g⁻¹]` (Wu et al. 2020: 102.8)"
    g_lma_slope::FT = FT(102.8)
end

"""
    leaf_sip_budget(leaf::LeafSIPProperties, optics; surface = nothing, μ_in = 1)

Evaluate the Wu et al. (2020) `Leaf-SIP` spectral-invariants leaf optical budget
(`surface = nothing`, matching published `SIP_Model.m` which omits cuticular
Fresnel reflection) or a surface-coupled `Leaf-SIP+Fresnel` budget when an
[`AbstractLeafSurface`](@ref) is supplied (`surface = ConeSurface(...)` or
`RoughDielectricSurface(...)`).

Note: Wu et al. (2020)'s default slopes (`alpha_lma = 600`, `f_lma_slope = 2765`,
`g_lma_slope = 102.8`) were empirically calibrated for bare `surface = nothing`;
adding an external `surface` restores physical blue/red cuticular reflectance
floors, and `alpha_lma`/`f_lma_slope`/`g_lma_slope` can be recalibrated when
fitting NIR/SWIR plateaus jointly with surface coupling.

Computes absorptance directly via `-expm1(-k_elem)` without `1 - R - T`
subtraction, preserving exact `R + T + A == 1` into the nonabsorbing limit.
"""
function leaf_sip_budget(
    leaf::LeafSIPProperties{FTL},
    optics::AbstractProspectProperties;
    surface::Union{Nothing,AbstractLeafSurface} = nothing,
    μ_in::Real = one(_leaf_base_float(promote_type(FTL, eltype(optics.nᵣ)))),
) where {FTL<:Real}
    _validate_leaf_contents((
        leaf.Ccab, leaf.Cchla, leaf.Cchlb, leaf.Ccar, leaf.Canth, leaf.Cbrown,
        leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
    ))
    _validate_optics_arrays(optics)

    BFT_L = _leaf_base_float(FTL)
    if prospect_primal(leaf.Ccab) > zero(BFT_L) &&
       (prospect_primal(leaf.Cchla) > zero(BFT_L) || prospect_primal(leaf.Cchlb) > zero(BFT_L))
        throw(ArgumentError("LeafSIPProperties cannot set both Ccab > 0 and Cchla/Cchlb > 0 (double-counts chlorophyll); set Ccab = 0 when using Cchla/Cchlb"))
    end
    if prospect_primal(leaf.Ccab) > zero(BFT_L) && optics isa MultiPigmentOpticalProperties
        throw(ArgumentError("LeafSIPProperties with Ccab > 0 cannot be paired with MultiPigmentOpticalProperties (ambiguous Kchla vs Kchlb); set Ccab = 0 and specify Cchla, Cchlb or use LeafOpticalProperties"))
    end
    if prospect_primal(leaf.Cm) > zero(BFT_L) &&
       (prospect_primal(leaf.Cprot) > zero(BFT_L) || prospect_primal(leaf.Ccbc) > zero(BFT_L))
        throw(ArgumentError("LeafSIPProperties cannot set both Cm > 0 and Cprot/Ccbc > 0 (double-counts LMA); set Cm = 0 when using Cprot + Ccbc"))
    end

    lma_eff = leaf.Cm + leaf.Cprot + leaf.Ccbc
    prospect_primal(lma_eff) > zero(prospect_primal(lma_eff)) ||
        throw(ArgumentError("LeafSIPProperties requires positive dry mass (Cm + Cprot + Ccbc > 0)"))
    prospect_primal(leaf.alpha_lma) > zero(prospect_primal(leaf.alpha_lma)) ||
        throw(ArgumentError("alpha_lma must be positive"))
    prospect_primal(leaf.f_lma_slope) > zero(prospect_primal(leaf.f_lma_slope)) ||
        throw(ArgumentError("f_lma_slope must be positive"))
    prospect_primal(leaf.g_lma_slope) >= zero(prospect_primal(leaf.g_lma_slope)) ||
        throw(ArgumentError("g_lma_slope must be nonnegative"))

    FTS = surface === nothing ? FTL : _surface_ft(surface)
    FT = float(promote_type(FTL, eltype(optics.nᵣ), FTS, typeof(μ_in)))
    BFT = _leaf_base_float(FT)
    μ_ft = FT(μ_in)
    zero(BFT) < prospect_primal(μ_ft) <= one(BFT) ||
        throw(ArgumentError("μ_in must be in (0, 1]"))

    f_lma = FT(leaf.f_lma_slope) * FT(lma_eff)
    g_lma = FT(leaf.g_lma_slope) * FT(lma_eff)
    # p = 1 - (1 - exp(-f_lma)) / f_lma = 1 + expm1(-f_lma) / f_lma
    p = one(FT) + expm1(-f_lma) / f_lma
    one_minus_p = -expm1(-f_lma) / f_lma
    two = one(FT) + one(FT)
    half = inv(two)
    q = one(FT) - two * exp(-g_lma)
    # Smooth |q| away from zero; at q == 0 use zero(FT)
    q_abs = prospect_primal(q) >= zero(BFT) ? q : -q

    nλ = length(optics.nᵣ)
    T_vec = Vector{FT}(undef, nλ)
    R_vec = Vector{FT}(undef, nλ)
    A_vec = Vector{FT}(undef, nλ)

    elem_scale = FT(lma_eff) * FT(leaf.alpha_lma)

    @inbounds for idx in 1:nλ
        k_chl = if optics isa MultiPigmentOpticalProperties
            FT(leaf.Cchla) * FT(optics.Kchla[idx]) + FT(leaf.Cchlb) * FT(optics.Kchlb[idx])
        else
            (FT(leaf.Ccab) + FT(leaf.Cchla) + FT(leaf.Cchlb)) * FT(optics.Kcab[idx])
        end
        k_tot = k_chl +
                FT(leaf.Ccar) * FT(optics.Kcar[idx]) +
                FT(leaf.Canth) * FT(optics.Kant[idx]) +
                FT(leaf.Cbrown) * FT(optics.Kb[idx]) +
                FT(leaf.Cw) * FT(optics.Kw[idx]) +
                FT(leaf.Cm) * FT(optics.Km[idx]) +
                FT(leaf.Cprot) * FT(optics.Kp[idx]) +
                FT(leaf.Ccbc) * FT(optics.Kcbc[idx])

        k_elem = k_tot / elem_scale
        loss0 = -expm1(-k_elem)
        w0 = one(FT) - loss0
        # 1 - p * w0 = (1 - p) + p * loss0 (strictly positive, cancellation-free)
        denom = one_minus_p + p * loss0
        w_int = w0 * one_minus_p / denom
        a_int = loss0 / denom
        asym = half * q * denom / (one(FT) - q_abs * p * w0)
        r_int = w_int * (half + asym)
        t_int = w_int * (half - asym)

        if surface === nothing
            T_vec[idx] = t_int
            R_vec[idx] = r_int
            A_vec[idx] = a_int
        else
            nr = FT(optics.nᵣ[idx])
            s_top = surface_interface_optics(surface, μ_ft, nr)
            s_hem = surface_interface_optics(surface, FT(90), nr; hemispherical = true)
            # Couple external dielectric interfaces I_top ⋆ M_sip ⋆ I_bot without air gaps
            b = _couple_external_interfaces_to_medium(
                s_top.r_ext, s_top.t_ext, s_top.r_int_diffuse, s_top.t_int_diffuse,
                r_int, t_int, a_int, r_int, t_int, a_int,
                s_hem.r_int_diffuse, s_hem.t_int_diffuse,
            )
            T_vec[idx] = b.transmission
            R_vec[idx] = b.reflectance
            A_vec[idx] = b.absorptance
        end
    end

    return (; transmission = T_vec, reflectance = R_vec, absorptance = A_vec)
end
