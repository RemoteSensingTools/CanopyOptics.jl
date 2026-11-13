"""
    LeafMultiPigmentProperties{FT<:Real} <: AbstractLeafProperties

Leaf biochemical and structural state with separate chlorophyll *a* (`Cchla`)
and chlorophyll *b* (`Cchlb`) pools alongside carotenoids (`Ccar`), anthocyanins
(`Canth`), brown pigments (`Cbrown`), water (`Cw`), bulk dry matter (`Cm`),
protein (`Cprot`), and carbon-based constituents (`Ccbc`).

When constructed with keywords:
- Omitting `Cm` selects **PROSPECT-PRO dry-matter mode** (`Cm = 0`, defaulting `Cprot = 0.001`, `Ccbc = 0.009`).
- Supplying `Cm` explicitly (including `Cm = 0` or `ForwardDiff.Dual`) selects **bulk `Cm` mode**,
  defaulting `Cprot = 0, Ccbc = 0` without double-counting dry mass.
- Setting both `Cm > 0` and `Cprot > 0` (or `Ccbc > 0`) throws an `ArgumentError`.
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
    "Anthocyanin content `[µg cm⁻²]` (mass basis matching `OPTI_2021` and `MultiPigmentOpticalProperties`)"
    Canth::FT
    "Brown pigments content in arbitrary units"
    Cbrown::FT
    "Equivalent water thickness `[cm]`"
    Cw::FT
    "Bulk dry matter content `[g cm⁻²]` (0 when `Cprot + Ccbc` are used)"
    Cm::FT
    "Protein content `[g cm⁻²]`"
    Cprot::FT
    "Carbon-based constituents content `[g cm⁻²]`"
    Ccbc::FT

    function LeafMultiPigmentProperties{FT}(
        N, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
    ) where {FT<:Real}
        return new{FT}(
            FT(N), FT(Cchla), FT(Cchlb), FT(Ccar), FT(Canth),
            FT(Cbrown), FT(Cw), FT(Cm), FT(Cprot), FT(Ccbc),
        )
    end
end

function LeafMultiPigmentProperties(
    N::Real, Cchla::Real, Cchlb::Real, Ccar::Real, Canth::Real,
    Cbrown::Real, Cw::Real, Cm::Real, Cprot::Real, Ccbc::Real,
)
    FT = _canopy_parameter_type(
        N, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
    )
    return LeafMultiPigmentProperties{FT}(
        N, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
    )
end

function LeafMultiPigmentProperties{FT}(;
    N = _UNSET_LEAF_KW,
    Cchla = _UNSET_LEAF_KW,
    Cchlb = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    n_v = N === _UNSET_LEAF_KW ? FT(BFT(1.4)) : FT(N)
    ca_v = Cchla === _UNSET_LEAF_KW ? FT(BFT(30.0)) : FT(Cchla)
    cb_v = Cchlb === _UNSET_LEAF_KW ? FT(BFT(10.0)) : FT(Cchlb)
    car_v = Ccar === _UNSET_LEAF_KW ? FT(BFT(10.0)) : FT(Ccar)
    ant_v = Canth === _UNSET_LEAF_KW ? FT(BFT(0.5)) : FT(Canth)
    brn_v = Cbrown === _UNSET_LEAF_KW ? zero(FT) : FT(Cbrown)
    cw_v = Cw === _UNSET_LEAF_KW ? FT(BFT(0.012)) : FT(Cw)
    cm_v, cprot_v, ccbc_v = _resolve_pro_dry_matter(FT, Cm, Cprot, Ccbc, "LeafMultiPigmentProperties")
    return LeafMultiPigmentProperties{FT}(
        n_v, ca_v, cb_v, car_v, ant_v, brn_v, cw_v, cm_v, cprot_v, ccbc_v,
    )
end

function LeafMultiPigmentProperties(;
    N = _UNSET_LEAF_KW,
    Cchla = _UNSET_LEAF_KW,
    Cchlb = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
)
    T_raw = promote_type(
        _kw_type(N), _kw_type(Cchla), _kw_type(Cchlb), _kw_type(Ccar), _kw_type(Canth),
        _kw_type(Cbrown), _kw_type(Cw), _kw_type(Cm), _kw_type(Cprot), _kw_type(Ccbc),
    )
    FT = T_raw === Union{} ? Float64 : float(T_raw)
    return LeafMultiPigmentProperties{FT}(;
        N, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
    )
end

"""
    MultiPigmentOpticalProperties(λ_bnds = (399.5:2500.5)u"nm"; basis = :calibrated_in_vivo, FT = ...)
    MultiPigmentOpticalProperties{FT}(λ_bnds = (399.5:2500.5)u"nm"; basis = :calibrated_in_vivo)
    MultiPigmentOpticalProperties(base::LeafOpticalProperties; Kchla = ..., Kchlb = ..., ...)
    MultiPigmentOpticalProperties{FT<:Real, FT2} <: AbstractProspectProperties

Specific absorption cross-section spectra and refractive index `nᵣ` for leaf
models that resolve chlorophyll *a* (`Kchla`) and chlorophyll *b* (`Kchlb`)
separately (`[cm² µg⁻¹]`) across spectral bin edges `λ_bnds` or a custom grid `λ`.

Supported `basis` options when loading from bundled data:
- `:calibrated_in_vivo` (default): empirical in-vivo `Kchla, Kchlb` calibrated on
  `ANGERS` (`400–800 nm`, `data/multi_pigment_chlab_in_vivo.csv`), spliced to
  `OPTI_2021` for `λ > 800 nm`.
- `:prospect_mp_2017`: Zhang et al. (2017, *Sci. Rep.* 7:6429) pseudo-Voigt
  `Kchla, Kchlb, Kcar` on `400–800 nm`, spliced to `OPTI_2021` for `λ > 800 nm`.
- `:prospect_mp_plus_2022`: Zhang et al. (2022, *Sensors* 22:3025) pseudo-Voigt
  `Kchla, Kchlb, Kcar, Kant` on `400–800 nm`, spliced to `OPTI_2021` for `λ > 800 nm`.
- `:prospect_pro_pooled`: baseline control with `Kchla == Kchlb == Kcab` from `OPTI_2021`.
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
    "Specific absorption coefficient of anthocyanins `[cm² μg⁻¹]`"
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
    ) where {FT<:Real,FT2}
        nλ = length(nᵣ)
        nλ > 0 || throw(ArgumentError("MultiPigmentOpticalProperties requires non-empty spectra"))
        all(v -> length(v) == nλ, (Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)) ||
            throw(DimensionMismatch("All optical coefficient vectors must have length $(nλ)"))
        length(λ) == nλ ||
            throw(DimensionMismatch("Wavelength grid length $(length(λ)) must match optical vectors $(nλ)"))
        conv(v) = v isa Vector{FT} ? v : copyto!(Vector{FT}(undef, nλ), v)
        return new{FT,FT2}(
            λ, conv(nᵣ), conv(Kchla), conv(Kchlb), conv(Kcar),
            conv(Kant), conv(Kb), conv(Kw), conv(Km), conv(Kp), conv(Kcbc),
        )
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
        λ, nᵣ, Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc,
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
        base.λ, nᵣ, Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc,
    )
end

function Base.:(==)(a::MultiPigmentOpticalProperties, b::MultiPigmentOpticalProperties)
    return a.λ == b.λ && a.nᵣ == b.nᵣ && a.Kchla == b.Kchla && a.Kchlb == b.Kchlb &&
           a.Kcar == b.Kcar && a.Kant == b.Kant && a.Kb == b.Kb && a.Kw == b.Kw &&
           a.Km == b.Km && a.Kp == b.Kp && a.Kcbc == b.Kcbc
end

function Base.hash(a::MultiPigmentOpticalProperties, h::UInt)
    h = hash(:MultiPigmentOpticalProperties, h)
    for f in (a.λ, a.nᵣ, a.Kchla, a.Kchlb, a.Kcar, a.Kant, a.Kb, a.Kw, a.Km, a.Kp, a.Kcbc)
        h = hash(f, h)
    end
    return h
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
`carotenoid_transfer > 0`). Uses `complementary_excitation` to preserve
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

@doc raw"""
    prospect_budget(leaf, optics, surface::AbstractLeafSurface;
                    μ_in = one(FT), n_surface = nothing)

Compute a generalized PROSPECT plate budget accepting any [`AbstractLeafSurface`](@ref)
([`ConeSurface`](@ref) or [`RoughDielectricSurface`](@ref)), optional directional
incidence cosine ``\mu_{\mathrm{in}} \in (0, 1]``, and optional outer cuticular
refractive-index spectrum `n_surface` (e.g., [`load_wax_refractive_index`](@ref)),
for either [`LeafProspectProProperties`](@ref) or [`LeafMultiPigmentProperties`](@ref).

First-surface entry uses `surface_interface_optics(surface, μ_in, nr; n_surface)` and
symmetric diffuse plate interfaces use ``t_{\mathrm{int,diffuse}} = t_{\mathrm{hem,ext}} / n_r^2``
across all ``N \ge 1`` plates (smoothly supporting ``1 \le N < 2``), reducing
identically to `prospect_budget(leaf, optics)` when `surface == ConeSurface(40, 1)`
and `n_surface === nothing`. For distinct outer-epidermis vs internal aqueous
cell-wall interfaces, see [`DorsiventralLeafProperties`](@ref). Returns
`(; transmission, reflectance, absorptance, reflectance_volume)`.
"""
function prospect_budget(
    leaf::Union{LeafProspectProProperties{FTL},LeafMultiPigmentProperties{FTL}},
    optics::AbstractProspectProperties,
    surface::AbstractLeafSurface{FTS};
    μ_in::Real = one(_leaf_base_float(promote_type(FTL, eltype(optics.nᵣ), FTS))),
    n_surface::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {FTL<:Real,FTS<:Real}
    n_val = prospect_primal(leaf.N)
    isfinite(n_val) && n_val >= one(n_val) ||
        throw(ArgumentError("PROSPECT requires finite N >= 1"))
    if leaf isa LeafProspectProProperties
        _validate_leaf_contents((
            leaf.Ccab, leaf.Ccar, leaf.Canth, leaf.Cbrown,
            leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
        ))
        BFT_L = _leaf_base_float(FTL)
        if prospect_primal(leaf.Cm) > zero(BFT_L) &&
           (prospect_primal(leaf.Cprot) > zero(BFT_L) || prospect_primal(leaf.Ccbc) > zero(BFT_L))
            throw(ArgumentError("LeafProspectProProperties cannot set both Cm > 0 and Cprot/Ccbc > 0 (double-counts LMA)"))
        end
    else
        _validate_leaf_contents((
            leaf.Cchla, leaf.Cchlb, leaf.Ccar, leaf.Canth, leaf.Cbrown,
            leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
        ))
        BFT_L = _leaf_base_float(FTL)
        if prospect_primal(leaf.Cm) > zero(BFT_L) &&
           (prospect_primal(leaf.Cprot) > zero(BFT_L) || prospect_primal(leaf.Ccbc) > zero(BFT_L))
            throw(ArgumentError("LeafMultiPigmentProperties cannot set both Cm > 0 and Cprot/Ccbc > 0 (double-counts LMA)"))
        end
    end
    _validate_optics_arrays(optics)
    nλ = length(optics.nᵣ)
    if n_surface !== nothing
        length(n_surface) == nλ ||
            throw(DimensionMismatch("n_surface length $(length(n_surface)) must match optics length $(nλ)"))
        all(x -> (v = prospect_primal(x); isfinite(v) && v > one(v)), n_surface) ||
            throw(ArgumentError("n_surface refractive indices must be finite and > 1"))
    end

    K_total = _total_optical_depth_vector(leaf, optics)
    FTW = n_surface === nothing ? eltype(optics.nᵣ) : eltype(n_surface)
    FT = float(promote_type(eltype(K_total), typeof(leaf.N), eltype(optics.nᵣ), FTS, typeof(μ_in), FTW))
    BFT = _leaf_base_float(FT)
    μ_val = prospect_primal(μ_in)
    isfinite(μ_val) && zero(BFT) < μ_val <= one(BFT) ||
        throw(ArgumentError("μ_in must be finite and in (0, 1]"))

    T_vec = Vector{FT}(undef, nλ)
    R_vec = Vector{FT}(undef, nλ)
    A_vec = Vector{FT}(undef, nλ)
    Rvol_vec = Vector{FT}(undef, nλ)

    N_ft = FT(leaf.N)
    μ_ft = FT(μ_in)
    prep_surf = _prepare_surface_cheb(FT, surface, μ_ft, false)
    ns_off = n_surface === nothing ? 0 : firstindex(n_surface) - 1

    @inbounds for idx in 1:nλ
        nr = FT(optics.nᵣ[idx])
        ns_idx = n_surface === nothing ? nr : FT(n_surface[ns_off + idx])
        k_unit = FT(K_total[idx]) / N_ft
        surf_opt = _eval_prepared_surface(surface, prep_surf, μ_ft, nr, ns_idx, false)
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
    ::MultiPigmentOpticalProperties;
    kwargs...,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

function prospect(
    ::LeafProspectProProperties,
    ::MultiPigmentOpticalProperties;
    kwargs...,
)
    throw(ArgumentError(
        "LeafProspectProProperties (single Ccab pool) cannot be paired with MultiPigmentOpticalProperties (separate Kchla, Kchlb); use LeafMultiPigmentProperties or LeafOpticalProperties",
    ))
end

function prospect_budget(
    leaf::LeafMultiPigmentProperties{FTL},
    optics::AbstractProspectProperties;
    n_surface = nothing,
) where {FTL<:Real}
    FT = float(promote_type(FTL, eltype(optics.nᵣ)))
    b = prospect_budget(leaf, optics, ConeSurface{FT}(FT(40), one(FT)); n_surface = n_surface)
    return (;
        transmission = b.transmission,
        reflectance = b.reflectance,
        absorptance = b.absorptance,
    )
end

function prospect(
    leaf::LeafMultiPigmentProperties,
    optics::AbstractProspectProperties;
    n_surface = nothing,
)
    b = prospect_budget(leaf, optics; n_surface = n_surface)
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

When constructed with keywords, supplying `Cprot` or `Ccbc` without `Cm` defaults
`Cm = 0`, and supplying `Cchla` or `Cchlb` without `Ccab` defaults `Ccab = 0`,
avoiding double-counting dry mass or chlorophyll.
"""
struct LeafSIPProperties{FT<:Real} <: AbstractLeafProperties
    "Chlorophyll a+b content `[µg cm⁻²]` (0 when `Cchla/Cchlb` are used)"
    Ccab::FT
    "Optional separate Chlorophyll a content `[µg cm⁻²]` (requires `Ccab == 0`)"
    Cchla::FT
    "Optional separate Chlorophyll b content `[µg cm⁻²]` (requires `Ccab == 0`)"
    Cchlb::FT
    "Carotenoid content `[µg cm⁻²]`"
    Ccar::FT
    "Anthocyanin content `[µg cm⁻²]`"
    Canth::FT
    "Brown pigments content"
    Cbrown::FT
    "Equivalent water thickness `[cm]`"
    Cw::FT
    "Bulk dry matter content `[g cm⁻²]` (0 when `Cprot + Ccbc` are used)"
    Cm::FT
    "Protein content `[g cm⁻²]` (requires `Cm == 0`)"
    Cprot::FT
    "Carbon-based constituents content `[g cm⁻²]` (requires `Cm == 0`)"
    Ccbc::FT
    "Scattering element scaling factor `[cm² g⁻¹]` (Wu et al. 2020: 600)"
    alpha_lma::FT
    "Recollision probability LMA slope `[cm² g⁻¹]` (Wu et al. 2020: 2765.0)"
    f_lma_slope::FT
    "Asymmetry parameter LMA slope `[cm² g⁻¹]` (Wu et al. 2020: 102.8)"
    g_lma_slope::FT

    function LeafSIPProperties{FT}(
        Ccab, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw,
        Cm, Cprot, Ccbc, alpha_lma, f_lma_slope, g_lma_slope,
    ) where {FT<:Real}
        return new{FT}(
            FT(Ccab), FT(Cchla), FT(Cchlb), FT(Ccar), FT(Canth), FT(Cbrown), FT(Cw),
            FT(Cm), FT(Cprot), FT(Ccbc), FT(alpha_lma), FT(f_lma_slope), FT(g_lma_slope),
        )
    end
end

function LeafSIPProperties(
    Ccab::Real, Cchla::Real, Cchlb::Real, Ccar::Real, Canth::Real, Cbrown::Real, Cw::Real,
    Cm::Real, Cprot::Real, Ccbc::Real, alpha_lma::Real, f_lma_slope::Real, g_lma_slope::Real,
)
    FT = _canopy_parameter_type(
        Ccab, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw,
        Cm, Cprot, Ccbc, alpha_lma, f_lma_slope, g_lma_slope,
    )
    return LeafSIPProperties{FT}(
        Ccab, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw,
        Cm, Cprot, Ccbc, alpha_lma, f_lma_slope, g_lma_slope,
    )
end

function LeafSIPProperties{FT}(;
    Ccab = _UNSET_LEAF_KW,
    Cchla = _UNSET_LEAF_KW,
    Cchlb = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
    alpha_lma = _UNSET_LEAF_KW,
    f_lma_slope = _UNSET_LEAF_KW,
    g_lma_slope = _UNSET_LEAF_KW,
) where {FT<:Real}
    BFT = _leaf_base_float(FT)
    if Ccab === _UNSET_LEAF_KW && (Cchla !== _UNSET_LEAF_KW || Cchlb !== _UNSET_LEAF_KW)
        cab_v = zero(FT)
        ca_v = Cchla === _UNSET_LEAF_KW ? zero(FT) : FT(Cchla)
        cb_v = Cchlb === _UNSET_LEAF_KW ? zero(FT) : FT(Cchlb)
    else
        cab_v = Ccab === _UNSET_LEAF_KW ? FT(BFT(40.0)) : FT(Ccab)
        ca_v = Cchla === _UNSET_LEAF_KW ? zero(FT) : FT(Cchla)
        cb_v = Cchlb === _UNSET_LEAF_KW ? zero(FT) : FT(Cchlb)
    end
    if prospect_primal(cab_v) > zero(BFT) &&
       (prospect_primal(ca_v) > zero(BFT) || prospect_primal(cb_v) > zero(BFT))
        throw(ArgumentError("LeafSIPProperties cannot set both Ccab > 0 and Cchla/Cchlb > 0 (double-counts chlorophyll); omit Ccab when using Cchla/Cchlb"))
    end

    if Cm === _UNSET_LEAF_KW && (Cprot !== _UNSET_LEAF_KW || Ccbc !== _UNSET_LEAF_KW)
        cm_v = zero(FT)
        cprot_v = Cprot === _UNSET_LEAF_KW ? zero(FT) : FT(Cprot)
        ccbc_v = Ccbc === _UNSET_LEAF_KW ? zero(FT) : FT(Ccbc)
    else
        cm_v = Cm === _UNSET_LEAF_KW ? FT(BFT(0.008)) : FT(Cm)
        cprot_v = Cprot === _UNSET_LEAF_KW ? zero(FT) : FT(Cprot)
        ccbc_v = Ccbc === _UNSET_LEAF_KW ? zero(FT) : FT(Ccbc)
    end
    if prospect_primal(cm_v) > zero(BFT) &&
       (prospect_primal(cprot_v) > zero(BFT) || prospect_primal(ccbc_v) > zero(BFT))
        throw(ArgumentError("LeafSIPProperties cannot set both Cm > 0 and Cprot/Ccbc > 0 (double-counts LMA); omit Cm when using Cprot + Ccbc"))
    end

    car_v = Ccar === _UNSET_LEAF_KW ? FT(BFT(10.0)) : FT(Ccar)
    ant_v = Canth === _UNSET_LEAF_KW ? FT(BFT(0.5)) : FT(Canth)
    brn_v = Cbrown === _UNSET_LEAF_KW ? zero(FT) : FT(Cbrown)
    cw_v = Cw === _UNSET_LEAF_KW ? FT(BFT(0.012)) : FT(Cw)
    alp_v = alpha_lma === _UNSET_LEAF_KW ? FT(BFT(600.0)) : FT(alpha_lma)
    flm_v = f_lma_slope === _UNSET_LEAF_KW ? FT(BFT(2765.0)) : FT(f_lma_slope)
    glm_v = g_lma_slope === _UNSET_LEAF_KW ? FT(BFT(102.8)) : FT(g_lma_slope)

    return LeafSIPProperties{FT}(
        cab_v, ca_v, cb_v, car_v, ant_v, brn_v, cw_v,
        cm_v, cprot_v, ccbc_v, alp_v, flm_v, glm_v,
    )
end

function LeafSIPProperties(;
    Ccab = _UNSET_LEAF_KW,
    Cchla = _UNSET_LEAF_KW,
    Cchlb = _UNSET_LEAF_KW,
    Ccar = _UNSET_LEAF_KW,
    Canth = _UNSET_LEAF_KW,
    Cbrown = _UNSET_LEAF_KW,
    Cw = _UNSET_LEAF_KW,
    Cm = _UNSET_LEAF_KW,
    Cprot = _UNSET_LEAF_KW,
    Ccbc = _UNSET_LEAF_KW,
    alpha_lma = _UNSET_LEAF_KW,
    f_lma_slope = _UNSET_LEAF_KW,
    g_lma_slope = _UNSET_LEAF_KW,
)
    T_raw = promote_type(
        _kw_type(Ccab), _kw_type(Cchla), _kw_type(Cchlb), _kw_type(Ccar),
        _kw_type(Canth), _kw_type(Cbrown), _kw_type(Cw), _kw_type(Cm),
        _kw_type(Cprot), _kw_type(Ccbc), _kw_type(alpha_lma),
        _kw_type(f_lma_slope), _kw_type(g_lma_slope),
    )
    FT = T_raw === Union{} ? Float64 : float(T_raw)
    return LeafSIPProperties{FT}(;
        Ccab, Cchla, Cchlb, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc,
        alpha_lma, f_lma_slope, g_lma_slope,
    )
end

@doc raw"""
    leaf_sip_budget(leaf::LeafSIPProperties, optics; surface = nothing, μ_in = 1, n_surface = nothing)

Evaluate the Wu et al. (2020) `Leaf-SIP` spectral-invariants leaf optical budget
(`surface = nothing`, matching published `SIP_Model.m` which omits cuticular
Fresnel reflection) or a surface-coupled `Leaf-SIP+Fresnel` budget when an
[`AbstractLeafSurface`](@ref) is supplied (`surface = ConeSurface(...)` or
`RoughDielectricSurface(...)`, with optional `n_surface` wax dispersion).

Wu et al. (2020)'s default slopes (`alpha_lma = 600`, `f_lma_slope = 2765`,
`g_lma_slope = 102.8`) were empirically calibrated for bare `surface = nothing`;
adding an external `surface` restores physical blue/red cuticular reflectance
floors, and `alpha_lma`, `f_lma_slope`, and `g_lma_slope` can be recalibrated when
fitting NIR/SWIR plateaus jointly with surface coupling.

Computes absorptance directly via ``-\operatorname{expm1}(-k_{\mathrm{elem}})``
without ``1 - R - T`` subtraction, preserving ``R + T + A = 1`` into the
nonabsorbing limit.
"""
function leaf_sip_budget(
    leaf::LeafSIPProperties{FTL},
    optics::AbstractProspectProperties;
    surface::Union{Nothing,AbstractLeafSurface} = nothing,
    μ_in::Real = one(_leaf_base_float(promote_type(FTL, eltype(optics.nᵣ)))),
    n_surface::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {FTL<:Real}
    _validate_leaf_contents((
        leaf.Ccab, leaf.Cchla, leaf.Cchlb, leaf.Ccar, leaf.Canth, leaf.Cbrown,
        leaf.Cw, leaf.Cm, leaf.Cprot, leaf.Ccbc,
    ))
    _validate_optics_arrays(optics)
    nλ = length(optics.nᵣ)
    if n_surface !== nothing
        length(n_surface) == nλ ||
            throw(DimensionMismatch("n_surface length $(length(n_surface)) must match optics length $(nλ)"))
        all(x -> (v = prospect_primal(x); isfinite(v) && v > one(v)), n_surface) ||
            throw(ArgumentError("n_surface refractive indices must be finite and > 1"))
    end

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
    FTW = n_surface === nothing ? eltype(optics.nᵣ) : eltype(n_surface)
    FT = float(promote_type(FTL, eltype(optics.nᵣ), FTS, typeof(μ_in), FTW))
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

    T_vec = Vector{FT}(undef, nλ)
    R_vec = Vector{FT}(undef, nλ)
    A_vec = Vector{FT}(undef, nλ)

    elem_scale = FT(lma_eff) * FT(leaf.alpha_lma)
    prep_surf_top = surface === nothing ? nothing : _prepare_surface_cheb(FT, surface, μ_ft, false)
    ns_off = n_surface === nothing ? 0 : firstindex(n_surface) - 1

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
            ns_idx = n_surface === nothing ? nr : FT(n_surface[ns_off + idx])
            s_top = _eval_prepared_surface(surface, prep_surf_top, μ_ft, nr, ns_idx, false)
            # Note: s_top.r_int_diffuse and s_top.t_int_diffuse already use t_hem_ext / nr^2!
            b = _couple_external_interfaces_to_medium(
                s_top.r_ext, s_top.t_ext, s_top.r_int_diffuse, s_top.t_int_diffuse,
                r_int, t_int, a_int, r_int, t_int, a_int,
                s_top.r_int_diffuse, s_top.t_int_diffuse,
            )
            T_vec[idx] = b.transmission
            R_vec[idx] = b.reflectance
            A_vec[idx] = b.absorptance
        end
    end

    return (; transmission = T_vec, reflectance = R_vec, absorptance = A_vec)
end
