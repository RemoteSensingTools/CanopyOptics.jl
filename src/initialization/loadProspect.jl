@inline function _default_spectral_float(λ_bnds)
    ET = eltype(λ_bnds)
    return ET === Any ? Float64 : float(typeof(ustrip(zero(ET))))
end

function LeafOpticalProperties{FT}(
    λ_bnds::AbstractVector = (399.5:2500.5)u"nm",
) where {FT<:AbstractFloat}
    N = length(λ_bnds) - 1
    N >= 1 || throw(ArgumentError("λ_bnds must contain at least 2 boundary points"))
    λ_ref = unit2nm(λ_bnds)

    KS = readdlm(OPTI_2021, '\t', FT)
    λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc = ntuple(_ -> zeros(FT, N), 10)
    vars = (λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)
    λ = KS[:, 1] .* u"nm"
    i0_ref = firstindex(λ_ref)
    i0_bnds = firstindex(λ_bnds)
    @inbounds for i in 1:N
        r1 = λ_ref[i0_ref + i - 1]
        r2 = λ_ref[i0_ref + i]
        start = min(r1, r2)
        stop  = max(r1, r2)
        ind_all = findall((λ .≥ start) .& (λ .< stop))
        if isempty(ind_all)
            @warn "Some λ_bnds bins are empty, better coarsen the input grid $(λ_bnds[i0_bnds + i - 1]) -> $(λ_bnds[i0_bnds + i])"
        end
        for j in eachindex(vars)
            vars[j][i] = mean(view(KS, ind_all, j))
        end
    end
    ou = unit(first(λ_bnds))
    out_unit = ou == NoUnits ? u"nm" : ou
    return LeafOpticalProperties(
        uconvSpectral(λ_out .* u"nm", out_unit),
        nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc,
    )
end

function LeafOpticalProperties(
    λ_bnds::AbstractVector = (399.5:2500.5)u"nm";
    FT::Type{<:AbstractFloat} = _default_spectral_float(λ_bnds),
)
    return LeafOpticalProperties{FT}(λ_bnds)
end

"""
    createLeafOpticalStruct(λ_bnds = (399.5:2500.5)u"nm")

Deprecated constructor for [`LeafOpticalProperties`](@ref); call
`LeafOpticalProperties(λ_bnds)` instead.
"""
function createLeafOpticalStruct(λ_bnds::AbstractVector = (399.5:2500.5)u"nm")
    Base.depwarn(
        "`createLeafOpticalStruct` is deprecated, use `LeafOpticalProperties` instead.",
        :createLeafOpticalStruct,
    )
    return LeafOpticalProperties(λ_bnds)
end

"""
    unit2nm(λ_bnds)

Convert a spectral grid `λ_bnds` with `Unitful` spectral units (or unitless
values assumed to be in `nm`) to nanometers (`u"nm"`).
"""
function unit2nm(λ_bnds)
    try
        return uconvSpectral(λ_bnds, u"nm")
    catch
        @warn "Assuming unit of [nm] here for optical struct"
        return λ_bnds .* u"nm"
    end
end

"""
    uconvSpectral(x, out_unit)

Convert spectral array `x` to `out_unit` using `UnitfulEquivalences.Spectral()`,
returning `x` unchanged when `unit(first(x)) == out_unit`.
"""
function uconvSpectral(in_arr, out_unit)
    return unit(first(in_arr)) == out_unit ? in_arr : uconvert.(out_unit, in_arr, Spectral())
end

const _PACKAGE_DATA_DIR = normpath(joinpath(@__DIR__, "..", "..", "data"))

@inline function _eval_pseudo_voigt_spectrum!(
    out::AbstractVector{FT},
    λ_nm::AbstractVector{FT},
    idx_range::AbstractUnitRange{Int},
    spec_dict::AbstractDict,
) where {FT<:AbstractFloat}
    K_v = FT.(spec_dict["K_v"])
    K_h = FT.(spec_dict["K_h"])
    K_w = FT.(spec_dict["K_w"])
    K_p = FT.(spec_dict["K_p"])
    np = length(K_v)
    (length(K_h) == np && length(K_w) == np && length(K_p) == np) ||
        throw(DimensionMismatch("Pseudo-Voigt parameter vectors must have equal length"))
    ln2 = log(FT(2))
    half = inv(FT(2))
    @inbounds for i in idx_range
        w = λ_nm[i]
        acc = zero(FT)
        for j in 1:np
            x = (w - K_p[j]) / (half * K_w[j])
            x2 = x * x
            g_part = exp(-ln2 * x2)
            l_part = inv(one(FT) + x2)
            acc += K_h[j] * (K_v[j] * g_part + (one(FT) - K_v[j]) * l_part)
        end
        out[i] = acc
    end
    return out
end

function _bin_1nm_columns(
    λ_1nm::AbstractVector{FT},
    cols::Tuple,
    λ_bnds,
    ::Type{FT},
) where {FT<:AbstractFloat}
    N = length(λ_bnds) - 1
    N >= 1 || throw(ArgumentError("λ_bnds must contain at least 2 boundary points"))
    λ_ref = unit2nm(λ_bnds)
    λ_q = λ_1nm .* u"nm"
    outs = ntuple(_ -> zeros(FT, N), length(cols))
    i0_ref = firstindex(λ_ref)
    i0_bnds = firstindex(λ_bnds)
    @inbounds for i in 1:N
        r1 = λ_ref[i0_ref + i - 1]
        r2 = λ_ref[i0_ref + i]
        start = min(r1, r2)
        stop = max(r1, r2)
        ind_all = findall((λ_q .>= start) .& (λ_q .< stop))
        if isempty(ind_all)
            @warn "Some λ_bnds bins are empty, better coarsen the input grid $(λ_bnds[i0_bnds + i - 1]) -> $(λ_bnds[i0_bnds + i])"
        end
        for j in eachindex(cols)
            outs[j][i] = mean(view(cols[j], ind_all))
        end
    end
    return outs
end

function MultiPigmentOpticalProperties{FT}(
    λ_bnds::AbstractVector = (399.5:2500.5)u"nm";
    basis::Symbol = :calibrated_in_vivo,
) where {FT<:AbstractFloat}
    length(λ_bnds) - 1 >= 1 || throw(ArgumentError("λ_bnds must contain at least 2 boundary points"))
    basis in (:calibrated_in_vivo, :prospect_mp_2017, :prospect_mp_plus_2022, :prospect_pro_pooled) ||
        throw(ArgumentError(
            "Unsupported multi-pigment basis :$(basis); expected :calibrated_in_vivo, :prospect_mp_2017, :prospect_mp_plus_2022, or :prospect_pro_pooled",
        ))
    if basis === :prospect_pro_pooled
        return MultiPigmentOpticalProperties(LeafOpticalProperties{FT}(λ_bnds))
    end
    KS = readdlm(OPTI_2021, '\t', FT)
    λ_1nm = KS[:, 1]
    nᵣ_1nm = KS[:, 2]
    Kcab_1nm = KS[:, 3]
    Kcar_1nm = copy(KS[:, 4])
    Kant_1nm = copy(KS[:, 5])
    Kb_1nm = KS[:, 6]
    Kw_1nm = KS[:, 7]
    Km_1nm = KS[:, 8]
    Kp_1nm = KS[:, 9]
    Kcbc_1nm = KS[:, 10]

    Kchla_1nm = copy(Kcab_1nm)
    Kchlb_1nm = copy(Kcab_1nm)

    if basis === :calibrated_in_vivo
        csv_path = joinpath(_PACKAGE_DATA_DIR, "multi_pigment_chlab_in_vivo.csv")
        tab, _ = readdlm(csv_path, ',', FT; comments = true, comment_char = '#', header = true)
        size(tab, 1) == 401 || throw(DimensionMismatch("Expected 401 rows (400:800 nm) in $csv_path"))
        @inbounds for i in 1:401
            Kchla_1nm[i] = tab[i, 2]
            Kchlb_1nm[i] = tab[i, 3]
        end
    else
        voigt_cfg = TOML.parsefile(joinpath(_PACKAGE_DATA_DIR, "prospect_mp_voigt.toml"))
        cutoff_nm = FT(voigt_cfg["cutoff_nm"])
        idx_vnir = 1:findlast(w -> w <= cutoff_nm, λ_1nm)
        key = basis === :prospect_mp_2017 ? "prospect_mp_2017" : "prospect_mp_plus_2022"
        fam = voigt_cfg[key]
        _eval_pseudo_voigt_spectrum!(Kchla_1nm, λ_1nm, idx_vnir, fam["chla"])
        _eval_pseudo_voigt_spectrum!(Kchlb_1nm, λ_1nm, idx_vnir, fam["chlb"])
        _eval_pseudo_voigt_spectrum!(Kcar_1nm, λ_1nm, idx_vnir, fam["car"])
        if haskey(fam, "anth")
            _eval_pseudo_voigt_spectrum!(Kant_1nm, λ_1nm, idx_vnir, fam["anth"])
        end
    end

    λ_out, nᵣ, Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc = _bin_1nm_columns(
        λ_1nm,
        (λ_1nm, nᵣ_1nm, Kchla_1nm, Kchlb_1nm, Kcar_1nm, Kant_1nm, Kb_1nm, Kw_1nm, Km_1nm, Kp_1nm, Kcbc_1nm),
        λ_bnds,
        FT,
    )
    ou = unit(first(λ_bnds))
    out_unit = ou == NoUnits ? u"nm" : ou
    return MultiPigmentOpticalProperties(
        uconvSpectral(λ_out .* u"nm", out_unit),
        nᵣ, Kchla, Kchlb, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc,
    )
end

function MultiPigmentOpticalProperties(
    λ_bnds::AbstractVector = (399.5:2500.5)u"nm";
    basis::Symbol = :calibrated_in_vivo,
    FT::Type{<:AbstractFloat} = _default_spectral_float(λ_bnds),
)
    return MultiPigmentOpticalProperties{FT}(λ_bnds; basis = basis)
end

"""
    createMultiPigmentOpticalStruct(λ_bnds = (399.5:2500.5)u"nm"; basis = :calibrated_in_vivo, FT = ...)

Deprecated constructor for [`MultiPigmentOpticalProperties`](@ref); call
`MultiPigmentOpticalProperties(λ_bnds; basis = ...)` instead.
"""
function createMultiPigmentOpticalStruct(
    λ_bnds::AbstractVector = (399.5:2500.5)u"nm";
    basis::Symbol = :calibrated_in_vivo,
    FT::Type{<:AbstractFloat} = _default_spectral_float(λ_bnds),
)
    Base.depwarn(
        "`createMultiPigmentOpticalStruct` is deprecated, use `MultiPigmentOpticalProperties` instead.",
        :createMultiPigmentOpticalStruct,
    )
    return MultiPigmentOpticalProperties{FT}(λ_bnds; basis = basis)
end

@doc raw"""
    load_wax_refractive_index(λ_bnds = (399.5:2500.5)u"nm"; FT = ...)

Load the Stuckens et al. (2009, *Remote Sens. Environ.* 113(12), 2560–2573;
DOI `10.1016/j.rse.2009.07.014`) epicuticular leaf wax refractive-index dispersion
spectrum ``n_{\mathrm{wax}}(\lambda)`` (`data/epicuticular_wax_stuckens2009.csv`)
averaged onto `λ_bnds`.

Unlike aqueous internal leaf refractive index ``n_r(\lambda)`` (which drops from
``\sim 1.42`` in the VNIR to ``\sim 1.27`` near ``2450\;\mathrm{nm}`` due to
liquid-water dispersion), epicuticular wax exhibits weak hydrocarbon dispersion
(``n_{\mathrm{wax}} \approx 1.48`` at ``400\;\mathrm{nm}`` to ``1.40`` at
``2500\;\mathrm{nm}``), sustaining cuticular surface specular reflectance in the
SWIR when passed via `n_surface` to [`prospect_budget`](@ref),
[`dorsiventral_budget`](@ref), or [`prepare_leaf_bsdf`](@ref).
"""
function load_wax_refractive_index(
    λ_bnds = (399.5:2500.5)u"nm";
    FT::Type{<:AbstractFloat} = _default_spectral_float(λ_bnds),
)
    csv_path = joinpath(_PACKAGE_DATA_DIR, "epicuticular_wax_stuckens2009.csv")
    tab, _ = readdlm(csv_path, ',', FT; comments = true, comment_char = '#', header = true)
    size(tab, 1) == 2101 || throw(DimensionMismatch("Expected 2101 rows (400:2500 nm) in $csv_path"))
    λ_1nm = tab[:, 1]
    nwax_1nm = tab[:, 2]
    (nwax_out,) = _bin_1nm_columns(λ_1nm, (nwax_1nm,), λ_bnds, FT)
    return nwax_out
end
