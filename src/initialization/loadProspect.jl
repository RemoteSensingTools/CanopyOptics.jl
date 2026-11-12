"""
    createLeafOpticalStruct(λ_bnds)
    createLeafOpticalStruct()

Load the bundled PROSPECT-PRO leaf optical database (`400–2500 nm` at `1 nm`
resolution, Féret et al., 2021) and bin-average refractive index and specific
absorption cross-sections over spectral bin edges `λ_bnds`, returning a
[`LeafOpticalProperties`](@ref) table.

`λ_bnds` may carry any `Unitful` spectral unit (`u"nm"`, `u"μm"`, `u"cm^-1"`)
or be unitless (interpreted as `nm` with a warning). Calling
`createLeafOpticalStruct()` uses the native `(399.5:2500.5)u"nm"` `1 nm` grid.

# Examples
```jldoctest
julia> using Unitful

julia> opti = createLeafOpticalStruct((400.0:100.0:2400.0)u"nm");

julia> length(opti.nᵣ)
20
```
"""
function createLeafOpticalStruct(λ_bnds)
    N = length(λ_bnds) - 1
    N >= 1 || throw(ArgumentError("λ_bnds must contain at least 2 boundary points"))
    FT = float(eltype(ustrip(first(λ_bnds))))
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

function createLeafOpticalStruct()
    λ_bnds = (399.5:2500.5)u"nm"
    return createLeafOpticalStruct(λ_bnds)
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
    λ_ref = unit2nm(λ_bnds)
    N = length(λ_bnds) - 1
    N >= 1 || throw(ArgumentError("λ_bnds must contain at least 2 boundary points"))
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

"""
    createMultiPigmentOpticalStruct(λ_bnds = (399.5:2500.5)u"nm"; basis::Symbol = :calibrated_in_vivo)

Load a [`MultiPigmentOpticalProperties`](@ref) table with separate chlorophyll *a*
(`Kchla`) and chlorophyll *b* (`Kchlb`) specific absorption cross-sections `[cm² µg⁻¹]`
averaged onto `λ_bnds`, using consistent mass-based pigment units (`[cm² µg⁻¹]` for
`Kchla`, `Kchlb`, `Kcar`, and `Kant`) matching Féret et al. (2017, 2021) `OPTI_2021`.

Supported `basis` options:
- `:calibrated_in_vivo` (default): Empirical in-vivo `Kchla, Kchlb` calibrated on
  `ANGERS` (`400–800 nm`, inner 5-fold species-grouped cross-validation at `8°`
  integrating-sphere incidence; `data/multi_pigment_chlab_in_vivo.csv`), spliced to
  `OPTI_2021` for `λ > 800 nm`.
- `:prospect_mp_2017`: Published Zhang et al. (2017, *Sci. Rep.* 7:6429) pseudo-Voigt
  `Kchla, Kchlb, Kcar` on `400–800 nm` (`data/prospect_mp_voigt.toml`), spliced to
  `OPTI_2021` for `λ > 800 nm` so Lorentzian wings do not leak absorption into NIR/SWIR.
- `:prospect_mp_plus_2022`: Published Zhang et al. (2022, *Sensors* 22:3025) pseudo-Voigt
  `Kchla, Kchlb, Kcar, Kant` on `400–800 nm` (`data/prospect_mp_voigt.toml`), spliced to
  `OPTI_2021` for `λ > 800 nm`.
- `:prospect_pro_pooled`: Baseline control with `Kchla == Kchlb == Kcab` from `OPTI_2021`.
"""
function createMultiPigmentOpticalStruct(
    λ_bnds = (399.5:2500.5)u"nm";
    basis::Symbol = :calibrated_in_vivo,
    FT::Type{<:AbstractFloat} = float(eltype(ustrip(first(λ_bnds)))),
)
    basis in (:calibrated_in_vivo, :prospect_mp_2017, :prospect_mp_plus_2022, :prospect_pro_pooled) ||
        throw(ArgumentError(
            "Unsupported multi-pigment basis :$(basis); expected :calibrated_in_vivo, :prospect_mp_2017, :prospect_mp_plus_2022, or :prospect_pro_pooled",
        ))
    if basis === :prospect_pro_pooled
        ou_p = unit(first(λ_bnds))
        bnds_ft = ou_p == NoUnits ? FT.(λ_bnds) : FT.(ustrip.(λ_bnds)) .* ou_p
        return MultiPigmentOpticalProperties(createLeafOpticalStruct(bnds_ft))
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

"""
    load_wax_refractive_index(λ_bnds = (399.5:2500.5)u"nm")

Load the Stuckens et al. (2009, *Remote Sens. Environ.* 113(12), 2560–2573;
DOI `10.1016/j.rse.2009.07.014`) epicuticular
leaf wax refractive-index dispersion spectrum `n_wax(λ)` (`data/epicuticular_wax_stuckens2009.csv`)
averaged onto `λ_bnds`.

Unlike aqueous internal leaf refractive index `nᵣ(λ)` (which drops from `~1.42` in
VNIR to `~1.27` near `2450 nm` due to liquid-water dispersion), epicuticular wax
exhibits weak hydrocarbon dispersion (`n_wax ≈ 1.48` at `400 nm` to `1.40` at `2500 nm`),
sustaining cuticular surface specular reflectance in the SWIR when passed via
`n_surface` to [`prospect_budget`](@ref), [`dorsiventral_budget`](@ref), or
[`prepare_leaf_bsdf`](@ref).
"""
function load_wax_refractive_index(
    λ_bnds = (399.5:2500.5)u"nm";
    FT::Type{<:AbstractFloat} = float(eltype(ustrip(first(λ_bnds)))),
)
    csv_path = joinpath(_PACKAGE_DATA_DIR, "epicuticular_wax_stuckens2009.csv")
    tab, _ = readdlm(csv_path, ',', FT; comments = true, comment_char = '#', header = true)
    size(tab, 1) == 2101 || throw(DimensionMismatch("Expected 2101 rows (400:2500 nm) in $csv_path"))
    λ_1nm = tab[:, 1]
    nwax_1nm = tab[:, 2]
    (nwax_out,) = _bin_1nm_columns(λ_1nm, (nwax_1nm,), λ_bnds, FT)
    return nwax_out
end
