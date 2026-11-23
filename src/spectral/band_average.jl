# Physical constants (SI 2019 exact values) for wavenumber Planck integration.
# With wavenumber ν in cm⁻¹, M_ν(ν, T) = π B_ν(ν, T) = c₁ ν³ / expm1(c₂ ν / T)
# gives hemispheric blackbody spectral exitance in W m⁻² (cm⁻¹)⁻¹, whose integral
# over [0, ∞) equals σ T⁴ in W m⁻².
const _PLANCK_H_JS = 6.62607015e-34
const _LIGHT_SPEED_MS = 299792458.0
const _BOLTZMANN_JK = 1.380649e-23
const _PLANCK_C1_FLUX_CM = π * 2.0e8 * _PLANCK_H_JS * _LIGHT_SPEED_MS^2
const _PLANCK_C2_CM_K = 100.0 * _PLANCK_H_JS * _LIGHT_SPEED_MS / _BOLTZMANN_JK
const _PLANCK_MIN_TEMP_K = 1.0
const _PLANCK_MAX_TEMP_K = 1.0e8

const _DEFAULT_SOLAR_CSV = joinpath(_PACKAGE_DATA_DIR, "reference_solar_global.csv")

# Positive nodes and weights of 8-point Gauss-Legendre quadrature on [-1, 1] (17 digits).
const _GL8_POS_NODES = (
    0.18343464249564980,
    0.52553240991632899,
    0.79666647741362674,
    0.96028985649753623,
)
const _GL8_POS_WEIGHTS = (
    0.36268378337836198,
    0.31370664587788729,
    0.22238103445337447,
    0.10122853629037626,
)

# Native PROSPECT-PRO 1-nm wavelength bin limits [nm] and corresponding wavenumber
# coverage limits [cm⁻¹] for dataSpec_PRO.csv (400:1:2500 nm).
const _PROSPECT_MIN_NM_EDGE = 399.5
const _PROSPECT_MAX_NM_EDGE = 2500.5
const _PROSPECT_MIN_NU_CM = 1.0e7 / _PROSPECT_MAX_NM_EDGE
const _PROSPECT_MAX_NU_CM = 1.0e7 / _PROSPECT_MIN_NM_EDGE

const _NATIVE_PROSPECT_LOCK = ReentrantLock()
const _NATIVE_PROSPECT_F64 =
    Ref{Union{Nothing,LeafOpticalProperties{Float64,Vector{typeof(1.0u"nm")}}}}(nothing)
const _NATIVE_PROSPECT_F32 =
    Ref{Union{Nothing,LeafOpticalProperties{Float32,Vector{typeof(1.0f0u"nm")}}}}(nothing)

function _native_prospect_optics(::Type{Float64})
    return @lock _NATIVE_PROSPECT_LOCK begin
        opt = _NATIVE_PROSPECT_F64[]
        if opt === nothing
            opt = LeafOpticalProperties{Float64}()
            _NATIVE_PROSPECT_F64[] = opt
        end
        opt
    end
end

function _native_prospect_optics(::Type{Float32})
    return @lock _NATIVE_PROSPECT_LOCK begin
        opt = _NATIVE_PROSPECT_F32[]
        if opt === nothing
            opt = LeafOpticalProperties{Float32}()
            _NATIVE_PROSPECT_F32[] = opt
        end
        opt
    end
end

_native_prospect_optics(::Type{FT}) where {FT<:AbstractFloat} = LeafOpticalProperties{FT}()

@inline _canon_zero(x::FT) where {FT} = iszero(x) ? zero(FT) : x

function _copy_canonical_vec(::Type{FT}, src::AbstractVector) where {FT}
    n = length(src)
    dst = Vector{FT}(undef, n)
    i0 = firstindex(src)
    @inbounds for k in 1:n
        dst[k] = _canon_zero(FT(src[i0 + k - 1]))
    end
    return dst
end

@noinline _throw_edges_length() =
    throw(ArgumentError("edges must contain at least 2 boundary points"))
@noinline _throw_edges_invalid() =
    throw(ArgumentError("edges must be finite, nonnegative, and strictly ascending in cm^-1"))
@noinline _throw_values_dim(nv::Int, ne::Int) =
    throw(DimensionMismatch("values length ($nv) must equal length(edges) - 1 ($(ne - 1))"))
@noinline _throw_values_nonfinite() =
    throw(ArgumentError("spectral bin values must be finite"))
@noinline _throw_bands_shape(sz) =
    throw(DimensionMismatch("bands must have size (2, nbnd), got $sz"))
@noinline _throw_dest_dim(nd::Int, nb::Int) =
    throw(DimensionMismatch("destination length ($nd) must match number of bands ($nb)"))
@noinline _throw_invalid_band(b::Int, ν_lo, ν_hi) =
    throw(ArgumentError("band $b ([$ν_lo, $ν_hi] cm^-1) must have finite nonnegative limits with ν_lo < ν_hi"))
@noinline _throw_uncovered_band(b::Int, ν_lo, ν_hi, e_min, e_max) =
    throw(ArgumentError("band $b ([$ν_lo, $ν_hi] cm^-1) is not fully covered by spectral grid [$e_min, $e_max] cm^-1"))
@noinline _throw_uncovered_weight_band(ν_lo, ν_hi, w_min, w_max) =
    throw(ArgumentError("interval [$ν_lo, $ν_hi] cm^-1 is not fully covered by TabulatedWeight grid [$w_min, $w_max] cm^-1"))
@noinline _throw_uncovered_grid(ν_lo, ν_hi, m_lo, m_hi) =
    throw(ArgumentError("grid [$ν_lo, $ν_hi] cm^-1 is not fully covered by leaf model domain [$m_lo, $m_hi] cm^-1"))
@noinline _throw_uncovered_bin(i::Int, e_lo, e_hi) =
    throw(ArgumentError("CompositeLeaf does not cover bin $i ([$e_lo, $e_hi] cm^-1)"))
@noinline _throw_aliased_dest() =
    throw(ArgumentError("dest must not alias input arrays or weight storage"))
@noinline _throw_outside_parent_band(b::Int, ν_lo, ν_hi, par_lo, par_hi) =
    throw(ArgumentError("band $b ([$ν_lo, $ν_hi] cm^-1) must lie within parent_band [$par_lo, $par_hi] cm^-1"))

@inline _has_same_elements_parent(::Base.ReshapedArray) = true
@inline _has_same_elements_parent(::Base.PermutedDimsArrays.PermutedDimsArray) = true
@inline _has_same_elements_parent(::Base.ReinterpretArray) = true
@inline _has_same_elements_parent(::SubArray) = false
@inline function _has_same_elements_parent(a::AbstractArray)
    p = parent(a)
    return p !== a && length(p) == length(a) && Base.dataids(p) == Base.dataids(a)
end

@inline _unwrap_alias_bijection(a::AbstractArray) =
    _has_same_elements_parent(a) ? _unwrap_alias_bijection(parent(a)) : a

@inline function _strided_byte_span(s::StridedArray)
    elsz = Base.elsize(typeof(s))
    p0 = UInt(pointer(s))
    sz = size(s)
    st = strides(s)
    min_off = 0
    max_off = 0
    @inbounds for d in eachindex(sz, st)
        span_d = (sz[d] - 1) * st[d]
        if span_d >= 0
            max_off += span_d
        else
            min_off += span_d
        end
    end
    lo = p0 - UInt(-min_off) * UInt(elsz)
    hi = p0 + UInt(max_off) * UInt(elsz) + UInt(elsz - 1)
    return lo, hi, elsz
end

@inline function _strided_lattice_disjoint(
    a::StridedVector,
    lo_a::UInt,
    b::StridedArray,
    lo_b::UInt,
    elsz::Int,
)
    sa = abs(stride(a, 1))
    sa <= 1 && return false
    ue = UInt(elsz)
    base = min(lo_a, lo_b)
    delta = max(lo_a, lo_b) - base
    (delta % ue) == 0 || return false
    if ndims(b) == 1
        sb = abs(stride(b, 1))
        g = gcd(sa, sb)
        g <= 1 && return false
        ug = UInt(g) * ue
        qa = Int(div((lo_a - base) % ug, ue))
        qb = Int(div((lo_b - base) % ug, ue))
        return qa != qb
    elseif ndims(b) == 2 && stride(b, 1) == 1
        d1 = size(b, 1)
        sb2 = abs(stride(b, 2))
        g = gcd(sa, sb2)
        g <= d1 && return false
        ug = UInt(g) * ue
        qa = Int(div((lo_a - base) % ug, ue))
        qb = Int(div((lo_b - base) % ug, ue))
        return mod(qa - qb, g) >= d1
    end
    return false
end

@inline function _arrays_might_alias(dest::AbstractArray, src::AbstractArray)
    a = _unwrap_alias_bijection(dest)
    b = _unwrap_alias_bijection(src)
    Base.mightalias(a, b) || return false
    (isempty(a) || isempty(b)) && return false
    if a isa StridedArray && b isa StridedArray
        ea = Base.elsize(typeof(a))
        eb = Base.elsize(typeof(b))
        if ea > 0 && eb > 0
            lo_a, hi_a, _ = _strided_byte_span(a)
            lo_b, hi_b, _ = _strided_byte_span(b)
            (hi_a < lo_b || hi_b < lo_a) && return false
            if ea == eb && a isa StridedVector && _strided_lattice_disjoint(a, lo_a, b, lo_b, ea)
                return false
            end
        end
    end
    if b isa SubArray && isbits(b.indices) && !_arrays_might_alias(a, parent(b))
        return false
    end
    if a isa SubArray && isbits(a.indices) && !_arrays_might_alias(parent(a), b)
        return false
    end
    return true
end

@inline function _validate_wavenumber_edges(edges::AbstractVector)
    n_edges = length(edges)
    n_edges >= 2 || _throw_edges_length()
    i0 = firstindex(edges)
    @inbounds begin
        prev = prospect_primal(edges[i0])
        (isfinite(prev) && prev >= zero(prev)) || _throw_edges_invalid()
        for k in 1:(n_edges - 1)
            cur = prospect_primal(edges[i0 + k])
            (isfinite(cur) && cur > prev) || _throw_edges_invalid()
            prev = cur
        end
    end
    return n_edges - 1
end

@inline function _validate_bin_values(values::AbstractVector, n_bins::Int)
    length(values) == n_bins || _throw_values_dim(length(values), n_bins + 1)
    @inbounds for i in eachindex(values)
        isfinite(prospect_primal(values[i])) || _throw_values_nonfinite()
    end
    return nothing
end

@inline function _validate_bands_shape(bands::AbstractMatrix)
    size(bands, 1) == 2 || _throw_bands_shape(size(bands))
    return size(bands, 2)
end

# Return 0-based bin offset k in 0:(length(edges)-2) of the first bin [edges[i0+k], edges[i0+k+1]]
# whose primal upper edge strictly exceeds primal(ν_lo).
@inline function _first_overlapping_bin(edges::AbstractVector, ν_lo::Real)
    i0 = firstindex(edges)
    p_lo = prospect_primal(ν_lo)
    lo = 0
    hi = length(edges) - 2
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if prospect_primal(edges[i0 + mid + 1]) <= p_lo
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo
end

# Primal-tie-preserving interval intersection [max(ν_lo, e_lo), min(ν_hi, e_hi)] so
# ForwardDiff.Dual tangents on band limits ν_lo / ν_hi survive at exact grid knots.
@inline function _clamp_overlap_lo(ν_lo, e_lo, ::Type{AccFT}) where {AccFT}
    return prospect_primal(ν_lo) >= prospect_primal(e_lo) ? AccFT(ν_lo) : AccFT(e_lo)
end

@inline function _clamp_overlap_hi(ν_hi, e_hi, ::Type{AccFT}) where {AccFT}
    return prospect_primal(ν_hi) <= prospect_primal(e_hi) ? AccFT(ν_hi) : AccFT(e_hi)
end

"""
    AbstractSpectralWeight

Supertype for wavenumber spectral weighting functions ``w(\\nu)`` (per `cm⁻¹`) used
by [`band_average`](@ref), [`band_average!`](@ref), [`subband_weights`](@ref),
[`subband_weights!`](@ref), [`CanopyOptics.band_weight`](@ref), and
[`CanopyOptics.band_weight!`](@ref).

Concrete subtypes provided by `CanopyOptics` are [`FlatWeight`](@ref),
[`PlanckWeight`](@ref), and [`TabulatedWeight`](@ref) (including
[`ReferenceSolar`](@ref)). Custom subtypes may implement
[`CanopyOptics.weight_interval_integral`](@ref).
"""
abstract type AbstractSpectralWeight end

"""
    FlatWeight()

Uniform spectral weight ``w(\\nu) = 1`` per `cm⁻¹`, weighting fine bins purely by
their wavenumber overlap ``\\Delta\\nu`` with each band.

# Examples
```jldoctest
julia> FlatWeight() isa AbstractSpectralWeight
true
```
"""
struct FlatWeight <: AbstractSpectralWeight end

@doc raw"""
    PlanckWeight(T)

Blackbody spectral exitance weight ``M_\nu(\nu, T) = \pi B_\nu(\nu, T)`` in
``\mathrm{W\,m^{-2}\,(cm^{-1})^{-1}}`` at temperature `T` in `K`:

```math
M_\nu(\nu, T) = \frac{2\pi \times 10^8 h c^2\,\nu^3}{\exp(100\,h c\,\nu / (k_B T)) - 1},
```

integrated in coordinates relative to each overlap's lower wavenumber edge via
composite 8-point Gauss–Legendre quadrature. Temperature `T` must be finite and satisfy
`1 <= T <= 10^8` K.

# Examples
```jldoctest
julia> PlanckWeight(300.0)
PlanckWeight{Float64}(300.0)
```
"""
struct PlanckWeight{FT<:Real} <: AbstractSpectralWeight
    T::FT

    function PlanckWeight{FT}(T::Real) where {FT<:Real}
        t_val = _canon_zero(FT(T))
        tp = prospect_primal(t_val)
        (isfinite(tp) && _PLANCK_MIN_TEMP_K <= tp <= _PLANCK_MAX_TEMP_K) || throw(ArgumentError(
            "PlanckWeight temperature T must be finite and in [$_PLANCK_MIN_TEMP_K, $_PLANCK_MAX_TEMP_K] K, got $T",
        ))
        return new{FT}(t_val)
    end
end

PlanckWeight(T::Real) = PlanckWeight{float(typeof(T))}(float(T))

Base.:(==)(a::PlanckWeight, b::PlanckWeight) = a.T == b.T
Base.hash(a::PlanckWeight, h::UInt) = hash(_canon_zero(a.T), hash(:PlanckWeight, h))

"""
    TabulatedWeight(edges, w)

Piecewise-constant spectral weight density `w` (per `cm⁻¹`) defined on ascending
wavenumber bin edges `edges` in `cm⁻¹` (`length(edges) == length(w) + 1`).

After conversion to element type `FT`, `edges` must contain at least two finite,
nonnegative, strictly ascending points, and `w` must contain finite nonnegative bin
densities. When [`band_average`](@ref) or [`subband_weights`](@ref) evaluates a parent
band not fully covered by `[first(edges), last(edges)]` (or a band over which the
integrated weight vanishes), evaluation falls back explicitly to [`FlatWeight`](@ref)
over that band as a finite placeholder so `0/0` never occurs.

# Examples
```jldoctest
julia> tw = TabulatedWeight([10000.0, 15000.0, 20000.0], [1.2, 0.8]);

julia> length(tw.w)
2
```
"""
struct TabulatedWeight{FT<:Real} <: AbstractSpectralWeight
    edges::Vector{FT}
    w::Vector{FT}

    function TabulatedWeight{FT}(edges::AbstractVector, w::AbstractVector) where {FT<:Real}
        ne = length(edges)
        ne >= 2 || _throw_edges_length()
        n_bins = ne - 1
        length(w) == n_bins ||
            throw(DimensionMismatch("weight vector length ($(length(w))) must equal length(edges) - 1 ($n_bins)"))
        e_vec = _copy_canonical_vec(FT, edges)
        w_vec = _copy_canonical_vec(FT, w)
        _validate_wavenumber_edges(e_vec)
        @inbounds for i in eachindex(w_vec)
            wp = prospect_primal(w_vec[i])
            (isfinite(wp) && wp >= zero(wp)) ||
                throw(ArgumentError("TabulatedWeight densities must be finite and nonnegative"))
        end
        return new{FT}(e_vec, w_vec)
    end
end

function TabulatedWeight(edges::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    FT = float(promote_type(eltype(edges), eltype(w)))
    return TabulatedWeight{FT}(edges, w)
end

Base.:(==)(a::TabulatedWeight, b::TabulatedWeight) = a.edges == b.edges && a.w == b.w
function Base.hash(a::TabulatedWeight, h::UInt)
    return hash(a.w, hash(a.edges, hash(:TabulatedWeight, h)))
end

@doc raw"""
    ReferenceSolar(grid, irradiance; grid_unit = :nm) -> TabulatedWeight
    ReferenceSolar(path = ...) -> TabulatedWeight{Float64}

Construct a [`TabulatedWeight`](@ref) (per `cm⁻¹`) from a clear-sky global
(direct + diffuse) surface solar irradiance spectrum.

# Grid and `|dλ/dν|` Jacobian Conventions
- `grid_unit = :nm` (default): `irradiance` is spectral irradiance per nanometer
  (``\mathrm{W\,m^{-2}\,nm^{-1}}``).
  - Point samples (`length(grid) == length(irradiance) >= 2`): `grid` holds strictly
    ascending wavelengths ``\lambda_k`` in `nm`. Interior wavelength bin edges sit at
    midpoints ``\frac{1}{2}(\lambda_k + \lambda_{k+1})`` and outer edges extend by half
    the adjacent spacing; each bin density is converted to per `cm⁻¹` via the exact
    bin-averaged Jacobian ``|d\lambda / d\nu|``:
    ``\bar{E}_{\nu,i} = E_{\lambda,k}\,\lambda_{\mathrm{lo}}\lambda_{\mathrm{hi}} / 10^7``.
    Whole-bin integrated flux ``\bar{E}_{\nu,i}\,\Delta\nu_i = E_{\lambda,k}\,\Delta\lambda_k``
    is conserved; a band boundary cutting through the interior of a bin integrates the
    piecewise-constant per-`cm⁻¹` density ``\bar{E}_{\nu,i}``.
  - Bin edges (`length(grid) == length(irradiance) + 1`): `grid` holds strictly
    ascending wavelength bin edges in `nm`, converted with the same per-bin Jacobian.
- `grid_unit = :cm_inv`: `irradiance` is already per `cm⁻¹`, on strictly ascending
  wavenumber point samples (`length(grid) == length(irradiance) >= 2`) or bin edges
  (`length(grid) == length(irradiance) + 1`).

# CSV Table Format and Zero-Irradiance UV Convention
Calling `ReferenceSolar()` with no arguments loads `data/reference_solar_global.csv` if
that file exists, or throws an `ArgumentError` when no default table is installed; pass
`ReferenceSolar(path)` to load a CSV file from an explicit path. Expected CSV format:
- UTF-8 comma-separated values (`,`), optional `#` comment lines, optional single
  column-header row with non-numeric text in both columns (e.g.
  `wavelength_nm,irradiance_W_m2_nm`), and at least two numeric rows of
  `(wavelength_nm, global_irradiance_W_m2_nm)` with strictly ascending positive
  wavelengths in `nm` and finite nonnegative global irradiance in ``\mathrm{W\,m^{-2}\,nm^{-1}}``.
- Include an explicit zero-irradiance row below `280 nm` (e.g. `195.0,0.0` paired with
  `280.0,0.0`) so UV bands below the stratospheric ozone cutoff are covered by the
  table: partially illuminated bands (such as `263–345 nm`) then weight strictly by
  their illuminated `> 280 nm` portion, whereas fully dark bands (such as `200–263 nm`)
  have zero physical solar flux (`CanopyOptics.band_weight == 0.0`, and share `0.0` in
  [`subband_weights`](@ref) inside an illuminated parent band).
- On a band whose total integrated solar weight is zero (`0/0` solar-weighted mean
  undefined), [`band_average`](@ref) explicitly returns the [`FlatWeight`](@ref)
  wavenumber average as a finite placeholder optical property.

# Examples
```jldoctest
julia> tw = ReferenceSolar(400.0:100.0:700.0, [1.1, 1.5, 1.4, 1.2]);

julia> tw isa TabulatedWeight{Float64}
true
```
"""
function ReferenceSolar(path::AbstractString = _DEFAULT_SOLAR_CSV)
    if !isfile(path)
        if path == _DEFAULT_SOLAR_CSV
            throw(ArgumentError(
                "ReferenceSolar() found no default solar table at \"$path\"; " *
                "place a two-column CSV (wavelength_nm, irradiance_W_m2_nm) there, " *
                "pass ReferenceSolar(csv_path), or supply arrays via " *
                "ReferenceSolar(grid, irradiance; grid_unit = :nm)",
            ))
        end
        throw(ArgumentError("ReferenceSolar CSV file not found at \"$path\""))
    end
    return _load_reference_solar_csv(path)
end

function _load_reference_solar_csv(path::AbstractString)
    λ_nm = Float64[]
    irr_nm = Float64[]
    saw_first_row = false
    line_no = 0
    for raw_line in eachline(path)
        line_no += 1
        line_str = line_no == 1 ? chopprefix(raw_line, "\ufeff") : raw_line
        s = strip(line_str)
        (isempty(s) || startswith(s, '#')) && continue
        idx_hash = findfirst('#', s)
        if idx_hash !== nothing
            s = strip(SubString(s, 1, prevind(s, idx_hash)))
            isempty(s) && continue
        end
        parts = split(s, ',')
        length(parts) == 2 || throw(ArgumentError(
            "expected 2 comma-separated columns (wavelength_nm, irradiance_W_m2_nm) on line $line_no ($path:$line_no)",
        ))
        col1 = strip(parts[1])
        col2 = strip(parts[2])
        (isempty(col1) || isempty(col2)) && throw(ArgumentError(
            "empty column entry on line $line_no ($path:$line_no)",
        ))
        v1 = tryparse(Float64, col1)
        v2 = tryparse(Float64, col2)
        if v1 === nothing && v2 === nothing
            if !saw_first_row
                saw_first_row = true
                continue
            end
            throw(ArgumentError("unexpected non-numeric row on line $line_no ($path:$line_no)"))
        elseif v1 === nothing || v2 === nothing
            throw(ArgumentError("malformed numeric entry on line $line_no ($path:$line_no)"))
        end
        (isfinite(v1) && v1 > 0.0 && isfinite(v2) && v2 >= 0.0) || throw(ArgumentError(
            "non-finite or out-of-range solar CSV value on line $line_no ($path:$line_no)",
        ))
        if !isempty(λ_nm) && v1 <= last(λ_nm)
            throw(ArgumentError(
                "wavelength must be strictly ascending on line $line_no ($path:$line_no)",
            ))
        end
        saw_first_row = true
        push!(λ_nm, v1)
        push!(irr_nm, v2)
    end
    length(λ_nm) >= 2 || throw(ArgumentError(
        "ReferenceSolar CSV at \"$path\" must contain at least 2 numeric rows",
    ))
    return ReferenceSolar(λ_nm, irr_nm; grid_unit = :nm)
end

function ReferenceSolar(
    grid::AbstractVector{<:Real},
    irradiance::AbstractVector{<:Real};
    grid_unit::Symbol = :nm,
)
    grid_unit in (:nm, :cm_inv) ||
        throw(ArgumentError("grid_unit must be :nm or :cm_inv, got :$grid_unit"))
    ng = length(grid)
    ni = length(irradiance)
    ni >= 1 || throw(ArgumentError("irradiance must be non-empty"))
    (ng == ni || ng == ni + 1) || throw(DimensionMismatch(
        "grid length ($ng) must equal length(irradiance) ($ni) for point samples " *
        "or length(irradiance) + 1 ($(ni + 1)) for bin edges",
    ))
    FT = float(promote_type(eltype(grid), eltype(irradiance)))
    BFT = _leaf_base_float(FT)

    if grid_unit === :nm
        if ng == ni
            ni >= 2 || throw(ArgumentError("point-sampled wavelength grid requires at least 2 points"))
            λ_edges = _point_grid_to_edges(grid, FT, true)
        else
            _validate_positive_ascending(grid)
            λ_edges = copyto!(Vector{FT}(undef, ng), grid)
        end
        edges_nu = Vector{FT}(undef, ni + 1)
        w_nu = Vector{FT}(undef, ni)
        i0_i = firstindex(irradiance)
        inv_scale = inv(BFT(1.0e7))
        @inbounds for m in 1:(ni + 1)
            edges_nu[m] = BFT(1.0e7) / λ_edges[ni + 2 - m]
        end
        @inbounds for i in 1:ni
            k = ni + 1 - i
            e_val = FT(irradiance[i0_i + k - 1])
            ep = prospect_primal(e_val)
            (isfinite(ep) && ep >= zero(ep)) ||
                throw(ArgumentError("solar irradiance values must be finite and nonnegative"))
            λ_lo = λ_edges[k]
            λ_hi = λ_edges[k + 1]
            w_nu[i] = e_val * (λ_lo * λ_hi * inv_scale)
        end
        return TabulatedWeight{FT}(edges_nu, w_nu)
    else
        if ng == ni
            ni >= 2 || throw(ArgumentError("point-sampled wavenumber grid requires at least 2 points"))
            edges_nu = _point_grid_to_edges(grid, FT, false)
            return TabulatedWeight{FT}(edges_nu, irradiance)
        else
            return TabulatedWeight{FT}(grid, irradiance)
        end
    end
end

function _validate_positive_ascending(grid::AbstractVector)
    n = length(grid)
    n >= 2 || throw(ArgumentError("grid must contain at least 2 points"))
    i0 = firstindex(grid)
    @inbounds begin
        prev = prospect_primal(grid[i0])
        (isfinite(prev) && prev > zero(prev)) ||
            throw(ArgumentError("wavelength grid values must be finite and strictly positive"))
        for k in 1:(n - 1)
            cur = prospect_primal(grid[i0 + k])
            (isfinite(cur) && cur > prev) ||
                throw(ArgumentError("wavelength grid values must be finite and strictly ascending"))
            prev = cur
        end
    end
    return nothing
end

function _point_grid_to_edges(grid::AbstractVector, ::Type{FT}, strictly_positive::Bool) where {FT<:Real}
    n = length(grid)
    i0 = firstindex(grid)
    BFT = _leaf_base_float(FT)
    half = inv(BFT(2))
    edges = Vector{FT}(undef, n + 1)
    @inbounds begin
        g1 = FT(grid[i0])
        g2 = FT(grid[i0 + 1])
        p1 = prospect_primal(g1)
        p2 = prospect_primal(g2)
        (isfinite(p1) && isfinite(p2) && p2 > p1) ||
            throw(ArgumentError("point grid must be finite and strictly ascending"))
        e_first = g1 - half * (g2 - g1)
        ep1 = prospect_primal(e_first)
        if strictly_positive
            ep1 > zero(ep1) || throw(ArgumentError("implied lower wavelength bin edge ($ep1 nm) must be positive"))
        else
            ep1 >= zero(ep1) || throw(ArgumentError("implied lower wavenumber bin edge ($ep1 cm^-1) must be nonnegative"))
        end
        edges[1] = e_first
        edges[2] = half * (g1 + g2)
        prev_g = g2
        prev_p = p2
        for k in 3:n
            gk = FT(grid[i0 + k - 1])
            pk = prospect_primal(gk)
            (isfinite(pk) && pk > prev_p) ||
                throw(ArgumentError("point grid must be finite and strictly ascending"))
            edges[k] = half * (prev_g + gk)
            prev_g = gk
            prev_p = pk
        end
        gn_1 = FT(grid[i0 + n - 2])
        edges[n + 1] = prev_g + half * (prev_g - gn_1)
    end
    return edges
end

@doc raw"""
    wavelength_to_wavenumber(wavelengths_nm) -> edges
    wavelength_to_wavenumber(wavelengths_nm, values; density = false) -> (edges, bin_values)

Convert point samples on a regular ascending wavelength grid `wavelengths_nm` in `nm`
(such as PROSPECT's `400:1:2500`) to ascending wavenumber bin edges `edges` in `cm⁻¹`
(`length(edges) == length(wavelengths_nm) + 1`) and reordered wavenumber bin means
`bin_values`.

Each sample ``\lambda_k`` on a regular grid with step ``\Delta\lambda`` spans the
wavelength bin ``[\lambda_{\mathrm{lo}}, \lambda_{\mathrm{hi}}] = [\lambda_k - \Delta\lambda/2,\, \lambda_k + \Delta\lambda/2]``,
mapped to ``[\nu_{\mathrm{lo}}, \nu_{\mathrm{hi}}] = [10^7/\lambda_{\mathrm{hi}},\, 10^7/\lambda_{\mathrm{lo}}]``
in reverse index order (`i = N + 1 - k`).

# Jacobian Convention (`density`)
- `density = false` (default): `values` is an intensive dimensionless or per-state
  quantity (such as reflectance ``\rho(\lambda)`` or transmittance ``\tau(\lambda)``),
  so `bin_values[i] = values[N + 1 - i]`. Apply physical spectral weights on the native
  fine grid before coarsening, because resampling intensive bin means in wavenumber
  preserves the flat-in-wavenumber (`FlatWeight`) integral.
- `density = true`: `values` is a spectral density per nanometer (such as spectral
  irradiance ``E_\lambda`` in ``\mathrm{W\,m^{-2}\,nm^{-1}}``). Converting
  ``E_\lambda\,|d\lambda| = E_\nu\,d\nu`` over ``[\lambda_{\mathrm{lo}}, \lambda_{\mathrm{hi}}]``
  applies the exact bin-averaged Jacobian ``|d\lambda / d\nu|`` (`10^7 / \nu^2`):

```math
\bar{E}_{\nu,i} = E_{\lambda,k}\,\frac{\lambda_{\mathrm{lo}}\,\lambda_{\mathrm{hi}}}{10^7},
```

  conserving whole-bin flux ``\bar{E}_{\nu,i}(\nu_{\mathrm{hi}} - \nu_{\mathrm{lo}}) = E_{\lambda,k}\Delta\lambda``
  to floating-point precision. A band boundary that cuts through the interior of a bin
  integrates the constant per-`cm⁻¹` density ``\bar{E}_{\nu,i}`` (for example, a single
  uniform `400–800 nm` bin at `1 W m⁻² nm⁻¹` yields `133.33 W m⁻²` over `600–800 nm`,
  versus `200 W m⁻²` in the piecewise-constant per-`nm` representation).

# Examples
```jldoctest
julia> edges, rho_nu = wavelength_to_wavenumber(400.0:500.0:1400.0, [0.05, 0.45, 0.30]);

julia> length(edges) == 4 && rho_nu == [0.30, 0.45, 0.05]
true
```

See also [`wavenumber_to_wavelength`](@ref), [`ReferenceSolar`](@ref).
"""
function wavelength_to_wavenumber(wavelengths_nm::AbstractVector{<:Real})
    FT = float(eltype(wavelengths_nm))
    return _regular_nm_to_nu_edges(wavelengths_nm, FT)[1]
end

function wavelength_to_wavenumber(
    wavelengths_nm::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    density::Bool = false,
)
    n = length(wavelengths_nm)
    length(values) == n ||
        throw(DimensionMismatch("values length ($(length(values))) must match wavelengths_nm length ($n)"))
    FT = float(promote_type(eltype(wavelengths_nm), eltype(values)))
    BFT = _leaf_base_float(FT)
    edges, λ1, dλ = _regular_nm_to_nu_edges(wavelengths_nm, FT)
    bin_values = Vector{FT}(undef, n)
    i0_v = firstindex(values)
    half = inv(BFT(2))
    inv_scale = inv(BFT(1.0e7))
    @inbounds for i in 1:n
        k = n + 1 - i
        v = FT(values[i0_v + k - 1])
        isfinite(prospect_primal(v)) || _throw_values_nonfinite()
        if density
            λ_lo = λ1 + (BFT(k - 1) - half) * dλ
            λ_hi = λ1 + (BFT(k) - half) * dλ
            bin_values[i] = v * (λ_lo * λ_hi * inv_scale)
        else
            bin_values[i] = v
        end
    end
    return edges, bin_values
end

function _regular_nm_to_nu_edges(wavelengths_nm::AbstractVector, ::Type{FT}) where {FT<:Real}
    n = length(wavelengths_nm)
    n >= 2 || throw(ArgumentError("wavelengths_nm must contain at least 2 points to define regular bin spacing"))
    BFT = _leaf_base_float(FT)
    i0 = firstindex(wavelengths_nm)
    λ_first = FT(wavelengths_nm[i0])
    λ_last = FT(wavelengths_nm[i0 + n - 1])
    p_first = prospect_primal(λ_first)
    p_last = prospect_primal(λ_last)
    (isfinite(p_first) && isfinite(p_last) && p_first > zero(p_first) && p_last > p_first) ||
        throw(ArgumentError("wavelengths_nm must be finite, positive, and strictly ascending"))

    dλ = (λ_last - λ_first) / BFT(n - 1)
    dp = prospect_primal(dλ)
    half = inv(BFT(2))
    prospect_primal(λ_first - half * dλ) > zero(dp) ||
        throw(ArgumentError("lower bin edge λ_1 - Δλ/2 must be strictly positive"))

    @inbounds for k in 1:(n - 1)
        a = prospect_primal(FT(wavelengths_nm[i0 + k - 1]))
        b = prospect_primal(FT(wavelengths_nm[i0 + k]))
        (isfinite(b) && b > a) ||
            throw(ArgumentError("wavelengths_nm must be finite, positive, and strictly ascending"))
        tol = BFT(16) * eps(BFT) * max(abs(a), abs(b)) + BFT(16) * eps(dp)
        abs((b - a) - dp) <= tol ||
            throw(ArgumentError("wavelengths_nm must be a regular uniformly spaced grid in nm"))
    end

    edges = Vector{FT}(undef, n + 1)
    @inbounds for i in 1:(n + 1)
        m = n + 1 - i
        λ_edge = λ_first + (BFT(m) - half) * dλ
        edges[i] = BFT(1.0e7) / λ_edge
    end
    _validate_wavenumber_edges(edges)
    return edges, λ_first, dλ
end

"""
    wavenumber_to_wavelength(edges) -> wavelengths_nm
    wavenumber_to_wavelength(edges, bin_values; density = false) -> (wavelengths_nm, values)

Convert ascending positive wavenumber bin edges `edges` in `cm⁻¹` (and optional bin
means `bin_values`) back to ascending wavelength bin midpoints `wavelengths_nm` in
`nm` (and reordered wavelength values `values`), inverting
[`wavelength_to_wavenumber`](@ref).

When `density = true`, per-`cm⁻¹` spectral densities are converted back to per-`nm`
densities via ``E_{\\lambda,k} = \\bar{E}_{\\nu,i}\\,10^7 / (\\lambda_{\\mathrm{lo}}\\,\\lambda_{\\mathrm{hi}})``.

# Examples
```jldoctest
julia> edges, e_nu = wavelength_to_wavenumber(400.0:200.0:1000.0, [1.0, 1.5, 1.2, 0.9]; density = true);

julia> λ_rt, e_rt = wavenumber_to_wavelength(edges, e_nu; density = true);

julia> λ_rt ≈ [400.0, 600.0, 800.0, 1000.0] && e_rt ≈ [1.0, 1.5, 1.2, 0.9]
true
```
"""
function wavenumber_to_wavelength(edges::AbstractVector{<:Real})
    n_bins = _validate_wavenumber_edges(edges)
    i0 = firstindex(edges)
    prospect_primal(edges[i0]) > zero(prospect_primal(edges[i0])) ||
        throw(ArgumentError("edges must be strictly positive in cm^-1 to convert to wavelength"))
    FT = float(eltype(edges))
    BFT = _leaf_base_float(FT)
    half = inv(BFT(2))
    wavelengths_nm = Vector{FT}(undef, n_bins)
    @inbounds for k in 1:n_bins
        λ_lo = BFT(1.0e7) / FT(edges[i0 + n_bins + 1 - k])
        λ_hi = BFT(1.0e7) / FT(edges[i0 + n_bins - k])
        wavelengths_nm[k] = half * (λ_lo + λ_hi)
    end
    return wavelengths_nm
end

function wavenumber_to_wavelength(
    edges::AbstractVector{<:Real},
    bin_values::AbstractVector{<:Real};
    density::Bool = false,
)
    n_bins = _validate_wavenumber_edges(edges)
    _validate_bin_values(bin_values, n_bins)
    i0_e = firstindex(edges)
    prospect_primal(edges[i0_e]) > zero(prospect_primal(edges[i0_e])) ||
        throw(ArgumentError("edges must be strictly positive in cm^-1 to convert to wavelength"))
    FT = float(promote_type(eltype(edges), eltype(bin_values)))
    BFT = _leaf_base_float(FT)
    half = inv(BFT(2))
    wavelengths_nm = Vector{FT}(undef, n_bins)
    values = Vector{FT}(undef, n_bins)
    i0_v = firstindex(bin_values)
    @inbounds for k in 1:n_bins
        i = n_bins + 1 - k
        λ_lo = BFT(1.0e7) / FT(edges[i0_e + i])
        λ_hi = BFT(1.0e7) / FT(edges[i0_e + i - 1])
        wavelengths_nm[k] = half * (λ_lo + λ_hi)
        v = FT(bin_values[i0_v + i - 1])
        values[k] = density ? v * (BFT(1.0e7) / (λ_lo * λ_hi)) : v
    end
    return wavelengths_nm, values
end

_weight_eltype(::FlatWeight) = Union{}
_weight_eltype(::PlanckWeight{FT}) where {FT} = FT
_weight_eltype(::TabulatedWeight{FT}) where {FT} = FT
_weight_eltype(::AbstractSpectralWeight) = Float64

@inline _weight_might_alias(::AbstractArray, ::FlatWeight) = false
@inline _weight_might_alias(::AbstractArray, ::PlanckWeight) = false
@inline _weight_might_alias(dest::AbstractArray, tw::TabulatedWeight) =
    _arrays_might_alias(dest, tw.edges) || _arrays_might_alias(dest, tw.w)

@inline _any_might_alias(dest::AbstractArray, x::AbstractArray) = _arrays_might_alias(dest, x)
@inline _any_might_alias(::AbstractArray, _) = false

@generated function _weight_might_alias(dest::AbstractArray, w::W) where {W<:AbstractSpectralWeight}
    exprs = Expr[]
    for i in 1:fieldcount(W)
        if !(fieldtype(W, i) <: Number)
            push!(exprs, :(_any_might_alias(dest, getfield(w, $i))))
        end
    end
    isempty(exprs) && return :(false)
    return foldr((e1, e2) -> :($e1 || $e2), exprs)
end

# Evaluate M_ν(a_w + rel, T) cleanly without cancellation near ν = 0, overflow at large
# x = c₂ ν / T, or catastrophic subtraction in shifted Wien-tail exponents.
@inline function _planck_exitance_rel(a_w, rel, inv_T, ν_shift, use_shift::Bool)
    c1 = _PLANCK_C1_FLUX_CM
    c2 = _PLANCK_C2_CM_K
    ν = a_w + rel
    x = (c2 * ν) * inv_T
    if !use_shift
        if iszero(prospect_primal(x))
            # Exact Rayleigh-Jeans limit: c₁ ν³ / (c₂ ν / T) = (c₁ / c₂) T ν².
            return (c1 / c2) * ((ν * ν) / inv_T)
        elseif prospect_primal(x) <= 36.0
            return c1 * (ν / expm1(x)) * (ν * ν)
        else
            neg_ex = exp(-x)
            return c1 * (ν * ν * ν) * neg_ex / (one(neg_ex) - neg_ex)
        end
    else
        # Multiplies by exp(c₂ ν_shift / T), which is identical across all bins of a band
        # and cancels in normalized ratios while keeping x_rel = c₂ (ν - ν_shift) / T O(1).
        dν_shift = (a_w - ν_shift) + rel
        x_rel = (c2 * dν_shift) * inv_T
        neg_ex = exp(-x)
        return c1 * (ν * ν * ν) * exp(-x_rel) / (one(neg_ex) - neg_ex)
    end
end

# Integrate M_ν(ν, T) (optionally scaled by exp(c₂ ν_shift / T) when use_shift = true)
# over [a, b] in coordinates relative to the lower edge `a` using composite 8-point
# Gauss-Legendre quadrature. Truncates negligible far-Wien tails (Δx > 80 beyond
# max(x_a, 36), relative contribution < 1e-31 << eps(Float64)) directly in relative span
# so `span` never cancels to zero and panel count stays bounded (<= 78).
@inline function _integrate_planck_interval(a, b, T, ν_shift = 0.0, use_shift::Bool = false)
    WFT = promote_type(typeof(a), typeof(b), typeof(T), Float64)
    a_w = WFT(a)
    b_w = WFT(b)
    raw_span = b_w - a_w
    prospect_primal(raw_span) <= 0 && return zero(WFT)
    inv_T = inv(WFT(T))
    c2 = _PLANCK_C2_CM_K
    p_inv_T = prospect_primal(inv_T)
    p_a = prospect_primal(a_w)
    x_a = c2 * p_a * p_inv_T
    x_from_ref = use_shift ? c2 * (p_a - ν_shift) * p_inv_T : x_a
    x_from_ref > 750.0 && return zero(WFT)

    dx_cap = max(0.0, 36.0 - x_a) + 80.0
    max_span = (dx_cap / c2) / inv_T
    span = prospect_primal(raw_span) > prospect_primal(max_span) ? max_span : raw_span
    dx_dimless = c2 * prospect_primal(span) * p_inv_T
    n_panels = max(1, ceil(Int, dx_dimless / 1.5))
    panel_w = span / n_panels
    half_p = 0.5 * panel_w
    total = zero(WFT)
    @inbounds for p in 0:(n_panels - 1)
        mid_rel = (p + 0.5) * panel_w
        acc = zero(WFT)
        for q in 1:4
            off = half_p * _GL8_POS_NODES[q]
            acc += _GL8_POS_WEIGHTS[q] * (
                _planck_exitance_rel(a_w, mid_rel - off, inv_T, ν_shift, use_shift) +
                _planck_exitance_rel(a_w, mid_rel + off, inv_T, ν_shift, use_shift)
            )
        end
        total += half_p * acc
    end
    return total
end

"""
    CanopyOptics.weight_interval_integral(weight::AbstractSpectralWeight, a, b)

Integrate spectral weight density `weight` (per `cm⁻¹`) over wavenumber interval
`[a, b]` (`0 <= a <= b` in `cm⁻¹`), preserving exact additivity across contiguous
sub-intervals.

For [`TabulatedWeight`](@ref), `[a, b]` must be covered by `[first(weight.edges), last(weight.edges)]`
(throwing `ArgumentError` if uncovered) and returns `0` when the tabulated density
is zero across `[a, b]`. For normalized sub-band shares that share [`band_average`](@ref)'s
[`FlatWeight`](@ref) fallback and Wien-tail scaling, see [`subband_weights`](@ref).

This qualified extension hook is not exported; custom [`AbstractSpectralWeight`](@ref)
subtypes may specialize `CanopyOptics.weight_interval_integral` to plug into
[`band_average!`](@ref), [`CanopyOptics.band_weight!`](@ref), and
[`subband_weights!`](@ref).
"""
@inline weight_interval_integral(::FlatWeight, a::Real, b::Real) = b - a

@inline function weight_interval_integral(w::PlanckWeight, a::Real, b::Real)
    return _integrate_planck_interval(a, b, w.T, 0.0, false)
end

function weight_interval_integral(tw::TabulatedWeight, a::Real, b::Real)
    AccFT = promote_type(typeof(a), typeof(b), eltype(tw.edges), eltype(tw.w), Float64)
    p_a = prospect_primal(a)
    p_b = prospect_primal(b)
    t_edges = tw.edges
    t_min = prospect_primal(first(t_edges))
    t_max = prospect_primal(last(t_edges))
    (p_a >= t_min && p_b <= t_max) || _throw_uncovered_weight_band(p_a, p_b, t_min, t_max)
    p_a == p_b && return zero(AccFT)
    return _tabulated_covered_integral(tw, a, b, AccFT)
end

@inline function _tabulated_covered_integral(tw::TabulatedWeight, a, b, ::Type{AccFT}) where {AccFT}
    t_edges = tw.edges
    t_w = tw.w
    nt = length(t_w)
    j = _first_overlapping_bin(t_edges, a)
    p_b = prospect_primal(b)
    acc = zero(AccFT)
    @inbounds while j < nt
        lo = _clamp_overlap_lo(a, t_edges[j + 1], AccFT)
        hi = _clamp_overlap_hi(b, t_edges[j + 2], AccFT)
        dν = hi - lo
        if prospect_primal(dν) > 0
            acc += AccFT(t_w[j + 1]) * dν
        end
        prospect_primal(t_edges[j + 2]) >= p_b && break
        j += 1
    end
    return acc
end

@inline function _band_average_flat(
    values::AbstractVector,
    edges::AbstractVector,
    n_bins::Int,
    ν_lo,
    ν_hi,
    ::Type{AccFT},
) where {AccFT}
    i0_v = firstindex(values)
    i0_e = firstindex(edges)
    p_hi = prospect_primal(ν_hi)
    k = _first_overlapping_bin(edges, ν_lo)
    num = zero(AccFT)
    den = zero(AccFT)
    @inbounds while k < n_bins
        a = _clamp_overlap_lo(ν_lo, edges[i0_e + k], AccFT)
        b = _clamp_overlap_hi(ν_hi, edges[i0_e + k + 1], AccFT)
        dν = b - a
        if prospect_primal(dν) > 0
            num += AccFT(values[i0_v + k]) * dν
            den += dν
        end
        prospect_primal(edges[i0_e + k + 1]) >= p_hi && break
        k += 1
    end
    return num / den
end

@inline function _band_average_planck(
    values::AbstractVector,
    edges::AbstractVector,
    n_bins::Int,
    ν_lo,
    ν_hi,
    pw::PlanckWeight,
    ::Type{AccFT},
) where {AccFT}
    i0_v = firstindex(values)
    i0_e = firstindex(edges)
    p_lo = prospect_primal(ν_lo)
    p_hi = prospect_primal(ν_hi)
    # Shift exponent relative to ν_lo when x_lo > 36 so Wien-tail bands never underflow to 0/0.
    x_lo = _PLANCK_C2_CM_K * p_lo / prospect_primal(pw.T)
    use_shift = x_lo > 36.0
    ν_shift = use_shift ? p_lo : zero(p_lo)
    k = _first_overlapping_bin(edges, ν_lo)
    num = zero(AccFT)
    den = zero(AccFT)
    @inbounds while k < n_bins
        a = _clamp_overlap_lo(ν_lo, edges[i0_e + k], AccFT)
        b = _clamp_overlap_hi(ν_hi, edges[i0_e + k + 1], AccFT)
        if prospect_primal(b) > prospect_primal(a)
            dw = _integrate_planck_interval(a, b, pw.T, ν_shift, use_shift)
            num += AccFT(values[i0_v + k]) * dw
            den += dw
        end
        prospect_primal(edges[i0_e + k + 1]) >= p_hi && break
        k += 1
    end
    return prospect_primal(den) > 0 ? (num / den) : _band_average_flat(values, edges, n_bins, ν_lo, ν_hi, AccFT)
end

@inline function _band_average_tabulated(
    values::AbstractVector,
    edges::AbstractVector,
    n_bins::Int,
    ν_lo,
    ν_hi,
    tw::TabulatedWeight,
    ::Type{AccFT},
) where {AccFT}
    t_edges = tw.edges
    t_w = tw.w
    p_lo = prospect_primal(ν_lo)
    p_hi = prospect_primal(ν_hi)
    # Explicit fallback to FlatWeight when the weight table does not cover the band.
    if p_lo < prospect_primal(first(t_edges)) || p_hi > prospect_primal(last(t_edges))
        return _band_average_flat(values, edges, n_bins, ν_lo, ν_hi, AccFT)
    end
    i0_v = firstindex(values)
    i0_e = firstindex(edges)
    m_bins = length(t_w)
    i = _first_overlapping_bin(edges, ν_lo)
    j = _first_overlapping_bin(t_edges, ν_lo)
    cur = AccFT(ν_lo)
    hi_lim = AccFT(ν_hi)
    num = zero(AccFT)
    den = zero(AccFT)
    @inbounds while prospect_primal(cur) < p_hi
        e_next = AccFT(edges[i0_e + i + 1])
        t_next = AccFT(t_edges[j + 2])
        pe = prospect_primal(e_next)
        pt = prospect_primal(t_next)
        # Prefer hi_lim on primal ties at the upper band boundary so Dual tangents on ν_hi survive.
        if p_hi <= pe && p_hi <= pt
            nxt = hi_lim
        elseif pe <= pt
            nxt = e_next
        else
            nxt = t_next
        end
        dν = nxt - cur
        if prospect_primal(dν) > 0
            dw = AccFT(t_w[j + 1]) * dν
            num += AccFT(values[i0_v + i]) * dw
            den += dw
        end
        cur = nxt
        p_cur = prospect_primal(cur)
        if pe <= p_cur && i + 1 < n_bins
            i += 1
        end
        if pt <= p_cur && j + 1 < m_bins
            j += 1
        end
    end
    # Also fall back to FlatWeight if tabulated density vanishes everywhere across the band.
    return prospect_primal(den) > 0 ? (num / den) : _band_average_flat(values, edges, n_bins, ν_lo, ν_hi, AccFT)
end

@inline function _band_average_generic(
    values::AbstractVector,
    edges::AbstractVector,
    n_bins::Int,
    ν_lo,
    ν_hi,
    weight::AbstractSpectralWeight,
    ::Type{AccFT},
) where {AccFT}
    i0_v = firstindex(values)
    i0_e = firstindex(edges)
    p_hi = prospect_primal(ν_hi)
    k = _first_overlapping_bin(edges, ν_lo)
    num = zero(AccFT)
    den = zero(AccFT)
    @inbounds while k < n_bins
        a = _clamp_overlap_lo(ν_lo, edges[i0_e + k], AccFT)
        b = _clamp_overlap_hi(ν_hi, edges[i0_e + k + 1], AccFT)
        if prospect_primal(b) > prospect_primal(a)
            dw = AccFT(weight_interval_integral(weight, a, b))
            num += AccFT(values[i0_v + k]) * dw
            den += dw
        end
        prospect_primal(edges[i0_e + k + 1]) >= p_hi && break
        k += 1
    end
    return prospect_primal(den) > 0 ? (num / den) : _band_average_flat(values, edges, n_bins, ν_lo, ν_hi, AccFT)
end

# Resolve once for `parent_band = [par_lo, par_hi]` whether `weight` falls back to
# `FlatWeight` (`use_flat == true`), plus the parent Wien-tail reference wavenumber for
# `PlanckWeight`. Shared by `band_average!(...; parent_band)` and `subband_weights!` so
# mixed-coverage sub-band means and shares always use the same effective weight.
@inline function _resolve_parent_weight(::FlatWeight, par_lo, par_hi, ::Type{AccFT}) where {AccFT}
    p_lo = prospect_primal(par_lo)
    return (true, zero(p_lo), false, zero(AccFT))
end

@inline function _resolve_parent_weight(pw::PlanckWeight, par_lo, par_hi, ::Type{AccFT}) where {AccFT}
    p_par_lo = prospect_primal(par_lo)
    x_lo = _PLANCK_C2_CM_K * p_par_lo / prospect_primal(pw.T)
    use_shift = x_lo > 36.0
    ν_shift = use_shift ? p_par_lo : zero(p_par_lo)
    tot = AccFT(_integrate_planck_interval(par_lo, par_hi, pw.T, ν_shift, use_shift))
    return (prospect_primal(tot) <= 0, ν_shift, use_shift, tot)
end

@inline function _resolve_parent_weight(tw::TabulatedWeight, par_lo, par_hi, ::Type{AccFT}) where {AccFT}
    t_edges = tw.edges
    p_par_lo = prospect_primal(par_lo)
    p_par_hi = prospect_primal(par_hi)
    covered = p_par_lo >= prospect_primal(first(t_edges)) && p_par_hi <= prospect_primal(last(t_edges))
    if !covered
        return (true, zero(p_par_lo), false, zero(AccFT))
    end
    tot = _tabulated_covered_integral(tw, par_lo, par_hi, AccFT)
    return (prospect_primal(tot) <= 0, zero(p_par_lo), false, tot)
end

@inline function _resolve_parent_weight(weight::AbstractSpectralWeight, par_lo, par_hi, ::Type{AccFT}) where {AccFT}
    p_par_lo = prospect_primal(par_lo)
    tot = AccFT(weight_interval_integral(weight, par_lo, par_hi))
    return (prospect_primal(tot) <= 0, zero(p_par_lo), false, tot)
end

@inline _parent_band_eltype(::Nothing) = Union{}
@inline _parent_band_eltype(pb::AbstractArray{<:Real}) = eltype(pb)
@inline _parent_band_eltype(pb::Tuple{Real,Real}) = promote_type(typeof(pb[1]), typeof(pb[2]))
@inline _parent_band_eltype(_) = Union{}

@inline _validate_parent_band(::Nothing, p_e_min, p_e_max, weight, ::Type{AccFT}) where {AccFT} =
    (false, p_e_min, p_e_max, false)

@inline function _validate_parent_band(parent_band, p_e_min, p_e_max, weight, ::Type{AccFT}) where {AccFT}
    par_lo, par_hi = _extract_parent_band_limits(parent_band)
    p_par_lo = prospect_primal(par_lo)
    p_par_hi = prospect_primal(par_hi)
    (isfinite(p_par_lo) && isfinite(p_par_hi) && p_par_lo >= zero(p_par_lo) && p_par_hi > p_par_lo) ||
        _throw_invalid_band(1, p_par_lo, p_par_hi)
    (p_par_lo >= p_e_min && p_par_hi <= p_e_max) ||
        _throw_uncovered_band(1, p_par_lo, p_par_hi, p_e_min, p_e_max)
    parent_flat, _, _, _ = _resolve_parent_weight(weight, par_lo, par_hi, AccFT)
    return (true, p_par_lo, p_par_hi, parent_flat)
end

@doc raw"""
    band_average(values, edges, bands, weight = FlatWeight(); parent_band = nothing)

Compute the flux-weighted mean of fine-grid spectral bin means `values` across each
band in `bands` by **exact wavenumber overlap** of the fine bins with every band:

```math
\bar{v}_b = \frac{\int_{\nu_{\mathrm{lo},b}}^{\nu_{\mathrm{hi},b}} v(\nu)\,w(\nu)\,d\nu}{\int_{\nu_{\mathrm{lo},b}}^{\nu_{\mathrm{hi},b}} w(\nu)\,d\nu}.
```

- `values`: `N` piecewise-constant bin means on `[edges[i], edges[i + 1]]`.
- `edges`: `N + 1` finite, nonnegative, strictly ascending wavenumber bin edges in `cm⁻¹`.
- `bands`: `(2, nbnd)` matrix of band limits `[ν_lo; ν_hi]` in `cm⁻¹` with `0 <= ν_lo < ν_hi`.
- `weight`: [`AbstractSpectralWeight`](@ref) (default [`FlatWeight`](@ref)).
- `parent_band`: optional parent interval (`(2, 1)` matrix or 2-element vector/tuple
  `[ν_lo, ν_hi]`). When `bands` holds sub-bands of `parent_band`, passing `parent_band`
  resolves the effective weight (including [`TabulatedWeight`](@ref) coverage/zero-flux
  fallback to `FlatWeight`) once for `parent_band` and uses it across all sub-bands so
  `sum(subband_weights(sub_bands, parent_band, weight) .* band_average(values, edges, sub_bands, weight; parent_band = parent_band))`
  always recombines to `only(band_average(values, edges, parent_band, weight))`.

Throws an `ArgumentError` naming the band if any band is invalid or not fully
covered by `[first(edges), last(edges)]` (never returning silent zeros). If `weight`
is a [`TabulatedWeight`](@ref) whose grid does not cover a parent band (or integrates to
zero over that band), that band falls back explicitly to [`FlatWeight`](@ref) as a
placeholder optical mean.

See [`band_average!`](@ref) for the in-place variant (zero heap allocation when
arrays and built-in `weight` share a matching `Float32` or `Float64` element type),
[`subband_weights`](@ref) for normalized sub-band partition shares, and
[`CanopyOptics.band_weight`](@ref) for physical band-integrated weights.

# Examples
```jldoctest
julia> edges = [1000.0, 2000.0, 3000.0];

julia> bands = [1500.0 1000.0; 2500.0 3000.0];

julia> band_average([0.2, 0.6], edges, bands, FlatWeight())
2-element Vector{Float64}:
 0.4
 0.4
```
"""
function band_average(
    values::AbstractVector{<:Real},
    edges::AbstractVector{<:Real},
    bands::AbstractMatrix{<:Real},
    weight::AbstractSpectralWeight = FlatWeight();
    parent_band = nothing,
)
    PFT = _parent_band_eltype(parent_band)
    OutFT = float(promote_type(eltype(values), eltype(edges), eltype(bands), _weight_eltype(weight), PFT))
    dest = similar(bands, OutFT, axes(bands, 2))
    return band_average!(dest, values, edges, bands, weight; parent_band = parent_band)
end

"""
    band_average!(dest, values, edges, bands, weight = FlatWeight(); parent_band = nothing) -> dest

Overwrite `dest` in place with the exact-overlap band means of `values` over `bands`
weighted by `weight`, and return `dest` (zero heap allocation when `dest`, input
arrays, and built-in `weight` share a matching `Float32` or `Float64` element type;
`FlatWeight()` matches both).

`dest` must not alias `values`, `edges`, `bands`, `parent_band`, or `weight`'s
storage (throws `ArgumentError` if aliased; disjoint views of a shared buffer are
permitted). Throws `DimensionMismatch` when array dimensions disagree and
`ArgumentError` naming the band if any band in `bands` is invalid, lies outside
`parent_band`, or is not fully covered by `edges`. See [`band_average`](@ref).
"""
function band_average!(
    dest::AbstractVector,
    values::AbstractVector{<:Real},
    edges::AbstractVector{<:Real},
    bands::AbstractMatrix{<:Real},
    weight::AbstractSpectralWeight = FlatWeight();
    parent_band = nothing,
)
    (_arrays_might_alias(dest, values) ||
     _arrays_might_alias(dest, edges) ||
     _arrays_might_alias(dest, bands) ||
     _weight_might_alias(dest, weight) ||
     (parent_band isa AbstractArray && _arrays_might_alias(dest, parent_band))) && _throw_aliased_dest()

    n_bins = _validate_wavenumber_edges(edges)
    _validate_bin_values(values, n_bins)
    nbnd = _validate_bands_shape(bands)
    length(dest) == nbnd || _throw_dest_dim(length(dest), nbnd)

    e_min = first(edges)
    e_max = last(edges)
    p_e_min = prospect_primal(e_min)
    p_e_max = prospect_primal(e_max)
    r0 = firstindex(bands, 1)
    c0 = firstindex(bands, 2)
    d0 = firstindex(dest)
    DestFT = eltype(dest)
    PFT = _parent_band_eltype(parent_band)
    AccFT = promote_type(DestFT, eltype(values), eltype(edges), eltype(bands), _weight_eltype(weight), PFT, Float64)

    has_parent, p_par_lo, p_par_hi, parent_flat =
        _validate_parent_band(parent_band, p_e_min, p_e_max, weight, AccFT)

    @inbounds for b in 1:nbnd
        ν_lo = bands[r0, c0 + b - 1]
        ν_hi = bands[r0 + 1, c0 + b - 1]
        p_lo = prospect_primal(ν_lo)
        p_hi = prospect_primal(ν_hi)
        (isfinite(p_lo) && isfinite(p_hi) && p_lo >= zero(p_lo) && p_hi > p_lo) ||
            _throw_invalid_band(b, p_lo, p_hi)
        (p_lo >= p_e_min && p_hi <= p_e_max) ||
            _throw_uncovered_band(b, p_lo, p_hi, p_e_min, p_e_max)
        if has_parent && (p_lo < p_par_lo || p_hi > p_par_hi)
            _throw_outside_parent_band(b, p_lo, p_hi, p_par_lo, p_par_hi)
        end

        val = if parent_flat || weight isa FlatWeight
            _band_average_flat(values, edges, n_bins, ν_lo, ν_hi, AccFT)
        elseif weight isa PlanckWeight
            _band_average_planck(values, edges, n_bins, ν_lo, ν_hi, weight, AccFT)
        elseif weight isa TabulatedWeight
            _band_average_tabulated(values, edges, n_bins, ν_lo, ν_hi, weight, AccFT)
        else
            _band_average_generic(values, edges, n_bins, ν_lo, ν_hi, weight, AccFT)
        end
        dest[d0 + b - 1] = DestFT(val)
    end
    return dest
end

"""
    CanopyOptics.band_weight(bands, weight = FlatWeight())

Integrate `weight` (per `cm⁻¹`) across each band column `[ν_lo, ν_hi]` of the
`(2, nbnd)` matrix `bands`, returning the additive physical band integrals
``W_b = \\int_{\\nu_{\\mathrm{lo},b}}^{\\nu_{\\mathrm{hi},b}} w(\\nu)\\,d\\nu``
(for example ``\\mathrm{W\\,m^{-2}}`` for [`PlanckWeight`](@ref) or solar irradiance
[`TabulatedWeight`](@ref), and ``\\Delta\\nu`` in `cm⁻¹` for [`FlatWeight`](@ref)).

For `TabulatedWeight`, every band in `bands` must be covered by
`[first(weight.edges), last(weight.edges)]` (throwing `ArgumentError` otherwise), and
zero-density regions integrate to `0` so `W(B) = ∑ W(B_k)` holds identically across
contiguous sub-bands. For dimensionless sub-band flux shares `w_j` (`∑ w_j = 1`)
that include `band_average`'s [`FlatWeight`](@ref) fallback and Wien-tail scaling,
use [`subband_weights`](@ref).

This function is not exported; call `CanopyOptics.band_weight` (or in-place
[`CanopyOptics.band_weight!`](@ref), which runs with zero heap allocation when `dest`,
`bands`, and built-in `weight` share a matching `Float32` or `Float64` element type).

# Examples
```jldoctest
julia> CanopyOptics.band_weight([1000.0 2000.0; 2000.0 3500.0], FlatWeight())
2-element Vector{Float64}:
 1000.0
 1500.0
```
"""
function band_weight(
    bands::AbstractMatrix{<:Real},
    weight::AbstractSpectralWeight = FlatWeight(),
)
    OutFT = float(promote_type(eltype(bands), _weight_eltype(weight)))
    dest = similar(bands, OutFT, axes(bands, 2))
    return band_weight!(dest, bands, weight)
end

"""
    CanopyOptics.band_weight!(dest, bands, weight = FlatWeight()) -> dest

Overwrite `dest` in place with the band-integrated weights of `bands` under
`weight`, and return `dest` (zero heap allocation when `dest`, `bands`, and built-in
`weight` share a matching `Float32` or `Float64` element type; `dest` must not alias
`bands` or `weight`'s storage). See [`CanopyOptics.band_weight`](@ref).
"""
function band_weight!(
    dest::AbstractVector,
    bands::AbstractMatrix{<:Real},
    weight::AbstractSpectralWeight = FlatWeight(),
)
    (_arrays_might_alias(dest, bands) || _weight_might_alias(dest, weight)) && _throw_aliased_dest()
    nbnd = _validate_bands_shape(bands)
    length(dest) == nbnd || _throw_dest_dim(length(dest), nbnd)
    r0 = firstindex(bands, 1)
    c0 = firstindex(bands, 2)
    d0 = firstindex(dest)
    DestFT = eltype(dest)
    @inbounds for b in 1:nbnd
        ν_lo = bands[r0, c0 + b - 1]
        ν_hi = bands[r0 + 1, c0 + b - 1]
        p_lo = prospect_primal(ν_lo)
        p_hi = prospect_primal(ν_hi)
        (isfinite(p_lo) && isfinite(p_hi) && p_lo >= zero(p_lo) && p_hi > p_lo) ||
            _throw_invalid_band(b, p_lo, p_hi)
        dest[d0 + b - 1] = DestFT(weight_interval_integral(weight, ν_lo, ν_hi))
    end
    return dest
end

"""
    subband_weights(sub_bands, parent_band, weight = FlatWeight())

Compute normalized dimensionless sub-band flux shares `w_j` (`0 <= w_j <= 1`,
`sum(w_j) ≈ 1`) for a contiguous partition `sub_bands` (`(2, m)` matrix) of a single
`parent_band` (`(2, 1)` matrix or 2-element vector/tuple `[ν_lo, ν_hi]`).

`subband_weights` resolves the effective weight once on `parent_band`, matching
`band_average(v, edges, sub_bands, weight; parent_band = parent_band)`:
- If `weight` is a [`TabulatedWeight`](@ref) that does not cover `parent_band`, or
  whose integral across the entire `parent_band` is zero, all `m` sub-bands fall back
  jointly to `FlatWeight` (`Δν_j / Δν_parent`); if `parent_band` is covered with
  positive integrated weight while an individual sub-band has zero density (e.g. a
  dark `< 280 nm` UV sub-interval), that sub-band receives share `0`.
- If `weight` is a [`PlanckWeight`](@ref), all sub-bands share `parent_band`'s
  relative Wien-tail exponent shift so relative shares never underflow to `0/0` in
  `Float32` or `Float64`.

Consequently,
`sum(subband_weights(sub_bands, parent_band, weight) .* band_average(v, edges, sub_bands, weight; parent_band = parent_band))`
reproduces `only(band_average(v, edges, parent_band, weight))` across all weights,
coverage regimes, and precisions. See [`subband_weights!`](@ref) for the in-place
variant (zero heap allocation when arrays and built-in `weight` share a matching
`Float32` or `Float64` element type).

# Examples
```jldoctest
julia> sub = [12850.0 14285.7; 14285.7 16000.0];

julia> w = subband_weights(sub, [12850.0; 16000.0;;], FlatWeight());

julia> isapprox(sum(w), 1.0)
true
```
"""
function subband_weights(
    sub_bands::AbstractMatrix{<:Real},
    parent_band,
    weight::AbstractSpectralWeight = FlatWeight(),
)
    p_lo, p_hi = _extract_parent_band_limits(parent_band)
    OutFT = float(promote_type(eltype(sub_bands), typeof(p_lo), typeof(p_hi), _weight_eltype(weight)))
    dest = similar(sub_bands, OutFT, axes(sub_bands, 2))
    return subband_weights!(dest, sub_bands, parent_band, weight)
end

"""
    subband_weights!(dest, sub_bands, parent_band, weight = FlatWeight()) -> dest

Overwrite `dest` in place with the normalized sub-band shares of `sub_bands` within
`parent_band` under `weight`, and return `dest` (zero heap allocation when `dest`,
`sub_bands`, `parent_band`, and built-in `weight` share a matching `Float32` or
`Float64` element type; `dest` must not alias `sub_bands`, `parent_band`, or
`weight`'s storage). See [`subband_weights`](@ref).
"""
function subband_weights!(
    dest::AbstractVector,
    sub_bands::AbstractMatrix{<:Real},
    parent_band,
    weight::AbstractSpectralWeight = FlatWeight(),
)
    (_arrays_might_alias(dest, sub_bands) ||
     _weight_might_alias(dest, weight) ||
     (parent_band isa AbstractArray && _arrays_might_alias(dest, parent_band))) &&
        _throw_aliased_dest()

    m = _validate_bands_shape(sub_bands)
    m >= 1 || throw(ArgumentError("sub_bands must contain at least 1 sub-band"))
    length(dest) == m || _throw_dest_dim(length(dest), m)

    par_lo, par_hi = _extract_parent_band_limits(parent_band)
    p_par_lo = prospect_primal(par_lo)
    p_par_hi = prospect_primal(par_hi)
    (isfinite(p_par_lo) && isfinite(p_par_hi) && p_par_lo >= zero(p_par_lo) && p_par_hi > p_par_lo) ||
        _throw_invalid_band(1, p_par_lo, p_par_hi)

    r0 = firstindex(sub_bands, 1)
    c0 = firstindex(sub_bands, 2)
    d0 = firstindex(dest)
    @inbounds begin
        prospect_primal(sub_bands[r0, c0]) == p_par_lo ||
            throw(ArgumentError("sub_bands must start at parent_band lower limit ($p_par_lo cm^-1)"))
        for j in 1:m
            s_lo = sub_bands[r0, c0 + j - 1]
            s_hi = sub_bands[r0 + 1, c0 + j - 1]
            ps_lo = prospect_primal(s_lo)
            ps_hi = prospect_primal(s_hi)
            (isfinite(ps_lo) && isfinite(ps_hi) && ps_hi > ps_lo) || _throw_invalid_band(j, ps_lo, ps_hi)
            if j > 1
                ps_lo == prospect_primal(sub_bands[r0 + 1, c0 + j - 2]) ||
                    throw(ArgumentError("sub_bands must be contiguous without gaps or overlaps"))
            end
        end
        prospect_primal(sub_bands[r0 + 1, c0 + m - 1]) == p_par_hi ||
            throw(ArgumentError("sub_bands must end at parent_band upper limit ($p_par_hi cm^-1)"))
    end

    DestFT = eltype(dest)
    AccFT = promote_type(DestFT, eltype(sub_bands), typeof(par_lo), typeof(par_hi), _weight_eltype(weight), Float64)
    parent_span = AccFT(par_hi) - AccFT(par_lo)
    use_flat, ν_shift, use_shift, tot_parent = _resolve_parent_weight(weight, par_lo, par_hi, AccFT)

    if use_flat
        @inbounds for j in 1:m
            dν = AccFT(sub_bands[r0 + 1, c0 + j - 1]) - AccFT(sub_bands[r0, c0 + j - 1])
            dest[d0 + j - 1] = DestFT(dν / parent_span)
        end
        return dest
    elseif weight isa PlanckWeight
        tot = zero(AccFT)
        @inbounds for j in 1:m
            s_lo = sub_bands[r0, c0 + j - 1]
            s_hi = sub_bands[r0 + 1, c0 + j - 1]
            tot += _integrate_planck_interval(s_lo, s_hi, weight.T, ν_shift, use_shift)
        end
        if prospect_primal(tot) > 0
            @inbounds for j in 1:m
                s_lo = sub_bands[r0, c0 + j - 1]
                s_hi = sub_bands[r0 + 1, c0 + j - 1]
                wj = _integrate_planck_interval(s_lo, s_hi, weight.T, ν_shift, use_shift)
                dest[d0 + j - 1] = DestFT(wj / tot)
            end
        else
            @inbounds for j in 1:m
                dν = AccFT(sub_bands[r0 + 1, c0 + j - 1]) - AccFT(sub_bands[r0, c0 + j - 1])
                dest[d0 + j - 1] = DestFT(dν / parent_span)
            end
        end
        return dest
    elseif weight isa TabulatedWeight
        @inbounds for j in 1:m
            s_lo = sub_bands[r0, c0 + j - 1]
            s_hi = sub_bands[r0 + 1, c0 + j - 1]
            wj = _tabulated_covered_integral(weight, s_lo, s_hi, AccFT)
            dest[d0 + j - 1] = DestFT(wj / tot_parent)
        end
        return dest
    else
        tot = zero(AccFT)
        @inbounds for j in 1:m
            s_lo = sub_bands[r0, c0 + j - 1]
            s_hi = sub_bands[r0 + 1, c0 + j - 1]
            tot += AccFT(weight_interval_integral(weight, s_lo, s_hi))
        end
        if prospect_primal(tot) > 0
            @inbounds for j in 1:m
                s_lo = sub_bands[r0, c0 + j - 1]
                s_hi = sub_bands[r0 + 1, c0 + j - 1]
                wj = AccFT(weight_interval_integral(weight, s_lo, s_hi))
                dest[d0 + j - 1] = DestFT(wj / tot)
            end
        else
            @inbounds for j in 1:m
                dν = AccFT(sub_bands[r0 + 1, c0 + j - 1]) - AccFT(sub_bands[r0, c0 + j - 1])
                dest[d0 + j - 1] = DestFT(dν / parent_span)
            end
        end
        return dest
    end
end

@inline function _extract_parent_band_limits(parent_band::AbstractMatrix{<:Real})
    size(parent_band) == (2, 1) ||
        throw(DimensionMismatch("parent_band matrix must have size (2, 1), got $(size(parent_band))"))
    r0 = firstindex(parent_band, 1)
    c0 = firstindex(parent_band, 2)
    return @inbounds parent_band[r0, c0], parent_band[r0 + 1, c0]
end

@inline function _extract_parent_band_limits(parent_band::AbstractVector{<:Real})
    length(parent_band) == 2 ||
        throw(DimensionMismatch("parent_band vector must have length 2, got $(length(parent_band))"))
    i0 = firstindex(parent_band)
    return @inbounds parent_band[i0], parent_band[i0 + 1]
end

@inline function _extract_parent_band_limits(parent_band::Tuple{Real,Real})
    return parent_band[1], parent_band[2]
end

@noinline function _extract_parent_band_limits(parent_band)
    throw(ArgumentError(
        "parent_band must be a (2, 1) matrix, 2-element vector, or 2-tuple of real limits, got $(typeof(parent_band))",
    ))
end

# Fast linear two-pointer resampling of fine bin means `src_vals` on `src_edges`
# onto covered target grid `dst_edges` (both strictly ascending).
function _resample_wavenumber_bins!(
    dst_vals::AbstractVector{OutFT},
    src_vals::AbstractVector,
    src_edges::AbstractVector,
    dst_edges::AbstractVector,
) where {OutFT}
    n_dst = length(dst_edges) - 1
    n_src = length(src_vals)
    i0_se = firstindex(src_edges)
    i0_sv = firstindex(src_vals)
    i0_de = firstindex(dst_edges)
    i0_dv = firstindex(dst_vals)
    AccFT = promote_type(OutFT, eltype(src_vals), eltype(src_edges), eltype(dst_edges), Float64)

    s = _first_overlapping_bin(src_edges, dst_edges[i0_de])
    @inbounds for d in 0:(n_dst - 1)
        ν_lo = dst_edges[i0_de + d]
        ν_hi = dst_edges[i0_de + d + 1]
        p_lo = prospect_primal(ν_lo)
        p_hi = prospect_primal(ν_hi)
        while s + 1 < n_src && prospect_primal(src_edges[i0_se + s + 1]) <= p_lo
            s += 1
        end
        k = s
        num = zero(AccFT)
        den = zero(AccFT)
        while k < n_src
            a = _clamp_overlap_lo(ν_lo, src_edges[i0_se + k], AccFT)
            b = _clamp_overlap_hi(ν_hi, src_edges[i0_se + k + 1], AccFT)
            dν = b - a
            if prospect_primal(dν) > 0
                num += AccFT(src_vals[i0_sv + k]) * dν
                den += dν
            end
            prospect_primal(src_edges[i0_se + k + 1]) >= p_hi && break
            k += 1
        end
        s = k > 0 && prospect_primal(src_edges[i0_se + k]) > p_hi ? k - 1 : k
        s = min(s, n_src - 1)
        dst_vals[i0_dv + d] = OutFT(num / den)
    end
    return dst_vals
end

@inline function _validate_leaf_rt_pair(rho, tau, ::Type{BFT}) where {BFT<:AbstractFloat}
    rp = prospect_primal(rho)
    tp = prospect_primal(tau)
    tol = BFT(16) * eps(BFT)
    (isfinite(rp) && isfinite(tp) && rp >= zero(rp) && tp >= zero(tp) &&
     rp <= one(rp) + tol && tp <= one(tp) + tol && rp + tp <= one(rp) + tol) ||
        throw(ArgumentError(
            "leaf reflectance and transmittance must be finite, nonnegative, and satisfy rho + tau <= 1 (got rho = $rp, tau = $tp)",
        ))
    return nothing
end

"""
    TabulatedLeaf(edges, rho, tau)

Measured or precomputed leaf reflectance `rho` and transmittance `tau` bin means on
ascending wavenumber bin edges `edges` in `cm⁻¹` (`length(edges) == length(rho) + 1 == length(tau) + 1`).

All bin values must be finite, nonnegative, and satisfy `rho[i] + tau[i] <= 1`
(within floating-point rounding). Evaluate on any covered wavenumber grid via
[`leaf_spectrum`](@ref).

# Examples
```jldoctest
julia> leaf = TabulatedLeaf([800.0, 2500.0, 4000.0], [0.03, 0.04], [0.0, 0.0]);

julia> leaf_spectrum(leaf, [1000.0, 2000.0, 3500.0])
([0.03, 0.03666666666666667], [0.0, 0.0])
```
"""
struct TabulatedLeaf{FT<:Real}
    edges::Vector{FT}
    rho::Vector{FT}
    tau::Vector{FT}

    function TabulatedLeaf{FT}(
        edges::AbstractVector,
        rho::AbstractVector,
        tau::AbstractVector,
    ) where {FT<:Real}
        ne = length(edges)
        ne >= 2 || _throw_edges_length()
        n_bins = ne - 1
        (length(rho) == n_bins && length(tau) == n_bins) || throw(DimensionMismatch(
            "rho ($(length(rho))) and tau ($(length(tau))) lengths must equal length(edges) - 1 ($n_bins)",
        ))
        e_vec = _copy_canonical_vec(FT, edges)
        r_vec = _copy_canonical_vec(FT, rho)
        t_vec = _copy_canonical_vec(FT, tau)
        _validate_wavenumber_edges(e_vec)
        BFT = _leaf_base_float(FT)
        @inbounds for k in 1:n_bins
            _validate_leaf_rt_pair(r_vec[k], t_vec[k], BFT)
        end
        return new{FT}(e_vec, r_vec, t_vec)
    end
end

function TabulatedLeaf(
    edges::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    tau::AbstractVector{<:Real},
)
    FT = float(promote_type(eltype(edges), eltype(rho), eltype(tau)))
    return TabulatedLeaf{FT}(edges, rho, tau)
end

function Base.:(==)(a::TabulatedLeaf, b::TabulatedLeaf)
    return a.edges == b.edges && a.rho == b.rho && a.tau == b.tau
end

function Base.hash(a::TabulatedLeaf, h::UInt)
    return hash(a.tau, hash(a.rho, hash(a.edges, hash(:TabulatedLeaf, h))))
end

"""
    ConstantLeaf(rho, tau; ν_min = ..., ν_max = ...)

Wavenumber-independent leaf reflectance `rho` and transmittance `tau` over
`[ν_min, ν_max]` in `cm⁻¹` (defaulting to `[0, ∞)` when `ν_min` and `ν_max` are
omitted).

`rho` and `tau` must be finite, nonnegative, and satisfy `rho + tau <= 1`. Pass to
[`leaf_spectrum`](@ref) directly or as a regional/fallback component inside
[`CompositeLeaf`](@ref).

# Examples
```jldoctest
julia> cl = ConstantLeaf(0.03, 0.0);

julia> leaf_spectrum(cl, [100.0, 500.0, 1000.0])
([0.03, 0.03], [0.0, 0.0])
```
"""
struct ConstantLeaf{FT<:Real}
    rho::FT
    tau::FT
    ν_min::FT
    ν_max::FT

    function ConstantLeaf{FT}(rho::Real, tau::Real, ν_min::Real, ν_max::Real) where {FT<:Real}
        BFT = _leaf_base_float(FT)
        r_v = _canon_zero(FT(rho))
        t_v = _canon_zero(FT(tau))
        lo_v = _canon_zero(FT(ν_min))
        hi_v = _canon_zero(FT(ν_max))
        _validate_leaf_rt_pair(r_v, t_v, BFT)
        p_lo = prospect_primal(lo_v)
        p_hi = prospect_primal(hi_v)
        (isfinite(p_lo) && p_lo >= zero(p_lo) && !isnan(p_hi) && p_hi > p_lo) ||
            throw(ArgumentError("ConstantLeaf domain [ν_min, ν_max] must satisfy 0 <= ν_min < ν_max"))
        return new{FT}(r_v, t_v, lo_v, hi_v)
    end
end

function ConstantLeaf(
    rho::Real,
    tau::Real;
    ν_min = _UNSET_LEAF_KW,
    ν_max = _UNSET_LEAF_KW,
)
    FT = float(promote_type(typeof(rho), typeof(tau), _kw_type(ν_min), _kw_type(ν_max)))
    BFT = _leaf_base_float(FT)
    lo = ν_min === _UNSET_LEAF_KW ? zero(FT) : FT(ν_min)
    hi = ν_max === _UNSET_LEAF_KW ? FT(BFT(Inf)) : FT(ν_max)
    return ConstantLeaf{FT}(FT(rho), FT(tau), lo, hi)
end

function Base.:(==)(a::ConstantLeaf, b::ConstantLeaf)
    return a.rho == b.rho && a.tau == b.tau && a.ν_min == b.ν_min && a.ν_max == b.ν_max
end

function Base.hash(a::ConstantLeaf, h::UInt)
    h = hash(:ConstantLeaf, h)
    h = hash(_canon_zero(a.rho), h)
    h = hash(_canon_zero(a.tau), h)
    h = hash(_canon_zero(a.ν_min), h)
    return hash(_canon_zero(a.ν_max), h)
end

"""
    CompositeLeaf(models...)

Ordered composition of leaf spectrum models in which **the first model covering each
wavenumber bin wins**.

Nested `CompositeLeaf` instances are flattened on construction. During
[`leaf_spectrum(comp, edges)`](@ref), bin `i` (`[edges[i], edges[i + 1]]`) is
assigned to the earliest model in `models` whose domain covers the entire bin
`[edges[i], edges[i + 1]]`; if any bin of `edges` remains uncovered after all models
have been checked, an `ArgumentError` naming the uncovered bin is thrown. When splicing adjacent spectral regions, bound each regional [`ConstantLeaf`](@ref) to
its physical wavenumber range and include the regional domain boundaries in `edges` so no
bin straddles a boundary (for example, extending `400 nm` PROSPECT optics onto
`ν >= ν_pro_hi`, `ConstantLeaf(0.04, 0.0; ν_min = 2500.0, ν_max = ν_pro_lo)` across
`2.5–4.0 µm`, and `ConstantLeaf(0.03, 0.0; ν_min = 0.0, ν_max = 2500.0)` beyond
`4.0 µm`, with `[2500.0, ν_pro_lo, ν_pro_hi]` inserted into `edges`).

# Examples
```jldoctest
julia> comp = CompositeLeaf(
           TabulatedLeaf([10000.0, 20000.0], [0.45], [0.40]),
           ConstantLeaf(0.03, 0.0),
       );

julia> leaf_spectrum(comp, [2000.0, 10000.0, 20000.0, 30000.0])
([0.03, 0.45, 0.03], [0.0, 0.4, 0.0])
```
"""
struct CompositeLeaf{FT,T<:Tuple}
    models::T
end

_leaf_model_eltype(::LeafProspectProProperties{FT}) where {FT} = FT
_leaf_model_eltype(::TabulatedLeaf{FT}) where {FT} = FT
_leaf_model_eltype(::ConstantLeaf{FT}) where {FT} = FT
_leaf_model_eltype(::CompositeLeaf{FT}) where {FT} = FT
_leaf_model_eltype(::Any) = Union{}

@inline _flatten_leaf_model(c::CompositeLeaf) = c.models
@inline _flatten_leaf_model(m) = (m,)

@inline _flatten_leaf_models() = ()
@inline _flatten_leaf_models(m, rest...) =
    (_flatten_leaf_model(m)..., _flatten_leaf_models(rest...)...)

@inline _promote_leaf_models_eltype(::Tuple{}) = Union{}
@inline _promote_leaf_models_eltype(t::Tuple) =
    promote_type(_leaf_model_eltype(first(t)), _promote_leaf_models_eltype(Base.tail(t)))

function CompositeLeaf(models...)
    isempty(models) && throw(ArgumentError("CompositeLeaf requires at least one leaf model"))
    flat = _flatten_leaf_models(models...)
    T_raw = _promote_leaf_models_eltype(flat)
    FT = T_raw === Union{} ? Union{} : float(T_raw)
    return CompositeLeaf{FT,typeof(flat)}(flat)
end

Base.:(==)(a::CompositeLeaf, b::CompositeLeaf) = a.models == b.models
Base.hash(a::CompositeLeaf, h::UInt) = hash(a.models, hash(:CompositeLeaf, h))

"""
    CanopyOptics.spectral_domain(model) -> (ν_min, ν_max)

Return the contiguous wavenumber coverage interval `[ν_min, ν_max]` in `cm⁻¹` of an
interval-domain leaf model (`LeafProspectProProperties`, [`TabulatedLeaf`](@ref), or
[`ConstantLeaf`](@ref)), or `nothing` when coverage must be probed per bin.

This qualified extension hook is not exported. Custom leaf models with a contiguous
wavenumber domain can specialize `CanopyOptics.spectral_domain(model)` so
[`CompositeLeaf`](@ref) evaluates them in a single vectorized slice call; models with
non-contiguous domains that leave `spectral_domain(model) === nothing` may signal an
uncovered bin in `leaf_spectrum(model, bin_edges)` by throwing
[`CanopyOptics.SpectralCoverageError`](@ref) (or returning `nothing`).
"""
spectral_domain(::LeafProspectProProperties) = (_PROSPECT_MIN_NU_CM, _PROSPECT_MAX_NU_CM)
spectral_domain(m::TabulatedLeaf) = (first(m.edges), last(m.edges))
spectral_domain(m::ConstantLeaf) = (m.ν_min, m.ν_max)
spectral_domain(::Any) = nothing

"""
    CanopyOptics.SpectralCoverageError(msg)

Exception thrown by custom [`leaf_spectrum`](@ref) implementations (when
[`CanopyOptics.spectral_domain`](@ref) is `nothing`) to report that a queried bin lies
outside the model's spectral coverage.

[`CompositeLeaf`](@ref) catches `SpectralCoverageError` (and `nothing` returns) during
per-bin probing and delegates the uncovered bin to subsequent component models, while
allowing `ArgumentError` and all other exceptions from invalid parameters or coefficients
to propagate untouched.
"""
struct SpectralCoverageError <: Exception
    msg::String
end

Base.showerror(io::IO, e::SpectralCoverageError) = print(io, "SpectralCoverageError: ", e.msg)

"""
    leaf_spectrum(model, edges) -> (rho, tau)
    leaf_spectrum(leaf::LeafProspectProProperties, edges; optics = nothing) -> (rho, tau)

Evaluate bin-mean leaf hemispherical reflectance `rho` and transmittance `tau` on
ascending wavenumber bin edges `edges` in `cm⁻¹` (`length(rho) == length(tau) == length(edges) - 1`).

`model` must fully cover every bin of `edges`; otherwise an `ArgumentError` is
thrown (never returning silent zeros). Built-in methods support:
- [`LeafProspectProProperties`](@ref): evaluates [`prospect_budget`](@ref) on the
  bundled `1 nm` PROSPECT-PRO optical table (`399.5–2500.5 nm`,
  `3999.20–25031.29 cm⁻¹`) and averages onto `edges` by exact wavenumber overlap
  (accepting an optional `optics` keyword to supply a custom regular-nm
  [`LeafOpticalProperties`](@ref) table).
- [`TabulatedLeaf`](@ref): resamples measured `(rho, tau)` spectra onto `edges` by
  exact wavenumber overlap.
- [`ConstantLeaf`](@ref): fills `(rho, tau)` across `[ν_min, ν_max]`.
- [`CompositeLeaf`](@ref): selects the first component model covering each bin.

For PROSPECT pigment absorption shares per bin, see [`leaf_pigment_shares`](@ref).

# Examples
```jldoctest
julia> edges = wavelength_to_wavenumber(400.0:1.0:2500.0);

julia> rho, tau = leaf_spectrum(LeafProspectProProperties(), edges);

julia> length(rho) == 2101 && all(0 .< rho .+ tau .< 1)
true
```
"""
function leaf_spectrum end

function leaf_spectrum(
    leaf::LeafProspectProProperties{FT_L},
    edges::AbstractVector{<:Real};
    optics::Union{Nothing,LeafOpticalProperties} = nothing,
) where {FT_L}
    n_bins = _validate_wavenumber_edges(edges)
    OutFT = float(promote_type(FT_L, eltype(edges)))
    BFT = _leaf_base_float(OutFT)

    opt, ν_native = _prepare_prospect_grid_and_optics(edges, optics, BFT)
    budget = prospect_budget(leaf, opt)

    # Reverse native 1-nm wavelength arrays into ascending wavenumber order.
    r_nat = reverse(budget.reflectance)
    t_nat = reverse(budget.transmission)

    rho = Vector{OutFT}(undef, n_bins)
    tau = Vector{OutFT}(undef, n_bins)
    if !(eltype(edges) <: ForwardDiff.Dual) && length(edges) == length(ν_native) && edges == ν_native
        copyto!(rho, r_nat)
        copyto!(tau, t_nat)
    else
        _resample_wavenumber_bins!(rho, r_nat, ν_native, edges)
        _resample_wavenumber_bins!(tau, t_nat, ν_native, edges)
    end
    return rho, tau
end

function _prepare_prospect_grid_and_optics(
    edges::AbstractVector,
    optics::Union{Nothing,LeafOpticalProperties},
    ::Type{BFT},
) where {BFT<:AbstractFloat}
    e_lo = prospect_primal(first(edges))
    e_hi = prospect_primal(last(edges))
    if optics === nothing
        (e_lo >= _PROSPECT_MIN_NU_CM && e_hi <= _PROSPECT_MAX_NU_CM) ||
            _throw_uncovered_grid(e_lo, e_hi, _PROSPECT_MIN_NU_CM, _PROSPECT_MAX_NU_CM)
        opt = _native_prospect_optics(BFT)
        ν_native = wavelength_to_wavenumber(BFT(400):one(BFT):BFT(2500))
        if BFT(e_lo) < ν_native[1]
            ν_native[1] = BFT(e_lo)
        end
        if BFT(e_hi) > ν_native[end]
            ν_native[end] = BFT(e_hi)
        end
        return opt, ν_native
    else
        λ_nm = ustrip.(unit2nm(optics.λ))
        ν_native = wavelength_to_wavenumber(BFT.(λ_nm))
        m_lo = prospect_primal(first(ν_native))
        m_hi = prospect_primal(last(ν_native))
        (e_lo >= m_lo && e_hi <= m_hi) || _throw_uncovered_grid(e_lo, e_hi, m_lo, m_hi)
        return optics, ν_native
    end
end

function leaf_spectrum(model::TabulatedLeaf{FT_M}, edges::AbstractVector{<:Real}) where {FT_M}
    n_bins = _validate_wavenumber_edges(edges)
    e_lo = prospect_primal(first(edges))
    e_hi = prospect_primal(last(edges))
    m_lo = prospect_primal(first(model.edges))
    m_hi = prospect_primal(last(model.edges))
    (e_lo >= m_lo && e_hi <= m_hi) || _throw_uncovered_grid(e_lo, e_hi, m_lo, m_hi)

    OutFT = float(promote_type(FT_M, eltype(edges)))
    rho = Vector{OutFT}(undef, n_bins)
    tau = Vector{OutFT}(undef, n_bins)
    if !(eltype(edges) <: ForwardDiff.Dual) && !(FT_M <: ForwardDiff.Dual) &&
       length(edges) == length(model.edges) && edges == model.edges
        copyto!(rho, model.rho)
        copyto!(tau, model.tau)
    else
        _resample_wavenumber_bins!(rho, model.rho, model.edges, edges)
        _resample_wavenumber_bins!(tau, model.tau, model.edges, edges)
    end
    return rho, tau
end

function leaf_spectrum(model::ConstantLeaf{FT_M}, edges::AbstractVector{<:Real}) where {FT_M}
    n_bins = _validate_wavenumber_edges(edges)
    e_lo = prospect_primal(first(edges))
    e_hi = prospect_primal(last(edges))
    m_lo = prospect_primal(model.ν_min)
    m_hi = prospect_primal(model.ν_max)
    (e_lo >= m_lo && e_hi <= m_hi) || _throw_uncovered_grid(e_lo, e_hi, m_lo, m_hi)

    OutFT = float(promote_type(FT_M, eltype(edges)))
    rho = fill(OutFT(model.rho), n_bins)
    tau = fill(OutFT(model.tau), n_bins)
    return rho, tau
end

function leaf_spectrum(comp::CompositeLeaf{FT_C}, edges::AbstractVector{<:Real}) where {FT_C}
    n_bins = _validate_wavenumber_edges(edges)
    OutFT = float(promote_type(FT_C, eltype(edges)))
    rho = Vector{OutFT}(undef, n_bins)
    tau = Vector{OutFT}(undef, n_bins)
    covered = falses(n_bins)
    _apply_composite_models!(rho, tau, covered, comp.models, edges)
    i0_e = firstindex(edges)
    @inbounds for i in 1:n_bins
        covered[i] || _throw_uncovered_bin(
            i,
            prospect_primal(edges[i0_e + i - 1]),
            prospect_primal(edges[i0_e + i]),
        )
    end
    return rho, tau
end

@inline function _apply_composite_models!(
    ::Vector{OutFT},
    ::Vector{OutFT},
    ::BitVector,
    ::Tuple{},
    ::AbstractVector,
) where {OutFT}
    return nothing
end

function _apply_composite_models!(
    rho::Vector{OutFT},
    tau::Vector{OutFT},
    covered::BitVector,
    models::Tuple{Any,Vararg{Any}},
    edges::AbstractVector,
) where {OutFT}
    all(covered) && return nothing
    model = first(models)
    n_bins = length(rho)
    i0_e = firstindex(edges)
    dom = spectral_domain(model)
    if dom !== nothing
        m_lo = prospect_primal(dom[1])
        m_hi = prospect_primal(dom[2])
        first_win = 0
        last_win = 0
        @inbounds for i in 1:n_bins
            if !covered[i]
                e_lo = prospect_primal(edges[i0_e + i - 1])
                e_hi = prospect_primal(edges[i0_e + i])
                if e_lo >= m_lo && e_hi <= m_hi
                    first_win == 0 && (first_win = i)
                    last_win = i
                end
            end
        end
        if first_win > 0
            n_sub = last_win - first_win + 1
            sub_edges = @view edges[(i0_e + first_win - 1):(i0_e + last_win)]
            r_sub, t_sub = leaf_spectrum(model, sub_edges)
            (length(r_sub) == n_sub && length(t_sub) == n_sub) || throw(DimensionMismatch(
                "leaf_spectrum returned ($(length(r_sub)), $(length(t_sub))) bins for $n_sub-bin slice",
            ))
            i0_r = firstindex(r_sub)
            i0_t = firstindex(t_sub)
            @inbounds for i in first_win:last_win
                if !covered[i]
                    e_lo = prospect_primal(edges[i0_e + i - 1])
                    e_hi = prospect_primal(edges[i0_e + i])
                    if e_lo >= m_lo && e_hi <= m_hi
                        idx = i - first_win
                        rho[i] = OutFT(r_sub[i0_r + idx])
                        tau[i] = OutFT(t_sub[i0_t + idx])
                        covered[i] = true
                    end
                end
            end
        end
    else
        # Fallback for custom leaf models that do not specialize `spectral_domain`:
        # probe uncovered bins individually.
        bin_buf = Vector{eltype(edges)}(undef, 2)
        @inbounds for i in 1:n_bins
            if !covered[i]
                bin_buf[1] = edges[i0_e + i - 1]
                bin_buf[2] = edges[i0_e + i]
                ok, r_i, t_i = _try_single_bin_spectrum(model, bin_buf, OutFT)
                if ok
                    rho[i] = r_i
                    tau[i] = t_i
                    covered[i] = true
                end
            end
        end
    end
    return _apply_composite_models!(rho, tau, covered, Base.tail(models), edges)
end

function _try_single_bin_spectrum(model, bin_edges, ::Type{OutFT}) where {OutFT}
    try
        res = leaf_spectrum(model, bin_edges)
        res === nothing && return false, zero(OutFT), zero(OutFT)
        r, t = res
        (length(r) == 1 && length(t) == 1) || throw(DimensionMismatch(
            "leaf_spectrum must return 1-element vectors for a 1-bin grid, got ($(length(r)), $(length(t)))",
        ))
        return true, OutFT(first(r)), OutFT(first(t))
    catch err
        err isa SpectralCoverageError && return false, zero(OutFT), zero(OutFT)
        rethrow()
    end
end

@doc raw"""
    leaf_pigment_shares(leaf::LeafProspectProProperties, edges; optics = nothing) -> (; chlorophyll, carotenoids)

Compute per-bin chlorophyll and carotenoid shares of absorbed light on ascending
wavenumber bin edges `edges` in `cm⁻¹` for a [`LeafProspectProProperties`](@ref)
leaf, returning `(; chlorophyll, carotenoids)` with `length(edges) - 1` elements
each.

`edges` must be fully covered by the PROSPECT-PRO domain (`399.5–2500.5 nm`,
`3999.20–25031.29 cm⁻¹`), or an `ArgumentError` is thrown. At each native `1 nm`
PROSPECT bin, directly computed material absorptance ``\alpha`` from
[`prospect_budget`](@ref) is partitioned across [`prospect_constituent_coefficients`](@ref)
using [`HomogeneousAbsorption`](@ref) and [`partition_absorption`](@ref), with
chlorophyll, carotenoids, anthocyanins, brown material, water, bulk dry matter, protein,
and carbon constituents included in the denominator ``\sum_j C_j K_j(\lambda)``.

# Zero-Absorption and Coarse-Band Aggregation Conventions
- Native `1 nm` bins with zero total absorption coefficient (`AbsorptionNoMaterial`)
  return `0.0` for both pigment shares (rather than `NaN`). Multi-bin intervals in
  `edges` with positive integrated absorptance weight native shares by ``\alpha(\nu)\,d\nu``;
  if integrated absorptance underflows to zero while constituent coefficients remain
  nonzero, resampling falls back to the flat-in-wavenumber mean of ``C_i K_i / \sum_j C_j K_j``.
- For incident-spectrum-weighted coarse bands, obtain native-grid absorptance from
  `prospect_budget(leaf, LeafOpticalProperties()).absorptance` (rather than `1 - rho - tau`,
  which can round to zero in weak-absorption bins) via [`wavelength_to_wavenumber`](@ref),
  and divide `band_average(α .* shares.chlorophyll, edges, bands, w)` by
  `band_average(α, edges, bands, w)` (returning `0` when band-mean absorptance is zero).
  For photon-count shares, replace `w` in both calls with a photon-flux weight integrating
  ``w(\nu) / (100\,h c\,\nu)`` across each bin.

# Examples
```jldoctest
julia> edges = wavelength_to_wavenumber(400.0:1.0:700.0);

julia> shares = leaf_pigment_shares(LeafProspectProProperties(), edges);

julia> all(0 .<= shares.chlorophyll .+ shares.carotenoids .<= 1)
true
```
"""
function leaf_pigment_shares(
    leaf::LeafProspectProProperties{FT_L},
    edges::AbstractVector{<:Real};
    optics::Union{Nothing,LeafOpticalProperties} = nothing,
) where {FT_L}
    n_bins = _validate_wavenumber_edges(edges)
    OutFT = float(promote_type(FT_L, eltype(edges)))
    BFT = _leaf_base_float(OutFT)

    opt, ν_native = _prepare_prospect_grid_and_optics(edges, optics, BFT)
    budget = prospect_budget(leaf, opt)
    products = prospect_constituent_coefficients(leaf, opt)

    n_nat = length(budget.absorptance)
    PartFT = eltype(budget.absorptance)
    f_chl_nat = Vector{PartFT}(undef, n_nat)
    f_car_nat = Vector{PartFT}(undef, n_nat)
    α_nat = Vector{PartFT}(undef, n_nat)
    α_chl_nat = Vector{PartFT}(undef, n_nat)
    α_car_nat = Vector{PartFT}(undef, n_nat)

    @inbounds for k in 1:n_nat
        rev_k = n_nat + 1 - k
        q_abs = budget.absorptance[rev_k]
        coeffs_k = (
            chlorophyll = PartFT(products.chlorophyll[rev_k]),
            carotenoids = PartFT(products.carotenoids[rev_k]),
            anthocyanins = PartFT(products.anthocyanins[rev_k]),
            brown_material = PartFT(products.brown_material[rev_k]),
            water = PartFT(products.water[rev_k]),
            dry_matter = PartFT(products.dry_matter[rev_k]),
            protein = PartFT(products.protein[rev_k]),
            carbon_constituents = PartFT(products.carbon_constituents[rev_k]),
        )
        f_chl, f_car, q_chl, q_car = _native_pigment_partition(q_abs, coeffs_k)
        f_chl_nat[k] = f_chl
        f_car_nat[k] = f_car
        α_nat[k] = q_abs
        α_chl_nat[k] = q_chl
        α_car_nat[k] = q_car
    end

    chl_out = Vector{OutFT}(undef, n_bins)
    car_out = Vector{OutFT}(undef, n_bins)
    if !(eltype(edges) <: ForwardDiff.Dual) && length(edges) == length(ν_native) && edges == ν_native
        copyto!(chl_out, f_chl_nat)
        copyto!(car_out, f_car_nat)
    else
        _resample_absorption_shares!(
            chl_out, car_out,
            α_nat, α_chl_nat, α_car_nat, f_chl_nat, f_car_nat,
            ν_native, edges,
        )
    end
    return (chlorophyll = chl_out, carotenoids = car_out)
end

@inline function _native_pigment_partition(
    q_abs::FT,
    coeffs::NamedTuple{Names,NTuple{N,FT}},
) where {FT<:AbstractFloat,Names,N}
    res = partition_absorption(HomogeneousAbsorption(), q_abs, coeffs)
    if res.status === AbsorptionNoMaterial
        return zero(FT), zero(FT), zero(FT), zero(FT)
    end
    return (
        res.fractions.chlorophyll,
        res.fractions.carotenoids,
        res.constituents.chlorophyll,
        res.constituents.carotenoids,
    )
end

@inline function _native_pigment_partition(
    q_abs::FT,
    coeffs::NamedTuple{Names,NTuple{N,FT}},
) where {FT<:Real,Names,N}
    scale = maximum(prospect_primal, values(coeffs))
    if iszero(scale)
        return zero(FT), zero(FT), zero(FT), zero(FT)
    end
    scaled = map(k -> k / scale, coeffs)
    denom = sum(values(scaled))
    f_chl = scaled.chlorophyll / denom
    f_car = scaled.carotenoids / denom
    return f_chl, f_car, q_abs * f_chl, q_abs * f_car
end

function _resample_absorption_shares!(
    chl_out::Vector{OutFT},
    car_out::Vector{OutFT},
    α_nat::Vector,
    α_chl_nat::Vector,
    α_car_nat::Vector,
    f_chl_nat::Vector,
    f_car_nat::Vector,
    src_edges::Vector,
    dst_edges::AbstractVector,
) where {OutFT}
    n_dst = length(dst_edges) - 1
    n_src = length(α_nat)
    i0_de = firstindex(dst_edges)
    AccFT = promote_type(OutFT, eltype(α_nat), eltype(src_edges), eltype(dst_edges), Float64)

    s = _first_overlapping_bin(src_edges, dst_edges[i0_de])
    @inbounds for d in 0:(n_dst - 1)
        ν_lo = dst_edges[i0_de + d]
        ν_hi = dst_edges[i0_de + d + 1]
        p_lo = prospect_primal(ν_lo)
        p_hi = prospect_primal(ν_hi)
        while s + 1 < n_src && prospect_primal(src_edges[s + 2]) <= p_lo
            s += 1
        end
        k = s
        acc_α = zero(AccFT)
        acc_α_chl = zero(AccFT)
        acc_α_car = zero(AccFT)
        acc_f_chl = zero(AccFT)
        acc_f_car = zero(AccFT)
        acc_dν = zero(AccFT)
        while k < n_src
            a = _clamp_overlap_lo(ν_lo, src_edges[k + 1], AccFT)
            b = _clamp_overlap_hi(ν_hi, src_edges[k + 2], AccFT)
            dν = b - a
            if prospect_primal(dν) > 0
                acc_α += AccFT(α_nat[k + 1]) * dν
                acc_α_chl += AccFT(α_chl_nat[k + 1]) * dν
                acc_α_car += AccFT(α_car_nat[k + 1]) * dν
                acc_f_chl += AccFT(f_chl_nat[k + 1]) * dν
                acc_f_car += AccFT(f_car_nat[k + 1]) * dν
                acc_dν += dν
            end
            prospect_primal(src_edges[k + 2]) >= p_hi && break
            k += 1
        end
        s = k > 0 && prospect_primal(src_edges[k + 1]) > p_hi ? k - 1 : k
        s = min(s, n_src - 1)
        if prospect_primal(acc_α) > 0
            chl_out[d + 1] = OutFT(acc_α_chl / acc_α)
            car_out[d + 1] = OutFT(acc_α_car / acc_α)
        else
            chl_out[d + 1] = OutFT(acc_f_chl / acc_dν)
            car_out[d + 1] = OutFT(acc_f_car / acc_dν)
        end
    end
    return chl_out, car_out
end
