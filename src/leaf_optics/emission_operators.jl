const EMISSION_FATE_NAMES = (:front, :back, :reabsorbed)

"""
    emission_operators(probabilities, fates) -> (front, back, reabsorbed)

Multiply conditional emission probabilities by optical escape fates on one
(emission, excitation) grid, giving the operators
`M_F = p .* F`, `M_B = p .* B` and `M_A = p .* A` of one source.
`probabilities[m, e]` holds, column by column, the output of
[`emission_probabilities`](@ref) (any spectrum or support rule) or an all-zero
dark column. Rows are emission bins and columns excitation bins.

Dark-column contract: a zero (dark) column is permitted only in excitation bins
where this source produces no photons. Photons produced in a dark column have no
outcome and are silently lost; [`emission_closure_residual`](@ref) cannot detect
this because it does not see production. For an excitation bin with production
and empty emission support, expand the emission grid or treat it as an error.
Check `all(iszero, produced[iszero.(vec(sum(probabilities; dims = 1)))])` before
applying the operators.

`fates` is either a NamedTuple `(front = F, back = B, reabsorbed = A)` of three
matrices (keys in any order), or a matrix of `(front, back, reabsorbed)` NamedTuples
as returned by broadcasting [`leaf_emission_fates`](@ref) over an emission-escape
column and an excitation-birth row. All arrays have the size and floating-point
element type of `probabilities`; mismatches throw on the host. Arrays may live on
any backend that supports broadcasting; no clipping, renormalization or yield is
applied. Call once per source: sources may differ in probabilities and in fates.
Multiplying an operator by already-produced photons per excitation bin gives
hemispheric face totals in photons; see [`emission_closure_residual`](@ref).
"""
function emission_operators(
    probabilities::AbstractMatrix{FT},
    fates::NamedTuple{Names,<:Tuple{AbstractMatrix,AbstractMatrix,AbstractMatrix}},
) where {FT<:AbstractFloat,Names}
    all(name -> name in Names, EMISSION_FATE_NAMES) ||
        throw(ArgumentError("fates must be named front, back and reabsorbed"))
    F, B, A = fates.front, fates.back, fates.reabsorbed
    all(x -> eltype(x) === FT, (F, B, A)) ||
        throw(ArgumentError("fates must have the floating-point type of the probabilities"))
    all(x -> size(x) == size(probabilities), (F, B, A)) ||
        throw(DimensionMismatch("fates and probabilities must share the (emission, excitation) grid"))
    return (front = probabilities .* F, back = probabilities .* B, reabsorbed = probabilities .* A)
end

function emission_operators(
    probabilities::AbstractMatrix{FT},
    fates::AbstractMatrix{NamedTuple{EMISSION_FATE_NAMES,NTuple{3,FT}}},
) where {FT<:AbstractFloat}
    size(fates) == size(probabilities) ||
        throw(DimensionMismatch("fates and probabilities must share the (emission, excitation) grid"))
    return (
        front = broadcast((p, x) -> p * x.front, probabilities, fates),
        back = broadcast((p, x) -> p * x.back, probabilities, fates),
        reabsorbed = broadcast((p, x) -> p * x.reabsorbed, probabilities, fates),
    )
end

emission_operators(::AbstractMatrix{<:AbstractFloat}, ::Any) = throw(ArgumentError(
    "fates must be three matrices named front, back and reabsorbed, or a matrix of " *
    "(front, back, reabsorbed) NamedTuples, with the probabilities' floating-point type",
))
emission_operators(::Any, ::Any) =
    throw(ArgumentError("emission probabilities must be a floating-point matrix"))

"""
    emission_closure_residual(operators, probabilities) -> vector

Photon-closure diagnostic per excitation column:
`r[e] = Σ_m (M_F + M_B + M_A - p)[m, e]` for `operators` from
[`emission_operators`](@ref). Differences are formed per element before the
emission sum, so `r[e]` is a roundoff-level residual: the operators of a column sum
to the column sum of `p` (one for a conditional spectrum, zero for a dark column).
The caller chooses a tolerance; nothing is clipped or renormalized.

`r` checks the fates weighted by `p`, not production: a dark column has `r[e] = 0`
even if photons are produced there, which would be lost. Dark columns are allowed
only where production is zero (see [`emission_operators`](@ref)).
"""
function emission_closure_residual(
    operators::NamedTuple{Names,<:Tuple{AbstractMatrix,AbstractMatrix,AbstractMatrix}},
    probabilities::AbstractMatrix{FT},
) where {FT<:AbstractFloat,Names}
    all(name -> name in Names, EMISSION_FATE_NAMES) ||
        throw(ArgumentError("operators must be named front, back and reabsorbed"))
    F, B, A = operators.front, operators.back, operators.reabsorbed
    all(x -> eltype(x) === FT, (F, B, A)) ||
        throw(ArgumentError("operators must have the floating-point type of the probabilities"))
    all(x -> size(x) == size(probabilities), (F, B, A)) ||
        throw(DimensionMismatch("operators and probabilities must share the (emission, excitation) grid"))
    return vec(sum(F .+ B .+ A .- probabilities; dims = 1))
end

"""
    band_signal(response, operator, produced_photons) -> FT

Linear band observation `response' * operator * produced_photons`, in photons for
photon-basis weights.
`operator[m, e]` is one face operator of one source (from
[`emission_operators`](@ref)); `produced_photons[e]` holds photons already produced
by that source in each excitation bin (absorption, allocation and yield applied
once, by the caller). Production must be zero in dark (all-zero) operator columns;
photons produced there are lost (see [`emission_operators`](@ref)).

`response[m]` are dimensionless photon-basis weights on the emission grid, for
example a detector quantum efficiency times a filter transmission. Emission-bin
measures are already inside the probabilities and are not applied again. A power
responsivity (A/W) is not a photon weight: if the caller folds the photon energy
`hc/λ_m` of the emission wavelength into the weights, the result is in energy
units. That conversion is the caller's and always uses emission wavelengths.

This is a linear functional of hemispheric face totals. It knows nothing of
instrument geometry or fluorometers and gives no directional radiance.
For several sources, call once per source and sum.
All three arrays share one floating-point type and a matching grid; mismatches
throw. Device arrays are supported by their matrix-vector product and `dot`.
"""
function band_signal(
    response::AbstractVector{FT},
    operator::AbstractMatrix{FT},
    produced_photons::AbstractVector{FT},
) where {FT<:AbstractFloat}
    size(operator) == (length(response), length(produced_photons)) || throw(DimensionMismatch(
        "operator size $(size(operator)) does not match response length $(length(response)) " *
        "and produced-photon length $(length(produced_photons))",
    ))
    return dot(response, operator * produced_photons)
end

band_signal(::AbstractVector, ::AbstractMatrix, ::AbstractVector) = throw(ArgumentError(
    "response, operator and produced photons must share one floating-point element type",
))
