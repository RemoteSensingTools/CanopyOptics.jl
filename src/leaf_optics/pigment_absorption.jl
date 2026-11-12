"""
    AbstractAbsorptionPartition

Physical rule for assigning spectral material absorption to constituents.
"""
abstract type AbstractAbsorptionPartition end

"""
    HomogeneousAbsorption()

Assign absorption in one wavelength bin in proportion to constituent absorption
coefficients. All constituents must have the same proportions along every optical
path. This does not describe spatially segregated pigments in layered mesophyll.
"""
struct HomogeneousAbsorption <: AbstractAbsorptionPartition end

"""
    AbsorptionPartitionStatus

Status code reporting the availability and compatibility of a constituent
absorption partition calculation.
"""
@enum AbsorptionPartitionStatus::UInt8 begin
    AbsorptionPartitionValid
    AbsorptionNoMaterial
    AbsorptionInvalidInput
    AbsorptionIncompatibleInput
end

@doc """
    AbsorptionPartitionValid

Valid constituent absorption partition (including zero illumination of absorbing material).
""" AbsorptionPartitionValid

@doc """
    AbsorptionNoMaterial

Zero material absorption coefficient paired with zero total absorbed photons.
""" AbsorptionNoMaterial

@doc """
    AbsorptionInvalidInput

Negative or nonfinite total absorbed photons or constituent absorption coefficients.
""" AbsorptionInvalidInput

@doc """
    AbsorptionIncompatibleInput

Positive total absorbed photons supplied with an empty or all-zero constituent set.
""" AbsorptionIncompatibleInput

"""
    AbsorptionPartitionResult(total_absorbed, constituents, fractions, status)

Spectral absorbed photons and their named constituent contributions, in the same
photon-flux units and area convention as the caller's input. `constituents` is a
concrete NamedTuple, so adding an absorber does not change a fixed field layout.
`fractions` retains unitless probabilities before photon multiplication, so
selection does not lose pigment identity when constituent photons underflow.
Fractions are undefined (NaN) for absent absorption coefficients; absorbed photons
are then zero. Invalid results contain NaN payloads; `absorption_valid` checks status.
This result owns no wavelength grid, device storage or physiological state.
"""
struct AbsorptionPartitionResult{FT,C<:NamedTuple}
    total_absorbed::FT
    constituents::C
    fractions::C
    status::AbsorptionPartitionStatus
end

"""
    absorption_valid(result::AbsorptionPartitionResult) -> Bool

Return `true` when `result` represents a valid absorption partition, including
material with zero absorption coefficient (`AbsorptionNoMaterial`).
"""
@inline absorption_valid(result::AbsorptionPartitionResult) =
    result.status == AbsorptionPartitionValid || result.status == AbsorptionNoMaterial

@inline function invalid_absorption_partition(absorbed, coefficients, status)
    unavailable = oftype(absorbed, NaN)
    return AbsorptionPartitionResult(
        unavailable,
        map(_ -> unavailable, coefficients),
        map(_ -> unavailable, coefficients),
        status,
    )
end

"""
    partition_absorption(HomogeneousAbsorption(), absorbed_photons, coefficients)

Partition the supplied total absorbed photons **at one wavelength/bin** using
named constituent concentration × absorption-coefficient products in a common
optical-depth convention. Inputs must share FT. No PROSPECT database or leaf
composition is chosen here; their assembly belongs to the host.

For each constituent, Q_i = Q_abs k_i / sum(k_i). Scale the nonnegative coefficient
products before summing to avoid overflow. This defines constituent probabilities;
it never normalizes or alters the supplied optical absorption budget. Obtain Q_abs
independently from the optical solver (e.g. `prospect_budget`).

Negative/nonfinite inputs return `AbsorptionInvalidInput`. An empty/all-zero
coefficient set with nonzero Q_abs returns `AbsorptionIncompatibleInput`; with
Q_abs=0 it returns `AbsorptionNoMaterial`. Zero illumination of absorbing material
is instead `AbsorptionPartitionValid`. All paths are scalar and nonthrowing for
same-FT numerical inputs. Mixed-precision tuples have no implicit conversion.
"""
@inline function partition_absorption(
    ::HomogeneousAbsorption,
    absorbed::FT,
    coefficients::NamedTuple{Names,<:NTuple{N,FT}},
) where {FT<:AbstractFloat,Names,N}
    valid(x) = isfinite(x) && x >= zero(x)
    if !valid(absorbed) || !all(valid, values(coefficients))
        return invalid_absorption_partition(absorbed, coefficients, AbsorptionInvalidInput)
    end
    scale = N == 0 ? zero(FT) : maximum(values(coefficients))
    if iszero(scale)
        iszero(absorbed) || return invalid_absorption_partition(
            absorbed,
            coefficients,
            AbsorptionIncompatibleInput,
        )
        return AbsorptionPartitionResult(
            absorbed,
            map(zero, coefficients),
            map(_ -> oftype(absorbed, NaN), coefficients),
            AbsorptionNoMaterial,
        )
    end
    scaled = map(k -> k / scale, coefficients)
    denominator = sum(values(scaled))
    fractions = map(k -> k / denominator, scaled)
    photons = map(f -> absorbed * f, fractions)
    return AbsorptionPartitionResult(absorbed, photons, fractions, AbsorptionPartitionValid)
end
