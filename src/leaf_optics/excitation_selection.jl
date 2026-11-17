"""
    AbstractAntennaSelection

Rule for transferring absorbed constituent photons into an antenna excitation pool.
"""
abstract type AbstractAntennaSelection end

"""
    ChlorophyllOnly()

Transfer only chlorophyll absorption into the antenna pool (requires `:chlorophyll`).
"""
struct ChlorophyllOnly <: AbstractAntennaSelection end

"""
    CarotenoidAntenna(transfer_fraction)

Chlorophyll absorption plus an explicit fraction of carotenoid-absorbed photons
enters the antenna pool. Requires `:chlorophyll` and `:carotenoids` contributions.
The finite fraction in [0,1] is a photon-transfer probability, not an energy
efficiency, PSII allocation or physiological nonphotochemical-quenching rate.
"""
struct CarotenoidAntenna{FT<:AbstractFloat} <: AbstractAntennaSelection
    transfer_fraction::FT
    function CarotenoidAntenna{FT}(fraction::FT) where {FT<:AbstractFloat}
        isfinite(fraction) && zero(FT) <= fraction <= one(FT) ||
            throw(ArgumentError("carotenoid transfer fraction must be finite and in [0,1]"))
        new{FT}(fraction)
    end
end
CarotenoidAntenna(fraction::FT) where {FT<:AbstractFloat} = CarotenoidAntenna{FT}(fraction)

"""
    antenna_excitation(selection, partition)

Selected antenna photons in the partition's photon/area units, before photosystem
allocation. Invalid partition status, a missing required pigment name or an invalid
selection returns a typed NaN. Unselected absorption remains in `partition`; it
is not removed from material energy accounting or assigned to a physiological KN.

Custom `AbstractAntennaSelection` types implement `select_absorbed_photons` for a
valid partition, returning a scalar of the same FT. Neither a fixed absorber list
nor a wavelength database is imposed on custom selections.
"""
@inline function antenna_excitation(
    selection::AbstractAntennaSelection,
    partition::AbsorptionPartitionResult,
)
    absorption_valid(partition) || return oftype(partition.total_absorbed, NaN)
    selected = select_absorbed_photons(selection, partition)
    isfinite(selected) && zero(selected) <= selected <= partition.total_absorbed ||
        return oftype(partition.total_absorbed, NaN)
    return selected
end

"""
    select_absorbed_photons(selection, partition)

Extension hook for an antenna selection applied to a valid absorption partition.
Implementations return selected photons with the input FT, or typed NaN for missing
required constituents. Normally call `antenna_excitation`, which checks status and
the returned photon bound. Built-in selections accumulate the smaller selected or
unselected probability before multiplying by photons, using the complementary
probability identity near total selection. This retains zero-pigment and all-selected
limits even if constituent photon amounts have underflowed.
This equivalent probability identity improves rounding without clipping, changing
constituent diagnostics or repairing an optical flux budget.
"""
function select_absorbed_photons end

@inline function complementary_excitation(partition, selected, unselected)
    total = partition.total_absorbed
    return selected <= unselected ? total * selected : total - total * unselected
end

@inline function select_absorbed_photons(
    ::ChlorophyllOnly,
    partition::AbsorptionPartitionResult{FT},
) where {FT}
    fractions = partition.fractions
    haskey(fractions, :chlorophyll) || return FT(NaN)
    partition.status == AbsorptionNoMaterial && return zero(FT)
    selected = fractions.chlorophyll
    unselected = sum(map(keys(fractions), values(fractions)) do name, value
        name === :chlorophyll ? zero(value) : value
    end)
    return complementary_excitation(partition, selected, unselected)
end

@inline function select_absorbed_photons(
    selection::CarotenoidAntenna{FT},
    partition::AbsorptionPartitionResult{FT},
) where {FT}
    fractions = partition.fractions
    haskey(fractions, :chlorophyll) && haskey(fractions, :carotenoids) || return FT(NaN)
    partition.status == AbsorptionNoMaterial && return zero(FT)
    fraction = selection.transfer_fraction
    selected = fractions.chlorophyll + fraction * fractions.carotenoids
    unselected = sum(
        map(keys(fractions), values(fractions)) do name, value
            name === :chlorophyll ? zero(value) :
            (name === :carotenoids ? (one(fraction)-fraction)*value : value)
        end,
    )
    return complementary_excitation(partition, selected, unselected)
end

"""
    photosystem_excitation(antenna_photons, allocation_fraction)

Allocate antenna photons to one declared photosystem, with an explicit finite
fraction in [0,1]. Values refer to the same wavelength bin and photon/area units;
the caller chooses and integrates the spectral domain. A PSII allocation is applied
once here, not again in a leaf model receiving the resulting PSII excitation.
Invalid numerical inputs return typed NaN; no empirical fraction is supplied.
"""
@inline function photosystem_excitation(antenna::FT, fraction::FT) where {FT<:AbstractFloat}
    isfinite(antenna) &&
    antenna >= zero(FT) &&
    isfinite(fraction) &&
    zero(FT) <= fraction <= one(FT) || return FT(NaN)
    return fraction * antenna
end

"""
    allocate_excitation(photons, fractions) -> (sinks, unallocated)

Allocate `photons` jointly to any number of named sinks with explicit fractions,
returning `(sinks = NamedTuple, unallocated = FT)`. `sinks` has the keys and order
of `fractions`; each sink receives `fraction * photons`. Zero, one, two or more
sinks are allowed. Values refer to one wavelength bin and photon/area unit.

The function does not know what `photons` are; the caller's fractions state the
convention. Fractions `f` of antenna excitation apply to `antenna_excitation`
output. Fractions `a = h*f` of leaf-absorbed photons, with `h` the antenna share of
absorbed photons, apply to total absorption. Applying `a` after antenna selection
counts the transfer loss twice; applying `f` to total absorption omits it.

`photons` and every fraction share one floating-point type; mixed precision is a
`MethodError`, as for [`photosystem_excitation`](@ref). Valid input has finite
`photons ≥ 0`, finite fractions in [0,1] and a fraction sum at most
`1 + N*eps(FT)` for `N` sinks (left-to-right summation of `N` rounded fractions).
Otherwise every sink and `unallocated` is typed NaN. With no sinks, valid photons
are returned unallocated. `unallocated = max(1 - Σf, 0) * photons`: never negative
and never above `photons`; within the allowance a sum above one leaves zero.
A single sink reproduces `photosystem_excitation(photons, fraction)` exactly.
The scalar call is type stable and allocation free for any `N`; broadcast it over
CPU or GPU arrays. No empirical fraction, yield or rate constant is supplied.
"""
@inline function allocate_excitation(
    photons::FT,
    fractions::NamedTuple{Names,<:NTuple{N,FT}},
) where {FT<:AbstractFloat,Names,N}
    # Indexed left-to-right loop and Val(N) tuples: type stable and allocation free
    # for any N (Base `map`/`values` over tuples longer than 32 allocate).
    valid = isfinite(photons) && photons >= zero(FT)
    total = zero(FT)
    for i in 1:N
        f = fractions[i]
        valid &= isfinite(f) && zero(FT) <= f <= one(FT)
        total += f
    end
    valid &= total <= one(FT) + N * eps(FT)
    if !valid
        unavailable = FT(NaN)
        return (
            sinks = NamedTuple{Names}(ntuple(_ -> unavailable, Val(N))),
            unallocated = unavailable,
        )
    end
    sinks = NamedTuple{Names}(ntuple(i -> fractions[i] * photons, Val(N)))
    return (sinks = sinks, unallocated = max(one(FT) - total, zero(FT)) * photons)
end
