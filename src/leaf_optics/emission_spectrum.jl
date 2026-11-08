"""
    TabulatedPhotonSpectrum(wavelength_nm, relative_photon_mass)

Host-prepared discrete photon measure on finite, positive, strictly increasing
wavelengths in nm. Both vectors must have the same floating-point element type.
Masses are finite and nonnegative with at least one positive entry; they need not
sum to one. For a tabulated photon *density*, the caller must first multiply by
the corresponding quadrature measures exactly once.

The constructor copies its inputs and retains raw masses so a rare supported tail
is not lost by premature normalization. The vectors remain mutable Julia arrays;
`emission_probabilities` revalidates them before preparation. This object contains
neither absorptance nor fluorescence yield, and is not a device kernel argument.
"""
struct TabulatedPhotonSpectrum{FT<:AbstractFloat}
    wavelength_nm::Vector{FT}
    relative_photon_mass::Vector{FT}

    function TabulatedPhotonSpectrum(wavelength_nm::AbstractVector{FT},
                                    relative_photon_mass::AbstractVector{FT}) where {FT<:AbstractFloat}
        wavelengths = collect(wavelength_nm)
        masses = collect(relative_photon_mass)
        validate_photon_spectrum(wavelengths, masses)
        new{FT}(wavelengths, masses)
    end
end

TabulatedPhotonSpectrum(::AbstractVector, ::AbstractVector) =
    throw(ArgumentError("Wavelengths and photon masses must have the same floating-point element type"))

function validate_photon_spectrum(wavelengths, masses)
    isempty(wavelengths) && throw(ArgumentError("A photon spectrum cannot be empty"))
    length(wavelengths) == length(masses) || throw(DimensionMismatch("Wavelength and photon-mass lengths differ"))
    all(x -> isfinite(x) && x > zero(x), wavelengths) ||
        throw(ArgumentError("Wavelengths must be finite and positive"))
    all(i -> wavelengths[i] > wavelengths[i-1], 2:length(wavelengths)) ||
        throw(ArgumentError("Wavelengths must be strictly increasing"))
    all(x -> isfinite(x) && x >= zero(x), masses) ||
        throw(ArgumentError("Photon masses must be finite and nonnegative"))
    any(x -> x > zero(x), masses) || throw(ArgumentError("A photon spectrum must have positive mass"))
    return nothing
end

"""Excitation-dependent emission support; extend [`emission_allowed`](@ref)."""
abstract type AbstractEmissionSupport end

"""
    StrictStokes()

Permit only emission wavelengths strictly longer than excitation. Equality is
excluded. This explicit approximation omits thermal anti-Stokes excitation.
Conditioning changes the reference spectral shape while retaining the separately
prescribed physiological photon yield; removed weights do not establish additional
quenching or heat production.
"""
struct StrictStokes <: AbstractEmissionSupport end

"""
    emission_allowed(support, emission_nm, excitation_nm) -> Bool

Scalar extension hook for an `AbstractEmissionSupport`. Host preparation supplies
finite positive wavelengths of the same floating-point type. Implementations
return a Boolean and can be device compatible without depending on a canopy,
physiology model, or spectrum storage layout.
"""
function emission_allowed end

@inline emission_allowed(::StrictStokes, emission_nm::FT, excitation_nm::FT) where {FT<:AbstractFloat} =
    emission_nm > excitation_nm

"""
    emission_probabilities(spectrum, support, excitation_nm)

Prepare conditional photon probabilities on the spectrum's **original** grid,
with exact zeros outside the explicit support. Excitation must be finite,
positive, and have the same floating-point type as the table. Empty retained
support is an error, not an inactive physical source; callers skip genuinely
inactive sources before requesting a conditional spectrum.

Select raw relative masses before max-scaled normalization, preserving a small
retained tail even when excluded weights are vastly larger. The result is a host
vector in the table's precision that can be transferred to any backend. No yield,
absorptance, energy conversion, or spectral-bin width is applied here. Convert
photon outcomes to energy using their emission wavelengths, not excitation.
"""
function emission_probabilities(spectrum::TabulatedPhotonSpectrum{FT},
                                support::AbstractEmissionSupport,
                                excitation_nm::FT) where {FT<:AbstractFloat}
    isfinite(excitation_nm) && excitation_nm > zero(FT) ||
        throw(ArgumentError("Excitation wavelength must be finite and positive"))
    validate_photon_spectrum(spectrum.wavelength_nm, spectrum.relative_photon_mass)
    probabilities = similar(spectrum.relative_photon_mass)
    for i in eachindex(probabilities)
        allowed = emission_allowed(support, spectrum.wavelength_nm[i], excitation_nm)
        allowed isa Bool || throw(ArgumentError("An emission support rule must return Bool"))
        probabilities[i] = allowed ? spectrum.relative_photon_mass[i] : zero(FT)
    end
    scale = maximum(probabilities)
    scale > zero(FT) || throw(ArgumentError("No positive photon mass remains within emission support"))
    probabilities ./= scale
    total = sum(probabilities)
    isfinite(total) && total > zero(FT) || throw(ArgumentError("Unrepresentable photon normalization"))
    probabilities ./= total
    return probabilities
end

emission_probabilities(::TabulatedPhotonSpectrum, ::AbstractEmissionSupport, ::Any) =
    throw(ArgumentError("Excitation wavelength must match the spectrum's floating-point type"))
