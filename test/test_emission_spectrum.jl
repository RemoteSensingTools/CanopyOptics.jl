# A separate, bounded Stokes band demonstrates extension without editing the
# tabulated source, optical escape, or physiology implementations.
struct FiniteStokesBand{FT} <: AbstractEmissionSupport
    minimum_shift_nm::FT
    maximum_emission_nm::FT
end
CanopyOptics.emission_allowed(p::FiniteStokesBand{FT}, em::FT, ex::FT) where {FT} =
    em - ex >= p.minimum_shift_nm && em <= p.maximum_emission_nm

struct InvalidSupport <: AbstractEmissionSupport end
CanopyOptics.emission_allowed(::InvalidSupport, em, ex) = 1

@testset "Explicit conditional photon spectra" begin
    for FT in (Float32, Float64)
        wavelengths = FT[600, 650, 700, 750]
        masses = FT[1, 2, 3, 4]
        spectrum = TabulatedPhotonSpectrum(wavelengths, masses)
        @test spectrum isa TabulatedPhotonSpectrum{FT}
        wavelengths[1] = FT(599)
        masses[1] = FT(99)
        @test spectrum.wavelength_nm[1] == FT(600)
        @test spectrum.relative_photon_mass[1] == one(FT)
        for excitation in FT[400, 600, 625, 650, 700, 749]
            result = @inferred emission_probabilities(spectrum, StrictStokes(), excitation)
            @test result isa Vector{FT}
            @test length(result) == length(spectrum.wavelength_nm)
            reference = BigFloat.(spectrum.relative_photon_mass)
            reference[spectrum.wavelength_nm .<= excitation] .= 0
            reference ./= sum(reference)
            @test result ≈ FT.(reference) rtol=4eps(FT)
            @test all(iszero, result[spectrum.wavelength_nm .<= excitation])
            @test sum(result) ≈ one(FT) rtol=4eps(FT)
            # h*c cancels: emitted photon energy cannot exceed excitation energy.
            @test sum(result ./ spectrum.wavelength_nm) <= inv(excitation)*(one(FT)+4eps(FT))
        end
        alternative = emission_probabilities(spectrum, FiniteStokesBand(FT(25), FT(700)), FT(625))
        @test alternative ≈ FT[0, 2/5, 3/5, 0]
        @test_throws ArgumentError emission_probabilities(spectrum, InvalidSupport(), FT(500))
        @test_throws ArgumentError emission_probabilities(spectrum, StrictStokes(), FT(750))
        for invalid in FT[0, -1, NaN, Inf]
            @test_throws ArgumentError emission_probabilities(spectrum, StrictStokes(), invalid)
        end
        other = FT === Float32 ? Float64 : Float32
        @test_throws ArgumentError emission_probabilities(spectrum, StrictStokes(), other(500))
        @test_throws ArgumentError emission_probabilities(spectrum, StrictStokes(), 500)
        @test_throws ArgumentError TabulatedPhotonSpectrum(FT[600, 700], other[1, 2])
        @test_throws DimensionMismatch TabulatedPhotonSpectrum(FT[600, 700], FT[1])
        for wave in (FT[], FT[0], FT[-1], FT[Inf], FT[NaN], FT[700, 600], FT[600, 600])
            @test_throws ArgumentError TabulatedPhotonSpectrum(wave, ones(FT, length(wave)))
        end
        for mass in (FT[0, 0], FT[-1, 1], FT[NaN, 1], FT[Inf, 1])
            @test_throws ArgumentError TabulatedPhotonSpectrum(FT[600, 700], mass)
        end
        # Retain raw tiny support before normalizing; global normalization would
        # erase it. Even the smallest subnormal is a valid relative mass.
        tiny = nextfloat(zero(FT))
        extreme = TabulatedPhotonSpectrum(FT[600, 700, 750], FT[floatmax(FT), tiny, 2tiny])
        @test emission_probabilities(extreme, StrictStokes(), FT(650)) ≈ FT[0, 1/3, 2/3]
        huge = TabulatedPhotonSpectrum(FT[600, 700], fill(floatmax(FT), 2))
        @test emission_probabilities(huge, StrictStokes(), FT(500)) == fill(FT(1/2), 2)
        # Unequal spectral measures are applied by the caller, once.
        density = FT[1, 2, 1]
        measure = FT[1, 3, 2]
        integrated = TabulatedPhotonSpectrum(FT[650, 700, 750], density .* measure)
        @test emission_probabilities(integrated, StrictStokes(), FT(600)) ≈ FT[1/9, 6/9, 2/9]
        # Public mutable arrays are revalidated, including changes in length.
        for mutate! in (
            s -> (s.wavelength_nm[2] = s.wavelength_nm[1]),
            s -> (s.relative_photon_mass[1] = FT(-1)),
            s -> fill!(s.relative_photon_mass, zero(FT)),
            s -> pop!(s.relative_photon_mass),
        )
            s = TabulatedPhotonSpectrum(FT[600, 700], FT[1, 2])
            mutate!(s)
            @test_throws Union{ArgumentError,DimensionMismatch} emission_probabilities(s, StrictStokes(), FT(500))
        end
        for count in (1, 7, 2101)
            wave = FT.(collect(600:(599+count)))
            raw = FT.(exp2.(range(-100, 100; length=count == 1 ? 2 : count)))[1:count]
            s = TabulatedPhotonSpectrum(wave, raw)
            p = emission_probabilities(s, StrictStokes(), FT(500))
            ref = BigFloat.(raw); ref ./= sum(ref)
            @test p ≈ FT.(ref) rtol=8eps(FT)
            @test sum(p) ≈ one(FT) rtol=8eps(FT)
        end
    end
end
