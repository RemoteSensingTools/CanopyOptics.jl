using CanopyOptics, Test, LinearAlgebra

function check_emission_spectrum_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        wavelength = FT[640, 660, 680, 700, 720, 740, 760]
        spectrum = TabulatedPhotonSpectrum(wavelength, FT[1, 3, 7, 9, 6, 3, 1])
        excitation = FT[450, 650, 700]
        conditional = hcat([emission_probabilities(spectrum, StrictStokes(), ex) for ex in excitation]...)
        escape = [RegularPlateEscape(b, FT(0.6)) for b in FT[3, 2, 1, 0.5, 0.1, 0.01, 0]]
        birth = [ExponentialBirthProfile(a, one(FT), FT(0.2)) for a in FT[4, 2, 1]]
        # An excitation distribution is not reapplied as absorptance or yield.
        # These are already-produced internal photons, including an inactive case.
        photons = FT[1 0 0 3; 0 0 2 4; 0 0 1 5]
        for reverse_face in (false, true)
            sources = reverse_face ? [ExponentialBirthProfile(p.excitation_depth, p.backward_fraction, p.forward_fraction) for p in birth] : birth
            fates = leaf_emission_fates.(reshape(escape, :, 1), reshape(sources, 1, :))
            operators = map(field -> conditional .* getproperty.(fates, field), (:front, :back, :reabsorbed))
            for order in (1:size(photons, 2), reverse(1:size(photons, 2)), 1:0)
                input = to_device(photons[:, order])
                output = map(op -> to_device(op)*input, operators)
                allowed = emission_allowed.(Ref(StrictStokes()), reshape(to_device(wavelength), :, 1), reshape(to_device(excitation), 1, :))
                synchronize()
                host = map(Array, output)
                @test Array(allowed) == (wavelength .> transpose(excitation))
                @test all(eltype(x) === FT for x in host)
                @test all(all(isfinite, x) && all(>=(zero(FT)), x) for x in host)
                for (j, original) in enumerate(order)
                    produced = sum(photons[:, original])
                    outcomes = [sum(x[:, j]) for x in host]
                    @test sum(outcomes) ≈ produced rtol=32eps(FT) atol=zero(FT)
                    emitted_energy = sum(sum(x[:, j] ./ wavelength) for x in host)
                    excitation_energy = sum(photons[:, original] ./ excitation)
                    @test emitted_energy <= excitation_energy*(one(FT)+32eps(FT))
                    for k in eachindex(operators), m in eachindex(wavelength)
                        # Independent high-precision application of the prepared
                        # matrix: accumulation only, separate from GPU BLAS.
                        reference = sum(BigFloat(operators[k][m,e])*BigFloat(photons[e,original]) for e in eachindex(excitation))
                        @test host[k][m,j] ≈ FT(reference) rtol=16eps(FT) atol=zero(FT)
                    end
                end
                @test all(size(x) == (length(wavelength), length(order)) for x in host)
            end
        end
    end
end
