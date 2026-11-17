using CanopyOptics, Test, LinearAlgebra

# Batch application of joint allocation, per-source operators and band signals on a
# backend. `to_device` moves host arrays; the host scalar/matrix path is the reference.
function check_emission_operator_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        tiny = nextfloat(zero(FT))
        agrees(a, b) = isnan(b) ? isnan(a) : isapprox(a, b; rtol = 8eps(FT), atol = tiny)

        # Joint allocation over a spectrum of photons, including invalid entries,
        # with shared fractions (Ref) and with per-bin fractions (array).
        photons = FT[0, 1, 1000, tiny, floatmax(FT), -1, NaN, Inf, 2, 3]
        shared = (psii = FT(0.5), psi = FT(0.4), other = FT(0.05))
        per_bin = [
            (psii = FT(a), psi = FT(b), other = FT(c)) for (a, b, c) in (
                (0.5, 0.4, 0.05), (0.6, 0.4, 0), (1, 0, 0), (0, 0, 0), (0.3, 0.3, 0.3),
                (0.5, 0.4, 0.05), (0.5, 0.4, 0.05), (0.5, 0.4, 0.05), (0.7, 0.4, 0), (-0.1, 0.5, 0),
            )
        ]
        for fractions in (Ref(shared), per_bin)
            host = allocate_excitation.(photons, fractions)
            for order in (1:length(photons), reverse(1:length(photons)), 1:0)
                f = fractions isa Ref ? fractions : to_device(per_bin[order])
                device = allocate_excitation.(to_device(photons[order]), f)
                synchronize()
                result = Array(device)
                @test eltype(result) === eltype(host)
                @test length(result) == length(order)
                for (r, h) in zip(result, host[order])
                    @test keys(r.sinks) == keys(h.sinks)
                    @test agrees(r.unallocated, h.unallocated)
                    @test all(agrees(x, y) for (x, y) in zip(values(r.sinks), values(h.sinks)))
                end
            end
        end
        empty = Array(allocate_excitation.(to_device(FT[1, 2]), Ref(NamedTuple())))
        @test [x.unallocated for x in empty] == FT[1, 2]

        # Two sources with different probabilities and fates on the device.
        emission = FT[660, 680, 700, 720, 740, 760]
        excitation = FT[450, 650, 700, 770]
        spectra = (
            psii = TabulatedPhotonSpectrum(emission, FT[2, 8, 5, 3, 2, 1]),
            psi = TabulatedPhotonSpectrum(emission, FT[1, 2, 4, 6, 5, 3]),
        )
        probabilities = map(spectra) do spectrum
            reduce(hcat, [ex < last(emission) ? emission_probabilities(spectrum, StrictStokes(), ex) :
                          zeros(FT, length(emission)) for ex in excitation])
        end
        escape = (
            psii = [RegularPlateEscape(b, FT(0.6)) for b in FT[3, 2, 1, 0.5, 0.1, 0]],
            psi = [RegularPlateEscape(b, FT(0.8)) for b in FT[2, 1, 0.6, 0.3, 0.05, 0]],
        )
        birth = (
            psii = [ExponentialBirthProfile(a, one(FT), FT(0.2)) for a in FT[4, 2, 1, 0]],
            psi = [ExponentialBirthProfile(a, FT(0.3), one(FT)) for a in FT[3, 1, 0.5, 0]],
        )
        fates = map((e, b) -> leaf_emission_fates.(reshape(e, :, 1), reshape(b, 1, :)), escape, birth)
        host_ops = map(emission_operators, probabilities, fates)
        # Production is zero in the dark (last) excitation column: the contract.
        produced = (psii = FT[3, 2, 1, 0], psi = FT[0.5, 1, 0.25, 0])
        response = FT[0, 1, 1, 0.5, 0, 0]
        faces = (front = :front, back = :back, reabsorbed = :reabsorbed)
        for name in keys(spectra)
            p = to_device(probabilities[name])
            q = to_device(produced[name])
            # Caller recipe: dark columns from probability column sums, no production there.
            dark = Array(iszero.(vec(sum(p; dims = 1))))
            @test dark == [false, false, false, true]
            @test all(iszero, produced[name][dark])
            # Fates as a device matrix of NamedTuples and as device matrices.
            device_fates = to_device(fates[name])
            split_fates = map(k -> to_device(getproperty.(fates[name], k)), faces)
            for f in (device_fates, split_fates)
                ops = emission_operators(p, f)
                residual = emission_closure_residual(ops, p)
                synchronize()
                for k in keys(ops)
                    @test Array(getproperty(ops, k)) == getproperty(host_ops[name], k)
                end
                @test all(abs.(Array(residual)) .<= 16eps(FT))
                @test iszero(Array(residual)[end])
                # Closure against all produced photons of this source.
                outcomes = sum(sum(Array(getproperty(ops, k) * q)) for k in keys(ops))
                @test outcomes ≈ sum(produced[name]) rtol = 32eps(FT)
                # A perturbed entry (exact device addition of a one-hot matrix)
                # moves only its column's residual; in the dark column, exactly.
                for (m, e) in ((2, 4), (5, 1))
                    bump = zeros(FT, size(probabilities[name]))
                    bump[m, e] = FT(0.25)
                    bad = (front = ops.front, back = ops.back .+ to_device(bump), reabsorbed = ops.reabsorbed)
                    moved = Array(emission_closure_residual(bad, p))
                    base = Array(residual)
                    others = setdiff(1:length(base), e)
                    @test moved[others] == base[others]
                    @test moved[e] ≈ base[e] + FT(0.25) atol = 16eps(FT)
                    dark[e] && @test moved[e] === FT(0.25)
                end
                signals = map(faces) do k
                    signal = band_signal(to_device(response), getproperty(ops, k), q)
                    expected = band_signal(response, getproperty(host_ops[name], k), produced[name])
                    @test signal isa FT
                    @test signal ≈ expected rtol = 16eps(FT)
                    signal
                end
                # Matched permutations of excitation columns (probabilities, fates and
                # produced photons) leave band signals unchanged; zero excitation
                # columns give an exactly zero signal.
                select(x, cols) = x isa NamedTuple ? map(y -> y[:, cols], x) : x[:, cols]
                for cols in ([3, 1, 4, 2], [4, 3, 2, 1], Int[])
                    permuted = emission_operators(p[:, cols], select(f, cols))
                    synchronize()
                    for k in keys(faces)
                        signal = band_signal(to_device(response), getproperty(permuted, k), q[cols])
                        @test signal isa FT
                        if isempty(cols)
                            @test signal === zero(FT)
                            @test isempty(Array(emission_closure_residual(permuted, p[:, cols])))
                        else
                            @test signal ≈ signals[k] rtol = 16eps(FT)
                        end
                    end
                end
            end
        end
    end
end
