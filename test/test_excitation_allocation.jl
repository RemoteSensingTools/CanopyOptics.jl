# Joint allocation of photons to any number of named sinks. Payload status follows
# photosystem_excitation: invalid numerical input gives typed NaN everywhere.

# Function barrier: measure allocations of a compiled call (Julia 1.10 counts
# boxed values when @allocated runs at @testset scope).
function excitation_allocation_bytes(photons, fractions)
    allocate_excitation(photons, fractions)
    return @allocated allocate_excitation(photons, fractions)
end

all_nan(result) = isnan(result.unallocated) && all(isnan, values(result.sinks))

@testset "Joint excitation allocation to named sinks" begin
    for FT in (Float32, Float64)
        other = FT === Float32 ? Float64 : Float32
        photons = FT.([0, 0.25, 1, 1000])
        push!(photons, floatmax(FT), nextfloat(zero(FT)))

        # No sinks: everything stays unallocated.
        for q in photons
            empty = @inferred allocate_excitation(q, NamedTuple())
            @test empty.sinks === NamedTuple()
            @test empty.unallocated === q
        end
        for bad in FT.([-1, NaN, Inf, -Inf])
            @test isnan(allocate_excitation(bad, NamedTuple()).unallocated)
        end

        # One sink reproduces photosystem_excitation exactly, including NaN status.
        for q in photons, f in FT.([0, 0.3, 0.5, 1])
            one_sink = @inferred allocate_excitation(q, (psii = f,))
            @test keys(one_sink.sinks) == (:psii,)
            @test one_sink.sinks.psii === photosystem_excitation(q, f)
            @test one_sink.unallocated === (one(FT) - f) * q
            @test zero(FT) <= one_sink.unallocated <= q
        end
        for f in (FT(NaN), FT(Inf), -FT(Inf), -one(FT), -nextfloat(zero(FT)), nextfloat(one(FT)))
            @test all_nan(allocate_excitation(one(FT), (psii = f,)))
            @test isnan(photosystem_excitation(one(FT), f))
        end
        for q in FT.([-1, NaN, Inf, -Inf])
            @test all_nan(allocate_excitation(q, (psii = FT(0.5),)))
            @test isnan(photosystem_excitation(q, FT(0.5)))
        end

        # Three sinks, closure, bounds and exact per-sink products.
        fractions = (psii = FT(0.45), psi = FT(0.35), other = FT(0.1))
        for q in photons
            r = @inferred allocate_excitation(q, fractions)
            @test isbitstype(typeof(r))
            @test keys(r.sinks) == keys(fractions)
            @test all(x -> x isa FT, values(r.sinks))
            @test r.unallocated isa FT
            for name in keys(fractions)
                @test getproperty(r.sinks, name) === getproperty(fractions, name) * q
            end
            @test zero(FT) <= r.unallocated <= q
            # Closure in BigFloat (the FT sum can overflow at floatmax); subnormal
            # photons have only an absolute, few-ULP resolution.
            tiny = 4nextfloat(zero(FT))
            closure = sum(BigFloat.(values(r.sinks))) + BigFloat(r.unallocated)
            @test closure ≈ BigFloat(q) rtol = 4eps(FT) atol = tiny
            reference = (1 - sum(BigFloat.(values(fractions)))) * BigFloat(q)
            @test r.unallocated ≈ FT(reference) rtol = 8eps(FT) atol = tiny
        end

        # Reordered keys give reordered results with identical per-name values.
        reordered = (other = fractions.other, psii = fractions.psii, psi = fractions.psi)
        r1 = allocate_excitation(FT(7), fractions)
        r2 = allocate_excitation(FT(7), reordered)
        @test keys(r2.sinks) == (:other, :psii, :psi)
        @test all(name -> getproperty(r1.sinks, name) === getproperty(r2.sinks, name), keys(fractions))
        @test r1.unallocated ≈ r2.unallocated rtol = 2eps(FT)

        # Fractions summing exactly to one leave nothing unallocated.
        complete = (a = FT(0.5), b = FT(0.25), c = FT(0.25))
        r = allocate_excitation(FT(3), complete)
        @test r.unallocated === zero(FT)
        @test sum(values(r.sinks)) === FT(3)

        # Normalized fractions whose rounded sum lies slightly above or below one
        # are accepted; unallocated never becomes negative.
        for n in 2:9
            weights = ntuple(i -> FT(mod(7i^2 + 3i, 11) + 1), n)
            shares = weights ./ sum(weights)
            named = NamedTuple{ntuple(i -> Symbol(:sink_, i), n)}(shares)
            for q in FT.([1, 3, 1000])
                r = allocate_excitation(q, named)
                @test !isnan(r.unallocated)
                @test zero(FT) <= r.unallocated <= q * n * eps(FT)
                @test sum(values(r.sinks)) + r.unallocated ≈ q rtol = (n + 1) * eps(FT)
            end
        end
        # The documented allowance: a sum of 1 + N*eps(FT) passes, beyond it fails.
        @test allocate_excitation(one(FT), (a = one(FT), b = eps(FT))).unallocated === zero(FT)
        @test allocate_excitation(one(FT), (a = one(FT), b = 2eps(FT))).unallocated === zero(FT)
        @test all_nan(allocate_excitation(one(FT), (a = one(FT), b = 3eps(FT))))

        # Sums above one, negative and nonfinite fractions, invalid photons.
        for bad in (
            (a = FT(0.6), b = FT(0.5)),
            (a = one(FT), b = one(FT), c = zero(FT)),
            (a = FT(-0.1), b = FT(0.5)),
            (a = FT(NaN), b = FT(0.5)),
            (a = FT(Inf), b = FT(0.5)),
            (a = FT(0.5), b = -FT(Inf)),
            (a = floatmax(FT), b = zero(FT)),
        )
            r = @inferred allocate_excitation(one(FT), bad)
            @test all_nan(r)
            @test keys(r.sinks) == keys(bad)
            @test all(x -> x isa FT, values(r.sinks))
        end
        for q in FT.([-1, NaN, Inf, -Inf])
            @test all_nan(@inferred allocate_excitation(q, fractions))
        end

        # No implicit promotion: mixed precision or non-float values are rejected.
        @test_throws MethodError allocate_excitation(one(FT), (a = FT(0.5), b = other(0.5)))
        @test_throws MethodError allocate_excitation(one(other), fractions)
        @test_throws MethodError allocate_excitation(one(FT), (a = 1, b = 0))
        @test_throws MethodError allocate_excitation(1, (a = FT(0.5),))

        # Allocation free for zero, one, three and many sinks, valid or not.
        many = NamedTuple{ntuple(i -> Symbol(:s, i), 40)}(ntuple(_ -> FT(1 / 41), 40))
        r = @inferred allocate_excitation(FT(41), many)
        @test abs(r.unallocated - one(FT)) <= FT(41) * 40eps(FT)
        @test sum(values(r.sinks)) + r.unallocated ≈ FT(41) rtol = 41eps(FT)
        for f in (NamedTuple(), (psii = FT(0.5),), fractions, many, (a = FT(0.6), b = FT(0.5)))
            @test excitation_allocation_bytes(FT(2), f) == 0
        end
    end
end

@testset "Allocation conventions: antenna fractions f versus absorbed fractions a = h*f" begin
    for FT in (Float32, Float64)
        model = HomogeneousAbsorption()
        selector = CarotenoidAntenna(FT(0.6))
        # Small synthetic spectrum: absorbed photons and constituent products per bin.
        absorbed = FT.([120, 300, 450, 80, 10])
        coefficients = [
            (chlorophyll = FT(c), carotenoids = FT(a), water = FT(w)) for (c, a, w) in (
                (0.2, 0.6, 0.05),
                (0.5, 0.4, 0.02),
                (0.9, 0.1, 0.01),
                (0.7, 0.0, 0.05),
                (0.05, 0.0, 0.8),
            )
        ]
        # Fractions of antenna excitation; three named sinks.
        f = (psii = FT(0.48), psi = FT(0.42), other = FT(0.05))
        for (q, k) in zip(absorbed, coefficients)
            partition = partition_absorption(model, q, k)
            antenna = antenna_excitation(selector, partition)
            h = antenna / q
            a = map(x -> h * x, f)
            antenna_route = allocate_excitation(antenna, f)
            absorbed_route = allocate_excitation(q, a)
            for name in keys(f)
                @test getproperty(absorbed_route.sinks, name) ≈
                      getproperty(antenna_route.sinks, name) rtol = 8eps(FT)
            end
            # The absorbed-photon remainder also contains the non-antenna photons.
            @test absorbed_route.unallocated ≈ (q - antenna) + antenna_route.unallocated rtol = 16eps(FT)
            # Pitfalls: a after antenna selection counts the transfer loss twice;
            # f on total absorption omits it.
            double = allocate_excitation(antenna, a)
            omitted = allocate_excitation(q, f)
            @test double.sinks.psii ≈ h * antenna_route.sinks.psii rtol = 8eps(FT)
            @test omitted.sinks.psii ≈ antenna_route.sinks.psii / h rtol = 8eps(FT)
            h < one(FT) - sqrt(eps(FT)) && @test !(double.sinks.psii ≈ antenna_route.sinks.psii)
        end
    end
end
