# Independent species absorption from a mixture of optical path lengths.
# Integrate k_i exp(-sum(k)*s) by its convergent power series, then accumulate
# separate incoming path weights in BigFloat. No allocation-ratio formula is used.
function path_absorption_reference(coefficients, incident)
    setprecision(256) do
        k=map(BigFloat, coefficients)
        total=sum(values(k))
        lengths=BigFloat.((0.125, 0.75, 2.5))
        weights=BigFloat.((0.25, 0.5, 0.25))
        species=map(k) do coefficient
            flux=zero(BigFloat)
            for (length, weight) in zip(lengths, weights)
                term=length
                integral=term
                for n = 1:256
                    term *= -total*length/BigFloat(n+1)
                    integral += term
                end
                flux += BigFloat(incident)*weight*coefficient*integral
            end
            flux
        end
        absorbed=sum(
            BigFloat(incident)*weight*(-expm1(-total*length)) for
            (length, weight) in zip(lengths, weights)
        )
        (; absorbed, species)
    end
end

struct ExtraAntenna{FT} <: AbstractAntennaSelection
    fraction::FT
end
CanopyOptics.select_absorbed_photons(p::ExtraAntenna, r::AbsorptionPartitionResult) =
    p.fraction*r.constituents.extra

function partition_allocations(q, k, selection)
    partition=partition_absorption(HomogeneousAbsorption(), q, k)
    antenna_excitation(selection, partition)
    return @allocated antenna_excitation(
        selection,
        partition_absorption(HomogeneousAbsorption(), q, k),
    )
end

@testset "Homogeneous pigment partition and independent optical paths" begin
    for FT in (Float32, Float64)
        rtol=FT==Float32 ? FT(2e-6) : FT(2e-14)
        for k in (
                (
                    chlorophyll = FT(0.31),
                    carotenoids = FT(0.17),
                    water = FT(0.09),
                    extra = FT(0.23),
                ),
                (
                    chlorophyll = FT(0.8),
                    carotenoids = zero(FT),
                    water = zero(FT),
                    extra = zero(FT),
                ),
            ),
            incident in FT.([0, 0.01, 100, 5000])

            reference=path_absorption_reference(k, incident)
            q=FT(reference.absorbed)
            result=@inferred partition_absorption(HomogeneousAbsorption(), q, k)
            @test isbitstype(typeof(result))
            @test result.status==AbsorptionPartitionValid
            @test absorption_valid(result)
            @test result.total_absorbed===q
            @test isapprox(sum(values(result.fractions)), one(FT); rtol)
            @test all(x->x isa FT, values(result.fractions))
            @test all(x->x isa FT, values(result.constituents))
            @test isapprox(sum(values(result.constituents)), q; rtol)
            for name in keys(k)
                @test isapprox(
                    getproperty(result.constituents, name),
                    FT(getproperty(reference.species, name));
                    rtol,
                )
            end
            original=result.constituents
            chl=@inferred antenna_excitation(ChlorophyllOnly(), result)
            @test chl isa FT
            @test isapprox(chl, FT(reference.species.chlorophyll); rtol)
            for eta in FT.([0, 0.3, 1]), fraction in FT.([0, 0.5, 1])
                selector=CarotenoidAntenna(eta)
                antenna=@inferred antenna_excitation(selector, result)
                expected=reference.species.chlorophyll+BigFloat(eta)*reference.species.carotenoids
                @test antenna isa FT
                @test isapprox(antenna, FT(expected); rtol)
                @test zero(FT)<=antenna<=q
                psii=@inferred photosystem_excitation(antenna, fraction)
                @test psii isa FT
                @test zero(FT)<=psii<=antenna
                @test isapprox(psii, FT(BigFloat(fraction)*expected); rtol)
                @test result.constituents==original
            end
            custom=@inferred antenna_excitation(ExtraAntenna(FT(0.4)), result)
            @test custom≈FT(0.4)*result.constituents.extra
            @test partition_allocations(q, k, ChlorophyllOnly())==0
        end
    end
end

@testset "Pigment domains, zero material and photon-selection endpoints" begin
    for FT in (Float32, Float64)
        model=HomogeneousAbsorption()
        k=(chlorophyll = FT(0.7), carotenoids = FT(0.3), extra = zero(FT))
        zeros=map(zero, k)
        valid_result=partition_absorption(model, one(FT), k)
        for z in (zeros, NamedTuple())
            no_material=@inferred partition_absorption(model, zero(FT), z)
            @test no_material.status==AbsorptionNoMaterial
            @test absorption_valid(no_material)
            @test no_material.total_absorbed===zero(FT)
            @test all(iszero, values(no_material.constituents))
            @test all(isnan, values(no_material.fractions))
            impossible=@inferred partition_absorption(model, one(FT), z)
            @test impossible.status==AbsorptionIncompatibleInput
            @test !absorption_valid(impossible)
            @test isnan(impossible.total_absorbed)
            @test all(isnan, values(impossible.constituents))
            @test isnan(antenna_excitation(ChlorophyllOnly(), impossible))
        end
        # Pigment identity survives constituent-photon underflow.
        tiny=nextfloat(zero(FT))
        absent=(
            chlorophyll = zero(FT),
            carotenoids = zero(FT),
            extra_a = one(FT),
            extra_b = one(FT),
        )
        subnormal=partition_absorption(model, tiny, absent)
        @test all(iszero, values(subnormal.constituents))
        @test subnormal.fractions.extra_a===FT(0.5)
        @test subnormal.fractions.extra_b===FT(0.5)
        @test antenna_excitation(ChlorophyllOnly(), subnormal)===zero(FT)
        for eta in FT.([0, 0.5, 1])
            @test antenna_excitation(CarotenoidAntenna(eta), subnormal)===zero(FT)
        end
        balanced=partition_absorption(
            model,
            tiny,
            (chlorophyll = one(FT), carotenoids = one(FT)),
        )
        @test antenna_excitation(ChlorophyllOnly(), balanced)===tiny*FT(0.5)
        @test antenna_excitation(CarotenoidAntenna(one(FT)), balanced)===tiny
        @test partition_absorption(model, zero(FT), k).status==AbsorptionPartitionValid
        @test isnan(
            antenna_excitation(
                ChlorophyllOnly(),
                partition_absorption(model, zero(FT), NamedTuple()),
            ),
        )
        @test antenna_excitation(
            ChlorophyllOnly(),
            partition_absorption(model, FT(3), (chlorophyll = one(FT),)),
        )===FT(3)
        @test isnan(
            antenna_excitation(
                CarotenoidAntenna(one(FT)),
                partition_absorption(model, FT(3), (chlorophyll = one(FT),)),
            ),
        )
        @test isnan(
            antenna_excitation(
                ChlorophyllOnly(),
                partition_absorption(model, FT(3), (extra = one(FT),)),
            ),
        )
        for bad in FT.([-1, NaN, Inf, -Inf])
            invalid=@inferred partition_absorption(model, bad, k)
            @test invalid.status==AbsorptionInvalidInput
            @test isnan(antenna_excitation(ChlorophyllOnly(), invalid))
            for key in keys(k)
                invalid_k=merge(k, NamedTuple{(key,)}((bad,)))
                invalid=partition_absorption(model, one(FT), invalid_k)
                @test invalid.status==AbsorptionInvalidInput
                @test all(isnan, values(invalid.constituents))
            end
            @test_throws ArgumentError CarotenoidAntenna(bad)
            @test isnan(photosystem_excitation(bad, one(FT)))
            @test isnan(photosystem_excitation(one(FT), bad))
        end
        @test_throws ArgumentError CarotenoidAntenna(nextfloat(one(FT)))
        @test isnan(photosystem_excitation(one(FT), nextfloat(one(FT))))
        @test_throws MethodError partition_absorption(
            model,
            one(FT),
            (a = one(FT), b = one(FT==Float32 ? Float64 : Float32)),
        )
        # Finite constituent products whose sum overflows: the probabilities do not.
        huge=(chlorophyll = floatmax(FT), carotenoids = floatmax(FT), extra = zero(FT))
        result=partition_absorption(model, one(FT), huge)
        @test result.status==AbsorptionPartitionValid
        @test result.constituents.chlorophyll===FT(0.5)
        @test result.constituents.carotenoids===FT(0.5)
        for total in FT.([nextfloat(zero(FT)), 0.01, 1, 1000, floatmax(FT)]),
            a in FT.([nextfloat(zero(FT)), eps(FT), 0.1, 1, floatmax(FT)]),
            b in FT.([nextfloat(zero(FT)), eps(FT), 0.3, 1, floatmax(FT)]),
            eta in (zero(FT), FT(0.5), prevfloat(one(FT)), one(FT))

            coefficients=(chlorophyll = a, carotenoids = b, extra = zero(FT))
            result=partition_absorption(model, total, coefficients)
            value=antenna_excitation(CarotenoidAntenna(eta), result)
            @test isfinite(value)
            @test zero(FT)<=value<=total
            @test eta!=one(FT) || value===total
            @test photosystem_excitation(value, one(FT))===value
        end
    end
end
