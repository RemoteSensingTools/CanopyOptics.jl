struct BatchExtraAntenna{FT} <: AbstractAntennaSelection
    fraction::FT
end
CanopyOptics.select_absorbed_photons(p::BatchExtraAntenna, r::AbsorptionPartitionResult) =
    p.fraction*r.constituents.extra

function check_pigment_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        k=[
            (chlorophyll = FT(c), carotenoids = FT(a), water = FT(w), extra = FT(e)) for
            (c, a, w, e) in (
                (0.3, 0.1, 0.2, 0.4),
                (0.7, 0.3, 0, 0),
                (1, 0, 0, 0),
                (0, 0, 0, 0),
                (-1, 1, 0, 0),
                (NaN, 1, 0, 0),
                (Inf, 1, 0, 0),
            )
        ]
        photons=FT.([0, 1, 1000, 1, 1, 1, 1])
        append!(k, [k[1], k[1], k[1], k[4]])
        append!(photons, FT.([-1, NaN, Inf, 0]))
        append!(
            k,
            [
                (
                    chlorophyll = zero(FT),
                    carotenoids = zero(FT),
                    water = one(FT),
                    extra = one(FT),
                ),
                (
                    chlorophyll = one(FT),
                    carotenoids = one(FT),
                    water = zero(FT),
                    extra = zero(FT),
                ),
            ],
        )
        append!(photons, [nextfloat(zero(FT)), nextfloat(zero(FT))])
        model=HomogeneousAbsorption()
        expected=partition_absorption.(Ref(model), photons, k)
        actual_device=partition_absorption.(Ref(model), to_device(photons), to_device(k))
        synchronize()
        actual=Array(actual_device)
        @test eltype(actual)===eltype(expected)
        for (a, b) in zip(actual, expected)
            @test a.status===b.status
            @test absorption_valid(a)==absorption_valid(b)
            for (x, y) in zip(
                (a.total_absorbed, values(a.constituents)..., values(a.fractions)...),
                (b.total_absorbed, values(b.constituents)..., values(b.fractions)...),
            )
                @test isnan(x) ? isnan(y) :
                      isapprox(x, y; rtol = 8eps(FT), atol = nextfloat(zero(FT)))
            end
        end
        for selector in (
            ChlorophyllOnly(),
            CarotenoidAntenna(zero(FT)),
            CarotenoidAntenna(FT(0.4)),
            CarotenoidAntenna(one(FT)),
            BatchExtraAntenna(FT(0.5)),
        )
            reference=antenna_excitation.(Ref(selector), expected)
            selected=antenna_excitation.(Ref(selector), actual_device)
            synchronize()
            @test eltype(selected)===FT
            for (a, b) in zip(Array(selected), reference)
                @test isnan(a) ? isnan(b) :
                      isapprox(a, b; rtol = 8eps(FT), atol = nextfloat(zero(FT)))
            end
            for fraction in FT.([0, 0.5, 1, NaN, -1, 2])
                ps=Array(photosystem_excitation.(selected, fraction))
                for (a, b) in zip(ps, photosystem_excitation.(reference, fraction))
                    @test isnan(a) ? isnan(b) :
                          isapprox(a, b; rtol = 8eps(FT), atol = nextfloat(zero(FT)))
                end
            end
        end
        order=vcat(collect(2:length(k)), 1)
        reordered=Array(
            partition_absorption.(
                Ref(model),
                to_device(photons[order]),
                to_device(k[order]),
            ),
        )
        @test isequal(reordered, actual[order])
        empty_result=partition_absorption.(
            Ref(model),
            to_device(FT[]),
            to_device(eltype(k)[]),
        )
        @test isempty(Array(empty_result))
        @test eltype(empty_result)===eltype(actual_device)
        zero_material=Array(
            partition_absorption.(Ref(model), to_device(FT.([0, 1])), Ref(NamedTuple())),
        )
        @test zero_material[1].status===AbsorptionNoMaterial
        @test zero_material[2].status===AbsorptionIncompatibleInput
        @test all(
            isnan,
            Array(antenna_excitation.(Ref(ChlorophyllOnly()), to_device(zero_material))),
        )
        # Same spectral-bin extent, type and reduction on each backend. At eta=1
        # and photosystem fraction=1, the selected arrays must be exactly identical.
        spectrum=FT[FT(2)^power for power = -30:30]
        coefficients=fill(
            (
                chlorophyll = FT(0.73),
                carotenoids = FT(0.27),
                water = zero(FT),
                extra = zero(FT),
            ),
            length(spectrum),
        )
        q=to_device(spectrum)
        partition=partition_absorption.(Ref(model), q, to_device(coefficients))
        antenna=antenna_excitation.(Ref(CarotenoidAntenna(one(FT))), partition)
        psii=photosystem_excitation.(antenna, one(FT))
        synchronize()
        @test Array(antenna)==spectrum
        @test Array(psii)==spectrum
        @test sum(psii)===sum(q)
        partial=antenna_excitation.(Ref(CarotenoidAntenna(FT(0.3))), partition)
        @test zero(FT)<=sum(photosystem_excitation.(partial, FT(0.7)))<=sum(q)
    end
end
