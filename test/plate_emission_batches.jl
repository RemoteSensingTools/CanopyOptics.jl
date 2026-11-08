using CanopyOptics, Test

function plate_emission_batch_cases(::Type{FT}) where {FT}
    switch = sqrt(sqrt(eps(FT)))
    depths = FT.((
        0,
        1e-20,
        1e-7,
        0.01,
        1,
        prevfloat(FT(8)),
        8,
        nextfloat(FT(8)),
        30,
        700,
        1e20,
        floatmax(FT),
    ))
    bs = FT.((0, 1e-20, prevfloat(switch), switch, nextfloat(switch), 1, 30, floatmax(FT)))
    escapes = RegularPlateEscape{FT}[]
    births = ExponentialBirthProfile{FT}[]
    expected = NamedTuple{(:front, :back, :reabsorbed),Tuple{FT,FT,FT}}[]
    for a in depths,
        b in bs,
        t in FT.((1, 0.6, 1e-20, floatmin(FT))),
        weights in ((1, 0), (1, 0.3), (1, 1))

        f, g = FT.(weights)
        push!(escapes, RegularPlateEscape(b, t))
        push!(births, ExponentialBirthProfile(a, f, g))
        ref = plate_emission_reference(a, b, t, f, g)
        push!(expected, map(FT, ref))
    end
    tiny=nextfloat(zero(FT))
    for a in FT.((0, 1, 8, 100)),
        b in (zero(FT), tiny, FT(10)*tiny, FT(1000)*tiny, floatmin(FT), FT(1e-20))

        push!(escapes, RegularPlateEscape(b, tiny))
        push!(births, ExponentialBirthProfile(a, one(FT), zero(FT)))
        ref=plate_emission_reference(a, b, tiny, one(FT), zero(FT))
        push!(expected, map(FT, ref))
    end
    return escapes, births, expected
end

function check_plate_emission_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        escapes, births, expected=plate_emission_batch_cases(FT)
        tolerance=FT === Float32 ? 3e-5 : 2e-11
        for order in (eachindex(escapes), reverse(eachindex(escapes)), 1:0)
            e, p=to_device(escapes[order]), to_device(births[order])
            result=leaf_emission_fates.(e, p)
            uniform=leaf_emission_fates.(e, Ref(UniformBirthProfile()))
            synchronize()
            result, uniform=Array(result), Array(uniform)
            @test length(result)==length(order)
            for (i, index) in enumerate(order)
                ref=expected[index]
                @test all(isfinite, values(result[i]))
                @test isapprox(sum(values(result[i])), one(FT); atol = 16eps(FT), rtol = 0)
                for field in (:front, :back, :reabsorbed)
                    @test isapprox(
                        getproperty(result[i], field),
                        getproperty(ref, field);
                        rtol = tolerance,
                        atol = 8nextfloat(zero(FT)),
                    )
                end
                uniform_ref=plate_emission_reference(
                    zero(FT),
                    escapes[index].emission_depth,
                    escapes[index].interface_transmission,
                    one(FT),
                    zero(FT),
                )
                for field in (:front, :back, :reabsorbed)
                    @test isapprox(
                        getproperty(uniform[i], field),
                        FT(getproperty(uniform_ref, field));
                        rtol = tolerance,
                        atol = 8nextfloat(zero(FT)),
                    )
                end
            end
        end
    end
end
