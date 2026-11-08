using CanopyOptics, Test
using DelimitedFiles

include("plate_emission_reference.jl")

@testset "Regular plate emission construction and extension" begin
    @test_throws ArgumentError RegularPlateEscape(-1.0, 0.5)
    @test_throws ArgumentError RegularPlateEscape(Inf, 0.5)
    @test_throws ArgumentError RegularPlateEscape(0.0, 0.0)
    @test_throws ArgumentError RegularPlateEscape(0.0, 1.01)
    @test_throws ArgumentError ExponentialBirthProfile(1.0, 0.0, 0.0)
    @test_throws ArgumentError ExponentialBirthProfile(NaN, 1.0, 0.0)
    @test_throws ArgumentError ExponentialBirthProfile(1.0, -1.0, 0.0)
    @test_throws ArgumentError ExponentialBirthProfile(1.0, Inf, 1.0)
    @test_throws MethodError RegularPlateEscape(1.0f0, 0.5)
    @test_throws MethodError leaf_emission_fates(
        RegularPlateEscape(1.0f0, 0.5f0),
        ExponentialBirthProfile(1.0, 1.0, 0.0),
    )
    for FT in (Float32, Float64)
        tiny = floatmin(FT)
        p = ExponentialBirthProfile(FT(2), floatmax(FT), floatmax(FT))
        @test p.forward_fraction == p.backward_fraction == FT(0.5)
        p = ExponentialBirthProfile(FT(2), one(FT), tiny)
        @test p.backward_fraction == tiny
        @test isbitstype(typeof(p))
        e = RegularPlateEscape(FT(0.1), FT(0.6))
        @test @inferred(leaf_emission_fates(e, p)) isa
              NamedTuple{(:front, :back, :reabsorbed),Tuple{FT,FT,FT}}
        @test @inferred(leaf_emission_fates(e, UniformBirthProfile())) isa
              NamedTuple{(:front, :back, :reabsorbed),Tuple{FT,FT,FT}}
        for b in FT.((0, 1e-20, 0.01, 1, 100)), t in FT.((1, 0.6, 1e-20))
            e = RegularPlateEscape(b, t)
            uniform = leaf_emission_fates(e, UniformBirthProfile())
            zero_depth =
                leaf_emission_fates(e, ExponentialBirthProfile(zero(FT), FT(0.1), one(FT)))
            @test all(
                isapprox.(
                    values(uniform),
                    values(zero_depth);
                    atol = 8eps(FT),
                    rtol = 8eps(FT),
                ),
            )
        end
    end
end

plate_fate_allocation(e, b) = @allocated leaf_emission_fates(e, b)

@testset "Regular plate numerical switches and independent path quadrature" begin
    data, _ = readdlm(
        joinpath(@__DIR__, "reference", "plate_emission_reference.csv"),
        ',',
        Any;
        header = true,
    )
    for row in eachrow(data)
        FT = row[1] == "Float32" ? Float32 : Float64
        a, b, t, f, g = FT.(row[2:6])
        value =
            leaf_emission_fates(RegularPlateEscape(b, t), ExponentialBirthProfile(a, f, g))
        tolerance = FT === Float32 ? 3e-5 : 2e-11
        for i = 1:3
            @test isapprox(
                values(value)[i],
                FT(row[6+i]);
                rtol = tolerance,
                atol = 16floatmin(FT),
            )
        end
        @test row[10] < 1e-85
    end
    for FT in (Float32, Float64)
        switch = sqrt(sqrt(eps(FT)))
        depths = FT.((
            0,
            prevfloat(FT(1)),
            1,
            nextfloat(FT(1)),
            prevfloat(FT(8)),
            8,
            nextfloat(FT(8)),
        ))
        emission_depths = FT.((prevfloat(switch), switch, nextfloat(switch), 0.5, 1, 2))
        tolerance = FT === Float32 ? 3e-5 : 2e-11
        for a in depths, b in emission_depths, t in FT.((1, 0.6, 1e-5))
            value = leaf_emission_fates(
                RegularPlateEscape(b, t),
                ExponentialBirthProfile(a, one(FT), FT(0.3)),
            )
            ref = plate_emission_reference(a, b, t, one(FT), FT(0.3))
            for field in (:front, :back, :reabsorbed)
                @test isapprox(
                    getproperty(value, field),
                    FT(getproperty(ref, field));
                    rtol = tolerance,
                    atol = 16floatmin(FT),
                )
            end
        end
        for a in FT.((0.5, 1, 2, 8)),
            b in (prevfloat(a+one(FT)), a+one(FT), nextfloat(a+one(FT)))

            value=leaf_emission_fates(
                RegularPlateEscape(b, FT(0.6)),
                ExponentialBirthProfile(a, one(FT), zero(FT)),
            )
            ref=plate_emission_reference(a, b, FT(0.6), one(FT), zero(FT))
            @test all(
                isapprox.(values(value), FT.(values(ref)); rtol = tolerance, atol = 0),
            )
        end
        # Run the weak-loss and ordinary branches through a function barrier.
        for a in FT.((0, 8, 9)), b in FT.((1e-20, 0.1))
            e=RegularPlateEscape(b, FT(0.6))
            p=ExponentialBirthProfile(a, one(FT), FT(0.3))
            plate_fate_allocation(e, p)
            @test plate_fate_allocation(e, p) == 0
        end
    end
end

@testset "Regular plate independent high-precision balances" begin
    for FT in (Float32, Float64)
        # Finite depths include overflow in a+b and strongly localized sources.
        depths = FT.((0, 1e-20, 1e-7, 0.01, 1, 8, 9, 30, 700, 1e20, floatmax(FT)))
        transmission = FT.((1, 0.6, 0.01, 1e-20, floatmin(FT)))
        tolerance = FT === Float32 ? 3e-5 : 2e-11
        for a in depths,
            b in depths,
            t in transmission,
            weights in ((1, 0), (0, 1), (1, 0.3), (1, 1))

            f, g = FT.(weights)
            p = ExponentialBirthProfile(a, f, g)
            value = leaf_emission_fates(RegularPlateEscape(b, t), p)
            reference = plate_emission_reference(a, b, t, f, g)
            @test all(isfinite, values(value))
            @test all(x -> zero(FT) <= x <= one(FT)+8eps(FT), values(value))
            @test isapprox(sum(values(value)), one(FT); atol = 12eps(FT), rtol = 0)
            for field in (:front, :back, :reabsorbed)
                expected = FT(getproperty(reference, field))
                @test isapprox(
                    getproperty(value, field),
                    expected;
                    rtol = tolerance,
                    atol = 16floatmin(FT),
                )
            end
            mirrored = leaf_emission_fates(
                RegularPlateEscape(b, t),
                ExponentialBirthProfile(a, g, f),
            )
            @test value.front == mirrored.back
            @test value.back == mirrored.front
            @test value.reabsorbed == mirrored.reabsorbed
        end
        for a in FT.((0, 0.01, 1, 8, 100, 1e20)), t in FT.((1, 0.6, 0.01, 1e-8))
            # Resolve derivative in a regime where b/t is small.
            b = FT(1e-6)*eps(FT)*t
            value = leaf_emission_fates(
                RegularPlateEscape(b, t),
                ExponentialBirthProfile(a, one(FT), FT(0.3)),
            )
            @test isapprox(value.reabsorbed/b, (FT(2)-t)/(FT(2)*t); rtol = 16eps(FT))
            lossless = leaf_emission_fates(
                RegularPlateEscape(zero(FT), t),
                ExponentialBirthProfile(a, one(FT), zero(FT)),
            )
            @test lossless == (front = FT(0.5), back = FT(0.5), reabsorbed = zero(FT))
        end
    end
end
