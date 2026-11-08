using CanopyOptics, Test, DelimitedFiles

@testset "Diffuse regular plate preparation domains" begin
    for FT in (Float32, Float64)
        for k in FT.((-1, Inf, NaN))
            @test_throws ArgumentError diffuse_plate_optics(k, FT(1.4))
        end
        for nr in FT.((0, 1, Inf, NaN))
            @test_throws ArgumentError diffuse_plate_optics(one(FT), nr)
        end
        @test_throws DomainError diffuse_plate_optics(FT(1000), FT(1.4))
        @test_throws MethodError diffuse_plate_optics(FT(1), 1)
        clear=@inferred diffuse_plate_optics(zero(FT), FT(1.4))
        @test all(v->v isa FT, values(clear))
        @test clear.traversal_transmission===one(FT)
        @test clear.traversal_absorptance===zero(FT)
        @test iszero(clear.effective_depth)
        weak=diffuse_plate_optics(FT(1e-25), FT(1.4))
        @test weak.traversal_transmission===one(FT)
        @test weak.traversal_absorptance>zero(FT)
        @test weak.effective_depth>zero(FT)
        opaque=diffuse_plate_optics(FT(50), FT(1.4))
        @test opaque.traversal_absorptance===one(FT)
        @test opaque.traversal_transmission>zero(FT)
        @test isfinite(opaque.effective_depth)
        for k in FT.((0, 1e-25, 0.1, 1, 5, 50))
            plate=diffuse_plate_optics(k, FT(1.4))
            birth=ExponentialBirthProfile(
                plate.effective_depth,
                one(FT),
                (one(FT)-plate.interface_exit_transmission)*plate.traversal_transmission,
            )
            escape=RegularPlateEscape(
                plate.effective_depth,
                plate.interface_exit_transmission,
            )
            fates=@inferred leaf_emission_fates(escape, birth)
            @test all(v->v isa FT && isfinite(v) && zero(FT)<=v<=one(FT), values(fates))
            @test isapprox(sum(values(fates)), one(FT); rtol = 0, atol = 16eps(FT))
        end
    end
end

@testset "Independent angular optics and source depth integrals" begin
    rows=readdlm(
        joinpath(@__DIR__, "reference", "plate_preparation_reference.csv"),
        ',',
        String;
        skipstart = 1,
    )
    for row in eachrow(rows)
        FT=row[1]=="Float32" ? Float32 : Float64
        k, nr=(FT(parse(BigFloat, s)) for s in row[2:3])
        # SpecialFunctions can underflow E3 before the final factor of two;
        # preparation rejects that unresolved depth, even if its true tau rounds
        # to a positive subnormal. No epsilon transmission is manufactured.
        if iszero(CanopyOptics.prospect_internal_attenuation(k).transmission)
            @test_throws DomainError diffuse_plate_optics(k, nr)
            continue
        end
        plate=@inferred diffuse_plate_optics(k, nr)
        rtol=FT===Float32 ? FT(3e-5) : FT(2e-12)
        for (field, s) in zip(keys(plate), row[4:8])
            ref=FT(parse(BigFloat, s))
            @test isapprox(
                getproperty(plate, field),
                ref;
                rtol,
                atol = 8nextfloat(zero(FT)),
            )
        end
        tau, loss, t12, t21=plate.traversal_transmission,
        plate.traversal_absorptance,
        plate.interface_entry_transmission,
        plate.interface_exit_transmission
        r=one(FT)-t21
        for (entry, reference_column) in
            ((t12, 9), (FT(CanopyOptics.calctav(FT(40), nr)), 11))
            source=entry*loss/(t21+r*loss)
            expected=FT(parse(BigFloat, row[reference_column]))
            @test isapprox(source, expected; rtol, atol = 8nextfloat(zero(FT)))
            budget=CanopyOptics.prospect_spectral_budget(k, one(FT), entry, t12, nr)
            @test isapprox(
                source,
                budget.absorptance;
                rtol = 8eps(FT),
                atol = 8nextfloat(zero(FT)),
            )
        end
        @test isapprox(tau+loss, one(FT); rtol = 0, atol = 2eps(FT))
    end
end
