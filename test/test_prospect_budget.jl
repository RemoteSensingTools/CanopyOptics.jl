using DelimitedFiles

# General eigendecomposition and angular-integral references generated independently
# with mpmath; see reference/generate_prospect_reference.py. No runtime Python needed.
prospect_reference(name) =
    readdlm(joinpath(@__DIR__, "reference", name), ',', String; skipstart = 1)

parse_reference(FT, s) = FT(parse(BigFloat, s))

function repeated_plate_reference(t, alpha, m::Int)
    r=one(t)-t-alpha
    tm, rm=one(t), zero(t)
    for _ = 1:m
        tm, rm=t*tm/(1-r*rm), r+t*t*rm/(1-r*rm)
    end
    return (; transmission = tm, reflectance = rm, absorptance = one(t)-tm-rm)
end

@testset "PROSPECT independent spectral budgets" begin
    oracle=prospect_reference("prospect_reference.csv")
    for FT in (Float32, Float64)
        rtol=FT==Float32 ? FT(8e-6) : FT(2e-12)
        derivative_atol=FT==Float32 ? FT(2e-6) : FT(2e-12)
        derivative_rtol=FT==Float32 ? FT(2e-5) : FT(2e-9)
        for row in eachrow(oracle)
            k, n, talf, t12, nr=parse_reference.(FT, row[1:5])
            expected=parse_reference.(FT, row[6:8])
            result=@inferred CanopyOptics.prospect_spectral_budget(k, n, talf, t12, nr)
            @test all(x->x isa FT, values(result))
            @test all(x->isfinite(x)&&x>=0&&x<=1, values(result))
            @test abs(sum(values(result))-one(FT)) <= 12eps(FT)
            for (i, key) in enumerate(keys(result))
                # Relative absorption accuracy remains meaningful below eps(R/T).
                # The library E3 continued fraction is ~8 ppm low at k=5 in FP32;
                # fractional powers amplify transmission's relative error with N.
                value_rtol=FT==Float32 && key==:transmission ? FT(1e-4) : rtol
                @test isapprox(
                    getproperty(result, key),
                    expected[i];
                    rtol = value_rtol,
                    atol = FT(4)*nextfloat(zero(FT)),
                )
                dk=ForwardDiff.derivative(
                    x->getproperty(
                        CanopyOptics.prospect_spectral_budget(x, n, talf, t12, nr),
                        key,
                    ),
                    k,
                )
                dn=ForwardDiff.derivative(
                    x->getproperty(
                        CanopyOptics.prospect_spectral_budget(k, x, talf, t12, nr),
                        key,
                    ),
                    n,
                )
                @test dk isa FT
                @test dn isa FT
                @test isapprox(
                    dk,
                    parse_reference(FT, row[8+i]);
                    rtol = derivative_rtol,
                    atol = derivative_atol,
                )
                @test isapprox(
                    dn,
                    parse_reference(FT, row[11+i]);
                    rtol = derivative_rtol,
                    atol = derivative_atol,
                )
            end
        end
        # Internal E3/absorptance formulations meet at k=1, including first AD.
        for key in (:transmission, :absorptance)
            f(k) = getproperty(CanopyOptics.prospect_internal_attenuation(k), key)
            for k in (prevfloat(one(FT)), nextfloat(one(FT)))
                @test isapprox(f(k), f(one(FT)); rtol = 32eps(FT))
                @test isapprox(
                    ForwardDiff.derivative(f, k),
                    ForwardDiff.derivative(f, one(FT));
                    rtol = 64eps(FT),
                )
            end
        end
        zero_slab=CanopyOptics.prospect_internal_attenuation(zero(FT))
        @test zero_slab.transmission===one(FT)
        @test zero_slab.absorptance===zero(FT)
        @test ForwardDiff.derivative(
            x->CanopyOptics.prospect_internal_attenuation(x).absorptance,
            zero(FT),
        )===FT(2)
        @test ForwardDiff.derivative(
            x->CanopyOptics.prospect_internal_attenuation(x).transmission,
            zero(FT),
        )===FT(-2)
    end
end

@testset "PROSPECT extreme stacks and branch continuity" begin
    oracle=prospect_reference("prospect_stack_reference.csv")
    for FT in (Float32, Float64)
        rtol=FT==Float32 ? FT(1e-5) : FT(5e-12)
        for row in eachrow(oracle)
            t, alpha, m=parse_reference.(FT, row[1:3])
            # Cases outside this FT's representable input domain are tested in FP64.
            (!isfinite(m)||iszero(t)||(!iszero(parse(BigFloat, row[2]))&&iszero(alpha))) &&
                continue
            result=CanopyOptics.prospect_plate_stack(one(FT)-t-alpha, t, alpha, m)
            for (i, key) in enumerate(keys(result))
                @test isfinite(getproperty(result, key))
                @test isapprox(
                    getproperty(result, key),
                    parse_reference(FT, row[3+i]);
                    rtol,
                    atol = 4nextfloat(zero(FT)),
                )
            end
            @test abs(sum(values(result))-1)<=12eps(FT)
        end
        for m in (3, 1000), multiplier in FT.([0.999, 1, 1.001])
            t=FT(0.5)
            z=cbrt(eps(FT))/FT(m)^2*multiplier
            alpha=2z/(1+sqrt(1-4z))
            r=one(FT)-t-alpha
            result=CanopyOptics.prospect_plate_stack(r, t, alpha, FT(m))
            reference=setprecision(256) do
                repeated_plate_reference(BigFloat(t), BigFloat(alpha), m)
            end
            for key in keys(result)
                @test isapprox(
                    getproperty(result, key),
                    FT(getproperty(reference, key));
                    rtol,
                )
                derivative=ForwardDiff.derivative(alpha) do a
                    getproperty(
                        CanopyOptics.prospect_plate_stack(one(FT)-t-a, t, a, FT(m)),
                        key,
                    )
                end
                refderivative=setprecision(256) do
                    a=BigFloat(alpha)
                    h=a*big"1e-20"
                    (
                        getproperty(repeated_plate_reference(BigFloat(t), a+h, m), key) -
                        getproperty(repeated_plate_reference(BigFloat(t), a-h, m), key)
                    )/(2h)
                end
                @test isapprox(
                    derivative,
                    FT(refderivative);
                    rtol = FT==Float32 ? FT(2e-4) : FT(2e-9),
                )
            end
        end
        for m in FT.([0, 0.7, 1, 3])
            opaque=CanopyOptics.prospect_plate_stack(FT(0.2), zero(FT), FT(0.8), m)
            @test opaque == (
                iszero(m) ?
                (; transmission = one(FT), reflectance = zero(FT), absorptance = zero(FT)) :
                (; transmission = zero(FT), reflectance = FT(0.2), absorptance = FT(0.8))
            )
        end
    end
end

@testset "PROSPECT actual optical data and public API" begin
    native=createLeafOpticalStruct()
    for FT in (Float32, Float64)
        optics=LeafOpticalProperties(
            native.λ,
            (
                FT.(getproperty(native, key)) for
                key in (:nᵣ, :Kcab, :Kcar, :Kant, :Kb, :Kw, :Km, :Kp, :Kcbc)
            )...,
        )
        for n in FT.([1, 1.2, 1.7, 2, 3]), material in (false, true)
            zero_material=(;
                Ccab = zero(FT),
                Ccar = zero(FT),
                Canth = zero(FT),
                Cbrown = zero(FT),
                Cw = zero(FT),
                Cm = zero(FT),
                Cprot = zero(FT),
                Ccbc = zero(FT),
            )
            leaf=material ? LeafProspectProProperties{FT}(N = n) :
                 LeafProspectProProperties{FT}(; N = n, zero_material...)
            budget=@inferred prospect_budget(leaf, optics)
            t, r=prospect(leaf, optics)
            @test t==budget.transmission && r==budget.reflectance
            @test all(x->eltype(x)===FT, values(budget))
            @test all(x->all(y->isfinite(y)&&0<=y<=1, x), values(budget))
            @test maximum(
                abs.(budget.transmission .+ budget.reflectance .+ budget.absorptance .- 1),
            )<=12eps(FT)
            if !material
                @test all(iszero, budget.absorptance)
                # Analytic conservative plate law, independently of matrix series.
                t12=CanopyOptics.calctav.(FT(90), optics.nᵣ)
                talf=CanopyOptics.calctav.(FT(40), optics.nᵣ)
                r21=1 .- t12 ./ optics.nᵣ .^ 2
                t0=t12 ./ (1 .+ r21)
                ta0=talf ./ (1 .+ r21)
                expected=ta0 ./ (1 .+ (n-1) .* (1 .- t0))
                @test t≈expected rtol=12eps(FT)
                @test r≈(1 .- expected) rtol=12eps(FT)
            end
        end
        for row in eachrow(prospect_reference("prospect_fresnel_reference.csv"))
            nr, ta, t12=parse_reference.(FT, row)
            # Retained Fresnel closed form has looser FP32 accuracy than the stack.
            tolerance=FT==Float32 ? FT(3e-5) : FT(4e-13)
            @test CanopyOptics.calctav(FT(40), nr)≈ta rtol=tolerance
            @test CanopyOptics.calctav(FT(90), nr)≈t12 rtol=tolerance
        end
        for bad in FT.([NaN, Inf, -1, 0, 0.99])
            @test_throws ArgumentError prospect_budget(
                LeafProspectProProperties{FT}(N = bad),
                optics,
            )
        end
        for bad in FT.([NaN, Inf, -1])
            @test_throws ArgumentError prospect_budget(
                LeafProspectProProperties{FT}(Ccab = bad),
                optics,
            )
            @test_throws ArgumentError CanopyOptics.prospect_internal_attenuation(bad)
        end
        for bad in (FT(NaN), FT(Inf), one(FT), prevfloat(one(FT)), nextfloat(one(FT)))
            invalid=deepcopy(optics)
            invalid.nᵣ[1]=bad
            @test_throws ArgumentError prospect_budget(
                LeafProspectProProperties{FT}(),
                invalid,
            )
        end
        invalid=deepcopy(optics)
        invalid.Kcar[1]=-one(FT)
        @test_throws ArgumentError prospect_budget(LeafProspectProProperties{FT}(), invalid)
        # Zero material's public concentration derivative survives the host guards.
        f(c) = begin
            D=typeof(c)
            leaf=LeafProspectProProperties{D}(
                N = D(1.7),
                Ccab = c,
                Ccar = zero(c),
                Canth = zero(c),
                Cbrown = zero(c),
                Cw = zero(c),
                Cm = zero(c),
                Cprot = zero(c),
                Ccbc = zero(c),
            )
            sum(prospect_budget(leaf, optics).absorptance)
        end
        @test isfinite(ForwardDiff.derivative(f, zero(FT)))
        @test ForwardDiff.derivative(f, zero(FT))>0
    end
end
