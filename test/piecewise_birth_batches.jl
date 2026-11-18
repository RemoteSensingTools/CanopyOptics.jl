using CanopyOptics, Test

# Independent closed form for a layer-wise uniform birth profile: exact layer means
# of exp(-b z) and exp(-b (1-z)), then the two boundary balance equations solved
# explicitly. No production helper is called. At least 2048 bits, and enough to
# keep the b^2 term of exp(-b) next to one for subnormal b (1 - front - back).
function piecewise_plate_reference(weights, b, t)
    bits = 2048 + (iszero(b) ? 0 : 2 * max(0, -exponent(b)))
    setprecision(BigFloat, bits) do
        w = BigFloat.(collect(weights))
        w = w ./ sum(w)
        b, t = BigFloat(b), BigFloat(t)
        N = length(w)
        iszero(b) && return (front = BigFloat(1) / 2, back = BigFloat(1) / 2, reabsorbed = zero(b))
        first_front = zero(b)
        first_back = zero(b)
        for k = 1:N
            z0, z1 = BigFloat(k - 1) / N, BigFloat(k) / N
            first_front += w[k] * (exp(-b * z0) - exp(-b * z1)) / (2b * (z1 - z0))
            first_back += w[k] * (exp(-b * (1 - z1)) - exp(-b * (1 - z0))) / (2b * (z1 - z0))
        end
        q = (1 - t) * exp(-b)
        front = t * (first_front + q * first_back) / (1 - q * q)
        back = t * (first_back + q * first_front) / (1 - q * q)
        return (; front, back, reabsorbed = 1 - front - back)
    end
end

# Documented accuracy of the layered kernel, per fate, with u = nextfloat(0):
#   |value - reference| <= 8 (1 + b) eps |reference| + (N + 3) u
#                          + u / (2t)  (reabsorbed only, b < floatmin and 1 - t != 1).
# The (N + 3) u covers underflow of the at most N layer products per face and of
# the closure operations; u/(2t) is the subnormal product r (1 - exp(-b)) of the
# shared reflection closure divided by about t (the exponential path has it too).
function piecewise_tolerance(reference, b::FT, t::FT, N, field) where {FT}
    u = nextfloat(zero(FT))
    closure = field === :reabsorbed && b < floatmin(FT) && one(FT) - t != one(FT) ? u / (2t) : zero(FT)
    return 8 * (1 + b) * eps(FT) * abs(reference) + (N + 3) * u + closure
end
piecewise_close(value::FT, reference, b::FT, t::FT, N, field) where {FT} =
    abs(value - reference) <= piecewise_tolerance(reference, b, t, N, field)

function piecewise_birth_batch_cases(::Type{FT}) where {FT}
    escapes = RegularPlateEscape{FT}[]
    births = PiecewiseUniformBirthProfile{4,FT}[]
    weight_sets = (
        FT.((1, 1, 1, 1)),
        FT.((4, 3, 2, 1)),
        FT.((1, 2, 3, 4)),
        FT.((0, 1, 0, 0)),
        FT.((1, 0, 0, 1e-30)),
        FT.((floatmax(FT), 1, 0, floatmax(FT))),
    )
    for b in FT.((0, 1e-12, 1e-6, 0.01, 0.3, 1, 4, 30, 700)),
        t in FT.((1, 0.6, 1e-3, 1e-20)),
        w in weight_sets

        push!(escapes, RegularPlateEscape(b, t))
        push!(births, PiecewiseUniformBirthProfile(w))
    end
    # Nearly trapping interfaces down to the smallest subnormal, with b comparable
    # to t so that the escape probabilities stay substantial, and b >> t.
    tiny = nextfloat(zero(FT))
    for t in (tiny, 7tiny, 1000tiny, floatmin(FT) / 1024, floatmin(FT)),
        b in (t / 2, t, 3t, FT(1e-30), FT(0.3)),
        w in weight_sets[2:4]

        push!(escapes, RegularPlateEscape(b, t))
        push!(births, PiecewiseUniformBirthProfile(w))
    end
    return escapes, births
end

function check_piecewise_birth_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        escapes, births = piecewise_birth_batch_cases(FT)
        for order in (eachindex(escapes), reverse(eachindex(escapes)), 1:0)
            e, p = to_device(escapes[order]), to_device(births[order])
            result = leaf_emission_fates.(e, p)
            shared = leaf_emission_fates.(e, Ref(births[1]))
            synchronize()
            result, shared = Array(result), Array(shared)
            @test length(result) == length(order)
            for (i, index) in enumerate(order)
                escape = escapes[index]
                ref = piecewise_plate_reference(
                    births[index].layer_fractions,
                    escape.emission_depth,
                    escape.interface_transmission,
                )
                @test all(isfinite, values(result[i]))
                @test isapprox(sum(values(result[i])), one(FT); atol = 4eps(FT), rtol = 0)
                for field in (:front, :back, :reabsorbed)
                    @test piecewise_close(
                        getproperty(result[i], field),
                        getproperty(ref, field),
                        escape.emission_depth,
                        escape.interface_transmission,
                        4,
                        field,
                    )
                end
                # Same compiled kernel on host and device, up to fused operations.
                host = leaf_emission_fates(escape, births[1])
                @test all(isapprox.(values(shared[i]), values(host); rtol = 4eps(FT), atol = 0))
            end
        end
        # One thousand layers on the device, against the 2048-bit closed form.
        many = [PiecewiseUniformBirthProfile(ntuple(k -> FT(1 + (k % 7)), 1000)) for _ = 1:4]
        many_escapes = [RegularPlateEscape(b, FT(0.6)) for b in FT.((1e-12, 0.01, 1, 30))]
        result = Array(leaf_emission_fates.(to_device(many_escapes), to_device(many)))
        synchronize()
        for (i, escape) in enumerate(many_escapes)
            ref = piecewise_plate_reference(many[i].layer_fractions, escape.emission_depth, FT(0.6))
            @test isapprox(sum(values(result[i])), one(FT); atol = 4eps(FT), rtol = 0)
            for field in (:front, :back, :reabsorbed)
                @test piecewise_close(getproperty(result[i], field), getproperty(ref, field),
                                      escape.emission_depth, FT(0.6), 1000, field)
            end
        end
        # Layer fractions of exponential profiles broadcast over a device array.
        exponential = [
            ExponentialBirthProfile(a, one(FT), g) for
            a in FT.((0, 1e-7, 0.5, 1, 3, 30, 1e20)) for g in FT.((0, 0.2, 1))
        ]
        fractions = Array(birth_layer_fractions.(to_device(exponential), Ref(Val(4))))
        synchronize()
        for (k, p) in enumerate(exponential)
            @test fractions[k] isa NTuple{4,FT}
            @test all(isapprox.(fractions[k], birth_layer_fractions(p, Val(4)); rtol = 4eps(FT), atol = 0))
            @test isapprox(sum(fractions[k]), one(FT); atol = 4eps(FT), rtol = 0)
        end
    end
end
