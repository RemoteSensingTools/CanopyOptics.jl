# State-dependent (layer-wise uniform) birth profiles for the regular plate.
# All yields below are synthetic test inputs; the package owns none.
using CanopyOptics, Test

@isdefined(piecewise_plate_reference) || include("piecewise_birth_batches.jl")

const PIECEWISE_FATES = (:front, :back, :reabsorbed)

piecewise_fate_allocation(e, p) = @allocated leaf_emission_fates(e, p)
layer_fraction_allocation(p, n) = @allocated birth_layer_fractions(p, n)
layer_fraction_allocation(p, n, FT) = @allocated birth_layer_fractions(p, n, FT)

# Composite Simpson rule on [lo, hi] with `n` (even) intervals, Float64 or BigFloat.
function simpson(f, lo, hi, n)
    h = (hi - lo) / n
    total = f(lo) + f(hi)
    for i = 1:n-1
        total += (isodd(i) ? 4 : 2) * f(lo + i * h)
    end
    return total * h / 3
end
# One Richardson step on Simpson's rule (Boole's rule, error of order h^6).
boole(f, lo, hi, n) = (16 * simpson(f, lo, hi, 2n) - simpson(f, lo, hi, n)) / 15

# Normalized exponential birth density written out independently of the package.
function exponential_density(a, f, g)
    norm = iszero(a) ? (f + g) : (f + g) * (1 - exp(-a)) / a
    return z -> (f * exp(-a * z) + g * exp(-a * (1 - z))) / norm
end

fate_close(x, y; rtol, atol = 0) =
    all(isapprox(getproperty(x, k), getproperty(y, k); rtol, atol) for k in PIECEWISE_FATES)

@testset "Piecewise birth profile construction" begin
    @test_throws ArgumentError PiecewiseUniformBirthProfile(())
    @test_throws ArgumentError PiecewiseUniformBirthProfile(Float64[])
    @test_throws ArgumentError PiecewiseUniformBirthProfile{0,Float64}(())
    @test_throws ArgumentError PiecewiseUniformBirthProfile((1.0, -1e-300))
    @test_throws ArgumentError PiecewiseUniformBirthProfile((1.0, NaN))
    @test_throws ArgumentError PiecewiseUniformBirthProfile((1.0, Inf))
    @test_throws ArgumentError PiecewiseUniformBirthProfile((0.0, 0.0, 0.0))
    @test_throws ArgumentError PiecewiseUniformBirthProfile([0.0f0])
    @test_throws ArgumentError PiecewiseUniformBirthProfile([-1.0f0, 2.0f0])
    @test_throws MethodError PiecewiseUniformBirthProfile((1.0, 1.0f0))
    @test_throws MethodError PiecewiseUniformBirthProfile((1, 2))
    @test_throws MethodError PiecewiseUniformBirthProfile([1, 2])
    @test_throws MethodError leaf_emission_fates(
        RegularPlateEscape(1.0f0, 0.5f0),
        PiecewiseUniformBirthProfile((1.0, 2.0)),
    )
    # The two-zone escape has its own zone-resolved weights; no method here.
    @test_throws MethodError leaf_emission_fates(
        DorsiventralPlateEscape(1.0, 1.0, 0.6, 0.6),
        PiecewiseUniformBirthProfile((1.0, 2.0)),
    )
    for FT in (Float32, Float64)
        p = PiecewiseUniformBirthProfile(FT.((3, 1, 0, 4)))
        @test p isa PiecewiseUniformBirthProfile{4,FT}
        @test isbitstype(typeof(p))
        @test p.layer_fractions == FT.((0.375, 0.125, 0, 0.5))
        @test PiecewiseUniformBirthProfile(FT[3, 1, 0, 4]) === p
        @test PiecewiseUniformBirthProfile{4,FT}(FT.((3, 1, 0, 4))) === p
        huge = PiecewiseUniformBirthProfile((floatmax(FT), floatmax(FT), zero(FT)))
        @test huge.layer_fractions == FT.((0.5, 0.5, 0))
        tiny = PiecewiseUniformBirthProfile((nextfloat(zero(FT)), zero(FT)))
        @test tiny.layer_fractions == (one(FT), zero(FT))
        w = FT.(abs.(sin.(1:37)))
        q = PiecewiseUniformBirthProfile(w)
        @test all(isapprox.(q.layer_fractions, FT.(w ./ sum(Float64.(w))); rtol = 2eps(FT)))
        # Stored fractions sum to one within two roundings (exact sum), for any N and
        # with a dominant block of weights followed by a tiny one.
        for weights in (
            w,
            ones(FT, 1000),
            [ones(FT, 128); eps(FT) / 4],
            [FT(1 + (k % 7)) for k = 1:1000],
            [k == 500 ? floatmax(FT) : FT(k) for k = 1:1000],
        )
            fractions = PiecewiseUniformBirthProfile(weights).layer_fractions
            @test abs(sum(big.(fractions)) - 1) <= 2eps(FT)
        end
    end
end

@testset "Piecewise birth fates: closure, inference, allocation" begin
    for FT in (Float32, Float64)
        profiles = (
            PiecewiseUniformBirthProfile((one(FT),)),
            PiecewiseUniformBirthProfile(FT.((4, 3, 2, 1))),
            PiecewiseUniformBirthProfile(FT.(abs.(cos.(1:40)))),
            PiecewiseUniformBirthProfile(FT.((0, 0, 1, 0, 0, 0, 0))),
        )
        depths = FT.((0, 1e-30, 1e-12, 1e-4, 0.05, 1, 6, 80, 700, 1e20, floatmax(FT)))
        transmissions = FT.((1, 0.6, 1e-3, 1e-20, floatmin(FT)))
        for p in profiles, b in depths, t in transmissions
            r = leaf_emission_fates(RegularPlateEscape(b, t), p)
            @test all(isfinite, values(r))
            @test all(x -> zero(FT) <= x <= one(FT) + 4eps(FT), values(r))
            @test isapprox(sum(values(r)), one(FT); atol = 4eps(FT), rtol = 0)
        end
        e = RegularPlateEscape(FT(0.7), FT(0.6))
        for p in profiles
            @test @inferred(leaf_emission_fates(e, p)) isa
                  NamedTuple{PIECEWISE_FATES,Tuple{FT,FT,FT}}
            for b in FT.((0, 1e-12, 0.5, 3, 700))
                # Weak-loss series, exp and expm1 branches through a function barrier.
                escape = RegularPlateEscape(b, FT(0.6))
                piecewise_fate_allocation(escape, p)
                @test piecewise_fate_allocation(escape, p) == 0
            end
        end
    end
end

@testset "Piecewise birth fates: accuracy against a 2048-bit closed form" begin
    # Documented accuracy (piecewise_tolerance): relative 8(1+b) eps plus (N + 3)
    # subnormal units; the relative part grows like b eps from rounding of the
    # exponent b*j/N. The worst relative error of normal results is also bounded.
    for FT in (Float32, Float64)
        worst = 0.0
        for N in (1, 2, 3, 7, 16),
            b in (1e-12, 1e-9, 1e-6, 1e-4, 0.01, 0.1, 0.5, 1, 2, 10, 100, 300, 700),
            t in (1, 0.6, 0.1, 1e-3, 1e-6),
            pattern in 1:3

            w = FT.(abs.(sin.(7.3 .* (1:N) .* pattern .+ 1)))
            pattern == 2 && N > 1 && (w[end] = 0)
            pattern == 3 && (w = reverse(w))
            p = PiecewiseUniformBirthProfile(w)
            value = leaf_emission_fates(RegularPlateEscape(FT(b), FT(t)), p)
            ref = piecewise_plate_reference(p.layer_fractions, FT(b), FT(t))
            for k in PIECEWISE_FATES
                expected = getproperty(ref, k)
                @test piecewise_close(getproperty(value, k), expected, FT(b), FT(t), N, k)
                if expected >= floatmin(FT)
                    error = Float64(abs(getproperty(value, k) - expected) / expected) / eps(FT)
                    worst = max(worst, error / (1 + FT(b)))
                end
            end
        end
        @test worst <= 5
    end
end

@testset "Piecewise birth fates: reductions to uniform and mirror symmetry" begin
    for FT in (Float32, Float64)
        for b in FT.((0, 1e-20, 1e-6, 0.01, 1, 30, 700)), t in FT.((1, 0.6, 1e-20))
            e = RegularPlateEscape(b, t)
            uniform = leaf_emission_fates(e, UniformBirthProfile())
            # One layer: first flights coincide with the uniform path; the
            # cancellation-free first loss differs from 1 - survival at roundoff.
            single = leaf_emission_fates(e, PiecewiseUniformBirthProfile((FT(2.5),)))
            @test single.front == uniform.front
            @test single.back == uniform.back
            # The uniform path forms 1 - survival above its series switch (absolute eps).
            @test isapprox(single.reabsorbed, uniform.reabsorbed; rtol = 32eps(FT), atol = 4eps(FT))
            for N in (2, 3, 8, 33)
                flat = PiecewiseUniformBirthProfile(ntuple(_ -> FT(0.7), N))
                r = leaf_emission_fates(e, flat)
                @test fate_close(r, uniform; rtol = 8 * (1 + b) * eps(FT), atol = 4eps(FT))
            end
            # Mirror: exact for dyadic weights (identical normalization), else roundoff.
            p = PiecewiseUniformBirthProfile(FT.((1, 2, 0, 1)))
            m = PiecewiseUniformBirthProfile(FT.((1, 0, 2, 1)))
            r, s = leaf_emission_fates(e, p), leaf_emission_fates(e, m)
            @test r.front == s.back && r.back == s.front && r.reabsorbed == s.reabsorbed
            w = FT.(abs.(sin.(1:9)))
            r = leaf_emission_fates(e, PiecewiseUniformBirthProfile(w))
            s = leaf_emission_fates(e, PiecewiseUniformBirthProfile(reverse(w)))
            @test isapprox(r.front, s.back; rtol = 4eps(FT), atol = floatmin(FT))
            @test isapprox(r.back, s.front; rtol = 4eps(FT), atol = floatmin(FT))
            @test isapprox(r.reabsorbed, s.reabsorbed; rtol = 4eps(FT), atol = floatmin(FT))
        end
    end
end

@testset "Piecewise birth fates: limits" begin
    for FT in (Float32, Float64)
        # Cancellation-free 1 - L(x) on both sides of its series switch at x = 1.
        for x in FT.((0, 1e-30, 1e-8, 1e-3, 0.3, 0.99, 1, 1.01, 3, 50, 1e20)) ∪ (prevfloat(one(FT)), nextfloat(one(FT)))
            X = big(x)
            expected = iszero(X) ? zero(X) : 1 + expm1(-X) / X
            @test isapprox(CanopyOptics.emission_depth_deficit(x), expected; rtol = 4eps(FT), atol = 0)
        end
        p = PiecewiseUniformBirthProfile(FT.((5, 3, 1, 0.5, 0.2)))
        for t in FT.((1, 0.6, 0.01, 1e-8))
            # Lossless plate: (1/2, 1/2, 0) for any birth profile.
            @test leaf_emission_fates(RegularPlateEscape(zero(FT), t), p) ==
                  (front = FT(0.5), back = FT(0.5), reabsorbed = zero(FT))
            # Initial reabsorption slope (2-t)/(2t), independent of birth depth.
            b = FT(1e-6) * eps(FT) * t
            r = leaf_emission_fates(RegularPlateEscape(b, t), p)
            @test isapprox(r.reabsorbed / b, (FT(2) - t) / (FT(2) * t); rtol = 16eps(FT))
        end
        # Transparent interface: faces equal the first flights, no reflection.
        for b in FT.((1e-9, 0.2, 3))
            r = leaf_emission_fates(RegularPlateEscape(b, one(FT)), p)
            ref = piecewise_plate_reference(p.layer_fractions, b, one(FT))
            @test fate_close(r, map(FT, ref); rtol = 8 * (1 + b) * eps(FT))
        end
        # Opaque plate: only the face layers contribute, through the layer-mean
        # survival (1 - exp(-x))/x with x = b/5; deeper layers are below roundoff.
        w = p.layer_fractions
        for b in (FT(400), FT(700), FT(1e20), floatmax(FT))
            r = leaf_emission_fates(RegularPlateEscape(b, FT(0.6)), p)
            x = b / 5
            for (face, weight) in ((r.front, w[1]), (r.back, w[5]))
                @test isapprox(face, FT(0.3) * weight / x; rtol = 8eps(FT), atol = 16nextfloat(zero(FT)))
            end
            @test isapprox(r.reabsorbed, one(FT) - FT(0.3) * (w[1] + w[5]) / x; atol = 4eps(FT), rtol = 0)
        end
        # Nearly trapping interfaces down to the smallest subnormal. With b comparable
        # to t the escapes stay substantial: t/(2(t+b)) per face for b, t -> 0.
        tiny = nextfloat(zero(FT))
        for t in (FT(1e-20), floatmin(FT), floatmin(FT) / 1024, 1000tiny, 7tiny, tiny),
            b in unique((t / 2, t, 3t, FT(1e-30), FT(1e-12), FT(0.01), FT(2)))

            r = leaf_emission_fates(RegularPlateEscape(b, t), p)
            ref = piecewise_plate_reference(p.layer_fractions, b, t)
            @test isapprox(sum(values(r)), one(FT); atol = 4eps(FT), rtol = 0)
            for k in PIECEWISE_FATES
                @test piecewise_close(getproperty(r, k), getproperty(ref, k), b, t, 5, k)
            end
            if b == t
                @test isapprox(r.front, FT(0.25); rtol = 4eps(FT)) && isapprox(r.back, FT(0.25); rtol = 4eps(FT))
            end
        end
    end
end

@testset "Piecewise birth fates: subnormal results and underflowing layer products" begin
    # Regression (review): a deep source at large b gives a front escape just below
    # floatmin whose error is exponent rounding (relative), 5.8 subnormal units here.
    p = PiecewiseUniformBirthProfile(ntuple(k -> k >= 850 ? 1f0 : 0f0, 1000))
    r = leaf_emission_fates(RegularPlateEscape(100f0, 1f0), p)
    ref = piecewise_plate_reference(p.layer_fractions, 100f0, 1f0)
    @test r.front < floatmin(Float32)
    for k in PIECEWISE_FATES
        @test piecewise_close(getproperty(r, k), getproperty(ref, k), 100f0, 1f0, 1000, k)
    end
    for FT in (Float32, Float64)
        u = nextfloat(zero(FT))
        deep = FT === Float32 ? FT.((80, 95, 103, 150)) : FT.((708, 740, 745, 800))
        patterns(N) = (
            [k >= 0.85N ? one(FT) : zero(FT) for k = 1:N],      # sparse deep source
            [k == N ? one(FT) : zero(FT) for k = 1:N],          # last layer only
            [k <= 0.15N + 1 ? one(FT) : zero(FT) for k = 1:N],  # shallow: tiny back escape
            [FT(abs(sin(3.1k))) + (k == 1) for k = 1:N],
        )
        # Large b with subnormal or underflowing escapes, up to 1000 layers.
        for N in (10, 100, 1000), b in deep, t in FT.((1, 0.6, 1e-3)), w in patterns(N)
            q = PiecewiseUniformBirthProfile(w)
            r = leaf_emission_fates(RegularPlateEscape(b, t), q)
            ref = piecewise_plate_reference(q.layer_fractions, b, t)
            for k in PIECEWISE_FATES
                @test piecewise_close(getproperty(r, k), getproperty(ref, k), b, t, N, k)
            end
        end
        # Subnormal and near-floatmin b (first-flight loss b/2 branch), with t near
        # one, at eps scale (closure term u/(2t)) and subnormal.
        for N in (1, 10, 100), b in (2u, 7u, 2^20 * u, floatmin(FT), 100floatmin(FT), prevfloat(eps(FT)), eps(FT), nextfloat(eps(FT))),
            t in (one(FT), FT(1e-3), 4eps(FT), FT(1e-20), 3u),
            w in patterns(N)[[2, 4]]

            q = PiecewiseUniformBirthProfile(w)
            r = leaf_emission_fates(RegularPlateEscape(b, t), q)
            ref = piecewise_plate_reference(q.layer_fractions, b, t)
            @test isapprox(sum(values(r)), one(FT); atol = 4eps(FT), rtol = 0)
            for k in PIECEWISE_FATES
                @test piecewise_close(getproperty(r, k), getproperty(ref, k), b, t, N, k)
            end
        end
    end
end

@testset "Piecewise birth fates: many layers and dominant weights" begin
    # Compensated sums: accuracy and closure do not degrade with N (sequential sums
    # drifted by about N/10 roundings at N = 1000 in Float32).
    for FT in (Float32, Float64)
        profiles = (
            PiecewiseUniformBirthProfile(ntuple(_ -> one(FT), 1000)),
            PiecewiseUniformBirthProfile(ntuple(k -> FT(1 + (k % 7)), 1000)),
            PiecewiseUniformBirthProfile((ntuple(_ -> one(FT), 128)..., eps(FT) / 4)),
            PiecewiseUniformBirthProfile((ntuple(_ -> eps(FT) / 4, 127)..., one(FT))),
        )
        for (i, p) in enumerate(profiles),
            b in FT.((1e-12, 0.01, 1, 30, 700)),
            t in FT.((1, 0.6, 1e-3))

            e = RegularPlateEscape(b, t)
            r = leaf_emission_fates(e, p)
            ref = piecewise_plate_reference(p.layer_fractions, b, t)
            @test isapprox(sum(values(r)), one(FT); atol = 4eps(FT), rtol = 0)
            for k in PIECEWISE_FATES
                @test piecewise_close(getproperty(r, k), getproperty(ref, k), b, t, length(p.layer_fractions), k)
            end
            i == 1 && @test fate_close(r, leaf_emission_fates(e, UniformBirthProfile());
                                       rtol = 8 * (1 + b) * eps(FT), atol = 4eps(FT))
        end
        p = profiles[2]
        e = RegularPlateEscape(FT(0.5), FT(0.6))
        @test @inferred(leaf_emission_fates(e, p)) isa NamedTuple{PIECEWISE_FATES,Tuple{FT,FT,FT}}
        piecewise_fate_allocation(e, p)
        @test piecewise_fate_allocation(e, p) == 0
    end
end

@testset "Piecewise birth fates: a thin source layer against quadrature" begin
    # One lit layer among forty; first flights by Simpson quadrature of the
    # direct-path survival, reflections by the explicit two-face balance.
    N, k = 40, 13
    z0, z1 = (k - 1) / N, k / N
    for b in (0.05, 1.0, 5.0, 30.0), t in (1.0, 0.6, 0.05)
        first_front = simpson(z -> exp(-b * z), z0, z1, 2000) / (2 * (z1 - z0))
        first_back = simpson(z -> exp(-b * (1 - z)), z0, z1, 2000) / (2 * (z1 - z0))
        first_loss = simpson(z -> 1 - (exp(-b * z) + exp(-b * (1 - z))) / 2, z0, z1, 2000) / (z1 - z0)
        q = (1 - t) * exp(-b)
        front = t * (first_front + q * first_back) / (1 - q^2)
        back = t * (first_back + q * first_front) / (1 - q^2)
        # Reabsorbed: first-flight loss plus surviving photons lost on later passes.
        later = (first_front + first_back) - front - back
        expected = (; front, back, reabsorbed = first_loss + later)
        for FT in (Float32, Float64)
            w = ntuple(i -> i == k ? one(FT) : zero(FT), N)
            r = leaf_emission_fates(RegularPlateEscape(FT(b), FT(t)), PiecewiseUniformBirthProfile(w))
            tolerance = FT === Float32 ? 1e-5 * (1 + b) : 1e-12
            @test fate_close(r, expected; rtol = tolerance)
        end
    end
end

@testset "Layer fractions of the excitation profile" begin
    @test_throws MethodError birth_layer_fractions(ExponentialBirthProfile(1.0, 1.0, 0.0), Val(3), Float32)
    @test_throws MethodError birth_layer_fractions(UniformBirthProfile(), Val(3))
    @test_throws ArgumentError birth_layer_fractions(ExponentialBirthProfile(1.0, 1.0, 0.0), Val(0))
    @test_throws ArgumentError birth_layer_fractions(UniformBirthProfile(), Val(0), Float64)
    for FT in (Float32, Float64)
        @test @inferred(birth_layer_fractions(UniformBirthProfile(), Val(4), FT)) == ntuple(_ -> FT(0.25), 4)
        for a in FT.((0, 1e-30, 1e-7, 0.3, prevfloat(one(FT)), 1, 2, 9, 40, 700)),
            (f, g) in ((1, 0), (0, 1), (1, 0.3), (1, 1))

            p = ExponentialBirthProfile(a, FT(f), FT(g))
            P = @inferred birth_layer_fractions(p, Val(6))
            @test P isa NTuple{6,FT}
            @test birth_layer_fractions(p, Val(6), FT) === P
            @test all(>=(zero(FT)), P)
            @test isapprox(sum(P), one(FT); atol = 8eps(FT), rtol = 0)
            # Independent quadrature of the written-out density in BigFloat (a <= 40,
            # where 1000 Boole panels per layer resolve it), and the exact BigFloat
            # antiderivative for every depth.
            density = exponential_density(big(a), big(p.forward_fraction), big(p.backward_fraction))
            for k = 1:6
                lo, hi = big(k - 1) / 6, big(k) / 6
                A, F, G = big(a), big(p.forward_fraction), big(p.backward_fraction)
                exact = iszero(A) ? (hi - lo) :
                        (F * (exp(-A * lo) - exp(-A * hi)) + G * (exp(-A * (1 - hi)) - exp(-A * (1 - lo)))) /
                        ((F + G) * (1 - exp(-A)))
                exact >= floatmin(FT) / eps(FT) || continue
                @test isapprox(P[k], exact; rtol = 8 * (1 + a) * eps(FT))
                a <= 40 || continue
                @test isapprox(P[k], boole(density, lo, hi, 500); rtol = 8 * (1 + a) * eps(FT))
            end
        end
        # a -> 0 gives the uniform split; huge depth puts f in the first and g in the
        # last layer without overflow.
        @test all(isapprox.(birth_layer_fractions(ExponentialBirthProfile(FT(1e-30), one(FT), one(FT)), Val(5)),
                             FT(0.2); rtol = 2eps(FT)))
        for a in (FT(1e20), floatmax(FT))
            p = ExponentialBirthProfile(a, one(FT), FT(3))
            @test birth_layer_fractions(p, Val(3)) == (p.forward_fraction, zero(FT), p.backward_fraction)
        end
        p = ExponentialBirthProfile(FT(2), one(FT), FT(0.2))
        layer_fraction_allocation(p, Val(8))
        @test layer_fraction_allocation(p, Val(8)) == 0
        layer_fraction_allocation(UniformBirthProfile(), Val(8), FT)
        @test layer_fraction_allocation(UniformBirthProfile(), Val(8), FT) == 0
    end
end

# Second-order convergence of the layer-wise profile to the exponential one.
function piecewise_convergence(::Type{FT}, birth, escape, layers) where {FT}
    reference = leaf_emission_fates(escape, birth)
    map(layers) do N
        p = PiecewiseUniformBirthProfile(birth_layer_fractions(birth, Val(N)))
        r = leaf_emission_fates(escape, p)
        maximum(k -> abs(Float64(getproperty(r, k)) - Float64(getproperty(reference, k))), PIECEWISE_FATES)
    end
end

@testset "Piecewise profile converges to the exponential profile" begin
    for FT in (Float32, Float64)
        # From N = 8 on the sequence is asymptotic (N = 2, 4 give orders 1.8-1.9).
        layers = FT === Float64 ? (8, 16, 32, 64, 128, 256) : (8, 16, 32)
        for (a, f, g, b, t) in ((3.0, 1.0, 0.2, 2.0, 0.6), (8.0, 1.0, 0.0, 0.5, 0.3), (1.0, 0.4, 1.0, 5.0, 0.9))
            birth = ExponentialBirthProfile(FT(a), FT(f), FT(g))
            errors = piecewise_convergence(FT, birth, RegularPlateEscape(FT(b), FT(t)), layers)
            orders = [log2(errors[i] / errors[i+1]) for i = 1:length(errors)-1]
            @test all(diff(collect(errors)) .< 0)
            # Observed orders 1.95-2.00 in both precisions (midpoint-like layer error).
            @test all(o -> 1.9 <= o <= 2.1, orders)
        end
    end
end

@testset "Caller recipe: excitation fractions times layer yields" begin
    for FT in (Float32, Float64)
        N = 10
        excitation = ExponentialBirthProfile(FT(2.5), one(FT), FT(0.15))
        P = birth_layer_fractions(excitation, Val(N))
        escape = RegularPlateEscape(FT(0.8), FT(0.6))
        base = leaf_emission_fates(escape, excitation)
        unweighted = leaf_emission_fates(escape, PiecewiseUniformBirthProfile(P))

        # Depth-independent yield: same profile, total production equal to the yield.
        phi = FT(0.035)
        produced = map(x -> x * phi, P)
        @test isapprox(sum(produced), phi; rtol = 8eps(FT))
        constant = PiecewiseUniformBirthProfile(produced)
        @test all(isapprox.(constant.layer_fractions, PiecewiseUniformBirthProfile(P).layer_fractions; rtol = 8eps(FT)))
        @test fate_close(leaf_emission_fates(escape, constant), unweighted; rtol = 16eps(FT))

        # Yield falling linearly with depth (layer midpoints) moves production
        # toward the front: more front and less back escape. The reverse does the opposite.
        mid = ntuple(k -> (FT(k) - FT(0.5)) / N, N)
        falling = PiecewiseUniformBirthProfile(map((x, z) -> x * phi * (1 - FT(0.9) * z), P, mid))
        rising = PiecewiseUniformBirthProfile(map((x, z) -> x * phi * (FT(0.1) + FT(0.9) * z), P, mid))
        r_fall, r_rise = leaf_emission_fates(escape, falling), leaf_emission_fates(escape, rising)
        @test r_fall.front > unweighted.front > r_rise.front
        @test r_fall.back < unweighted.back < r_rise.back
        @test isapprox(unweighted.front, base.front; rtol = 2e-3)

        # A dark layer (zero yield) is allowed; total production excludes it.
        dark_yield = ntuple(k -> k in (4, 5) ? zero(FT) : phi, N)
        dark_production = map(*, P, dark_yield)
        @test isapprox(sum(dark_production), phi * (1 - P[4] - P[5]); rtol = 16eps(FT))
        dark = PiecewiseUniformBirthProfile(dark_production)
        @test iszero(dark.layer_fractions[4]) && iszero(dark.layer_fractions[5])
        @test isapprox(sum(values(leaf_emission_fates(escape, dark))), one(FT); atol = 8eps(FT))
    end
end

@testset "Piecewise profile broadcasting on the CPU" begin
    for FT in (Float32, Float64)
        escapes = [RegularPlateEscape(b, t) for b in FT.((0, 0.1, 1, 10)) for t in FT.((1, 0.5))]
        p = PiecewiseUniformBirthProfile(FT.((1, 3, 2)))
        shared = leaf_emission_fates.(escapes, Ref(p))
        @test shared == [leaf_emission_fates(e, p) for e in escapes]
        profiles = [PiecewiseUniformBirthProfile(FT.((1, k, 2))) for k = 1:length(escapes)]
        @test leaf_emission_fates.(escapes, profiles) ==
              [leaf_emission_fates(e, q) for (e, q) in zip(escapes, profiles)]
    end
end

# End to end: two sources with depth-varying allocation and yield, two lit excitation
# bins with different excitation depths and one entirely dark bin, through
# emission_operators and band_signal. The expected band signals are built from
# independent pieces only: spectral probabilities from photon densities times bin
# widths written out here, excitation layer fractions from the closed-form integral
# in BigFloat, and single-layer fates from the 2048-bit plate reference.
function end_to_end_reference_probabilities(density, width, emission, excitation)
    mass = [big(d) * big(w) * (em > excitation) for (d, w, em) in zip(density, width, emission)]
    total = sum(mass)
    return iszero(total) ? mass : mass ./ total
end

function end_to_end_reference_fractions(a, f, g, N)
    A, F, G = big(a), big(f), big(g)
    norm = iszero(A) ? F + G : (F + G) * (1 - exp(-A)) / A
    layer(lo, hi) = iszero(A) ? (F + G) * (hi - lo) :
                    (F * (exp(-A * lo) - exp(-A * hi)) + G * (exp(-A * (1 - hi)) - exp(-A * (1 - lo)))) / A
    return [layer(big(k - 1) / N, big(k) / N) / norm for k = 1:N]
end

@testset "Layered sources end to end through operators and band signals" begin
    for FT in (Float32, Float64)
        N = 6
        z = ntuple(k -> (FT(k) - FT(0.5)) / N, N)
        emission = FT[680, 700, 720, 740, 760]          # bin centres, nm
        width = FT[10, 20, 20, 20, 30]                  # unequal bin widths, nm
        density = (psii = FT[9, 6, 4, 3, 2], psi = FT[1, 3, 5, 6, 4])   # photons per nm
        spectra = map(d -> TabulatedPhotonSpectrum(emission, d .* width), density)
        escapes = [RegularPlateEscape(b, FT(0.6)) for b in FT[1.5, 0.9, 0.5, 0.3, 0.2]]
        excitation = FT[450, 650, 800]                 # 800 nm: no Stokes emission on the grid
        antenna = FT[300, 200, 0]                      # antenna excitation Q_e; zero in the dark bin
        depth, forward, backward = FT[4, 1.2, 1], one(FT), FT(0.25)
        profiles = [ExponentialBirthProfile(a, forward, backward) for a in depth]
        P = [birth_layer_fractions(p, Val(N)) for p in profiles]
        P_ref = [end_to_end_reference_fractions(a, forward, backward, N) for a in depth]
        for e in eachindex(depth)
            @test all(isapprox.(P[e], P_ref[e]; rtol = 8 * (1 + depth[e]) * eps(FT), atol = 0))
        end
        sources = (:psii, :psi)
        allocation = (psii = map(z -> FT(0.55) - FT(0.1) * z, z), psi = map(z -> FT(0.35) + FT(0.1) * z, z))
        yield = (
            psii = Base.setindex(map(z -> FT(0.03) * (1 - FT(0.7) * z), z), zero(FT), 4),  # dark layer
            psi = map(z -> FT(0.004) * (FT(0.5) + z), z),
        )
        response = FT[0.2, 1, 1, 0.5, 0.1]
        dark = [!any(em -> em > ex, emission) for ex in excitation]
        @test dark == [false, false, true]

        signals = map(sources) do s
            # Convention A: W = Q P alpha phi. Convention B conditions on the source
            # first: Q_s = Q sum(P alpha), P_s = P alpha / sum(P alpha), W = Q_s P_s phi.
            W = [ntuple(k -> antenna[e] * P[e][k] * allocation[s][k] * yield[s][k], N) for e in eachindex(excitation)]
            for e in eachindex(excitation)
                share = map(*, P[e], allocation[s])
                Qs = antenna[e] * sum(share)
                WB = map((x, y) -> Qs * (x / sum(share)) * y, share, yield[s])
                @test all(isapprox.(WB, W[e]; rtol = 8eps(FT), atol = 0))
                # Allocation also inside the yield counts it twice.
                twice = map((x, a, y) -> Qs * (x / sum(share)) * a * y, share, allocation[s], yield[s])
                antenna[e] > 0 && @test !isapprox(sum(twice), sum(W[e]); rtol = 0.1)
            end
            produced = [sum(w) for w in W]
            @test all(iszero, produced[dark])            # dark-column contract
            births = [iszero(q) ? PiecewiseUniformBirthProfile(ntuple(_ -> one(FT), N)) :
                      PiecewiseUniformBirthProfile(w) for (q, w) in zip(produced, W)]
            p = reduce(hcat, [d ? zeros(FT, length(emission)) :
                              emission_probabilities(spectra[s], StrictStokes(), ex)
                              for (d, ex) in zip(dark, excitation)])
            # Spectral probabilities against the independent density-times-width form:
            # lit columns sum to one, the dark column is zero.
            p_ref = [end_to_end_reference_probabilities(density[s], width, emission, ex) for ex in excitation]
            for e in eachindex(excitation)
                @test all(isapprox.(p[:, e], p_ref[e]; rtol = 4eps(FT), atol = 0))
                @test isapprox(sum(big.(p[:, e])), dark[e] ? 0 : 1; atol = 4eps(FT))
                @test abs(sum(p_ref[e]) - (dark[e] ? 0 : 1)) <= 1e-60
                dark[e] && @test all(iszero, p[:, e])
            end
            fates = leaf_emission_fates.(reshape(escapes, :, 1), reshape(births, 1, :))
            ops = emission_operators(p, fates)
            @test maximum(abs, emission_closure_residual(ops, p)) <= 8eps(FT)
            front = band_signal(response, ops.front, produced)
            back = band_signal(response, ops.back, produced)
            # The placeholder profile of the dark bin does not matter.
            other = [iszero(q) ? PiecewiseUniformBirthProfile(ntuple(k -> FT(k), N)) : b
                     for (q, b) in zip(produced, births)]
            ops2 = emission_operators(p, leaf_emission_fates.(reshape(escapes, :, 1), reshape(other, 1, :)))
            @test band_signal(response, ops2.front, produced) == front
            # Explicit accumulation over excitation bin, layer and emission bin of
            # independent production times independent single-layer fates.
            expected_front = expected_back = big(0)
            for e in eachindex(excitation), k = 1:N
                production = big(antenna[e]) * P_ref[e][k] * big(allocation[s][k]) * big(yield[s][k])
                for m in eachindex(emission)
                    onehot = ntuple(j -> j == k ? 1 : 0, N)
                    ref = piecewise_plate_reference(onehot, escapes[m].emission_depth, escapes[m].interface_transmission)
                    weight = production * big(response[m]) * p_ref[e][m]
                    expected_front += weight * ref.front
                    expected_back += weight * ref.back
                end
            end
            @test isapprox(front, expected_front; rtol = 32eps(FT))
            @test isapprox(back, expected_back; rtol = 32eps(FT))
            (; front, back, expected_front, expected_back)
        end
        total = sum(x -> x.front, signals)
        @test isapprox(total, sum(x -> x.expected_front, signals); rtol = 32eps(FT))
    end
end
