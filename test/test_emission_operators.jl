# Emission operators from prepared pieces and their linear band observation.
# All spectra, fates, fractions and yields here are synthetic test inputs.
using LinearAlgebra: dot

# Reuse the second support policy defined by the conditional-spectrum tests.
@isdefined(FiniteStokesBand) || include("test_emission_spectrum.jl")

# Conditional probabilities column by column; a column with no supported emission
# is dark (all zeros), which the caller supplies instead of a conditional spectrum.
function conditional_matrix(spectrum, support, excitation)
    FT = eltype(excitation)
    columns = map(excitation) do ex
        any(em -> emission_allowed(support, em, ex), spectrum.wavelength_nm) ?
        emission_probabilities(spectrum, support, ex) :
        zeros(FT, length(spectrum.wavelength_nm))
    end
    return reduce(hcat, columns)
end

# Two synthetic sources with different spectra, support rules and fates (a regular
# plate for one, a dorsiventral plate for the other), plus a third named source.
function synthetic_emission_sources(::Type{FT}) where {FT}
    emission = FT[660, 680, 700, 720, 740, 760, 780]
    excitation = FT[450, 640, 690, 765, 800]
    first = TabulatedPhotonSpectrum(emission, FT[2, 8, 5, 3, 2, 1, 0.5])
    second = TabulatedPhotonSpectrum(emission, FT[0.5, 1, 2, 4, 6, 5, 3])
    third = TabulatedPhotonSpectrum(emission, FT[1, 1, 1, 1, 1, 1, 1])
    regular = [RegularPlateEscape(b, FT(0.6)) for b in FT[3, 2, 1, 0.5, 0.1, 0.01, 0]]
    regular_birth = [ExponentialBirthProfile(a, one(FT), FT(0.2)) for a in FT[4, 2, 1, 0.5, 0]]
    layered = [
        DorsiventralPlateEscape(b, FT(0.7) * b, FT(0.55), FT(0.65), FT(0.9)) for
        b in FT[2.5, 1.5, 0.8, 0.3, 0.1, 0.02, 0]
    ]
    layered_birth = [
        DorsiventralBirthProfile(FT(0.7), FT(0.3), a, one(FT), FT(0.1), a / 2, one(FT), FT(0.3))
        for a in FT[3, 1.5, 0.8, 0.2, 0]
    ]
    return (; emission, excitation, first, second, third, regular, regular_birth, layered, layered_birth)
end

# Independent explicit accumulation: scalar fates per element, BigFloat sums.
function accumulated_faces(probabilities, escape, birth, produced)
    faces = map(_ -> zeros(BigFloat, size(probabilities, 1)), (front = 1, back = 2, reabsorbed = 3))
    for e in axes(probabilities, 2), m in axes(probabilities, 1)
        fate = leaf_emission_fates(escape[m], birth[e])
        weight = BigFloat(probabilities[m, e]) * BigFloat(produced[e])
        faces.front[m] += weight * BigFloat(fate.front)
        faces.back[m] += weight * BigFloat(fate.back)
        faces.reabsorbed[m] += weight * BigFloat(fate.reabsorbed)
    end
    return faces
end

@testset "Emission operators from prepared probabilities and fates" begin
    for FT in (Float32, Float64)
        other = FT === Float32 ? Float64 : Float32
        s = synthetic_emission_sources(FT)
        p1 = conditional_matrix(s.first, StrictStokes(), s.excitation)
        p2 = conditional_matrix(s.second, FiniteStokesBand(FT(20), FT(770)), s.excitation)
        dark1 = [all(iszero, c) for c in eachcol(p1)]
        dark2 = [all(iszero, c) for c in eachcol(p2)]
        @test dark1 == [false, false, false, false, true]
        @test dark2 == [false, false, false, true, true]

        # Source 1: matrix of NamedTuples, as returned by broadcasting fates.
        fates1 = leaf_emission_fates.(reshape(s.regular, :, 1), reshape(s.regular_birth, 1, :))
        # Source 2: NamedTuple of three matrices, keys deliberately reordered.
        fates2_nt = leaf_emission_fates.(reshape(s.layered, :, 1), reshape(s.layered_birth, 1, :))
        fates2 = (
            reabsorbed = map(x -> x.reabsorbed, fates2_nt),
            front = map(x -> x.front, fates2_nt),
            back = map(x -> x.back, fates2_nt),
        )
        ops1 = @inferred emission_operators(p1, fates1)
        ops2 = @inferred emission_operators(p2, fates2)
        @test keys(ops1) == (:front, :back, :reabsorbed)
        @test all(x -> x isa Matrix{FT} && size(x) == size(p1), ops1)
        # Both accepted fate forms give identical operators.
        @test emission_operators(p2, fates2_nt) == ops2
        @test emission_operators(p1, map(k -> getproperty.(fates1, k), (front = :front, back = :back, reabsorbed = :reabsorbed))) == ops1
        # Element definition M = p .* fate, checked with independent scalar loops.
        for (p, ops, escape, birth) in ((p1, ops1, s.regular, s.regular_birth), (p2, ops2, s.layered, s.layered_birth))
            for e in axes(p, 2), m in axes(p, 1)
                fate = leaf_emission_fates(escape[m], birth[e])
                @test ops.front[m, e] === p[m, e] * fate.front
                @test ops.back[m, e] === p[m, e] * fate.back
                @test ops.reabsorbed[m, e] === p[m, e] * fate.reabsorbed
            end
        end

        # Photon closure per excitation column; dark columns are exactly zero.
        for (p, ops, dark) in ((p1, ops1, dark1), (p2, ops2, dark2))
            residual = @inferred emission_closure_residual(ops, p)
            @test residual isa Vector{FT}
            @test length(residual) == size(p, 2)
            @test all(abs.(residual) .<= 16eps(FT))
            @test all(iszero, residual[dark])
            for e in axes(p, 2)
                column = sum(BigFloat, ops.front[:, e]) + sum(BigFloat, ops.back[:, e]) +
                         sum(BigFloat, ops.reabsorbed[:, e])
                @test column ≈ (dark[e] ? 0 : 1) atol = 16eps(FT)
                dark[e] && @test all(x -> iszero(x[:, e]), ops)
            end
            @test all(x -> all(>=(zero(FT)), x), ops)
        end

        # Failure case of the diagnostic: perturbing one operator entry by an exactly
        # representable amount moves only that column's residual. In a dark column the
        # residual equals the perturbation exactly.
        delta = FT(0.25)
        for (e, m) in ((5, 3), (2, 4))
            bad = map(copy, ops1)
            bad.back[m, e] += delta
            base = emission_closure_residual(ops1, p1)
            residual = emission_closure_residual(bad, p1)
            others = setdiff(axes(p1, 2), e)
            @test residual[others] == base[others]
            @test residual[e] ≈ base[e] + delta atol = 16eps(FT)
            dark1[e] && @test residual[e] === delta
        end

        # Two sources with produced photons from a joint allocation and synthetic
        # yields applied by the caller; compare with independent accumulation.
        absorbed = FT[100, 250, 400, 60, 30]
        fractions = (psii = FT(0.52), psi = FT(0.4))
        yields = (psii = FT(0.03), psi = FT(0.006))
        allocated = [allocate_excitation(q, fractions).sinks for q in absorbed]
        # Contract: a source produces no photons in its dark columns (no emission
        # support on this grid). Recipe: dark columns from the probability column sums.
        dark_of(p) = iszero.(vec(sum(p; dims = 1)))
        @test dark_of(p1) == dark1 && dark_of(p2) == dark2
        unmasked = (psii = yields.psii .* getproperty.(allocated, :psii),
                    psi = yields.psi .* getproperty.(allocated, :psi))
        produced = (psii = ifelse.(dark1, zero(FT), unmasked.psii),
                    psi = ifelse.(dark2, zero(FT), unmasked.psi))
        @test all(iszero, produced.psii[dark_of(p1)]) && all(iszero, produced.psi[dark_of(p2)])
        reference1 = accumulated_faces(p1, s.regular, s.regular_birth, produced.psii)
        reference2 = accumulated_faces(p2, s.layered, s.layered_birth, produced.psi)
        for face in (:front, :back, :reabsorbed)
            total = getproperty(ops1, face) * produced.psii + getproperty(ops2, face) * produced.psi
            reference = getproperty(reference1, face) .+ getproperty(reference2, face)
            @test total ≈ FT.(reference) rtol = 16eps(FT)
        end
        # Closure against ALL produced photons, for each source and combined, from the
        # independent accumulation and from the operators themselves.
        outcome(ref) = sum(sum(BigFloat, getproperty(ref, f)) for f in (:front, :back, :reabsorbed))
        applied(ops, q) = sum(sum(BigFloat, getproperty(ops, f) * q) for f in (:front, :back, :reabsorbed))
        @test outcome(reference1) ≈ sum(BigFloat, produced.psii) rtol = 16eps(FT)
        @test outcome(reference2) ≈ sum(BigFloat, produced.psi) rtol = 16eps(FT)
        @test applied(ops1, produced.psii) ≈ sum(BigFloat, produced.psii) rtol = 32eps(FT)
        @test applied(ops2, produced.psi) ≈ sum(BigFloat, produced.psi) rtol = 32eps(FT)
        @test outcome(reference1) + outcome(reference2) ≈
              sum(BigFloat, produced.psii) + sum(BigFloat, produced.psi) rtol = 16eps(FT)

        # Contract violation made visible: photons produced in a dark column vanish,
        # while the closure residual still reports zero for that column. Only the
        # caller's production check detects it; nothing is clipped or renormalized.
        lost = sum(BigFloat, unmasked.psii[dark1]) + sum(BigFloat, unmasked.psi[dark2])
        @test lost > 0
        @test !all(iszero, unmasked.psii[dark_of(p1)])
        @test applied(ops1, unmasked.psii) + applied(ops2, unmasked.psi) ≈
              sum(BigFloat, unmasked.psii) + sum(BigFloat, unmasked.psi) - lost rtol = 32eps(FT)
        @test all(iszero, emission_closure_residual(ops1, p1)[dark1])

        # Extension: three named sources, one call each, no code change.
        p3 = conditional_matrix(s.third, StrictStokes(), s.excitation)
        sources = (
            psii = (probabilities = p1, fates = fates1),
            psi = (probabilities = p2, fates = fates2),
            third = (probabilities = p3, fates = fates1),
        )
        operators = map(src -> emission_operators(src.probabilities, src.fates), sources)
        fractions3 = (psii = FT(0.5), psi = FT(0.35), third = FT(0.1))
        yields3 = (psii = FT(0.03), psi = FT(0.006), third = FT(0.01))
        allocated3 = [allocate_excitation(q, fractions3).sinks for q in absorbed]
        produced3 = map(NamedTuple{keys(sources)}(keys(sources))) do name
            q = getproperty(yields3, name) .* getproperty.(allocated3, name)
            ifelse.(dark_of(sources[name].probabilities), zero(FT), q)
        end
        for (src, op, q) in zip(sources, operators, produced3)
            @test applied(op, q) ≈ sum(BigFloat, q) rtol = 32eps(FT)
        end
        response = FT[0, 1, 1, 0.5, 0, 0, 0]
        signals = map((op, q) -> band_signal(response, op.front, q), operators, produced3)
        @test keys(signals) == (:psii, :psi, :third)
        @test all(x -> x isa FT, values(signals))
        combined = sum(response .* sum(op.front * q for (op, q) in zip(operators, produced3)))
        @test sum(values(signals)) ≈ combined rtol = 16eps(FT)

        # Rejected inputs: grid, element type, names and non-float probabilities.
        @test_throws DimensionMismatch emission_operators(p1[:, 1:4], fates1)
        @test_throws DimensionMismatch emission_operators(p1[1:6, :], fates2)
        @test_throws DimensionMismatch emission_operators(
            p2, (front = fates2.front, back = fates2.back[:, 1:4], reabsorbed = fates2.reabsorbed))
        @test_throws ArgumentError emission_operators(other.(p1), fates1)
        @test_throws ArgumentError emission_operators(other.(p2), fates2)
        @test_throws ArgumentError emission_operators(
            p2, (front = fates2.front, back = other.(fates2.back), reabsorbed = fates2.reabsorbed))
        @test_throws ArgumentError emission_operators(
            p2, (front = fates2.front, back = fates2.back, absorbed = fates2.reabsorbed))
        @test_throws ArgumentError emission_operators(p2, (front = fates2.front, back = fates2.back))
        @test_throws ArgumentError emission_operators(p2, fates2.front)
        @test_throws ArgumentError emission_operators(round.(Int, p1), fates1)
        @test_throws DimensionMismatch emission_closure_residual(ops1, p1[:, 1:4])
        @test_throws ArgumentError emission_closure_residual(ops1, other.(p1))
    end
end

@testset "Band signal: linear observation of hemispheric face totals" begin
    for FT in (Float32, Float64)
        other = FT === Float32 ? Float64 : Float32
        # Hand-computed reference with exactly representable values:
        # M*q = [2.5, 1.75, 4], response'*(M*q) = 2.5 + 0.875.
        response = FT[1, 0.5, 0]
        operator = FT[0.25 0.5; 0.125 0.375; 0.5 0.75]
        produced = FT[2, 4]
        @test (@inferred band_signal(response, operator, produced)) === FT(3.375)
        @test band_signal(zeros(FT, 3), operator, produced) === zero(FT)
        @test band_signal(response, zeros(FT, 3, 2), produced) === zero(FT)
        @test band_signal(response, operator, zeros(FT, 2)) === zero(FT)
        @test band_signal(FT[], zeros(FT, 0, 2), produced) === zero(FT)

        # Linearity in each argument.
        r2 = FT[0.2, 1.5, 3]
        M2 = FT[1 0; 0.5 2; 0.25 0.125]
        q2 = FT[0.75, 3]
        α, β = FT(1.5), FT(-0.25)
        s(r, M, q) = band_signal(r, M, q)
        @test s(α .* response .+ β .* r2, operator, produced) ≈
              α * s(response, operator, produced) + β * s(r2, operator, produced) rtol = 8eps(FT)
        @test s(response, α .* operator .+ β .* M2, produced) ≈
              α * s(response, operator, produced) + β * s(response, M2, produced) rtol = 8eps(FT)
        @test s(response, operator, α .* produced .+ β .* q2) ≈
              α * s(response, operator, produced) + β * s(response, operator, q2) rtol = 8eps(FT)

        # Mismatched grids and element types are rejected.
        @test_throws DimensionMismatch band_signal(FT[1, 0], operator, produced)
        @test_throws DimensionMismatch band_signal(response, operator, FT[1, 2, 3])
        @test_throws DimensionMismatch band_signal(response, transpose(operator), produced)
        @test_throws ArgumentError band_signal(other.(response), operator, produced)
        @test_throws ArgumentError band_signal(response, other.(operator), produced)
        @test_throws ArgumentError band_signal(response, operator, other.(produced))
        @test_throws ArgumentError band_signal([1, 0, 0], [1 0; 0 1; 1 1], [1, 1])

        # Equality with the existing single-source path composed by hand: partition,
        # antenna selection, photosystem allocation, a caller yield, and the operator
        # product written as in the conditional-spectrum batch tests.
        src = synthetic_emission_sources(FT)
        p = conditional_matrix(src.first, StrictStokes(), src.excitation)
        fates = leaf_emission_fates.(reshape(src.regular, :, 1), reshape(src.regular_birth, 1, :))
        hand = map(field -> p .* getproperty.(fates, field), (front = :front, back = :back, reabsorbed = :reabsorbed))
        ops = emission_operators(p, fates)
        @test ops == hand
        coefficients = (chlorophyll = FT(0.6), carotenoids = FT(0.3), water = FT(0.1))
        absorbed = FT[100, 250, 400, 60, 30]
        antenna = antenna_excitation.(Ref(CarotenoidAntenna(FT(0.7))),
                                      partition_absorption.(Ref(HomogeneousAbsorption()), absorbed, Ref(coefficients)))
        q2 = FT(0.03) .* photosystem_excitation.(antenna, FT(0.55))
        q2_joint = FT(0.03) .* [allocate_excitation(a, (psii = FT(0.55),)).sinks.psii for a in antenna]
        @test q2_joint == q2
        # Contract: no production in this source's dark excitation column (800 nm).
        dark = iszero.(vec(sum(p; dims = 1)))
        @test dark == [false, false, false, false, true]
        q2 = ifelse.(dark, zero(FT), q2)
        band = FT.(abs.(src.emission .- FT(685)) .<= FT(10))
        for face in (:front, :back, :reabsorbed)
            observed = band_signal(band, getproperty(ops, face), q2)
            @test observed === dot(band, getproperty(hand, face) * q2)
            reference = sum(BigFloat(band[m]) * BigFloat(getproperty(hand, face)[m, e]) * BigFloat(q2[e])
                            for m in eachindex(band), e in eachindex(q2))
            @test observed ≈ FT(reference) rtol = 16eps(FT)
        end
    end
end
