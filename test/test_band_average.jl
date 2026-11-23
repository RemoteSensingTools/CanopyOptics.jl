using CanopyOptics
using ForwardDiff
using Test

# Minimal custom non-1-based AbstractVector / AbstractMatrix wrappers so generic-axis
# contracts are tested directly without requiring OffsetArrays in standalone runs.
struct ShiftedVec{T,V<:AbstractVector{T}} <: AbstractVector{T}
    parent::V
    offset::Int
end
Base.size(v::ShiftedVec) = size(v.parent)
Base.axes(v::ShiftedVec) = (v.offset:(v.offset + length(v.parent) - 1),)
Base.IndexStyle(::Type{<:ShiftedVec}) = IndexLinear()
Base.parent(v::ShiftedVec) = v.parent
Base.dataids(v::ShiftedVec) = Base.dataids(v.parent)
Base.@propagate_inbounds Base.getindex(v::ShiftedVec, i::Int) = v.parent[i - v.offset + 1]
Base.@propagate_inbounds Base.setindex!(v::ShiftedVec, val, i::Int) = (v.parent[i - v.offset + 1] = val; v)
Base.similar(v::ShiftedVec, ::Type{S}, dims::Dims) where {S} = similar(v.parent, S, dims)
function Base.similar(v::ShiftedVec, ::Type{S}, dims::Tuple{AbstractUnitRange}) where {S}
    r = dims[1]
    return ShiftedVec(Vector{S}(undef, length(r)), first(r))
end

struct OffsetBandMat{T,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::M
end
Base.size(m::OffsetBandMat) = size(m.parent)
Base.axes(m::OffsetBandMat) = (0:1, 10:(9 + size(m.parent, 2)))
Base.parent(m::OffsetBandMat) = m.parent
Base.dataids(m::OffsetBandMat) = Base.dataids(m.parent)
Base.@propagate_inbounds Base.getindex(m::OffsetBandMat, i::Int, j::Int) = m.parent[i + 1, j - 9]
function Base.similar(m::OffsetBandMat, ::Type{S}, dims::Tuple{AbstractUnitRange}) where {S}
    r = dims[1]
    return ShiftedVec(Vector{S}(undef, length(r)), first(r))
end

# Custom weight and custom leaf types to verify open multiple-dispatch extensibility.
struct LinearNuWeight{FT<:Real} <: AbstractSpectralWeight
    c::FT
    d::FT
end
function CanopyOptics.weight_interval_integral(w::LinearNuWeight, a::Real, b::Real)
    return w.c * (b - a) + (w.d / 2) * (b^2 - a^2)
end

struct WindowConstantLeaf{FT<:Real}
    lo::FT
    hi::FT
    r::FT
    t::FT
end
function CanopyOptics.leaf_spectrum(m::WindowConstantLeaf{FT}, edges::AbstractVector{<:Real}) where {FT}
    (first(edges) >= m.lo && last(edges) <= m.hi) ||
        throw(CanopyOptics.SpectralCoverageError("outside WindowConstantLeaf"))
    OutFT = float(promote_type(FT, eltype(edges)))
    n = length(edges) - 1
    return fill(OutFT(m.r), n), fill(OutFT(m.t), n)
end

struct BadCoeffCustomLeaf end
function CanopyOptics.leaf_spectrum(::BadCoeffCustomLeaf, ::AbstractVector{<:Real})
    throw(ArgumentError("invalid internal leaf absorption coefficient"))
end

struct BadLengthLeaf end
CanopyOptics.spectral_domain(::BadLengthLeaf) = (1000.0, 5000.0)
CanopyOptics.leaf_spectrum(::BadLengthLeaf, edges::AbstractVector) = (zeros(length(edges) + 2), zeros(length(edges) + 2))

# Independent BigFloat composite Simpson quadrature (n_steps = 10_000 -> error < 1e-16).
function bigfloat_planck_quadrature(a_in, b_in, T_in; n_steps::Int = 10_000)
    return setprecision(256) do
        h = BigFloat("6.62607015e-34")
        c = BigFloat("299792458")
        kb = BigFloat("1.380649e-23")
        c1 = BigFloat(π) * 2 * BigFloat(10)^8 * h * c^2
        c2 = 100 * h * c / kb
        a = BigFloat(a_in)
        b = BigFloat(b_in)
        T = BigFloat(T_in)
        step = (b - a) / (2 * n_steps)
        planck(ν) = iszero(ν) ? zero(BigFloat) : c1 * ν^3 / expm1(c2 * ν / T)
        acc = planck(a) + planck(b)
        for k in 1:(2 * n_steps - 1)
            ν = a + k * step
            acc += (isodd(k) ? 4 : 2) * planck(ν)
        end
        acc * step / 3
    end
end

function check_band_average_zero_alloc!(dest, vals, edges, bands, wt; parent_band = nothing)
    band_average!(dest, vals, edges, bands, wt; parent_band = parent_band)
    return @allocated band_average!(dest, vals, edges, bands, wt; parent_band = parent_band)
end

function check_band_weight_zero_alloc!(dest, bands, wt)
    CanopyOptics.band_weight!(dest, bands, wt)
    return @allocated CanopyOptics.band_weight!(dest, bands, wt)
end

function check_subband_weights_zero_alloc!(dest, sub_bands, parent_band, wt)
    subband_weights!(dest, sub_bands, parent_band, wt)
    return @allocated subband_weights!(dest, sub_bands, parent_band, wt)
end

@testset "Spectral band aggregation and leaf spectrum utilities" begin
    @testset "Export surface vs qualified extension hooks" begin
        for sym in (
            :AbstractSpectralWeight, :FlatWeight, :PlanckWeight, :TabulatedWeight,
            :ReferenceSolar, :band_average, :band_average!,
            :subband_weights, :subband_weights!,
            :wavelength_to_wavenumber, :wavenumber_to_wavelength,
            :TabulatedLeaf, :ConstantLeaf, :CompositeLeaf,
            :leaf_spectrum, :leaf_pigment_shares,
        )
            @test Base.isexported(CanopyOptics, sym)
        end
        for qsym in (:band_weight, :band_weight!, :weight_interval_integral, :spectral_domain, :SpectralCoverageError)
            @test !Base.isexported(CanopyOptics, qsym)
            @test isdefined(CanopyOptics, qsym)
        end
    end

    @testset "Constant spectrum, exact overlap, and analytic linear cases" begin
        for FT in (Float32, Float64)
            rtol_exact = FT === Float32 ? FT(2e-6) : FT(1e-13)

            # 1. Constant spectrum -> constant across all weight types, including mid-bin band cuts
            edges_c = collect(FT(500):FT(25):FT(3000))
            val_const = FT(0.37)
            vals_c = fill(val_const, length(edges_c) - 1)
            bands_c = FT[512.5 980.0 500.0; 1437.5 2763.0 3000.0]
            tw_c = TabulatedWeight(
                collect(FT(400):FT(60):FT(3100)),
                FT[1 + 0.3 * sin(FT(k)) for k in 1:(length(FT(400):FT(60):FT(3100)) - 1)],
            )
            for wt in (FlatWeight(), PlanckWeight(FT(295)), tw_c)
                res = @inferred band_average(vals_c, edges_c, bands_c, wt)
                @test eltype(res) === FT
                @test all(x -> isapprox(x, val_const; rtol = rtol_exact), res)
            end

            # 2. Exact bin-overlap verification when a band slices partially into two fine bins
            edges_ov = FT[100, 200, 400]
            vals_ov = FT[0.2, 0.8]
            bands_ov = FT[150 125; 250 350]
            res_ov = @inferred band_average(vals_ov, edges_ov, bands_ov, FlatWeight())
            @test isapprox(res_ov[1], FT(0.5); rtol = rtol_exact)
            @test isapprox(res_ov[2], FT(0.6); rtol = rtol_exact)

            # 3. Linear spectrum v(ν) = a + b*ν with flat and linear weights against analytic integrals
            n_fine = 4000
            L = FT(1000)
            U = FT(3000)
            edges_lin = collect(range(L, U; length = n_fine + 1))
            mids_lin = FT(0.5) .* (edges_lin[1:end-1] .+ edges_lin[2:end])
            a_v, b_v = FT(0.15), FT(2.0e-4)
            vals_lin = a_v .+ b_v .* mids_lin

            bands_lin = FT[1000 1234.56 1500; 3000 2678.9 2500]
            res_flat = @inferred band_average(vals_lin, edges_lin, bands_lin, FlatWeight())
            for col in 1:size(bands_lin, 2)
                lo, hi = Float64(bands_lin[1, col]), Float64(bands_lin[2, col])
                analytic_flat = Float64(a_v) + 0.5 * Float64(b_v) * (lo + hi)
                @test isapprox(res_flat[col], FT(analytic_flat); rtol = FT === Float32 ? FT(5e-6) : FT(1e-7))
            end

            c_w, d_w = FT(0.5), FT(5.0e-4)
            w_lin_bins = c_w .+ d_w .* mids_lin
            tw_lin = @inferred TabulatedWeight(edges_lin, w_lin_bins)
            res_tw_lin = @inferred band_average(vals_lin, edges_lin, bands_lin, tw_lin)
            res_exact_lin = @inferred band_average(vals_lin, edges_lin, bands_lin, LinearNuWeight(c_w, d_w))
            for col in 1:size(bands_lin, 2)
                lo, hi = Float64(bands_lin[1, col]), Float64(bands_lin[2, col])
                a64, b64 = Float64(a_v), Float64(b_v)
                c64, d64 = Float64(c_w), Float64(d_w)
                num_an = a64 * c64 * (hi - lo) +
                         0.5 * (a64 * d64 + b64 * c64) * (hi^2 - lo^2) +
                         (b64 * d64 / 3.0) * (hi^3 - lo^3)
                den_an = c64 * (hi - lo) + 0.5 * d64 * (hi^2 - lo^2)
                analytic_weighted = FT(num_an / den_an)
                tol_lin = FT === Float32 ? FT(1e-5) : FT(2e-7)
                @test isapprox(res_tw_lin[col], analytic_weighted; rtol = tol_lin)
                @test isapprox(res_exact_lin[col], analytic_weighted; rtol = tol_lin)
            end
        end
    end

    @testset "Partition of a band into sub-bands reproduces parent band value (including mixed coverage)" begin
        for FT in (Float32, Float64)
            rtol = FT === Float32 ? FT(2e-6) : FT(1e-13)
            edges = collect(FT(800):FT(7):FT(3200))
            n = length(edges) - 1
            vals = FT[0.25 + 0.2 * cos(FT(0.003) * edges[i]) + 0.05 * sin(FT(0.011) * edges[i]) for i in 1:n]

            # Sub-band splits deliberately fall inside interior fine bins (non-multiples of 7)
            parent = reshape(FT[950, 2950], 2, 1)
            sub_bands = FT[
                950    1311.3  1987.6  2501.9
                1311.3 1987.6  2501.9  2950
            ]
            tw_edges = collect(FT(900):FT(19):FT(3100))
            tw = TabulatedWeight(tw_edges, FT[1.5 + sin(FT(0.002) * tw_edges[j]) for j in 1:(length(tw_edges) - 1)])

            for wt in (FlatWeight(), PlanckWeight(FT(300)), tw)
                v_parent = only(@inferred band_average(vals, edges, parent, wt))
                v_sub = @inferred band_average(vals, edges, sub_bands, wt; parent_band = parent)
                w_sub = @inferred CanopyOptics.band_weight(sub_bands, wt)
                w_parent = only(@inferred CanopyOptics.band_weight(parent, wt))
                shares = @inferred subband_weights(sub_bands, parent, wt)
                @test eltype(w_sub) === FT && eltype(shares) === FT
                @test isapprox(sum(w_sub), w_parent; rtol = rtol)
                @test isapprox(sum(shares), one(FT); rtol = rtol)
                @test isapprox(sum(w_sub .* v_sub) / sum(w_sub), v_parent; rtol = rtol)
                @test isapprox(sum(shares .* v_sub), v_parent; rtol = rtol)
            end

            # Covered band containing an ozone-blocked zero-flux sub-band:
            # band_weight returns 0 on the blocked sub-band (preserving exact additivity)
            # and subband_weights reproduces band_average!
            tw_ozone = TabulatedWeight(FT[900, 2200, 3100], FT[2.5, 0.0])
            v_par_oz = only(band_average(vals, edges, parent, tw_ozone))
            v_sub_oz = band_average(vals, edges, sub_bands, tw_ozone; parent_band = parent)
            w_sub_oz = CanopyOptics.band_weight(sub_bands, tw_ozone)
            w_par_oz = only(CanopyOptics.band_weight(parent, tw_ozone))
            sh_oz = @inferred subband_weights(sub_bands, parent, tw_ozone)
            @test iszero(w_sub_oz[4]) && iszero(sh_oz[4])
            @test isapprox(sum(w_sub_oz), w_par_oz; rtol = rtol)
            @test isapprox(sum(sh_oz .* v_sub_oz), v_par_oz; rtol = rtol)

            # Uncovered parent band in subband_weights falls back cleanly to FlatWeight shares
            tw_narrow = TabulatedWeight(FT[1500, 2000], FT[1.0])
            sh_fb = @inferred subband_weights(sub_bands, parent, tw_narrow)
            v_sub_fb = @inferred band_average(vals, edges, sub_bands, tw_narrow; parent_band = parent)
            @test isapprox(
                sum(sh_fb .* v_sub_fb),
                only(band_average(vals, edges, parent, tw_narrow));
                rtol = rtol,
            )
            @test_throws ArgumentError CanopyOptics.band_weight(parent, tw_narrow)

            # Mixed coverage regression 1 (code review): parent [1000, 1004], split at 1001,
            # unit-width fine bins [0, 1, 0, 0], table edges [1001, 1002, 1004], densities [9, 1]
            edges_mc1 = FT[1000, 1001, 1002, 1003, 1004]
            vals_mc1 = FT[0, 1, 0, 0]
            par_mc1 = FT[1000; 1004;;]
            sub_mc1 = FT[1000 1001; 1001 1004]
            tw_mc1 = TabulatedWeight(FT[1001, 1002, 1004], FT[9, 1])
            v_par_mc1 = only(@inferred band_average(vals_mc1, edges_mc1, par_mc1, tw_mc1))
            v_sub_mc1 = @inferred band_average(vals_mc1, edges_mc1, sub_mc1, tw_mc1; parent_band = par_mc1)
            sh_mc1 = @inferred subband_weights(sub_mc1, par_mc1, tw_mc1)
            @test isapprox(v_par_mc1, FT(0.25); rtol = rtol)
            @test isapprox(v_sub_mc1, FT[0, 1 // 3]; rtol = rtol)
            @test isapprox(sum(sh_mc1 .* v_sub_mc1), v_par_mc1; rtol = rtol)

            # Mixed coverage regression 2 (physics review): edges [0, 1, 2, 3], values [0, 0, 1],
            # tabulated densities [1, 3] on [1, 2, 3], parent [0, 3] split at 1
            edges_mc2 = FT[0, 1, 2, 3]
            vals_mc2 = FT[0, 0, 1]
            par_mc2 = FT[0; 3;;]
            sub_mc2 = FT[0 1; 1 3]
            tw_mc2 = TabulatedWeight(FT[1, 2, 3], FT[1, 3])
            v_par_mc2 = only(@inferred band_average(vals_mc2, edges_mc2, par_mc2, tw_mc2))
            v_sub_mc2 = @inferred band_average(vals_mc2, edges_mc2, sub_mc2, tw_mc2; parent_band = par_mc2)
            sh_mc2 = @inferred subband_weights(sub_mc2, par_mc2, tw_mc2)
            @test isapprox(v_par_mc2, FT(1 // 3); rtol = rtol)
            @test isapprox(v_sub_mc2, FT[0, 0.5]; rtol = rtol)
            @test isapprox(sum(sh_mc2 .* v_sub_mc2), v_par_mc2; rtol = rtol)

            # Mixed coverage regression 3: table built by ReferenceSolar(grid, irradiance)
            sol_mc = ReferenceSolar(FT[400, 500, 600], FT[0.1, 1.0, 9.0])
            # Implied λ edges [350, 450, 550, 650] nm -> ν edges ≈ [15384.6, 18181.8, 22222.2, 28571.4] cm⁻¹
            edges_mc3 = collect(FT(14000):FT(500):FT(28000))
            vals_mc3 = FT[0.1 + 0.7 * (k / length(edges_mc3))^2 for k in 1:(length(edges_mc3) - 1)]
            par_mc3 = FT[14000; 28000;;]
            sub_mc3 = FT[14000 16000; 16000 28000] # child 1 uncovered, child 2 covered by sol_mc
            v_par_mc3 = only(@inferred band_average(vals_mc3, edges_mc3, par_mc3, sol_mc))
            v_sub_mc3 = @inferred band_average(vals_mc3, edges_mc3, sub_mc3, sol_mc; parent_band = par_mc3)
            sh_mc3 = @inferred subband_weights(sub_mc3, par_mc3, sol_mc)
            @test isapprox(sum(sh_mc3 .* v_sub_mc3), v_par_mc3; rtol = rtol)
            # Without parent_band, child 2 uses sol_mc directly and differs from flat child 2
            @test band_average(vals_mc3, edges_mc3, sub_mc3, sol_mc)[2] != v_sub_mc3[2]

            # Sub-band outside parent_band throws ArgumentError
            @test_throws ArgumentError band_average(vals_mc1, edges_mc1, FT[1000 1001; 1001 1004], tw_mc1; parent_band = FT[1001, 1004])
        end

        # Float32 UV Wien-tail subband_weights does not underflow to 0/0
        edges_uv = Float32[26000, 38000, 50000]
        vals_uv = Float32[0.05, 0.09]
        par_uv = Float32[26000; 50000;;]
        sub_uv = Float32[26000 38000; 38000 50000]
        pw32 = PlanckWeight(300.0f0)
        sh_uv = @inferred subband_weights(sub_uv, par_uv, pw32)
        @test all(isfinite, sh_uv) && isapprox(sum(sh_uv), 1.0f0; rtol = 1.0f-6)
        @test isapprox(
            sum(sh_uv .* band_average(vals_uv, edges_uv, sub_uv, pw32; parent_band = par_uv)),
            only(band_average(vals_uv, edges_uv, par_uv, pw32));
            rtol = 1.0f-6,
        )

        # BigFloat partitions with parent_band (non-dyadic endpoints & sub-Float64 ulp spans)
        setprecision(256) do
            e_bf = [BigFloat("1000.1"), BigFloat("1350.7"), BigFloat("2000.3")]
            v_bf = [BigFloat("0.2"), BigFloat("0.8")]
            par_bf_mat = reshape([BigFloat("1000.1"), BigFloat("2000.3")], 2, 1)
            par_bf_vec = [BigFloat("1000.1"), BigFloat("2000.3")]
            par_bf_tup = (BigFloat("1000.1"), BigFloat("2000.3"))
            sub_bf = [BigFloat("1000.1") BigFloat("1500.9"); BigFloat("1500.9") BigFloat("2000.3")]
            tw_bf_part = TabulatedWeight([BigFloat("1500.9"), BigFloat("2000.3")], [BigFloat("3.0")])

            for wt_bf in (FlatWeight(), tw_bf_part)
                v_par_bf = only(@inferred band_average(v_bf, e_bf, par_bf_mat, wt_bf))
                sh_bf = @inferred subband_weights(sub_bf, par_bf_mat, wt_bf)
                @test eltype(sh_bf) === BigFloat
                for pb_form in (par_bf_mat, par_bf_vec, par_bf_tup)
                    v_sub_bf = @inferred band_average(v_bf, e_bf, sub_bf, wt_bf; parent_band = pb_form)
                    @test eltype(v_sub_bf) === BigFloat
                    @test isapprox(sum(sh_bf .* v_sub_bf), v_par_bf; rtol = BigFloat("1e-60"))
                end
            end

            pw_bf = PlanckWeight(BigFloat(300))
            v_par_pw = only(@inferred band_average(v_bf, e_bf, par_bf_mat, pw_bf))
            v_sub_pw = @inferred band_average(v_bf, e_bf, sub_bf, pw_bf; parent_band = par_bf_mat)
            sh_pw = @inferred subband_weights(sub_bf, par_bf_mat, pw_bf)
            @test eltype(v_sub_pw) === BigFloat && eltype(sh_pw) === BigFloat
            @test isapprox(sum(sh_pw .* v_sub_pw), v_par_pw; rtol = BigFloat("1e-14"))

            # Ultra-narrow BigFloat parent span (< eps(Float64)) that collapses if narrowed to Float64
            base_bf = BigFloat(1000)
            δ_bf = BigFloat("1e-25")
            e_tiny = [base_bf, base_bf + δ_bf, base_bf + 2 * δ_bf]
            par_tiny = reshape([base_bf, base_bf + 2 * δ_bf], 2, 1)
            sub_tiny = [base_bf (base_bf + δ_bf); (base_bf + δ_bf) (base_bf + 2 * δ_bf)]
            v_par_tiny = only(band_average(v_bf, e_tiny, par_tiny, FlatWeight()))
            v_sub_tiny = band_average(v_bf, e_tiny, sub_tiny, FlatWeight(); parent_band = par_tiny)
            sh_tiny = subband_weights(sub_tiny, par_tiny, FlatWeight())
            @test isapprox(sum(sh_tiny .* v_sub_tiny), v_par_tiny; rtol = BigFloat("1e-60"))

            # Sub-band protruding outside BigFloat parent_band by 1e-30 (< eps(Float64)) must throw
            @test_throws ArgumentError band_average(
                v_bf,
                [BigFloat(1000) - BigFloat("1e-30"), BigFloat(1500), BigFloat(2000)],
                [BigFloat(1000) - BigFloat("1e-30") BigFloat(1500); BigFloat(1500) BigFloat(2000)],
                FlatWeight();
                parent_band = [BigFloat(1000), BigFloat(2000)],
            )
        end
    end

    @testset "Partial-coverage, narrowing, aliasing, and invalid-input errors" begin
        edges = [1000.0, 1500.0, 2000.0]
        vals = [0.3, 0.5]

        bands_hi = [1000.0 1500.0; 1500.0 2001.0]
        err_hi = try
            band_average(vals, edges, bands_hi, FlatWeight())
            nothing
        catch e
            e
        end
        @test err_hi isa ArgumentError
        @test occursin("band 2", err_hi.msg)
        @test occursin("2001", err_hi.msg)

        bands_lo = [999.5 1200.0; 1500.0 1800.0]
        err_lo = try
            band_average(vals, edges, bands_lo, FlatWeight())
            nothing
        catch e
            e
        end
        @test err_lo isa ArgumentError
        @test occursin("band 1", err_lo.msg)

        @test_throws ArgumentError band_average(vals, edges, [2500.0; 3000.0;;], FlatWeight())
        @test_throws ArgumentError band_average(vals, edges, [1500.0; 1500.0;;], FlatWeight())
        @test_throws ArgumentError band_average(vals, edges, [1600.0; 1400.0;;], FlatWeight())
        @test_throws DimensionMismatch band_average(vals, edges, zeros(3, 2), FlatWeight())
        @test_throws DimensionMismatch band_average([0.3], edges, [1000.0; 1500.0;;], FlatWeight())
        @test_throws ArgumentError band_average(vals, [1000.0, 1000.0, 2000.0], [1000.0; 1500.0;;], FlatWeight())
        @test_throws ArgumentError PlanckWeight(0.0)
        @test_throws ArgumentError PlanckWeight(0.5)
        @test_throws ArgumentError PlanckWeight(1.0e-20)
        @test_throws ArgumentError PlanckWeight(1.0e9)
        @test_throws ArgumentError PlanckWeight(-10.0)
        @test_throws ArgumentError PlanckWeight(NaN)
        @test_throws ArgumentError PlanckWeight{Float32}(1.0e40)
        @test_throws ArgumentError TabulatedWeight([100.0, 200.0], [-0.1])

        # Typed table constructors validate AFTER narrowing to stored FT (Inf and collapsed edges)
        @test_throws ArgumentError TabulatedWeight{Float32}([1000.0, 2000.0], [1.0e40])
        @test_throws ArgumentError TabulatedWeight{Float32}([1000.0, 1.0e40], [1.0])
        @test_throws ArgumentError TabulatedWeight{Float32}([1000.0, 1000.0 + 1.0e-7, 2000.0], [1.0, 2.0])
        @test_throws ArgumentError TabulatedLeaf{Float32}([1000.0, 1000.0 + 1.0e-7, 2000.0], [0.2, 0.3], [0.1, 0.1])
        @test_throws ArgumentError TabulatedLeaf{Float32}([1000.0, 1.0e40], [0.2], [0.1])

        # Aliased dest against values, edges, bands, parent_band, and weight storage rejected before mutating
        alias_buf = [0.2, 0.8]
        @test_throws ArgumentError band_average!(
            alias_buf, alias_buf, [100.0, 200.0, 300.0], [200.0 100.0; 300.0 200.0], FlatWeight(),
        )
        @test alias_buf == [0.2, 0.8]
        @test_throws ArgumentError band_average!(
            @view(alias_buf[1:2]), alias_buf, [100.0, 200.0, 300.0], [200.0 100.0; 300.0 200.0], FlatWeight(),
        )
        @test alias_buf == [0.2, 0.8]
        alias_edges = [100.0, 200.0, 300.0]
        @test_throws ArgumentError band_average!(
            @view(alias_edges[1:2]), [0.2, 0.8], alias_edges, [200.0 100.0; 300.0 200.0], FlatWeight(),
        )
        @test alias_edges == [100.0, 200.0, 300.0]
        alias_bands = [100.0 200.0; 200.0 300.0]
        @test_throws ArgumentError band_average!(
            @view(alias_bands[1, :]), [0.2, 0.8], [100.0, 200.0, 300.0], alias_bands, FlatWeight(),
        )
        @test alias_bands == [100.0 200.0; 200.0 300.0]
        alias_par = [100.0, 300.0]
        @test_throws ArgumentError band_average!(
            alias_par, [0.2, 0.8], [100.0, 200.0, 300.0], alias_bands, FlatWeight(); parent_band = alias_par,
        )
        @test alias_par == [100.0, 300.0]

        # Weight storage aliasing (TabulatedWeight w and edges) rejected before mutating weight
        tw_alias = TabulatedWeight([100.0, 200.0, 300.0], [0.7, 0.64])
        @test_throws ArgumentError band_average!(
            tw_alias.w, [0.7, 0.4], [100.0, 200.0, 300.0], [100.0 100.0; 200.0 300.0], tw_alias,
        )
        @test tw_alias.w == [0.7, 0.64]
        @test_throws ArgumentError band_average!(
            @view(tw_alias.edges[1:2]), [0.7, 0.4], [100.0, 200.0, 300.0], [100.0 100.0; 200.0 300.0], tw_alias,
        )
        @test tw_alias.edges == [100.0, 200.0, 300.0]
        @test_throws ArgumentError CanopyOptics.band_weight!(tw_alias.w, alias_bands, tw_alias)
        @test tw_alias.w == [0.7, 0.64]
        @test_throws ArgumentError subband_weights!(tw_alias.w, alias_bands, [100.0, 300.0], tw_alias)
        @test tw_alias.w == [0.7, 0.64]

        @test_throws ArgumentError CanopyOptics.band_weight!(@view(alias_bands[1, :]), alias_bands, FlatWeight())
        @test_throws ArgumentError subband_weights!(@view(alias_bands[1, :]), alias_bands, [100.0, 300.0], FlatWeight())

        # Overlapping reshaped view rejected before mutation
        ov_buf = [0.1, 0.2, 100.0, 200.0, 200.0, 300.0]
        @test_throws ArgumentError band_average!(
            @view(ov_buf[2:3]), [0.4, 0.6], [100.0, 200.0, 300.0], @view(reshape(ov_buf, 2, 3)[:, 2:3]), FlatWeight(),
        )
        @test ov_buf == [0.1, 0.2, 100.0, 200.0, 200.0, 300.0]

        # Byte-misaligned reinterpreted strided view overlapping by 4 bytes rejected before mutation
        raw_bytes = reinterpret(UInt8, Float64[0.1, 0.2, 0.3, 0.4, 0.5])
        mis_a = @view reinterpret(Float64, @view(raw_bytes[1:32]))[1:2:3]
        mis_b = @view reinterpret(Float64, @view(raw_bytes[5:36]))[1:2:3]
        @test_throws ArgumentError band_average!(
            mis_a, mis_b, [100.0, 200.0, 300.0], [100.0 200.0; 200.0 300.0], FlatWeight(),
        )

        # Disjoint views of a shared buffer (including reshaped views and views of reshaped buffers) succeed
        shared = [
            0.0, 0.0,                         # 1:2   -> dest
            0.2, 0.8,                         # 3:4   -> vals
            100.0, 200.0, 300.0,              # 5:7   -> edges
            150.0, 250.0, 100.0, 300.0, 0.0,  # 8:11  -> bands (12 is padding for 2×6 reshape)
            100.0, 300.0,                     # 13:14 -> parent_band
        ]
        d_vw = @view shared[1:2]
        v_vw = @view shared[3:4]
        e_vw = @view shared[5:7]
        b_reshape_vw = reshape(@view(shared[8:11]), 2, 2)
        p_reshape_vw = reshape(@view(shared[13:14]), 2, 1)
        ref_disj = band_average([0.2, 0.8], [100.0, 200.0, 300.0], [150.0 100.0; 250.0 300.0], FlatWeight())
        @test band_average!(d_vw, v_vw, e_vw, b_reshape_vw, FlatWeight(); parent_band = p_reshape_vw) == ref_disj
        @test check_band_average_zero_alloc!(d_vw, v_vw, e_vw, b_reshape_vw, FlatWeight(); parent_band = p_reshape_vw) == 0

        shared_mat_buf = [0.0, 0.0, 150.0, 250.0, 100.0, 300.0]
        d_mat_vw = @view shared_mat_buf[1:2]
        b_mat_vw = @view reshape(shared_mat_buf, 2, 3)[:, 2:3]
        b_sub_res_vw = @view reshape(@view(shared_mat_buf[3:6]), 2, 2)[:, 1:2]
        @test band_average!(d_mat_vw, v_vw, e_vw, b_mat_vw, FlatWeight()) == ref_disj
        @test check_band_average_zero_alloc!(d_mat_vw, v_vw, e_vw, b_mat_vw, FlatWeight()) == 0
        fill!(d_mat_vw, 0.0)
        @test band_average!(d_mat_vw, v_vw, e_vw, b_sub_res_vw, FlatWeight()) == ref_disj

        # Disjoint interleaved row views of a reshaped 3×2 buffer
        row_buf = [0.0, 100.0, 200.0, 0.0, 200.0, 300.0]
        d_row = @view row_buf[1:3:4]
        b_rows = @view reshape(row_buf, 3, 2)[2:3, :]
        @test band_average!(d_row, [0.2, 0.8], [100.0, 200.0, 300.0], b_rows, FlatWeight()) == [0.2, 0.8]
        @test subband_weights!(d_row, b_rows, [100.0, 300.0], FlatWeight()) == [0.5, 0.5]
        @test CanopyOptics.band_weight!(d_row, b_rows, FlatWeight()) == [100.0, 100.0]
    end

    @testset "ReferenceSolar() default loader, synthetic CSV parser, and dark-UV band rules" begin
        if isfile(CanopyOptics._DEFAULT_SOLAR_CSV)
            sol_default = @inferred ReferenceSolar()
            @test sol_default isa TabulatedWeight{Float64}
            @test length(sol_default.w) >= 2 && length(sol_default.edges) == length(sol_default.w) + 1
        else
            err_def = try
                ReferenceSolar()
                nothing
            catch e
                e
            end
            @test err_def isa ArgumentError
            @test occursin("reference_solar_global.csv", err_def.msg)
        end

        mktempdir() do tmpdir
            missing_csv = joinpath(tmpdir, "missing_reference_solar.csv")
            err_missing = try
                ReferenceSolar(missing_csv)
                nothing
            catch e
                e
            end
            @test err_missing isa ArgumentError
            @test occursin("missing_reference_solar.csv", err_missing.msg)

            csv_path = joinpath(tmpdir, "synthetic_solar.csv")
            # Explicit zero-irradiance rows below 280 nm (190–280 nm ozone-blocked),
            # illuminated UV/VIS above 280 nm.
            write(
                csv_path,
                """
                # Synthetic clear-sky global surface solar spectrum for unit testing
                wavelength_nm,irradiance_W_m2_nm
                190.0, 0.0
                280.0, 0.0
                282.0, 0.05
                310.0, 0.55
                345.0, 0.95
                500.0, 1.80
                800.0, 1.10
                """,
            )
            sol_csv = @inferred ReferenceSolar(csv_path)
            @test sol_csv isa TabulatedWeight{Float64}

            # Fine leaf grid covering 195–750 nm
            edges_uv_vis = wavelength_to_wavenumber(195.0:0.5:750.0)
            n_uv = length(edges_uv_vis) - 1
            mids_nu = 0.5 .* (edges_uv_vis[1:end-1] .+ edges_uv_vis[2:end])
            # Step change at 280 nm (35714.29 cm⁻¹) so dark (< 280 nm) and illuminated (> 280 nm) parts differ
            vals_step = [1.0e7 / ν < 280.0 ? 0.90 : 0.10 for ν in mids_nu]
            vals_step_alt = [1.0e7 / ν < 280.0 ? 0.40 : 0.10 for ν in mids_nu]

            band_dark_200_263 = [1.0e7 / 263.0; 1.0e7 / 200.0;;]
            band_part_263_345 = [1.0e7 / 345.0; 1.0e7 / 263.0;;]
            band_par_200_345  = [1.0e7 / 345.0; 1.0e7 / 200.0;;]
            sub_200_345       = [1.0e7 / 345.0  1.0e7 / 263.0; 1.0e7 / 263.0  1.0e7 / 200.0]

            # 1. Fully dark 200–263 nm band: covered with 0.0 W m⁻² physical solar flux;
            #    band_average returns the explicit FlatWeight placeholder rather than NaN.
            @test iszero(only(CanopyOptics.band_weight(band_dark_200_263, sol_csv)))
            @test only(band_average(vals_step, edges_uv_vis, band_dark_200_263, sol_csv)) ==
                  only(band_average(vals_step, edges_uv_vis, band_dark_200_263, FlatWeight()))

            # 2. Partially illuminated 263–345 nm band: positive solar flux; retains illuminated
            #    > 280 nm portion and ignores dark < 280 nm values (up to midpoint edge 281 nm).
            @test only(CanopyOptics.band_weight(band_part_263_345, sol_csv)) > 0.0
            v_part_1 = only(band_average(vals_step, edges_uv_vis, band_part_263_345, sol_csv))
            v_part_2 = only(band_average(vals_step_alt, edges_uv_vis, band_part_263_345, sol_csv))
            # Check strictly below 279 nm (since 280 and 282 midpoint edge is 281 nm, below 281 nm is zero density)
            vals_below_245 = [1.0e7 / ν < 235.0 ? 0.95 : 0.10 for ν in mids_nu]
            @test isapprox(v_part_1, 0.10; atol = 1e-12)
            @test isapprox(v_part_2, 0.10; atol = 1e-12)
            @test v_part_1 != only(band_average(vals_step, edges_uv_vis, band_part_263_345, FlatWeight()))
            @test only(band_average(vals_below_245, edges_uv_vis, band_par_200_345, sol_csv)) ≈ 0.10

            # 3. Parent 200–345 nm split into [263–345 nm, 200–263 nm] in wavenumber:
            #    illuminated child gets share 1.0, dark child gets share 0.0, exact recombination.
            sh_uv_sol = @inferred subband_weights(sub_200_345, band_par_200_345, sol_csv)
            @test sh_uv_sol[1] ≈ 1.0 && iszero(sh_uv_sol[2])
            v_sub_sol = band_average(vals_step, edges_uv_vis, sub_200_345, sol_csv; parent_band = band_par_200_345)
            @test sum(sh_uv_sol .* v_sub_sol) ≈ only(band_average(vals_step, edges_uv_vis, band_par_200_345, sol_csv))

            # Headerless valid CSV loads cleanly
            no_hdr_csv = joinpath(tmpdir, "no_header.csv")
            write(no_hdr_csv, "300.0,0.2\n500.0,1.6\n700.0,1.1\n")
            @test ReferenceSolar(no_hdr_csv) == ReferenceSolar([300.0, 500.0, 700.0], [0.2, 1.6, 1.1])

            # Malformed CSV files throw ArgumentError naming the offending line (never silently treated as header)
            @test_throws ArgumentError ReferenceSolar(joinpath(tmpdir, "nonexistent.csv"))
            for (fname, content, expected_line) in (
                ("bad_first_col2.csv", "400,oops\n500.0,1.0\n600.0,1.2\n", "line 1"),
                ("bad_first_col1.csv", "oops,1.0\n500.0,1.0\n600.0,1.2\n", "line 1"),
                ("empty_col_hdr.csv", "wavelength_nm,\n400.0,1.0\n500.0,1.2\n", "line 1"),
                ("two_headers.csv", "wavelength_nm,irradiance\nunits_nm,W_m2_nm\n400.0,1.0\n500.0,1.2\n", "line 2"),
                ("bad_later_row.csv", "# comment\nwavelength_nm,irradiance\n400.0,1.0\n500.0,oops\n", "line 4"),
                ("bad_cols.csv", "400.0,1.0,2.0\n500.0,1.2,2.1\n", "line 1"),
                ("bad_nan.csv", "400.0,NaN\n500.0,1.2\n", "line 1"),
                ("bad_neg.csv", "400.0,-0.1\n500.0,1.2\n", "line 1"),
                ("bad_desc.csv", "500.0,1.0\n400.0,1.2\n", "line 2"),
            )
                fpath = joinpath(tmpdir, fname)
                write(fpath, content)
                err_csv = try
                    ReferenceSolar(fpath)
                    nothing
                catch e
                    e
                end
                @test err_csv isa ArgumentError
                @test occursin(expected_line, err_csv.msg)
            end
        end
    end

    @testset "TabulatedWeight fallback to FlatWeight on uncovered or zero-flux bands" begin
        edges = collect(1000.0:100.0:10000.0)
        vals = [0.1 + 0.8 * (i / length(edges)) for i in 1:(length(edges) - 1)]

        # Weight table covers [3000, 7000] cm⁻¹, with zero density on [5000, 6000].
        tw = TabulatedWeight([3000.0, 5000.0, 6000.0, 7000.0], [4.0, 0.0, 1.0])
        bands = [
            1000.0  2500.0  3000.0  5100.0  7500.0
            2000.0  4000.0  7000.0  5900.0  9000.0
        ]
        res_tw = @inferred band_average(vals, edges, bands, tw)
        res_flat = band_average(vals, edges, bands, FlatWeight())

        @test all(isfinite, res_tw)
        @test res_tw[1] == res_flat[1]
        @test res_tw[2] == res_flat[2]
        @test res_tw[4] == res_flat[4]
        @test res_tw[5] == res_flat[5]
        @test res_tw[3] != res_flat[3]
        @test res_tw[3] < res_flat[3]
    end

    @testset "PlanckWeight against BigFloat quadrature, Stefan-Boltzmann, and cold-T partition" begin
        h = 6.62607015e-34
        c = 299792458.0
        kb = 1.380649e-23
        sigma_sb = 2 * π^5 * kb^4 / (15 * h^3 * c^2)
        T_ref = 300.0
        w_tot = only(CanopyOptics.band_weight([0.0; 25000.0;;], PlanckWeight(T_ref)))
        @test isapprox(w_tot, sigma_sb * T_ref^4; rtol = 1e-13)

        edges_lw = [10.0, 250.0, 630.0, 980.0, 1480.0, 2250.0, 3250.0]
        vals_lw = [0.92, 0.95, 0.97, 0.96, 0.94, 0.91]
        bands_lw = [10.0 250.0 500.0; 3250.0 1480.0 2000.0]

        for T in (250.0, 300.0, 310.0)
            pw = PlanckWeight(T)
            got = @inferred band_average(vals_lw, edges_lw, bands_lw, pw)
            got_w = @inferred CanopyOptics.band_weight(bands_lw, pw)
            for b in 1:size(bands_lw, 2)
                lo, hi = bands_lw[1, b], bands_lw[2, b]
                w_bf = bigfloat_planck_quadrature(lo, hi, T)
                @test isapprox(got_w[b], Float64(w_bf); rtol = 1e-13)

                num_bf = zero(BigFloat)
                for i in 1:length(vals_lw)
                    a = max(lo, edges_lw[i])
                    bi = min(hi, edges_lw[i + 1])
                    if bi > a
                        num_bf += BigFloat(vals_lw[i]) * bigfloat_planck_quadrature(a, bi, T)
                    end
                end
                @test isapprox(got[b], Float64(num_bf / w_bf); rtol = 1e-13)
            end
        end

        # Extreme Wien-tail Float32 and minimum supported cryogenic temperature T = 1 K
        edges_uv32 = Float32[38000, 44000, 50000]
        vals_uv32 = Float32[0.05, 0.08]
        res_uv32 = @inferred band_average(vals_uv32, edges_uv32, Float32[38000; 50000;;], PlanckWeight(300.0f0))
        @test eltype(res_uv32) === Float32
        @test isfinite(only(res_uv32)) && 0.05f0 <= only(res_uv32) <= 0.08f0

        # Cold-T regression on unit-bin partition [1000, 1004] split at 1001 with [0, 1, 0, 0]
        edges_cold = [1000.0, 1001.0, 1002.0, 1003.0, 1004.0]
        vals_cold = [0.0, 1.0, 0.0, 0.0]
        par_cold = [1000.0; 1004.0;;]
        sub_cold = [1000.0 1001.0; 1001.0 1004.0]
        for T_cold in (1.0, 10.0, 50.0)
            pw_cold = PlanckWeight(T_cold)
            v_par_c = only(@inferred band_average(vals_cold, edges_cold, par_cold, pw_cold))
            v_sub_c = @inferred band_average(vals_cold, edges_cold, sub_cold, pw_cold; parent_band = par_cold)
            v_sub_no_par = @inferred band_average(vals_cold, edges_cold, sub_cold, pw_cold)
            sh_c = @inferred subband_weights(sub_cold, par_cold, pw_cold)
            @test isapprox(sum(sh_c .* v_sub_c), v_par_c; rtol = 1e-13)
            @test isapprox(sum(sh_c .* v_sub_no_par), v_par_c; rtol = 1e-13)
        end
    end

    @testset "Converter round trip, duplicate/descending rejection, and Jacobian flux conservation" begin
        for FT in (Float32, Float64)
            rtol = FT === Float32 ? FT(2e-6) : FT(1e-13)
            λ_nm = collect(FT(400):FT(1):FT(2500))
            rho_nm = FT[0.05 + 0.4 * (1 / (1 + exp(-(Float64(λ) - 710.0) / 15.0))) for λ in λ_nm]
            irr_nm = FT[1.8 * exp(-((Float64(λ) - 550.0) / 600.0)^2) + 0.05 for λ in λ_nm]

            edges_only = @inferred wavelength_to_wavenumber(λ_nm)
            edges_r, rho_nu = @inferred wavelength_to_wavenumber(λ_nm, rho_nm; density = false)
            edges_e, irr_nu = @inferred wavelength_to_wavenumber(λ_nm, irr_nm; density = true)
            @test edges_only == edges_r == edges_e
            @test eltype(edges_r) === FT && eltype(rho_nu) === FT && eltype(irr_nu) === FT
            @test issorted(edges_r) && length(edges_r) == length(λ_nm) + 1

            λ_rt_only = @inferred wavenumber_to_wavelength(edges_r)
            λ_rt1, rho_rt = @inferred wavenumber_to_wavelength(edges_r, rho_nu; density = false)
            λ_rt2, irr_rt = @inferred wavenumber_to_wavelength(edges_e, irr_nu; density = true)
            @test isapprox(λ_rt_only, λ_nm; rtol = rtol)
            @test isapprox(λ_rt1, λ_nm; rtol = rtol)
            @test rho_rt == rho_nm
            @test isapprox(λ_rt2, λ_nm; rtol = rtol)
            @test isapprox(irr_rt, irr_nm; rtol = rtol)

            flux_nm = sum(Float64.(irr_nm)) * 1.0
            flux_nu = sum(Float64.(irr_nu) .* diff(Float64.(edges_e)))
            @test isapprox(FT(flux_nu), FT(flux_nm); rtol = rtol)

            sol = @inferred ReferenceSolar(λ_nm, irr_nm)
            @test sol isa TabulatedWeight{FT}
            @test isapprox(sol.edges, edges_e; rtol = rtol)
            @test isapprox(sol.w, irr_nu; rtol = rtol)

            λ_nonunif = vcat(collect(FT(280):FT(0.5):FT(400)), collect(FT(401):FT(1):FT(700)), collect(FT(705):FT(5):FT(1000)))
            e_nonunif = fill(one(FT), length(λ_nonunif))
            sol_nu = @inferred ReferenceSolar(λ_nonunif, e_nonunif)
            @test sol_nu isa TabulatedWeight{FT}
            full_band = reshape(FT[first(sol_nu.edges), last(sol_nu.edges)], 2, 1)
            @test isapprox(only(CanopyOptics.band_weight(full_band, sol_nu)), FT(1002.5 - 279.75); rtol = rtol)

            @test_throws ArgumentError wavelength_to_wavenumber(FT[400, 401, 405], FT[0.1, 0.2, 0.3])
        end

        # Duplicate and descending wavelength samples within spacing tolerance must throw
        @test_throws ArgumentError wavelength_to_wavenumber(Float32[400, 400, 400.0002])
        @test_throws ArgumentError wavelength_to_wavenumber(Float32[400, 400, 400.0002], Float32[0.1, 0.2, 0.3])
        @test_throws ArgumentError wavelength_to_wavenumber(Float32[400, 399.9999, 400.0002])
        @test_throws ArgumentError wavelength_to_wavenumber(Float32[400, 399.9999, 400.0002], Float32[0.1, 0.2, 0.3])
        @test_throws ArgumentError wavelength_to_wavenumber([400.0, 400.0, 400.0 + 1e-13])

        # Documented 400–800 nm single-bin mid-bin cut behavior (whole-bin 400 vs 600–800 nm cut 400/3)
        sol_1bin = ReferenceSolar([400.0, 800.0], [1.0]; grid_unit = :nm)
        @test isapprox(only(CanopyOptics.band_weight([12500.0; 25000.0;;], sol_1bin)), 400.0; rtol = 1e-14)
        @test isapprox(only(CanopyOptics.band_weight([12500.0; 1.0e7 / 600.0;;], sol_1bin)), 400.0 / 3.0; rtol = 1e-14)
    end

    @testset "Flat-in-wavelength versus flat-in-wavenumber on PROSPECT red-edge band" begin
        for FT in (Float32, Float64)
            λ_nm = FT(400):one(FT):FT(2500)
            edges = @inferred wavelength_to_wavenumber(λ_nm)
            leaf = LeafProspectProProperties{FT}()
            rho, tau = @inferred leaf_spectrum(leaf, edges)
            @test eltype(rho) === FT && eltype(tau) === FT
            @test length(rho) == 2101 && length(tau) == 2101

            red_edge_band = reshape(FT[12850, 16000], 2, 1)
            rho_flat_nu = only(@inferred band_average(rho, edges, red_edge_band, FlatWeight()))
            tau_flat_nu = only(@inferred band_average(tau, edges, red_edge_band, FlatWeight()))

            w_flat_lambda = @inferred ReferenceSolar(λ_nm, ones(FT, length(λ_nm)))
            rho_flat_lam = only(@inferred band_average(rho, edges, red_edge_band, w_flat_lambda))
            tau_flat_lam = only(@inferred band_average(tau, edges, red_edge_band, w_flat_lambda))

            for v in (rho_flat_nu, tau_flat_nu, rho_flat_lam, tau_flat_lam)
                @test isfinite(v)
                @test zero(FT) < v < one(FT)
            end
            @test zero(FT) < rho_flat_nu + tau_flat_nu < one(FT)
            @test zero(FT) < rho_flat_lam + tau_flat_lam < one(FT)

            @test rho_flat_lam > rho_flat_nu + FT(0.015)
            @test tau_flat_lam > tau_flat_nu + FT(0.015)

            exact_nm_band = reshape(FT[FT(1e7) / FT(800.5), FT(1e7) / FT(599.5)], 2, 1)
            _, rho_back_nm = wavenumber_to_wavelength(edges, rho)
            arith_nm_mean = sum(rho_back_nm[201:401]) / FT(201)
            @test isapprox(
                only(band_average(rho, edges, exact_nm_band, w_flat_lambda)),
                arith_nm_mean;
                rtol = FT === Float32 ? FT(5e-6) : FT(1e-13),
            )
        end
    end

    @testset "Leaf models, bounded CompositeLeaf, custom probe safety, and PROSPECT pigment shares" begin
        for FT in (Float32, Float64)
            rtol = FT === Float32 ? FT(2e-5) : FT(1e-12)

            tab1 = @inferred TabulatedLeaf(FT[2000, 4000, 6000], FT[0.10, 0.20], FT[0.05, 0.15])
            tab2 = @inferred TabulatedLeaf(FT[1000, 5000, 8000], FT[0.35, 0.45], FT[0.25, 0.30])
            c_fallback = @inferred ConstantLeaf(FT(0.03), zero(FT))
            c_swir = @inferred ConstantLeaf(FT(0.04), zero(FT); ν_min = FT(2500), ν_max = FT(4500))

            @test_throws ArgumentError leaf_spectrum(tab1, FT[1500, 3000, 5000])
            @test_throws ArgumentError leaf_spectrum(c_swir, FT[2000, 3000])
            @test_throws ArgumentError ConstantLeaf(FT(0.7), FT(0.5))
            @test_throws ArgumentError TabulatedLeaf(FT[1000, 2000], FT[0.6], FT[0.6])

            # CompositeLeaf: first model covering each bin wins
            comp = @inferred CompositeLeaf(tab1, tab2, c_fallback)
            query_edges = FT[500, 1000, 2000, 4000, 5000, 6000, 8000, 9000]
            r_comp, t_comp = @inferred leaf_spectrum(comp, query_edges)
            @test eltype(r_comp) === FT && eltype(t_comp) === FT
            @test isapprox(r_comp, FT[0.03, 0.35, 0.10, 0.20, 0.20, 0.45, 0.03]; rtol = rtol)
            @test isapprox(t_comp, FT[0.00, 0.25, 0.05, 0.15, 0.15, 0.30, 0.00]; rtol = rtol)

            comp_no_fb = CompositeLeaf(CompositeLeaf(tab1), tab2)
            @test length(comp_no_fb.models) == 2
            @test_throws ArgumentError leaf_spectrum(comp_no_fb, query_edges)

            # Bounded 4-region broadband CompositeLeaf across 10–50000 cm⁻¹ with domain edges inserted
            leaf_pro = LeafProspectProProperties{FT}(Ccab = FT(40), Ccar = FT(10), Cw = FT(0.012))
            dom_pro = CanopyOptics.spectral_domain(leaf_pro)
            ν_pro_lo, ν_pro_hi = FT(dom_pro[1]), FT(dom_pro[2])
            edges_400 = wavelength_to_wavenumber(FT(400):one(FT):FT(401))
            r400, t400 = leaf_spectrum(leaf_pro, edges_400)
            uv_leaf   = @inferred ConstantLeaf(r400[end], t400[end]; ν_min = ν_pro_hi)
            mwir_leaf = @inferred ConstantLeaf(FT(0.04), zero(FT); ν_min = FT(2500), ν_max = ν_pro_lo)
            tir_leaf  = @inferred ConstantLeaf(FT(0.03), zero(FT); ν_min = zero(FT), ν_max = FT(2500))
            full_leaf = @inferred CompositeLeaf(leaf_pro, uv_leaf, mwir_leaf, tir_leaf)

            broad_edges = sort!(unique!(vcat(
                collect(FT(10):FT(25):FT(50000)),
                FT[2500, ν_pro_lo, ν_pro_hi],
            )))
            r_full, t_full = @inferred leaf_spectrum(full_leaf, broad_edges)
            @test length(r_full) == length(broad_edges) - 1
            @test all(isfinite, r_full) && all(isfinite, t_full)
            @test all(zero(FT) .<= r_full .+ t_full .<= one(FT))

            # Check each of the 4 physical domains explicitly
            idx_tir  = findfirst(i -> broad_edges[i + 1] <= FT(2500), 1:length(r_full))
            idx_mwir = findfirst(i -> broad_edges[i] >= FT(2500) && broad_edges[i + 1] <= ν_pro_lo, 1:length(r_full))
            idx_uv   = findfirst(i -> broad_edges[i] >= ν_pro_hi, 1:length(r_full))
            @test r_full[idx_tir] == FT(0.03) && iszero(t_full[idx_tir])
            @test r_full[idx_mwir] == FT(0.04) && iszero(t_full[idx_mwir])
            @test r_full[idx_uv] == r400[end] && t_full[idx_uv] == t400[end]
            # A bin straddling ν_pro_hi without boundary insertion must throw rather than leak SWIR 0.04
            @test_throws ArgumentError leaf_spectrum(full_leaf, FT[24900, 26000])

            fine_swir_edges = collect(FT(4000):one(FT):FT(4100))
            r_swir_fine, t_swir_fine = @inferred leaf_spectrum(leaf_pro, fine_swir_edges)
            @test all(isfinite, r_swir_fine) && all(isfinite, t_swir_fine)

            custom_comp = CompositeLeaf(WindowConstantLeaf(FT(2000), FT(4000), FT(0.22), FT(0.11)), c_fallback)
            rc, tc = leaf_spectrum(custom_comp, FT[1000, 2000, 4000, 5000])
            @test isapprox(rc, FT[0.03, 0.22, 0.03]; rtol = rtol)
            @test isapprox(tc, FT[0.00, 0.11, 0.00]; rtol = rtol)

            # Custom leaf probe must NOT swallow genuine ArgumentErrors from bad coefficients
            @test_throws ArgumentError leaf_spectrum(CompositeLeaf(BadCoeffCustomLeaf(), c_fallback), FT[1500, 2500])

            # Custom leaf returning wrong length throws DimensionMismatch rather than out-of-bounds read
            @test_throws DimensionMismatch leaf_spectrum(CompositeLeaf(BadLengthLeaf(), c_fallback), FT[1500, 2500, 3500])

            # PROSPECT pigment absorption shares (with independent chlorophyll & carotenoid assertions!)
            opt_native = LeafOpticalProperties{FT}()
            edges_pro, α_direct = wavelength_to_wavenumber(
                FT(400):one(FT):FT(2500),
                prospect_budget(leaf_pro, opt_native).absorptance,
            )
            shares = @inferred leaf_pigment_shares(leaf_pro, edges_pro)
            @test eltype(shares.chlorophyll) === FT && eltype(shares.carotenoids) === FT
            @test length(shares.chlorophyll) == 2101 && length(shares.carotenoids) == 2101
            @test all(isfinite, shares.chlorophyll) && all(isfinite, shares.carotenoids)
            @test all(zero(FT) .<= shares.chlorophyll .<= one(FT))
            @test all(zero(FT) .<= shares.carotenoids .<= one(FT))
            @test all(shares.chlorophyll .+ shares.carotenoids .<= one(FT) .+ FT(16) * eps(FT))
            @test iszero(shares.chlorophyll[1]) && iszero(shares.carotenoids[1])

            # Independent constituent-ratio check at 500 nm (rev index 2101 + 1 - 101 = 2001) and 675 nm
            prods = prospect_constituent_coefficients(leaf_pro, opt_native)
            k_nm_500 = 101 # 400 + 100 = 500 nm
            idx_nu_500 = 2101 + 1 - k_nm_500
            tot_k_500 = prods.chlorophyll[k_nm_500] + prods.carotenoids[k_nm_500] +
                        prods.anthocyanins[k_nm_500] + prods.brown_material[k_nm_500] +
                        prods.water[k_nm_500] + prods.dry_matter[k_nm_500] +
                        prods.protein[k_nm_500] + prods.carbon_constituents[k_nm_500]
            exp_chl_500 = prods.chlorophyll[k_nm_500] / tot_k_500
            exp_car_500 = prods.carotenoids[k_nm_500] / tot_k_500
            @test shares.carotenoids[idx_nu_500] > FT(0.5)
            @test isapprox(shares.chlorophyll[idx_nu_500], exp_chl_500; rtol = rtol)
            @test isapprox(shares.carotenoids[idx_nu_500], exp_car_500; rtol = rtol)

            idx_675 = findfirst(k -> edges_pro[k] <= FT(1e7 / 675) <= edges_pro[k + 1], 1:2101)
            @test shares.chlorophyll[idx_675] > FT(0.9)

            # Coarse-band pigment share = weighted mean of (α * f) / weighted mean of α,
            # verified against independent constituent optical products (not leaf_pigment_shares!)
            tot_k_rev = reverse(
                prods.chlorophyll .+ prods.carotenoids .+ prods.anthocyanins .+
                prods.brown_material .+ prods.water .+ prods.dry_matter .+
                prods.protein .+ prods.carbon_constituents,
            )
            f_chl_indep = ifelse.(tot_k_rev .> zero(FT), reverse(prods.chlorophyll) ./ tot_k_rev, zero(FT))
            f_car_indep = ifelse.(tot_k_rev .> zero(FT), reverse(prods.carotenoids) ./ tot_k_rev, zero(FT))

            coarse_par_edges = FT[14000, 17000, 20000, 25000]
            coarse_as_bands = vcat(coarse_par_edges[1:end-1]', coarse_par_edges[2:end]')
            shares_c = @inferred leaf_pigment_shares(leaf_pro, coarse_par_edges)
            α_band_flat = band_average(α_direct, edges_pro, coarse_as_bands, FlatWeight())
            f_chl_expected = band_average(α_direct .* f_chl_indep, edges_pro, coarse_as_bands, FlatWeight()) ./ α_band_flat
            f_car_expected = band_average(α_direct .* f_car_indep, edges_pro, coarse_as_bands, FlatWeight()) ./ α_band_flat
            @test all(shares_c.carotenoids[2:3] .> FT(0.15))
            @test isapprox(shares_c.chlorophyll, f_chl_expected; rtol = rtol)
            @test isapprox(shares_c.carotenoids, f_car_expected; rtol = rtol)

            # Zero-absorption convention: completely transparent leaf returns 0.0 shares (no NaN)
            leaf_zero = LeafProspectProProperties{FT}(
                Ccab = zero(FT), Ccar = zero(FT), Canth = zero(FT), Cbrown = zero(FT),
                Cw = zero(FT), Cm = zero(FT), Cprot = zero(FT), Ccbc = zero(FT),
            )
            sh_zero = @inferred leaf_pigment_shares(leaf_zero, coarse_par_edges)
            @test all(iszero, sh_zero.chlorophyll) && all(iszero, sh_zero.carotenoids)

            @test_throws ArgumentError leaf_spectrum(leaf_pro, FT[3000, 5000])
            @test_throws ArgumentError leaf_pigment_shares(leaf_pro, FT[3000, 5000])
        end

        # Weak-absorption Float64 test: 1 - rho - tau rounds to 0 while direct prospect_budget.absorptance
        # (~10⁻¹⁷) preserves exact pigment shares on both resampled and weighted coarse bands
        leaf_weak = LeafProspectProProperties(
            Ccab = 1.0e-16, Ccar = 1.0e-16, Canth = 0.0, Cbrown = 0.0,
            Cw = 0.0, Cm = 0.0, Cprot = 0.0, Ccbc = 0.0,
        )
        opt64 = LeafOpticalProperties{Float64}()
        bud_weak = prospect_budget(leaf_weak, opt64)
        @test 0.0 < bud_weak.absorptance[276] < 5.0e-17 # 675 nm
        @test iszero(1.0 - bud_weak.reflectance[276] - bud_weak.transmission[276])
        edges_w64, α_weak_nu = wavelength_to_wavenumber(400.0:1.0:2500.0, bud_weak.absorptance)
        sh_weak_fine = leaf_pigment_shares(leaf_weak, edges_w64)
        sh_weak_coarse = leaf_pigment_shares(leaf_weak, [14500.0, 15000.0, 20500.0])
        @test sh_weak_coarse.chlorophyll[1] > 0.99
        @test sh_weak_coarse.carotenoids[2] > 0.25
        sol_vis = ReferenceSolar(400.0:1.0:800.0, fill(1.5, 401))
        weak_bands = [14500.0 18500.0; 15000.0 20500.0]
        f_car_weak_sol = band_average(α_weak_nu .* sh_weak_fine.carotenoids, edges_w64, weak_bands, sol_vis) ./
                         band_average(α_weak_nu, edges_w64, weak_bands, sol_vis)
        @test all(isfinite, f_car_weak_sol) && f_car_weak_sol[2] > 0.35

        # Photon-count pigment share recipe check: exact 1/ν integral on broad bin, zero-absorptance guard, and 550 nm conversion
        n_avogadro = 6.02214076e23
        umol_per_joule_550 = 1.0e6 / (n_avogadro * 6.62607015e-34 * 299792458.0 / 550.0e-9)
        @test isapprox(umol_per_joule_550, 4.59764; rtol = 5e-4)

        photon_shares_recipe = (α_in, sh_in, edges_in, bands_in, w_in) -> begin
            h_planck, c_light = 6.62607015e-34, 299792458.0
            w_lo, w_hi = first(w_in.edges), last(w_in.edges)
            ν_u = sort!(unique!(vcat(
                w_in.edges,
                filter(ν -> w_lo < ν < w_hi, edges_in),
                filter(ν -> w_lo < ν < w_hi, vec(bands_in)),
            )))
            w_u = band_average(w_in.w, w_in.edges, vcat(ν_u[1:end-1]', ν_u[2:end]'), FlatWeight())
            q_u = [
                iszero(wk) ? zero(wk) :
                (wk / (100.0 * h_planck * c_light)) * log(ν_u[k + 1] / ν_u[k]) / (ν_u[k + 1] - ν_u[k])
                for (k, wk) in enumerate(w_u)
            ]
            w_photon     = TabulatedWeight(ν_u, q_u)
            α_photon     = band_average(α_in, edges_in, bands_in, w_photon)
            α_chl_photon = band_average(α_in .* sh_in.chlorophyll, edges_in, bands_in, w_photon)
            α_car_photon = band_average(α_in .* sh_in.carotenoids, edges_in, bands_in, w_photon)
            f_chl_photon = ifelse.(α_photon .> 0, α_chl_photon ./ α_photon, zero.(α_photon))
            f_car_photon = ifelse.(α_photon .> 0, α_car_photon ./ α_photon, zero.(α_photon))
            return f_chl_photon, f_car_photon
        end

        # Broad single-bin weight [12500, 25000] cm⁻¹ with step pigment share at midpoint 18750 cm⁻¹:
        # exact ∫ (1/ν) dν gives log(1.5)/log(2) ≈ 0.5849625 (not 0.500!)
        w_broad = TabulatedWeight([12500.0, 25000.0], [1.0])
        edges_broad = [12500.0, 18750.0, 25000.0]
        sh_broad = (chlorophyll = [1.0, 0.0], carotenoids = [0.0, 1.0])
        band_broad = [12500.0; 25000.0;;]
        f_chl_ph_b, f_car_ph_b = photon_shares_recipe([1.0, 1.0], sh_broad, edges_broad, band_broad, w_broad)
        @test isapprox(only(f_chl_ph_b), log(1.5) / log(2.0); rtol = 1e-13)
        @test isapprox(only(f_car_ph_b), log(4.0 / 3.0) / log(2.0); rtol = 1e-13)

        # Zero absorptance returns 0.0 photon shares (no NaN)
        f_chl_ph_z, f_car_ph_z = photon_shares_recipe([0.0, 0.0], (chlorophyll = [0.0, 0.0], carotenoids = [0.0, 0.0]), edges_broad, band_broad, w_broad)
        @test iszero(only(f_chl_ph_z)) && iszero(only(f_car_ph_z))

        # Realistic PROSPECT leaf across PAR 400–700 nm
        leaf_ph = LeafProspectProProperties()
        edges_ph, α_ph = wavelength_to_wavenumber(400.0:1.0:2500.0, prospect_budget(leaf_ph, opt64).absorptance)
        sh_ph = leaf_pigment_shares(leaf_ph, edges_ph)
        par_400_700 = [1.0e7 / 700.0; 1.0e7 / 400.0;;]
        f_chl_par, f_car_par = photon_shares_recipe(α_ph, sh_ph, edges_ph, par_400_700, sol_vis)
        @test 0.0 < only(f_chl_par) < 1.0 && 0.0 < only(f_car_par) < 1.0
        @test only(f_chl_par) + only(f_car_par) <= 1.0 + 16 * eps(Float64)
    end

    @testset "Allocation-free mutators, generic axes, hash invariants, and ForwardDiff" begin
        for FT in (Float32, Float64)
            edges = collect(FT(1000):FT(5):FT(5000))
            vals = FT[0.2 + 0.1 * sin(FT(0.001) * edges[i]) for i in 1:(length(edges) - 1)]
            bands = FT[1200 2000 3100 4200; 1900 3050 4150 4950]
            sub_b = FT[1500 2500 3500; 2500 3500 4500]
            par_b = FT[1500; 4500;;]
            dest = zeros(FT, size(bands, 2))
            w_dest = zeros(FT, 3)
            tw = TabulatedWeight(collect(FT(1500):FT(15):FT(4500)), fill(FT(1.3), length(FT(1500):FT(15):FT(4500)) - 1))

            for wt in (FlatWeight(), PlanckWeight(FT(300)), tw)
                ret = @inferred band_average!(dest, vals, edges, bands, wt)
                @test ret === dest
                @test check_band_average_zero_alloc!(dest, vals, edges, bands, wt) == 0

                ret_sub = @inferred band_average!(w_dest, vals, edges, sub_b, wt; parent_band = par_b)
                @test ret_sub === w_dest
                @test check_band_average_zero_alloc!(w_dest, vals, edges, sub_b, wt; parent_band = par_b) == 0

                sh_ret = @inferred subband_weights!(w_dest, sub_b, par_b, wt)
                @test sh_ret === w_dest
                @test check_subband_weights_zero_alloc!(w_dest, sub_b, par_b, wt) == 0

                bw_ret = @inferred CanopyOptics.band_weight!(w_dest, sub_b, wt)
                @test bw_ret === w_dest
                @test check_band_weight_zero_alloc!(w_dest, sub_b, wt) == 0
            end

            # Generic indexing: 0-based vectors, @view slices, and offset band matrix
            vals_z = ShiftedVec(vals, 0)
            edges_z = ShiftedVec(edges, -5)
            bands_o = OffsetBandMat(bands)
            res_ref = band_average(vals, edges, bands, tw)
            res_gen = band_average(vals_z, edges_z, bands_o, tw)
            @test axes(res_gen) == (axes(bands_o, 2),)
            @test collect(res_gen) == res_ref
            @test band_average(@view(vals[1:end]), @view(edges[1:end]), @view(bands[:, :]), PlanckWeight(FT(300))) ==
                  band_average(vals, edges, bands, PlanckWeight(FT(300)))
        end

        # Base.:(==) and Base.hash invariant (including -0.0 vs +0.0 and mixed FT)
        eq_pairs = [
            (FlatWeight(), FlatWeight()),
            (PlanckWeight(300.0), PlanckWeight(300.0f0)),
            (TabulatedWeight([0.0, 200.0], [0.0]), TabulatedWeight([-0.0, 200.0], [-0.0])),
            (TabulatedLeaf([0.0, 200.0], [0.2], [0.0]), TabulatedLeaf([-0.0, 200.0], [0.2], [-0.0])),
            (ConstantLeaf(0.04, 0.0), ConstantLeaf(0.04, -0.0)),
            (ConstantLeaf(0.25, 0.0), ConstantLeaf(0.25f0, -0.0f0)),
            (CompositeLeaf(ConstantLeaf(0.25, 0.0)), CompositeLeaf(ConstantLeaf(0.25f0, -0.0f0))),
        ]
        for (a, b) in eq_pairs
            @test a == b
            @test hash(a) == hash(b)
            @test length(Set([a, b])) == 1
        end

        # ForwardDiff differentiability through PROSPECT Cw -> leaf_spectrum -> band_average
        edges_swir = wavelength_to_wavenumber(1400.0:5.0:1800.0)
        band_swir = [5600.0; 7000.0;;]
        f_cw(cw) = begin
            r, _ = leaf_spectrum(LeafProspectProProperties(Cw = cw), edges_swir)
            only(band_average(r, edges_swir, band_swir, FlatWeight()))
        end
        d_cw = ForwardDiff.derivative(f_cw, 0.012)
        @test isfinite(d_cw) && d_cw < 0.0
        @test isapprox(d_cw, (f_cw(0.012 + 1e-6) - f_cw(0.012 - 1e-6)) / 2e-6; rtol = 1e-5)

        # ForwardDiff w.r.t. PlanckWeight temperature T (both thermal and Wien-shifted UV)
        edges_th = [200.0, 800.0, 1600.0]
        vals_th = [0.92, 0.98]
        f_T(T) = only(band_average(vals_th, edges_th, [200.0; 1600.0;;], PlanckWeight(T)))
        d_T = ForwardDiff.derivative(f_T, 300.0)
        @test isfinite(d_T) && d_T > 0.0
        @test isapprox(d_T, (f_T(300.0 + 1e-4) - f_T(300.0 - 1e-4)) / 2e-4; rtol = 1e-6)

        edges_wien_ad = [12000.0, 15000.0, 18000.0]
        vals_wien_ad = [0.10, 0.40]
        f_T_wien(T) = only(band_average(vals_wien_ad, edges_wien_ad, [12000.0; 18000.0;;], PlanckWeight(T)))
        d_T_wien = ForwardDiff.derivative(f_T_wien, 300.0)
        @test isfinite(d_T_wien) && d_T_wien > 0.0
        @test isapprox(d_T_wien, (f_T_wien(300.0 + 1e-3) - f_T_wien(300.0 - 1e-3)) / 2e-3; rtol = 1e-6)

        # ForwardDiff at exact TabulatedWeight & fine-bin knot boundaries: verify BOTH value AND non-zero tangent!
        tw_knot = TabulatedWeight([1000.0, 1400.0, 2000.0], [1.0, 3.0])
        edges_k = [1000.0, 1500.0, 2000.0]
        vals_k = [0.2, 0.8]
        f_hi_tw(hi) = only(band_average(vals_k, edges_k, [1200.0; hi;;], tw_knot))
        d_hi = ForwardDiff.Dual(2000.0, -1.0)
        got_dual = f_hi_tw(d_hi)
        @test ForwardDiff.value(got_dual) ≈ f_hi_tw(2000.0)
        # Analytic left-derivative at hi = 2000: (2.4*2000 - 1300*3.0) / 2000^2 = 900 / 4e6 = 2.25e-4,
        # so with Dual step -1.0, partial == -2.25e-4 (non-zero; fails if endpoint tangent is dropped!)
        @test isapprox(ForwardDiff.partials(got_dual, 1), -2.25e-4; rtol = 1e-13)
        @test isapprox(ForwardDiff.partials(got_dual, 1), (f_hi_tw(2000.0 - 1e-6) - f_hi_tw(2000.0)) / 1e-6; rtol = 1e-5)

        f_lo_tw(lo) = only(band_average(vals_k, edges_k, [lo; 1800.0;;], tw_knot))
        d_lo = ForwardDiff.Dual(1000.0, 1.0)
        got_dual_lo = f_lo_tw(d_lo)
        @test ForwardDiff.value(got_dual_lo) ≈ f_lo_tw(1000.0)
        @test ForwardDiff.partials(got_dual_lo, 1) > 0.0
        @test isapprox(ForwardDiff.partials(got_dual_lo, 1), (f_lo_tw(1000.0 + 1e-6) - f_lo_tw(1000.0)) / 1e-6; rtol = 1e-5)

        f_lo_knot(lo) = only(band_average(vals_k, edges_k, [lo; 1800.0;;], FlatWeight()))
        d_lo_knot = ForwardDiff.derivative(f_lo_knot, 1000.0)
        @test isapprox(d_lo_knot, (f_lo_knot(1000.0 + 1e-6) - f_lo_knot(1000.0)) / 1e-6; rtol = 1e-5)
    end
end
