using Aqua
using CanopyOptics
using Documenter
using ExplicitImports
using ForwardDiff
using LinearAlgebra
using OffsetArrays
using Test
using Unitful

@testset "Package QA, Generic Indexing, Stability & Equivalence" begin
    @testset "Aqua.jl hygiene" begin
        Aqua.test_all(CanopyOptics)
        ext = Base.get_extension(CanopyOptics, :CanopyOpticsCUDAExt)
        if ext !== nothing
            @test isempty(Test.detect_ambiguities(CanopyOptics, ext; recursive = true))
        end
    end

    if VERSION >= v"1.11"
        @testset "ExplicitImports.jl hygiene" begin
            @test check_no_implicit_imports(CanopyOptics) === nothing
            # LazyArtifacts must be bound in CanopyOptics module scope so @artifact_str
            # can trigger on-demand artifact downloads on clean depots.
            @test check_no_stale_explicit_imports(CanopyOptics; ignore = (:LazyArtifacts,)) === nothing
            @test check_all_explicit_imports_via_owners(CanopyOptics) === nothing
            @test check_all_qualified_accesses_via_owners(CanopyOptics) === nothing
            @test check_no_self_qualified_accesses(CanopyOptics) === nothing
        end

        @testset "Public & exported docstring completeness" begin
            @test isempty(Docs.undocumented_names(CanopyOptics))
        end
    end

    @testset "Documenter.jl doctests" begin
        DocMeta.setdocmeta!(CanopyOptics, :DocTestSetup, :(using CanopyOptics); recursive = true)
        doctest(CanopyOptics; manual = false)
    end

    @testset "Type stability (@inferred), Float32 dielectric, and AD" begin
        for (mat, T, f) in (
            (PureIce(), 260.0f0, 5.0f0),
            (LiquidPureWater(), 285.0f0, 5.0f0),
            (LiquidSaltWater(35.0f0), 285.0f0, 5.0f0),
            (SoilMW(0.4f0, 0.3f0, 0.2f0, 1.5f0), 290.0f0, 1.4f0),
            (LeafUlabyElRayes1987(0.6f0, 0.02f0), 295.0f0, 5.0f0),
        )
            ε32 = @inferred dielectric(mat, T, f)
            @test ε32 isa ComplexF32
            ε64 = @inferred dielectric(mat, Float64(T), Float64(f))
            @test ε64 isa ComplexF64
            @test ε32 ≈ ComplexF32(ε64) rtol = 5e-5

            # AD wrt T and f
            dT_re = ForwardDiff.derivative(t -> real(dielectric(mat, t, Float64(f))), Float64(T))
            df_im = ForwardDiff.derivative(ν -> imag(dielectric(mat, Float64(T), ν)), Float64(f))
            @test isfinite(dT_re) && isfinite(df_im)
        end

        # AD wrt material struct parameters
        @test isfinite(ForwardDiff.derivative(s -> imag(dielectric(LiquidSaltWater(s), 285.0, 1.4)), 35.0))
        @test isfinite(ForwardDiff.derivative(mv -> real(dielectric(SoilMW(0.4, 0.3, mv, 1.5), 290.0, 1.4)), 0.2))
        @test isfinite(ForwardDiff.derivative(mg -> real(dielectric(LeafUlabyElRayes1987(mg, 0.02), 295.0, 5.0)), 0.6))

        # Concrete LeafDistribution & rectangular Γ
        lad = @inferred planophile_leaves(Float64)
        @test isconcretetype(typeof(lad))
        @inferred ((l) -> l.LD)(lad)

        μ16, _ = CanopyOptics.gauleg(16, 0.0, 1.0)
        μ8, _  = CanopyOptics.gauleg(8, 0.0, 1.0)
        g16 = @inferred G(μ16, lad)
        @test eltype(g16) === Float64
        Γ_rect = @inferred CanopyOptics.compute_lambertian_Γ(μ16, μ8, 0.42, 0.40, lad)
        @test size(Γ_rect) == (16, 8)

        # Polar endpoint & cancellation guards
        @test CanopyOptics.H(1.0, 0.0) == 0.0
        @test CanopyOptics.H(-1.0, 0.0) == 0.0
        @test (@inferred CanopyOptics.H(1.0f0, 0.0f0)) === 0.0f0
        dv_clamp = DirectionCosine(1.0 + 1e-15, 0.25)
        @test CanopyOptics._unit_sin_from_cos(dv_clamp.μ) == 0.0
        @test isfinite(dot(dv_clamp, DirectionCosine(0.5, 0.8)))
        @test isfinite(ForwardDiff.derivative(μ -> dot(DirectionCosine(μ, 0.25), DirectionCosine(0.5, 0.8)), 0.75))

        # DirectionCosine constructors, direction_from_angles, and deprecations
        @test DirectionCosine(1, 0) === DirectionCosine{Float64}(1.0, 0.0)
        @test DirectionCosine(μ = 0.6f0, ϕ = 0.2f0) === DirectionCosine{Float32}(0.6f0, 0.2f0)
        @test DirectionCosine(θ = 0.5, ϕ = 1.1) === direction_from_angles(0.5, 1.1)
        @test_throws ArgumentError DirectionCosine()
        @test_throws ArgumentError DirectionCosine(μ = 0.5, θ = 0.5)
        @test CanopyOptics.scattering_model_lambertian(
            DirectionCosine(0.8, 0.0),
            DirectionCosine(-0.5, 0.4),
            DirectionCosine(1.0, 0.0),
            0.4,
            0.2,
        ) ≈ 0.4 * 0.5 / π
        @test (@test_deprecated dirVector(0.5, 1.1)) === direction_from_angles(0.5, 1.1)
        @test Base.isdeprecated(CanopyOptics, :dirVector_μ)
        @test getfield(CanopyOptics, :dirVector_μ)(0.8, 0.2) === DirectionCosine(0.8, 0.2)

        # Primary optical constructors, parametric {FT} forms, mixed-real leaf constructors, and deprecations
        bnds_short = (400.0:500.0:2400.0)u"nm"
        opt32 = LeafOpticalProperties{Float32}(bnds_short)
        @test opt32 isa LeafOpticalProperties{Float32}
        @test LeafOpticalProperties(bnds_short; FT = Float32) == opt32
        mp32 = MultiPigmentOpticalProperties{Float32}(bnds_short)
        @test mp32 isa MultiPigmentOpticalProperties{Float32}
        @test MultiPigmentOpticalProperties(bnds_short; FT = Float32) == mp32
        @test (@test_deprecated createLeafOpticalStruct(bnds_short)) == LeafOpticalProperties(bnds_short)
        @test (@test_deprecated createMultiPigmentOpticalStruct(bnds_short)) == MultiPigmentOpticalProperties(bnds_short)
        @test_throws ArgumentError LeafOpticalProperties(Float64[])
        @test_throws ArgumentError MultiPigmentOpticalProperties(Float64[])
        @test_throws ArgumentError load_wax_refractive_index(Float64[])
        @test LeafProspectProProperties(2, 40.0, 10, 0.5, 0, 0.012, 0, 0.001, 0.009) isa LeafProspectProProperties{Float64}
        @test LeafMultiPigmentProperties(2, 30.0, 10, 10, 0.5, 0, 0.012, 0, 0.001, 0.009) isa LeafMultiPigmentProperties{Float64}
        @test LeafSIPProperties(40, 0, 0, 10.0, 0.5, 0, 0.012, 0.008, 0, 0, 600, 2765.0, 102.8) isa LeafSIPProperties{Float64}

        # Mixed-precision DiscreteAngularNormalization{Float32} with Float64 / Dual BiLambertianCanopyScattering
        μ32, w32 = CanopyOptics.gauleg(6, 0.0f0, 1.0f0)
        norm32 = DiscreteAngularNormalization(μ32, w32; atol = 1.0f-5)
        lad32 = planophile_leaves(Float32)
        Zpp_mp, Zmp_mp = compute_Z_matrices(
            BiLambertianCanopyScattering(0.4, 0.3), μ32, lad32, 0:2;
            angular_normalization = norm32,
        )
        @test eltype(Zpp_mp) === Float64 && eltype(Zmp_mp) === Float64
        cols_mp = vec(transpose(Float64.(w32)) * (Zpp_mp[:, :, 1] + Zmp_mp[:, :, 1])) / 2.0
        @test isapprox(cols_mp, ones(Float64, length(μ32)); atol = 1e-13, rtol = 1e-13)
        f_mp(r) = sum(last(compute_Z_matrices(
            BiLambertianCanopyScattering(r, 0.3), μ32, lad32, 0;
            angular_normalization = norm32,
        )))
        dZ_mp = ForwardDiff.derivative(f_mp, 0.4)
        @test isfinite(dZ_mp) && dZ_mp > 0
        @test isapprox(dZ_mp, (f_mp(0.4 + 1e-6) - f_mp(0.4 - 1e-6)) / 2e-6; rtol = 1e-6)

        # Mixed-precision / ForwardDiff.Dual DirectionCosine pairs in compute_Γ*
        @test isfinite(ForwardDiff.derivative(
            μ -> CanopyOptics.compute_Γ(
                BiLambertianCanopyScattering(0.4, 0.3),
                DirectionCosine(μ, 0.0),
                DirectionCosine(0.5f0, 0.4f0),
                lad,
            ),
            0.7,
        ))
        @test isfinite(ForwardDiff.derivative(
            μ -> CanopyOptics.compute_Γ_isotropic(
                BiLambertianCanopyScattering(0.4, 0.3),
                DirectionCosine(μ, 0.0),
                DirectionCosine(0.5f0, 0.4f0),
            ),
            0.7,
        ))
        @test isfinite(ForwardDiff.derivative(
            μ -> compute_Γ_directional(
                DorsiventralBiLambertianCanopyScattering(0.4, 0.3, 0.45, 0.3),
                DirectionCosine(μ, 0.0),
                DirectionCosine(0.5f0, 0.4f0),
                lad,
            ),
            0.7,
        ))
    end

    @testset "Base.:(==) and Base.hash invariant on container structs" begin
        bnds = (400.0:500.0:2400.0)u"nm"
        μ_dn, w_dn = CanopyOptics.gauleg(6, 0.0, 1.0)
        pairs_to_check = [
            (planophile_leaves(Float64), planophile_leaves(Float64)),
            (LeafOpticalProperties(bnds), LeafOpticalProperties(bnds)),
            (MultiPigmentOpticalProperties(bnds), MultiPigmentOpticalProperties(bnds)),
            (
                DiscreteAngularNormalization(vcat(μ_dn, [0.5]), vcat(w_dn, [0.0]); atol = 1e-12),
                DiscreteAngularNormalization(vcat(μ_dn, [0.5]), vcat(w_dn, [-0.0]); atol = 1e-12),
            ),
            (
                prepare_leaf_bsdf(
                    CoupledLeafBSDF(
                        dorsiventral_from_traits(
                            LeafMultiPigmentProperties();
                            surface_adaxial = RoughDielectricSurface(0.12, 1.0),
                            surface_abaxial = RoughDielectricSurface(0.25, 1.0),
                        ),
                    ),
                    MultiPigmentOpticalProperties(bnds),
                ),
                prepare_leaf_bsdf(
                    CoupledLeafBSDF(
                        dorsiventral_from_traits(
                            LeafMultiPigmentProperties();
                            surface_adaxial = RoughDielectricSurface(0.12, 1.0),
                            surface_abaxial = RoughDielectricSurface(0.25, 1.0),
                        ),
                    ),
                    MultiPigmentOpticalProperties(bnds),
                ),
            ),
            (LUTWoodReflectance([400.0, 800.0], [0.1, 0.3]), LUTWoodReflectance([400.0, 800.0], [0.1, 0.3])),
            (PolynomialWoodReflectance([0.2, 0.05]), PolynomialWoodReflectance([0.2, 0.05])),
            (
                FourSAILGeometrySet(spherical_leaves(); sza_deg = [30.0, 45.0], vza_deg = [0.0, 20.0], raa_deg = 90.0),
                FourSAILGeometrySet(spherical_leaves(); sza_deg = [30.0, 45.0], vza_deg = [0.0, 20.0], raa_deg = 90.0),
            ),
            (
                foursail([0.1, 0.4], [0.05, 0.4], 0.15, FourSAILGeometry(spherical_leaves(); sza_deg = 30.0, vza_deg = 10.0, raa_deg = 45.0), 3.0),
                foursail([0.1, 0.4], [0.05, 0.4], 0.15, FourSAILGeometry(spherical_leaves(); sza_deg = 30.0, vza_deg = 10.0, raa_deg = 45.0), 3.0),
            ),
            (TabulatedPhotonSpectrum([650.0, 700.0, 740.0], [0.1, 0.6, 0.3]), TabulatedPhotonSpectrum([650.0, 700.0, 740.0], [0.1, 0.6, 0.3])),
        ]
        for (x, y) in pairs_to_check
            @test x !== y || x isa LeafDistribution
            @test x == y
            @test hash(x) == hash(y)
            @test length(Set([x, y])) == 1
        end
    end

    @testset "OffsetArrays and @view generic indexing contracts" begin
        bnds_vec = collect((400.0:200.0:2400.0)u"nm")
        bnds_off = OffsetArray(bnds_vec, -3)
        bnds_vw  = @view bnds_vec[1:end]

        opti_ref = LeafOpticalProperties(bnds_vec)
        @test LeafOpticalProperties(bnds_off) == opti_ref
        @test LeafOpticalProperties(bnds_vw) == opti_ref

        opti_mp = MultiPigmentOpticalProperties(bnds_vec)
        @test MultiPigmentOpticalProperties(bnds_off) == opti_mp
        @test MultiPigmentOpticalProperties(bnds_vw) == opti_mp

        μ_dn, w_dn = CanopyOptics.gauleg(6, 0.0, 1.0)
        norm_ref = DiscreteAngularNormalization(μ_dn, w_dn; atol = 1e-12)
        norm_off = DiscreteAngularNormalization(OffsetArray(μ_dn, -2), OffsetArray(w_dn, 5); atol = 1e-12)
        @test norm_off == norm_ref
        @test compute_Z_matrices(BiLambertianCanopyScattering(0.4, 0.3), OffsetArray(μ_dn, -4), spherical_leaves(), 0; angular_normalization = norm_ref) ==
              compute_Z_matrices(BiLambertianCanopyScattering(0.4, 0.3), μ_dn, spherical_leaves(), 0; angular_normalization = norm_ref)

        nwax = load_wax_refractive_index(bnds_vec)
        @test load_wax_refractive_index(bnds_off) == nwax
        nwax_off = OffsetArray(nwax, 10)
        nwax_vw  = @view nwax[1:end]

        leaf_mp = LeafMultiPigmentProperties(Cchla = 30.0, Cchlb = 10.0)
        bud_ref = prospect_budget(leaf_mp, opti_mp; n_surface = nwax)
        bud_off = prospect_budget(leaf_mp, opti_mp; n_surface = nwax_off)
        bud_vw  = prospect_budget(leaf_mp, opti_mp; n_surface = nwax_vw)
        @test bud_off.reflectance ≈ bud_ref.reflectance && bud_vw.reflectance ≈ bud_ref.reflectance

        sip = LeafSIPProperties(Cchla = 30.0, Cchlb = 10.0)
        @test leaf_sip_budget(sip, opti_mp; n_surface = nwax_off).reflectance ≈
              leaf_sip_budget(sip, opti_mp; n_surface = nwax).reflectance

        dors = dorsiventral_from_traits(
            leaf_mp;
            surface_adaxial = RoughDielectricSurface(0.12, 1.0),
            surface_abaxial = RoughDielectricSurface(0.25, 1.0),
        )
        dv_ref = dorsiventral_budget(dors, opti_mp; n_surface_adaxial = nwax, n_surface_abaxial = nwax)
        dv_off = dorsiventral_budget(dors, opti_mp; n_surface_adaxial = nwax_off, n_surface_abaxial = nwax_off)
        @test dv_off.adaxial.reflectance ≈ dv_ref.adaxial.reflectance &&
              dv_off.adaxial.transmission ≈ dv_ref.adaxial.transmission

        # Wood reflectance preserves OffsetArray axes
        lut = LUTWoodReflectance(OffsetArray([400.0, 800.0, 1600.0], -2), OffsetArray([0.1, 0.25, 0.4], 5))
        poly = PolynomialWoodReflectance(OffsetArray([0.15, 0.05], -1); x_offset = 400.0, x_scale = 1000.0)
        grid_off = OffsetArray([500.0, 1000.0, 1500.0], -4)
        r_lut_off = wood_reflectance(lut, grid_off)
        r_poly_off = wood_reflectance(poly, grid_off)
        @test axes(r_lut_off) == axes(grid_off)
        @test axes(r_poly_off) == axes(grid_off)
        @test collect(r_lut_off) ≈ wood_reflectance(lut, collect(grid_off))

        # Canopy G, bulk_G, Ψ, Γ, and Z matrices on OffsetVector & @view
        lad = plagiophile_leaves()
        μ_full, _ = CanopyOptics.gauleg(12, 0.0, 1.0)
        μ8 = collect(μ_full[1:8])
        μ8_off = OffsetArray(μ8, -5)
        μ8_vw  = @view μ_full[1:8]

        g_off = G(μ8_off, lad)
        @test axes(g_off) == axes(μ8_off)
        @test collect(g_off) ≈ G(μ8, lad)
        @test G(μ8_vw, lad) ≈ G(μ8, lad)

        tree = TreeCanopy(
            LAI = 3.0,
            SAI = 0.5,
            leaf_kwargs = (; clumping = EmpiricalDirectionalClumping(Ω₀ = 0.75)),
        )
        bg_off = bulk_G(tree, μ8_off; clumped = true)
        @test axes(bg_off) == axes(μ8_off)
        @test collect(bg_off) ≈ bulk_G(tree, μ8; clumped = true)

        μ4 = collect(μ_full[1:4])
        μ4_off = OffsetArray(μ4, 7)
        Ψ_off = CanopyOptics.compute_Ψ(μ8_off, μ4_off, 0.6)
        Ψ_ref = CanopyOptics.compute_Ψ(μ8, μ4, 0.6)
        @test axes(Ψ_off[1]) == (axes(μ8_off, 1), axes(μ4_off, 1))
        @test collect(Ψ_off[1]) ≈ Ψ_ref[1] && collect(Ψ_off[2]) ≈ Ψ_ref[2]
        Γ_off = CanopyOptics.compute_lambertian_Γ(μ8_off, μ4_off, 0.4, 0.35, lad)
        @test axes(Γ_off) == (axes(μ8_off, 1), axes(μ4_off, 1))
        @test collect(Γ_off) ≈ CanopyOptics.compute_lambertian_Γ(μ8, μ4, 0.4, 0.35, lad)

        bilam = BiLambertianCanopyScattering(0.4, 0.35)
        spec  = SpecularCanopyScattering(1.45, 0.2)
        dors_sc = DorsiventralBiLambertianCanopyScattering(0.38, 0.40, 0.44, 0.40)
        @test compute_Z_matrices(bilam, μ8_off, lad, 0:3)[1] ≈ compute_Z_matrices(bilam, μ8, lad, 0:3)[1]
        @test compute_Z_matrices(spec, μ8_off, lad, 1; npol = 4)[1] ≈ compute_Z_matrices(spec, μ8, lad, 1; npol = 4)[1]
        @test compute_Z_matrices(spec, μ8_vw, lad, 1; npol = 4)[2] ≈ compute_Z_matrices(spec, μ8, lad, 1; npol = 4)[2]
        @test compute_Z_matrices(tree, μ8_off, 0:2)[1] ≈ compute_Z_matrices(tree, μ8, 0:2)[1]
        @test compute_Z_matrices_directional(dors_sc, μ8_off, lad, 3).Zmp ≈
              compute_Z_matrices_directional(dors_sc, μ8, lad, 3).Zmp

        prep_bsdf = prepare_leaf_bsdf(CoupledLeafBSDF(dors; η_t = 0.2, γ_t = 3.0), opti_mp; n_surface = nwax_off)
        zc_ref = compute_Z_matrices_coupled_bsdf(prep_bsdf, 1, μ4, lad, 2; quadrature = CanopyQuadrature(n_leaf = 8, n_azimuth = 12))
        zc_off = compute_Z_matrices_coupled_bsdf(prep_bsdf, 1, μ4_off, lad, 2; quadrature = CanopyQuadrature(n_leaf = 8, n_azimuth = 12))
        @test zc_off.Zpp ≈ zc_ref.Zpp && zc_off.Zmp ≈ zc_ref.Zmp

        # 4SAIL supports OffsetArray angle inputs on FourSAILGeometrySet and enforces 1-based indexing on kernel arrays
        gs_off = FourSAILGeometrySet(lad; sza_deg = OffsetArray([25.0, 40.0], -2), vza_deg = OffsetArray([10.0, 30.0], 4), raa_deg = 45.0)
        gs_ref = FourSAILGeometrySet(lad; sza_deg = [25.0, 40.0], vza_deg = [10.0, 30.0], raa_deg = 45.0)
        @test gs_off == gs_ref
        geom1 = gs_ref[1]
        @test_throws ArgumentError foursail(OffsetArray([0.1, 0.4], -1), [0.05, 0.4], 0.15, geom1, 3.0)
        res1 = foursail([0.1, 0.4], [0.05, 0.4], 0.15, geom1, 3.0)
        @test_throws ArgumentError foursail!(res1, [0.1, 0.4], OffsetArray([0.05, 0.4], 2), 0.15, geom1, 3.0)
    end

    @testset "32-node Chebyshev rough-surface accuracy and ForwardDiff agreement" begin
        for σ in (0.02, 0.12, 0.35), μ in (1.0, 0.92, 0.65, 0.30, 0.05)
            surf = RoughDielectricSurface(σ, 1.0)
            cheb_dir = CanopyOptics._prepare_surface_cheb(Float64, surf, μ, false)
            cheb_hem = CanopyOptics._prepare_surface_cheb(Float64, surf, μ, true)
            for n_eff in (1.32, 1.45, 1.54)
                s_exact = surface_interface_optics(surf, μ, n_eff; hemispherical = false)
                s_cheb  = CanopyOptics._eval_prepared_surface(surf, cheb_dir, μ, n_eff, n_eff, false)
                @test isapprox(s_cheb.r_ext, s_exact.r_ext; atol = 1e-14, rtol = 1e-14)
                @test isapprox(s_cheb.t_ext, s_exact.t_ext; atol = 1e-14, rtol = 1e-14)
                @test isapprox(s_cheb.r_int_diffuse, s_exact.r_int_diffuse; atol = 1e-14, rtol = 1e-14)

                sh_exact = surface_interface_optics(surf, μ, n_eff; hemispherical = true)
                sh_cheb  = CanopyOptics._eval_prepared_surface(surf, cheb_hem, μ, n_eff, n_eff, true)
                @test isapprox(sh_cheb.r_ext, sh_exact.r_ext; atol = 1e-14, rtol = 1e-14)
            end
        end

        # Verify ForwardDiff gradients through dorsiventral_budget & leaf_bsdf_spectrum (including nadir μ = 1.0)
        bnds = (400.0:500.0:2400.0)u"nm"
        opti_mp = MultiPigmentOpticalProperties(bnds)
        f_dv = function (θ)
            lf = LeafMultiPigmentProperties(Cchla = θ[1], Cchlb = 10.0)
            dv = dorsiventral_from_traits(
                lf;
                N_tot = θ[2],
                surface_adaxial = RoughDielectricSurface(θ[3], θ[4]),
                surface_abaxial = RoughDielectricSurface(0.25, 1.0),
            )
            b = dorsiventral_budget(dv, opti_mp; μ_ad = θ[5], hemispherical = false)
            return sum(b.adaxial.reflectance) + sum(b.adaxial.transmission)
        end
        x0 = [30.0, 2.1, 0.14, 1.02, 1.0] # exact nadir μ_ad = 1.0!
        g_ad = ForwardDiff.gradient(f_dv, x0)
        @test all(isfinite, g_ad)
        h_fd = 1e-6
        for k in eachindex(x0)
            xp = copy(x0); xp[k] += h_fd
            xm = copy(x0); xm[k] -= h_fd
            # μ_ad ∈ (0, 1] is at the upper boundary 1.0, so use backward step for k == 5
            fd_k = k == 5 ? (f_dv(x0) - f_dv(xm)) / h_fd : (f_dv(xp) - f_dv(xm)) / (2 * h_fd)
            @test isapprox(g_ad[k], fd_k; rtol = 1e-4, atol = 1e-6)
        end

        # ForwardDiff through compute_Z_matrices_coupled_bsdf at symmetric [σ_ad, σ_ab] = [0.15, 0.15],
        # η_t = 0.25 > 0, and μ containing exact nadir 1.0 (verifies zero Dual partial aliasing & zero NaN)
        lad_z = planophile_leaves(Float64)
        μ_nadir = [0.5, 1.0]
        q_z = CanopyQuadrature(n_leaf = 6, n_azimuth = 8)
        f_zc = function (σ_pair)
            dv_s = dorsiventral_from_traits(
                LeafMultiPigmentProperties(Cchla = 30.0, Cchlb = 10.0);
                N_tot = 2.0,
                surface_adaxial = RoughDielectricSurface(σ_pair[1], 1.0),
                surface_abaxial = RoughDielectricSurface(σ_pair[2], 1.0),
            )
            p_s = prepare_leaf_bsdf(CoupledLeafBSDF(dv_s; η_t = 0.25, γ_t = 3.0), opti_mp)
            z_s = compute_Z_matrices_coupled_bsdf(p_s, 1, μ_nadir, lad_z, 2; quadrature = q_z)
            return sum(z_s.Zmp) + 0.7 * sum(z_s.Zpm) + sum(z_s.ω_down)
        end
        σ0 = [0.15, 0.15]
        g_zc = ForwardDiff.gradient(f_zc, σ0)
        @test all(isfinite, g_zc)
        for k in 1:2
            sp = copy(σ0); sp[k] += h_fd
            sm = copy(σ0); sm[k] -= h_fd
            fd_zc = (f_zc(sp) - f_zc(sm)) / (2 * h_fd)
            @test isapprox(g_zc[k], fd_zc; rtol = 1e-4, atol = 1e-6)
        end

        # ForwardDiff at exact nadir (μ = 1 / 0°) & hotspot (sza == vza, raa == 0) in Z, hotspot, and 4SAIL
        g_zb = ForwardDiff.gradient(
            rt -> begin
                Zpp, Zmp = compute_Z_matrices(BiLambertianCanopyScattering(rt[1], rt[2]), μ_nadir, lad_z, 0:2)
                sum(Zpp) + sum(Zmp)
            end,
            [0.42, 0.38],
        )
        @test all(isfinite, g_zb)

        @test isfinite(ForwardDiff.derivative(h -> hotspot_separation(1.0, 1.0, 0.0, h), 0.1))
        @test isfinite(ForwardDiff.derivative(h -> hotspot_separation(0.8, 0.8, 0.0, h), 0.1))
        g_sail = ForwardDiff.gradient(
            x -> begin
                gm = FourSAILGeometry(lad_z; sza_deg = x[1], vza_deg = x[2], raa_deg = 0.0)
                rs = foursail([0.1, 0.45], [0.05, 0.42], 0.15, gm, x[3]; hotspot = x[4])
                sum(rs.rsot) + sum(rs.rsdt)
            end,
            [0.0, 0.0, 3.0, 0.05],
        )
        @test all(isfinite, g_sail)
        g_sail_hs = ForwardDiff.gradient(
            x -> begin
                gm = FourSAILGeometry(lad_z; sza_deg = 30.0, vza_deg = 30.0, raa_deg = 0.0)
                rs = foursail([0.1, 0.45], [0.05, 0.42], 0.15, gm, x[1]; hotspot = x[2])
                sum(rs.rsot)
            end,
            [3.0, 0.05],
        )
        @test all(isfinite, g_sail_hs)

        # MultiPigmentOpticalProperties from OffsetArray spectral vectors, dry SoilMW mᵥ = 0, & mixed-FT LeafDistribution hash
        mp_off = MultiPigmentOpticalProperties(LeafOpticalProperties(bnds); Kchla = OffsetArray(opti_mp.Kchla, -3), nᵣ = OffsetArray(opti_mp.nᵣ, 4))
        @test mp_off.Kchla == opti_mp.Kchla && mp_off.nᵣ == opti_mp.nᵣ
        ε_dry = dielectric(SoilMW(0.4, 0.3, 0.0, 1.5), 290.0, 1.4)
        @test isfinite(real(ε_dry)) && iszero(imag(ε_dry))
        lad32 = LeafDistribution(CanopyOptics.Uniform(0.0f0, 1.0f0), Float32(0.5))
        lad64 = LeafDistribution(CanopyOptics.Uniform(-0.0, 1.0), 0.5)
        @test lad32 == lad64
        @test hash(lad32) == hash(lad64)
        @test CanopyOptics.K(0.0f0, Float32(π) / 2.0f0) === 1.0f0
        @test CanopyOptics.K(0.2f0, Float32(π) / 2.0f0) === 0.0f0

        # ForwardDiff w.r.t. stream cosines μ (including nadir 1.0 in [μ0, 1.0]) across G and Z kernels
        g_mu_fn = μ0 -> sum(G([μ0, 1.0], lad_z))
        d_g_mu = ForwardDiff.derivative(g_mu_fn, 0.5)
        @test isfinite(d_g_mu)
        @test isapprox(d_g_mu, (g_mu_fn(0.5 + h_fd) - g_mu_fn(0.5 - h_fd)) / (2 * h_fd); rtol = 1e-5, atol = 1e-7)

        for npol_spec in (1, 4)
            spec_mu_fn = μ0 -> begin
                Zpp_s, Zmp_s = compute_Z_matrices(SpecularCanopyScattering(1.45, 0.2), [μ0, 1.0], lad_z, 1; npol = npol_spec)
                sum(Zpp_s) + sum(Zmp_s)
            end
            d_spec_mu = ForwardDiff.derivative(spec_mu_fn, 0.5)
            @test isfinite(d_spec_mu)
            @test isapprox(d_spec_mu, (spec_mu_fn(0.5 + h_fd) - spec_mu_fn(0.5 - h_fd)) / (2 * h_fd); rtol = 1e-5, atol = 1e-7)
        end
    end
end
