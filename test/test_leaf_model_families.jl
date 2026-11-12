using Test
using CanopyOptics
using ForwardDiff
using LinearAlgebra
using Unitful

const _u_nm = u"nm"

@testset "Leaf Model Families & Coupled BSDF / Canopy Extensions" begin
    opti64_full = createLeafOpticalStruct(400.0:100.0:2400.0)
    # Compact Float32 optical fixture without implicit Float64 promotion
    opti32 = LeafOpticalProperties{Float32,Vector{Float32}}(
        Float32.(collect(450.0:200.0:2350.0)),
        Float32.(opti64_full.nᵣ[1:2:end]),
        Float32.(opti64_full.Kcab[1:2:end]),
        Float32.(opti64_full.Kcar[1:2:end]),
        Float32.(opti64_full.Kant[1:2:end]),
        Float32.(opti64_full.Kb[1:2:end]),
        Float32.(opti64_full.Kw[1:2:end]),
        Float32.(opti64_full.Km[1:2:end]),
        Float32.(opti64_full.Kp[1:2:end]),
        Float32.(opti64_full.Kcbc[1:2:end]),
    )

    @testset "1. Surface Optics: calctav_smooth, Microfacet Quadrature, Étendue & AD" begin
        for FT in (Float32, Float64)
            nr = FT(1.45)
            t0_exact = FT(4) * nr / ((nr + one(FT))^2)
            t0 = calctav_smooth(zero(FT), nr)
            @test t0 isa FT
            @test isapprox(t0, t0_exact; atol = 5 * eps(FT))

            t40 = calctav_smooth(FT(40), nr)
            t90 = calctav_smooth(FT(90), nr)
            @test t40 isa FT && t90 isa FT
            @test isapprox(t40, FT(CanopyOptics.calctav(FT(40), nr)); rtol = 10 * eps(FT))
            @test isapprox(t90, FT(CanopyOptics.calctav(FT(90), nr)); rtol = 10 * eps(FT))

            # ForwardDiff smoothness right at α = 0°, 8°, 40°, and 90°
            for α0 in (FT(0), FT(8), FT(40), FT(90))
                dα = ForwardDiff.derivative(a -> calctav_smooth(a, nr), α0)
                dn = ForwardDiff.derivative(n -> calctav_smooth(α0, n), nr)
                @test isfinite(dα) && isfinite(dn)
            end
            # At α = 0°, d(calctav)/dα == 0 by even symmetry in α
            @test abs(ForwardDiff.derivative(a -> calctav_smooth(a, nr), zero(FT))) <= 10 * eps(FT)

            # Smith G1 nadir behavior & ForwardDiff
            @test smith_beckmann_G1(one(FT), FT(0.2)) === one(FT)
            g_nadir = ForwardDiff.gradient(x -> smith_beckmann_G1(x[1], x[2]), [one(FT), FT(0.2)])
            @test all(isfinite, g_nadir)

            # Étendue invariance across n_scale != 1
            for surf in (ConeSurface{FT}(FT(8), FT(1.08)), RoughDielectricSurface{FT}(FT(0.18), FT(1.08)))
                s = surface_interface_optics(surf, FT(0.8), nr)
                @test s.r_ext isa FT && s.t_int_diffuse isa FT
                @test zero(FT) < s.r_ext < one(FT)
                @test isapprox(s.r_ext + s.t_ext, one(FT); atol = 2 * eps(FT))
                @test isapprox(s.t_int_diffuse, s.t_hem_ext / (nr * nr); atol = 3 * eps(FT))
                @test zero(FT) < s.t_int_diffuse < inv(nr * nr) < one(FT)
            end
        end

        # Domain validation guards & Float32 keyword constructors
        @test ConeSurface(α_deg = 12.0f0) isa ConeSurface{Float32}
        @test RoughDielectricSurface(σ = 0.2f0) isa RoughDielectricSurface{Float32}
        @test calctav_smooth(40, 2) isa Float64
        @test_throws ArgumentError RoughDielectricSurface(0.005, 1.0)
        @test_throws ArgumentError ConeSurface(-1.0, 1.0)
        @test_throws ArgumentError calctav_smooth(10.0, 0.99)
        @test_throws ArgumentError surface_interface_optics(ConeSurface(40.0, 0.6), 1.0, 1.45)
        @test_throws ArgumentError beckmann_D(-0.1, 0.8, 0.2)
        @test_throws ArgumentError beckmann_D(0.1, 0.0, 0.2)
        @test_throws ArgumentError beckmann_D(0.1, 0.8, 0.0)
        @test smith_beckmann_G1(1.0 + 2 * eps(Float64), 0.2) === 1.0
        @test_throws ArgumentError smith_beckmann_G1(0.5, 0.0)
        @test_throws ArgumentError smith_beckmann_G1(1.1, 0.2)
        @test_throws ArgumentError microfacet_directional_reflectance(0.0, 1.45, 0.2)
        @test_throws ArgumentError microfacet_directional_reflectance(0.8, 1.0, 0.2)
        @test_throws ArgumentError microfacet_directional_reflectance(0.8, 1.45, -0.1)
        @test_throws ArgumentError microfacet_hemispherical_reflectance(0.95, 0.2)
        @test_throws ArgumentError microfacet_hemispherical_reflectance(1.45, 0.0)

        # High-order peak-clustered reference quadrature on (μ_o = 1 - ξ², ϕ_o = π - π η²) across all σ ∈ [0.02, 0.50]
        ξ96, wξ96 = CanopyOptics.gauleg(96, 0.0, 1.0)
        η128, wη128 = CanopyOptics.gauleg(128, 0.0, 1.0)
        function ref_hemi_Rsurf(μi, n_eff, σ)
            acc = 0.0
            for km in eachindex(ξ96)
                ξ = ξ96[km]
                μo = max(1e-12, 1.0 - ξ * ξ)
                wm = 2.0 * ξ * wξ96[km]
                for kp in eachindex(η128)
                    η = η128[kp]
                    ϕo = Float64(π) * (1.0 - η * η)
                    wp = 2.0 * (2.0 * Float64(π) * η * wη128[kp]) # factor 2 for [0, 2π] symmetry
                    M = surface_brdf_mueller(μi, 0.0, μo, ϕo, n_eff, σ, 1)
                    acc += wm * wp * M[1, 1] * μo
                end
            end
            return acc
        end

        for σ in (0.02, 0.05, 0.15, 0.30, 0.50), μi in (1.0, 0.7, 0.3, 0.1)
            r_fast = microfacet_directional_reflectance(μi, 1.45, σ)
            r_ref = ref_hemi_Rsurf(μi, 1.45, σ)
            @test 0.0 < r_fast < 1.0
            @test isapprox(r_fast, r_ref; rtol = 2.0e-3)
        end

        # ForwardDiff through microfacet_directional_reflectance at both μ_i = 1.0 and μ_i = 0.5
        g1 = ForwardDiff.gradient(x -> microfacet_directional_reflectance(x[1], x[2], x[3]), [1.0, 1.45, 0.15])
        g2 = ForwardDiff.gradient(x -> microfacet_directional_reflectance(x[1], x[2], x[3]), [0.5, 1.45, 0.15])
        fd_nadir = (microfacet_directional_reflectance(1.0, 1.45, 0.15) - microfacet_directional_reflectance(1.0 - 1e-5, 1.45, 0.15)) / 1e-5
        @test all(isfinite, g1) && all(isfinite, g2)
        @test g1[1] != 0.0
        @test isapprox(g1[1], fd_nadir; rtol = 2e-2)
    end

    @testset "2. Multi-Pigment Basis, Excitation Selection & Leaf-SIP" begin
        for (FT, opti) in ((Float32, opti32), (Float64, opti64_full))
            leaf_pro = LeafProspectProProperties{FT}(
                N = FT(1.6), Ccab = FT(42), Ccar = FT(9), Canth = FT(2),
                Cbrown = FT(0.1), Cw = FT(0.011), Cm = zero(FT), Cprot = FT(0.0012), Ccbc = FT(0.0085),
            )
            leaf_mp = LeafMultiPigmentProperties{FT}(
                N = FT(1.6), Cchla = FT(31.5), Cchlb = FT(10.5), Ccar = FT(9), Canth = FT(2),
                Cbrown = FT(0.1), Cw = FT(0.011), Cm = zero(FT), Cprot = FT(0.0012), Ccbc = FT(0.0085),
            )
            mp_opti = MultiPigmentOpticalProperties(opti)

            b_ref = prospect_budget(leaf_pro, opti)
            b_mp = prospect_budget(leaf_mp, mp_opti)
            @test eltype(b_mp.reflectance) === FT
            tol = FT === Float32 ? FT(5e-6) : FT(1e-14)
            @test maximum(abs.(b_mp.reflectance .- b_ref.reflectance)) <= tol
            @test maximum(abs.(b_mp.transmission .- b_ref.transmission)) <= tol
            @test maximum(abs.(b_mp.absorptance .- b_ref.absorptance)) <= tol

            # Generalized surface prospect_budget conservation (including nonabsorbing k = 0)
            leaf_zero = LeafMultiPigmentProperties{FT}(
                N = FT(1.4), Cchla = zero(FT), Cchlb = zero(FT), Ccar = zero(FT),
                Canth = zero(FT), Cbrown = zero(FT), Cw = zero(FT), Cm = zero(FT),
                Cprot = zero(FT), Ccbc = zero(FT),
            )
            mp_zero = MultiPigmentOpticalProperties(
                opti.λ, opti.nᵣ,
                zeros(FT, length(opti.nᵣ)), zeros(FT, length(opti.nᵣ)),
                zeros(FT, length(opti.nᵣ)), zeros(FT, length(opti.nᵣ)),
                zeros(FT, length(opti.nᵣ)), zeros(FT, length(opti.nᵣ)),
                zeros(FT, length(opti.nᵣ)), zeros(FT, length(opti.nᵣ)),
                zeros(FT, length(opti.nᵣ)),
            )
            for surf in (ConeSurface{FT}(zero(FT), one(FT)), RoughDielectricSurface{FT}(FT(0.15), FT(1.04)))
                bz = prospect_budget(leaf_zero, mp_zero, surf; μ_in = FT(0.75))
                @test all(iszero, bz.absorptance)
                @test maximum(abs.(bz.reflectance .+ bz.transmission .- one(FT))) <= tol

                ba = prospect_budget(leaf_mp, mp_opti, surf; μ_in = FT(0.75))
                @test maximum(abs.(ba.reflectance .+ ba.transmission .+ ba.absorptance .- one(FT))) <= tol
                @test all(ba.reflectance_volume .>= zero(FT))
            end

            # Constituent attribution + MultiChlorophyllAntenna
            depths = constituent_optical_depths(leaf_mp, mp_opti, 2)
            part = partition_absorption(HomogeneousAbsorption(), b_mp.absorptance[2], depths)
            @test absorption_valid(part)
            ant_full = antenna_excitation(MultiChlorophyllAntenna{FT}(one(FT), FT(0.35)), part)
            @test ant_full isa FT && zero(FT) <= ant_full <= b_mp.absorptance[2]

            depths0 = constituent_optical_depths(leaf_zero, mp_zero, 1)
            part0 = partition_absorption(HomogeneousAbsorption(), zero(FT), depths0)
            @test part0.status == AbsorptionNoMaterial
            @test antenna_excitation(MultiChlorophyllAntenna{FT}(), part0) === zero(FT)

            # Leaf-SIP vs Leaf-SIP+Fresnel + prospect tuple wrappers + MultiChlorophyllAntenna Float32
            @test MultiChlorophyllAntenna(0.8f0) isa MultiChlorophyllAntenna{Float32}
            @test MultiChlorophyllAntenna(chlb_transfer = 0.8f0) isa MultiChlorophyllAntenna{Float32}
            t_tup, r_tup = prospect(leaf_mp, mp_opti)
            @test t_tup == b_mp.transmission && r_tup == b_mp.reflectance
            t_surf, r_surf = prospect(leaf_mp, mp_opti, ConeSurface{FT}(FT(8), one(FT)); μ_in = FT(0.9))
            @test eltype(t_surf) === FT && eltype(r_surf) === FT

            # Guard: LeafProspectProProperties + MultiPigmentOpticalProperties throws clean ArgumentError
            @test_throws ArgumentError prospect_budget(leaf_pro, mp_opti)
            @test_throws ArgumentError prospect_budget(leaf_pro, mp_opti, ConeSurface{FT}(FT(40), one(FT)))
            @test_throws ArgumentError prospect(leaf_pro, mp_opti)
            @test_throws ArgumentError prospect_constituent_coefficients(leaf_pro, mp_opti)
            @test_throws ArgumentError constituent_optical_depths(leaf_pro, mp_opti, 1)

            sip = LeafSIPProperties{FT}(Ccab = FT(60), Ccar = FT(12), Cw = FT(0.015), Cm = FT(0.008))
            b_sip = leaf_sip_budget(sip, opti)
            b_sip_fres = leaf_sip_budget(sip, opti; surface = ConeSurface{FT}(FT(8), one(FT)))
            @test eltype(b_sip.reflectance) === FT
            @test maximum(abs.(b_sip.reflectance .+ b_sip.transmission .+ b_sip.absorptance .- one(FT))) <= tol
            @test maximum(abs.(b_sip_fres.reflectance .+ b_sip_fres.transmission .+ b_sip_fres.absorptance .- one(FT))) <= tol
            # Surface Fresnel restores blue/red surface floor above bare SIP
            @test b_sip_fres.reflectance[1] > b_sip.reflectance[1] + FT(0.02)

            # Valid Leaf-SIP separate Cchla + Cchlb + Cprot + Ccbc path (with Ccab = 0, Cm = 0)
            sip_mp = LeafSIPProperties{FT}(
                Ccab = zero(FT), Cchla = FT(45), Cchlb = FT(15), Ccar = FT(12),
                Cw = FT(0.015), Cm = zero(FT), Cprot = FT(0.001), Ccbc = FT(0.007),
            )
            b_sip_mp = leaf_sip_budget(sip_mp, mp_opti)
            @test eltype(b_sip_mp.reflectance) === FT
            @test maximum(abs.(b_sip_mp.reflectance .+ b_sip_mp.transmission .+ b_sip_mp.absorptance .- one(FT))) <= tol

            # Leaf-SIP double-counting / ambiguity guards
            @test_throws ArgumentError leaf_sip_budget(LeafSIPProperties{FT}(Ccab = FT(40), Cchla = FT(10)), opti)
            @test_throws ArgumentError leaf_sip_budget(LeafSIPProperties{FT}(Ccab = FT(40)), mp_opti)
            @test_throws ArgumentError leaf_sip_budget(LeafSIPProperties{FT}(Cm = FT(0.005), Cprot = FT(0.001)), opti)
        end
    end

    @testset "3. Dorsiventral Two-Zone & Single-Slab Exact Conservation, Reduction & Reciprocity" begin
        for (FT, opti) in ((Float32, opti32), (Float64, opti64_full))
            tol = FT === Float32 ? FT(5e-6) : FT(1e-13)

            # 3a. Exact reduction of 2-zone DorsiventralLeafProperties to prospect_budget
            N_pal, N_spo = FT(1.2), FT(1.5)
            N_tot = N_pal + N_spo
            f_pal = N_pal / N_tot
            f_spo = N_spo / N_tot
            leaf_whole = LeafProspectProProperties{FT}(
                N = N_tot, Ccab = FT(45), Ccar = FT(10), Canth = FT(1.5),
                Cbrown = FT(0.05), Cw = FT(0.014), Cm = FT(0.009), Cprot = zero(FT), Ccbc = zero(FT),
            )
            pal = LeafProspectProProperties{FT}(
                N = N_pal, Ccab = FT(45) * f_pal, Ccar = FT(10) * f_pal, Canth = FT(1.5) * f_pal,
                Cbrown = FT(0.05) * f_pal, Cw = FT(0.014) * f_pal, Cm = FT(0.009) * f_pal,
                Cprot = zero(FT), Ccbc = zero(FT),
            )
            spo = LeafProspectProProperties{FT}(
                N = N_spo, Ccab = FT(45) * f_spo, Ccar = FT(10) * f_spo, Canth = FT(1.5) * f_spo,
                Cbrown = FT(0.05) * f_spo, Cw = FT(0.014) * f_spo, Cm = FT(0.009) * f_spo,
                Cprot = zero(FT), Ccbc = zero(FT),
            )
            dors_sym = DorsiventralLeafProperties(
                pal, spo;
                surface_adaxial = ConeSurface{FT}(FT(40), one(FT)),
                surface_abaxial = ConeSurface{FT}(FT(40), one(FT)),
            )
            b_dors = dorsiventral_budget(dors_sym, opti)
            b_ref = prospect_budget(leaf_whole, opti)
            @test maximum(abs.(b_dors.adaxial.reflectance .- b_ref.reflectance)) <= tol
            @test maximum(abs.(b_dors.adaxial.transmission .- b_ref.transmission)) <= tol
            @test maximum(abs.(b_dors.adaxial.absorptance .- b_ref.absorptance)) <= tol

            # 3b. Single-slab N = 1 reduction to prospect_budget(N = 1)
            leaf_n1 = LeafProspectProProperties{FT}(
                N = one(FT), Ccab = FT(35), Ccar = FT(8), Canth = FT(1),
                Cbrown = zero(FT), Cw = FT(0.010), Cm = FT(0.007), Cprot = zero(FT), Ccbc = zero(FT),
            )
            pal_n1 = LeafProspectProProperties{FT}(
                N = one(FT), Ccab = FT(25), Ccar = FT(6), Canth = FT(1),
                Cbrown = zero(FT), Cw = FT(0.005), Cm = FT(0.004), Cprot = zero(FT), Ccbc = zero(FT),
            )
            spo_n1 = LeafProspectProProperties{FT}(
                N = one(FT), Ccab = FT(10), Ccar = FT(2), Canth = zero(FT),
                Cbrown = zero(FT), Cw = FT(0.005), Cm = FT(0.003), Cprot = zero(FT), Ccbc = zero(FT),
            )
            dors_n1 = DorsiventralLeafProperties(
                pal_n1, spo_n1;
                surface_adaxial = ConeSurface{FT}(FT(40), one(FT)),
                surface_abaxial = ConeSurface{FT}(FT(40), one(FT)),
                single_slab = true,
            )
            b_n1 = dorsiventral_budget(dors_n1, opti)
            b_ref_n1 = prospect_budget(leaf_n1, opti)
            @test maximum(abs.(b_n1.adaxial.reflectance .- b_ref_n1.reflectance)) <= tol
            @test maximum(abs.(b_n1.adaxial.transmission .- b_ref_n1.transmission)) <= tol
            @test maximum(abs.(b_n1.adaxial.absorptance .- b_ref_n1.absorptance)) <= tol

            # 3c. Strongly asymmetric leaf: conservation, face-flip invariance & hemispherical reciprocity!
            pal_asym = LeafProspectProProperties{FT}(
                N = FT(1.15), Ccab = FT(48), Ccar = FT(11), Canth = FT(3),
                Cbrown = zero(FT), Cw = FT(0.008), Cm = FT(0.006), Cprot = zero(FT), Ccbc = zero(FT),
            )
            spo_asym = LeafProspectProProperties{FT}(
                N = FT(2.25), Ccab = FT(12), Ccar = FT(3), Canth = zero(FT),
                Cbrown = zero(FT), Cw = FT(0.006), Cm = FT(0.004), Cprot = zero(FT), Ccbc = zero(FT),
            )
            s_ad = RoughDielectricSurface{FT}(FT(0.07), FT(1.05))
            s_ab = RoughDielectricSurface{FT}(FT(0.28), FT(0.97))

            d_fwd = DorsiventralLeafProperties(pal_asym, spo_asym; surface_adaxial = s_ad, surface_abaxial = s_ab)
            d_rev = DorsiventralLeafProperties(spo_asym, pal_asym; surface_adaxial = s_ab, surface_abaxial = s_ad)

            bf_dir = dorsiventral_budget(d_fwd, opti; μ_ad = FT(0.85), μ_ab = FT(0.85))
            br_dir = dorsiventral_budget(d_rev, opti; μ_ad = FT(0.85), μ_ab = FT(0.85))

            # Conservation on both faces
            for fb in (bf_dir.adaxial, bf_dir.abaxial)
                @test maximum(abs.(fb.reflectance .+ fb.transmission .+ fb.absorptance .- one(FT))) <= tol
                @test maximum(abs.(fb.absorptance_palisade .+ fb.absorptance_spongy .- fb.absorptance)) <= tol
            end
            # Face-flip symmetry
            @test maximum(abs.(bf_dir.adaxial.reflectance .- br_dir.abaxial.reflectance)) <= tol
            @test maximum(abs.(bf_dir.adaxial.transmission .- br_dir.abaxial.transmission)) <= tol
            @test maximum(abs.(bf_dir.adaxial.absorptance_palisade .- br_dir.abaxial.absorptance_spongy)) <= tol

            # Hemispherical (2π -> 2π) Helmholtz reciprocity: T_{ad->ab,hem} == T_{ab->ad,hem}
            bf_hem = dorsiventral_budget(d_fwd, opti; hemispherical = true)
            @test maximum(abs.(bf_hem.adaxial.transmission .- bf_hem.abaxial.transmission)) <= tol
            # Whereas R_ad != R_ab in visible due to palisade pigment concentration + spongy scattering!
            @test bf_hem.abaxial.reflectance[2] > bf_hem.adaxial.reflectance[2] + FT(0.02)
        end

        # ForwardDiff Jacobian wrt dorsiventral traits & rough surfaces
        function dors_vec(θ::AbstractVector{T}) where {T}
            p = LeafProspectProProperties{T}(N = θ[1], Ccab = θ[3], Cw = T(0.01), Cm = T(0.005))
            s = LeafProspectProProperties{T}(N = θ[2], Ccab = θ[4], Cw = T(0.01), Cm = T(0.005))
            d = DorsiventralLeafProperties(
                p, s;
                surface_adaxial = RoughDielectricSurface{T}(θ[5], θ[6]),
                surface_abaxial = RoughDielectricSurface{T}(T(0.20), one(T)),
            )
            b = dorsiventral_budget(d, opti64_full; μ_ad = T(1))
            return vcat(b.adaxial.reflectance, b.abaxial.reflectance, b.adaxial.transmission)
        end
        J_dors = ForwardDiff.jacobian(dors_vec, [1.2, 1.8, 40.0, 15.0, 0.12, 1.03])
        @test all(isfinite, J_dors)

        # 3d. AsymmetricSurfaceLeafProperties: N = 1, 1 < N < 2 ConeSurface & Rough guard, N >= 2, reciprocity & AD wrt N!
        @test AsymmetricSurfaceLeafProperties(
            LeafProspectProProperties{Float32}();
            surface_adaxial = ConeSurface(8.0f0),
            surface_abaxial = ConeSurface(40.0f0),
        ) isa AsymmetricSurfaceLeafProperties
        for N_cone in (1.0, 1.45, 2.3)
            lf_c = LeafProspectProProperties{Float64}(N = N_cone, Ccab = 40.0, Cw = 0.012, Cm = 0.007)
            asym_c = AsymmetricSurfaceLeafProperties(lf_c, ConeSurface(8.0, 1.0), ConeSurface(40.0, 1.0))
            bc = dorsiventral_budget(asym_c, opti64_full)
            @test maximum(abs.(bc.adaxial.reflectance .+ bc.adaxial.transmission .+ bc.adaxial.absorptance .- 1.0)) <= 1e-13
            @test maximum(abs.(bc.abaxial.reflectance .+ bc.abaxial.transmission .+ bc.abaxial.absorptance .- 1.0)) <= 1e-13
            bc_hem = dorsiventral_budget(asym_c, opti64_full; hemispherical = true)
            @test maximum(abs.(bc_hem.adaxial.transmission .- bc_hem.abaxial.transmission)) <= 1e-13
        end
        # Nonzero ForwardDiff derivative wrt N right at N = 1.0 and N = 2.0 for ConeSurface(..., 1)
        for N0 in (1.0, 1.5, 2.0)
            dN_cone = ForwardDiff.derivative(N0) do n_var
                lf = LeafProspectProProperties{typeof(n_var)}(N = n_var, Ccab = one(n_var) * 40, Cw = one(n_var) * 0.012, Cm = one(n_var) * 0.007)
                as = AsymmetricSurfaceLeafProperties(lf, ConeSurface(8.0, 1.0), ConeSurface(40.0, 1.0))
                sum(dorsiventral_budget(as, opti64_full).adaxial.reflectance)
            end
            @test isfinite(dN_cone) && dN_cone > 0.0
        end
        # RoughDielectricSurface AsymmetricSurfaceLeafProperties: continuous N >= 1 (including 1 < N < 2 convex bridge)
        for N_r in (1.0, 1.5, 2.4)
            asym_r = AsymmetricSurfaceLeafProperties(
                LeafProspectProProperties{Float64}(N = N_r, Ccab = 40.0, Cw = 0.012, Cm = 0.007),
                RoughDielectricSurface(0.08, 1.03), RoughDielectricSurface(0.24, 0.98),
            )
            br_hem = dorsiventral_budget(asym_r, opti64_full; hemispherical = true)
            @test maximum(abs.(br_hem.adaxial.transmission .- br_hem.abaxial.transmission)) <= 1e-13
            @test maximum(abs.(br_hem.adaxial.reflectance .+ br_hem.adaxial.transmission .+ br_hem.adaxial.absorptance .- 1.0)) <= 1e-13
        end
        g_ns_mid = ForwardDiff.gradient([1.0, 1.0]) do ns
            lf = LeafProspectProProperties{eltype(ns)}(N = eltype(ns)(1.5), Ccab = eltype(ns)(40), Cw = eltype(ns)(0.012), Cm = eltype(ns)(0.007))
            as = AsymmetricSurfaceLeafProperties(lf, ConeSurface{eltype(ns)}(eltype(ns)(8), ns[1]), ConeSurface{eltype(ns)}(eltype(ns)(40), ns[2]))
            bd = dorsiventral_budget(as, opti64_full)
            sum(bd.adaxial.reflectance) + 2 * sum(bd.abaxial.reflectance)
        end
        @test all(isfinite, g_ns_mid) && g_ns_mid[1] != 0.0 && g_ns_mid[2] != 0.0
        dN_rough2 = ForwardDiff.derivative(2.0) do n_var
            lf = LeafProspectProProperties{typeof(n_var)}(N = n_var, Ccab = one(n_var) * 40, Cw = one(n_var) * 0.012, Cm = one(n_var) * 0.007)
            as = AsymmetricSurfaceLeafProperties(lf, RoughDielectricSurface(0.08, 1.03), RoughDielectricSurface(0.24, 0.98))
            sum(dorsiventral_budget(as, opti64_full).adaxial.reflectance)
        end
        @test isfinite(dN_rough2) && dN_rough2 > 0.0
        # ForwardDiff gradient wrt [n_scale_ad, n_scale_ab] right at symmetric [1.0, 1.0] on ConeSurface(N = 2.2)
        g_ns_cone = ForwardDiff.gradient([1.0, 1.0]) do ns
            lf = LeafProspectProProperties{eltype(ns)}(N = eltype(ns)(2.2), Ccab = eltype(ns)(40), Cw = eltype(ns)(0.012), Cm = eltype(ns)(0.007))
            as = AsymmetricSurfaceLeafProperties(lf, ConeSurface{eltype(ns)}(eltype(ns)(8), ns[1]), ConeSurface{eltype(ns)}(eltype(ns)(40), ns[2]))
            bd = dorsiventral_budget(as, opti64_full)
            sum(bd.adaxial.reflectance) + 2 * sum(bd.abaxial.reflectance)
        end
        @test all(isfinite, g_ns_cone) && g_ns_cone[1] != 0.0 && g_ns_cone[2] != 0.0
    end

    @testset "4. CoupledLeafBSDF: Energy Conservation, Directional & Hovenier Reciprocity, Polarization" begin
        pal = LeafProspectProProperties{Float64}(N = 1.2, Ccab = 45.0, Ccar = 9.0, Cw = 0.012, Cm = 0.005)
        spo = LeafProspectProProperties{Float64}(N = 1.9, Ccab = 15.0, Ccar = 3.0, Cw = 0.010, Cm = 0.004)
        dors = DorsiventralLeafProperties(
            pal, spo;
            surface_adaxial = RoughDielectricSurface(0.08, 1.04),
            surface_abaxial = RoughDielectricSurface(0.24, 0.98),
        )

        μ28, wμ28 = CanopyOptics.gauleg(28, 0.0, 1.0)
        ϕ36, wϕ36 = CanopyOptics.gauleg(36, 0.0, Float64(π))

        for (η_t, γ_t) in ((0.0, 4.0), (0.65, 7.5))
            bsdf_mod = CoupledLeafBSDF(dors; η_t, γ_t)
            prep = prepare_leaf_bsdf(bsdf_mod, opti64_full)

            # 4a. Directional energy conservation across θ_i: ∫ f_r μ_o dω_o + ∫ f_t μ_o dω_o + A(μ_i) == 1
            for μi in (1.0, 0.65, 0.25)
                db = leaf_directional_budget(prep, μi; face = :adaxial)
                @test maximum(abs.(db.reflectance .+ db.transmission .+ db.absorptance .- 1.0)) <= 1e-13

                int_r = zeros(length(opti64_full.nᵣ))
                int_t = zeros(length(opti64_full.nᵣ))
                for km in eachindex(μ28)
                    μo = μ28[km]
                    wm = wμ28[km]
                    for kp in eachindex(ϕ36)
                        ϕo = ϕ36[kp]
                        wp = 2.0 * wϕ36[kp] # symmetry over [0, 2π]
                        sp = leaf_bsdf_spectrum(prep, μi, 0.0, μo, ϕo; face = :adaxial, npol = 1)
                        int_r .+= (wm * wp * μo) .* sp.brdf
                        int_t .+= (wm * wp * μo) .* sp.btdf
                    end
                end
                @test maximum(abs.(int_r .- db.reflectance)) <= 2.5e-3
                @test maximum(abs.(int_t .- db.transmission)) <= 2.5e-3
                @test maximum(abs.(int_r .+ int_t .+ db.absorptance .- 1.0)) <= 3.0e-3
            end

            # 4b. Exact Directional Helmholtz Reciprocity (BRDF, Hovenier 4x4 Mueller, and asymmetric BTDF!)
            P_hov = Diagonal([1.0, 1.0, -1.0, 1.0])
            for (μ1, ϕ1, μ2, ϕ2) in ((0.82, 0.2, 0.47, 1.9), (1.0, 0.0, 0.57, π), (0.35, -0.7, 0.91, 2.4))
                # Same-face BRDF & 4x4 Mueller reciprocity: ω_1 <-> ω_2
                s12 = leaf_bsdf_spectrum(prep, μ1, ϕ1, μ2, ϕ2; face = :adaxial, npol = 4)
                s21 = leaf_bsdf_spectrum(prep, μ2, ϕ2, μ1, ϕ1; face = :adaxial, npol = 4)
                @test maximum(abs.(s12.brdf .- s21.brdf)) <= 1e-13
                @test all(0.0 .<= s12.dolp .<= 1.0 + 1e-13)
                for idx in axes(s12.brdf_mueller, 3)
                    M12 = s12.brdf_mueller[:, :, idx]
                    M21 = s21.brdf_mueller[:, :, idx]
                    @test isapprox(M21, P_hov * transpose(M12) * P_hov; rtol = 1e-12, atol = 1e-14)
                end

                # Cross-face BTDF reciprocity:
                # Forward ray has source (μ1, ϕ1) on :adaxial (k_in = (-s1 cos ϕ1, -s1 sin ϕ1, -μ1))
                # and exits :abaxial toward (μ2, ϕ2) (k_out = (s2 cos ϕ2, s2 sin ϕ2, -μ2)).
                # Time-reversed ray enters :abaxial from source ω'_i = -k_out -> (μ2, ϕ2 + π) in -ẑ frame?
                # Wait: c_trans_dot = μ1*μ2 - s1*s2*cos(ϕ2 - ϕ1) is symmetric under (μ1, ϕ1) <-> (μ2, ϕ2)!
                s_ab_rev = leaf_bsdf_spectrum(prep, μ2, ϕ2, μ1, ϕ1; face = :abaxial, npol = 1)
                @test all(s12.btdf .>= 0.0)
                @test maximum(abs.(s12.btdf .- s_ab_rev.btdf)) <= 1e-13
            end

            # 4c. Physical 55° Brewster-like polarization signature: DoLP(650 nm, strong absorption) > DoLP(850 nm, NIR multiple scattering)
            s_55 = leaf_bsdf_spectrum(prep, cosd(55.0), 0.0, cosd(55.0), π; face = :adaxial, npol = 3)
            @test all(abs.(s_55.brdf_mueller[3, 1, :]) .<= 1e-14) # U == 0 in principal plane
            @test s_55.dolp[3] > s_55.dolp[5]
            # For a mesophytic rough cuticle (σ = 0.25) on whole leaf (N_tot ≈ 3.1), unpolarized NIR diffuse scattering strongly depolarizes 850 nm vs 650 nm
            dors_rough = DorsiventralLeafProperties(
                pal, spo;
                surface_adaxial = RoughDielectricSurface(0.25, 1.04),
                surface_abaxial = RoughDielectricSurface(0.25, 0.98),
            )
            prep_rough = prepare_leaf_bsdf(CoupledLeafBSDF(dors_rough; η_t, γ_t), opti64_full)
            s_55_rough = leaf_bsdf_spectrum(prep_rough, cosd(55.0), 0.0, cosd(55.0), π; face = :adaxial, npol = 3)
            @test s_55_rough.dolp[3] > s_55_rough.dolp[5] + 0.25

            # 4d. ForwardDiff gradient through CoupledLeafBSDF (including principal plane & nadir)
            bsdf_loss(x) = begin
                lf = LeafProspectProProperties{eltype(x)}(N = x[1], Ccab = x[2], Cw = eltype(x)(0.012), Cm = eltype(x)(0.006))
                sf = RoughDielectricSurface{eltype(x)}(x[3], x[4])
                md = CoupledLeafBSDF(lf, sf; η_t = x[5], γ_t = x[6])
                sp_nad = leaf_bsdf_spectrum(md, opti64_full, one(eltype(x)), zero(eltype(x)), eltype(x)(0.6), eltype(x)(π); npol = 3)
                sp_obl = leaf_bsdf_spectrum(md, opti64_full, eltype(x)(cosd(55)), zero(eltype(x)), eltype(x)(cosd(55)), eltype(x)(π); npol = 3)
                return sum(sp_nad.brdf) + sum(sp_nad.btdf) + sum(sp_obl.dolp)
            end
            g_bsdf = ForwardDiff.gradient(bsdf_loss, [1.45, 40.0, 0.14, 1.02, 0.35, 5.0])
            @test all(isfinite, g_bsdf)
        end

        # Float32 zero-promotion check on CoupledLeafBSDF + AsymmetricSurfaceLeafProperties BSDF + 5-arg mismatch guards
        bsdf32 = CoupledLeafBSDF(LeafProspectProProperties{Float32}(), RoughDielectricSurface{Float32}(0.15f0, 1.0f0); η_t = 0.3f0, γ_t = 4.0f0)
        sp32 = leaf_bsdf_spectrum(bsdf32, opti32, 0.8f0, 0.0f0, 0.7f0, Float32(π); npol = 3)
        @test eltype(sp32.brdf) === Float32
        @test eltype(sp32.btdf) === Float32
        @test eltype(sp32.brdf_mueller) === Float32
        @test eltype(sp32.dolp) === Float32

        asym_bsdf_leaf = AsymmetricSurfaceLeafProperties(
            LeafProspectProProperties{Float64}(N = 2.2, Ccab = 40.0, Cw = 0.012, Cm = 0.007),
            RoughDielectricSurface(0.08, 1.04), RoughDielectricSurface(0.26, 0.98),
        )
        bsdf_asym = CoupledLeafBSDF(asym_bsdf_leaf; η_t = 0.4, γ_t = 5.0)
        prep_asym = prepare_leaf_bsdf(bsdf_asym, opti64_full)
        s_asym_ad = leaf_bsdf_spectrum(prep_asym, 0.8, 0.1, 0.5, 2.0; face = :adaxial, npol = 3)
        s_asym_ab = leaf_bsdf_spectrum(prep_asym, 0.5, 2.0, 0.8, 0.1; face = :abaxial, npol = 3)
        @test maximum(abs.(s_asym_ad.btdf .- s_asym_ab.btdf)) <= 1e-13
        db_asym = leaf_directional_budget(prep_asym, 0.75; face = :abaxial)
        @test maximum(abs.(db_asym.reflectance .+ db_asym.transmission .+ db_asym.absorptance .- 1.0)) <= 1e-13

        # 5-arg CoupledLeafBSDF mismatch guards + surface_brdf_mueller guards
        @test_throws ArgumentError CoupledLeafBSDF(
            LeafProspectProProperties{Float64}(),
            RoughDielectricSurface(0.10, 1.0), RoughDielectricSurface(0.25, 1.0), 0.0, 4.0,
        )
        @test_throws ArgumentError CoupledLeafBSDF(
            dors, RoughDielectricSurface(0.15, 1.0), RoughDielectricSurface(0.24, 0.98), 0.0, 4.0,
        )
        @test_throws ArgumentError surface_brdf_mueller(0.8, 0.0, 0.7, π, 0.99, 0.2, 3)
        @test_throws ArgumentError surface_brdf_mueller(0.8, 0.0, 0.7, π, 1.45, 0.0, 3)

        # ForwardDiff gradient wrt [σ_ad, σ_ab] at symmetric [0.15, 0.15] preserves both adaxial & abaxial tangents!
        g_sym_sig = ForwardDiff.gradient([0.15, 0.15]) do sigs
            T = eltype(sigs)
            p = LeafProspectProProperties{T}(N = T(1.3), Ccab = T(35), Cw = T(0.012), Cm = T(0.006))
            dv = DorsiventralLeafProperties(p, p; surface_adaxial = RoughDielectricSurface{T}(sigs[1], one(T)), surface_abaxial = RoughDielectricSurface{T}(sigs[2], one(T)))
            sp = leaf_bsdf_spectrum(CoupledLeafBSDF(dv; η_t = T(0.35), γ_t = T(4.0)), opti64_full, T(0.8), zero(T), T(0.6), T(π); face = :adaxial, npol = 1)
            sum(sp.btdf)
        end
        @test all(isfinite, g_sym_sig) && g_sym_sig[1] != 0.0 && g_sym_sig[2] != 0.0
    end

    @testset "5. Dorsiventral & Coupled-BSDF Canopy Z Matrices" begin
        μ, w = CanopyOptics.gauleg(6, 0.0, 1.0)
        μ = collect(μ)
        quad = CanopyOptics.CanopyQuadrature(n_leaf = 36, n_azimuth = 32)

        # 5a. Symmetric reduction to BiLambertianCanopyScattering & asymmetric guard
        sym_dors = DorsiventralBiLambertianCanopyScattering(0.42, 0.08, 0.42, 0.08)
        ref_bi = BiLambertianCanopyScattering(0.42, 0.08)
        for LD in (spherical_leaves(), planophile_leaves2(), erectophile_leaves())
            Zpp_ref, Zmp_ref = compute_Z_matrices_aniso_analytic(ref_bi, μ, LD, 4; quadrature = quad)
            res_norm = compute_Z_matrices_directional(sym_dors, μ, LD, 4; quadrature = quad, normalize_by_ω = true)
            @test isapprox(res_norm.Zpp, Zpp_ref; rtol = 1e-13)
            @test isapprox(res_norm.Zmp, Zmp_ref; rtol = 1e-13)
            @test isapprox(res_norm.Zpm, Zmp_ref; rtol = 1e-13)
            @test isapprox(res_norm.Zmm, Zpp_ref; rtol = 1e-13)
        end

        asym_dors = DorsiventralBiLambertianCanopyScattering(0.06, 0.04, 0.18, 0.04)
        @test_throws ArgumentError compute_Z_matrices(asym_dors, μ, spherical_leaves(), 0)
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices_aniso(asym_dors, μ, spherical_leaves(), 0)
        @test_throws ArgumentError compute_Z_matrices_aniso_analytic(asym_dors, μ, spherical_leaves(), 2)
        @test_throws ArgumentError compute_Z_matrices_directional(asym_dors, [-0.1, 0.5], spherical_leaves(), 1)
        Z_aniso_sym = CanopyOptics.compute_Z_matrices_aniso(sym_dors, μ, spherical_leaves(), 1; quadrature = quad)
        Z_aniso_ref = CanopyOptics.compute_Z_matrices_aniso(ref_bi, μ, spherical_leaves(), 1; quadrature = quad)
        @test isapprox(Z_aniso_sym[1], Z_aniso_ref[1]; rtol = 1e-13)
        @test isapprox(Z_aniso_sym[2], Z_aniso_ref[2]; rtol = 1e-13)

        # 5b. Analytic closed-form vs direct 2D-quadrature compute_Γ_directional across m = 0:3
        LD_test = planophile_leaves2()
        res_dir = compute_Z_matrices_directional(asym_dors, μ, LD_test, 3; quadrature = quad, normalize_by_ω = false)
        G_vec = vec(CanopyOptics.G(μ, LD_test; quadrature = quad))
        quad_ref2d = CanopyOptics.CanopyQuadrature(n_leaf = 36, n_azimuth = 96)
        ϕ_az, w_az = CanopyOptics.gauleg(64, 0.0, 2π)

        for j in (2, 5), i in (1, 4), m in 0:3
            ff = m == 0 ? 2.0 : 4.0
            # Note: (1 / 2π) ∫_0^{2π} (2 Γ(Δϕ)) cos(m Δϕ) dΔϕ = (1 / π) ∫_0^{2π} Γ(Δϕ) cos(m Δϕ) dΔϕ,
            # multiplied by fourier_factor / G_j (no extra 1/2!). Because compute_Γ_directional integrates
            # piecewise-smooth max(0, d_in) with C^0 kinks numerically over ϕ_L, error is O(1/n_azimuth^2) ~ 1e-4.
            acc_pp = 0.0; acc_mp = 0.0; acc_pm = 0.0; acc_mm = 0.0
            for ia in eachindex(ϕ_az)
                dϕ = ϕ_az[ia]
                wt = w_az[ia] * cos(m * dϕ) / π
                acc_pp += wt * compute_Γ_directional(asym_dors, dirVector_μ(μ[j], 0.0), dirVector_μ(μ[i], dϕ), LD_test; quadrature = quad_ref2d)
                acc_mp += wt * compute_Γ_directional(asym_dors, dirVector_μ(μ[j], 0.0), dirVector_μ(-μ[i], dϕ), LD_test; quadrature = quad_ref2d)
                acc_pm += wt * compute_Γ_directional(asym_dors, dirVector_μ(-μ[j], 0.0), dirVector_μ(μ[i], dϕ), LD_test; quadrature = quad_ref2d)
                acc_mm += wt * compute_Γ_directional(asym_dors, dirVector_μ(-μ[j], 0.0), dirVector_μ(-μ[i], dϕ), LD_test; quadrature = quad_ref2d)
            end
            scale = ff / G_vec[j]
            @test isapprox(res_dir.Zpp[i, j, m + 1], acc_pp * scale; rtol = 5e-4, atol = 2e-5)
            @test isapprox(res_dir.Zmp[i, j, m + 1], acc_mp * scale; rtol = 5e-4, atol = 2e-5)
            @test isapprox(res_dir.Zpm[i, j, m + 1], acc_pm * scale; rtol = 5e-4, atol = 2e-5)
            @test isapprox(res_dir.Zmm[i, j, m + 1], acc_mm * scale; rtol = 5e-4, atol = 2e-5)
        end

        # 5c. Exact 4-block Helmholtz reciprocity when T_ad == T_ab (even with R_ad != R_ab!)
        for m in 0:3
            Wpp = res_dir.Zpp[:, :, m + 1] .* reshape(G_vec, 1, :)
            Wmm = res_dir.Zmm[:, :, m + 1] .* reshape(G_vec, 1, :)
            Wmp = res_dir.Zmp[:, :, m + 1] .* reshape(G_vec, 1, :)
            Wpm = res_dir.Zpm[:, :, m + 1] .* reshape(G_vec, 1, :)
            @test norm(Wmp .- transpose(Wmp)) <= 1e-12
            @test norm(Wpm .- transpose(Wpm)) <= 1e-12
            @test norm(Wpp .- transpose(Wmm)) <= 1e-12
        end

        # Column flux conservation at m = 0: 0.5 * sum(w .* (Zpp + Zmp)) ≈ ω_down
        col_dn = 0.5 .* vec(sum(w .* (res_dir.Zpp[:, :, 1] .+ res_dir.Zmp[:, :, 1]); dims = 1))
        col_up = 0.5 .* vec(sum(w .* (res_dir.Zpm[:, :, 1] .+ res_dir.Zmm[:, :, 1]); dims = 1))
        @test maximum(abs.(col_dn .- res_dir.ω_down)) <= 2e-3
        @test maximum(abs.(col_up .- res_dir.ω_up)) <= 2e-3
        # For planophile leaves, downward reflection uses mostly R_ad (0.06) while upward uses mostly R_ab (0.18),
        # with ω_up >= ω_down everywhere (approaching equality at horizontal grazing μ -> 0 and > +0.06 near nadir):
        @test all(res_dir.ω_up .>= res_dir.ω_down)
        @test res_dir.ω_up[end] > res_dir.ω_down[end] + 0.06

        # 5d. CoupledLeafBSDF canopy projection energy conservation & boundedness
        μ12, w12 = CanopyOptics.gauleg(12, 0.0, 1.0)
        μ12 = collect(μ12)
        bsdf_can = CoupledLeafBSDF(
            LeafProspectProProperties{Float64}(N = 1.5, Ccab = 40.0, Cw = 0.012, Cm = 0.006),
            RoughDielectricSurface(0.18, 1.0);
            η_t = 0.25, γ_t = 4.0,
        )
        prep_can = prepare_leaf_bsdf(bsdf_can, opti64_full)
        z_bsdf = compute_Z_matrices_coupled_bsdf(
            prep_can, 5, μ12, spherical_leaves(), 1;
            quadrature = CanopyQuadrature(n_leaf = 16, n_azimuth = 24),
        )
        flux_dn_bsdf = 0.5 .* vec(sum(w12 .* (z_bsdf.Zpp[:, :, 1] .+ z_bsdf.Zmp[:, :, 1]); dims = 1))
        @test all(0.0 .< z_bsdf.ω_down .<= 1.0)
        @test maximum(abs.(flux_dn_bsdf .- z_bsdf.ω_down)) <= 1.5e-2
        @test_throws ArgumentError compute_Z_matrices_coupled_bsdf(prep_can, 5, [0.0, 0.5], spherical_leaves(), 1)

        # 4-block canopy reciprocity on coupled BSDF (even with strongly asymmetric dorsiventral leaf!)
        dors_can = DorsiventralLeafProperties(
            LeafProspectProProperties{Float64}(N = 1.2, Ccab = 45.0, Cw = 0.010, Cm = 0.005),
            LeafProspectProProperties{Float64}(N = 1.9, Ccab = 15.0, Cw = 0.010, Cm = 0.005);
            surface_adaxial = RoughDielectricSurface(0.10, 1.04),
            surface_abaxial = RoughDielectricSurface(0.26, 0.98),
        )
        prep_dors_can = prepare_leaf_bsdf(CoupledLeafBSDF(dors_can; η_t = 0.45, γ_t = 5.0), opti64_full)
        z_dors_bsdf = compute_Z_matrices_coupled_bsdf(
            prep_dors_can, 2, μ12, planophile_leaves2(), 2;
            quadrature = CanopyQuadrature(n_leaf = 16, n_azimuth = 24),
        )
        G12 = vec(CanopyOptics.G(μ12, planophile_leaves2(); quadrature = CanopyQuadrature(n_leaf = 16, n_azimuth = 24)))
        for m in 0:2
            Wpp = z_dors_bsdf.Zpp[:, :, m + 1] .* reshape(G12, 1, :)
            Wmm = z_dors_bsdf.Zmm[:, :, m + 1] .* reshape(G12, 1, :)
            Wmp = z_dors_bsdf.Zmp[:, :, m + 1] .* reshape(G12, 1, :)
            Wpm = z_dors_bsdf.Zpm[:, :, m + 1] .* reshape(G12, 1, :)
            @test isapprox(Wmp, transpose(Wmp); rtol = 1e-12, atol = 1e-14)
            @test isapprox(Wpm, transpose(Wpm); rtol = 1e-12, atol = 1e-14)
            @test isapprox(Wpp, transpose(Wmm); rtol = 1e-12, atol = 1e-14)
        end

        # Conservative k = 0 leaf guarantees exact convex-combination stream albedo ω_down == ω_up == 1
        leaf_zero_can = LeafProspectProProperties{Float64}(
            N = 1.6, Ccab = 0.0, Ccar = 0.0, Canth = 0.0, Cbrown = 0.0,
            Cw = 0.0, Cm = 0.0, Cprot = 0.0, Ccbc = 0.0,
        )
        prep_zero_can = prepare_leaf_bsdf(
            CoupledLeafBSDF(leaf_zero_can, RoughDielectricSurface(0.18, 1.0); η_t = 0.3, γ_t = 4.0),
            opti64_full,
        )
        z_zero_can = compute_Z_matrices_coupled_bsdf(
            prep_zero_can, 1, μ12, planophile_leaves2(), 0;
            quadrature = CanopyQuadrature(n_leaf = 16, n_azimuth = 24),
        )
        @test all(0.0 .<= z_zero_can.ω_down .<= 1.0)
        @test all(0.0 .<= z_zero_can.ω_up .<= 1.0)
        @test all(isapprox.(z_zero_can.ω_down, 1.0; atol = 1e-13))
        @test all(isapprox.(z_zero_can.ω_up, 1.0; atol = 1e-13))
    end

    @testset "6. Dry-Matter Keyword Constructors (_UnsetLeafKw) & Short-Grid Optical Loader" begin
        # Default PRO mode: Cm = 0, Cprot = 0.001, Ccbc = 0.009
        lf_def = LeafProspectProProperties()
        @test lf_def isa LeafProspectProProperties{Float64}
        @test lf_def.Cm == 0.0 && lf_def.Cprot == 0.001 && lf_def.Ccbc == 0.009

        # Explicit PROSPECT-D Cm mode zeroes Cprot & Ccbc (preventing +0.010 g/cm² double counting!)
        lf_cm = LeafProspectProProperties(Cm = 0.006)
        @test lf_cm isa LeafProspectProProperties{Float64}
        @test lf_cm.Cm == 0.006 && lf_cm.Cprot == 0.0 && lf_cm.Ccbc == 0.0

        # Float32 preservation on untyped keyword constructors
        lf_cm32 = LeafProspectProProperties(Cm = 0.006f0)
        @test lf_cm32 isa LeafProspectProProperties{Float32}
        @test lf_cm32.Cprot === 0.0f0 && lf_cm32.Ccbc === 0.0f0

        mp_cm32 = LeafMultiPigmentProperties(Cm = 0.005f0)
        @test mp_cm32 isa LeafMultiPigmentProperties{Float32}
        @test mp_cm32.Cm === 0.005f0 && mp_cm32.Cprot === 0.0f0 && mp_cm32.Ccbc === 0.0f0

        sip_pro32 = LeafSIPProperties(Cprot = 0.001f0, Ccbc = 0.008f0)
        @test sip_pro32 isa LeafSIPProperties{Float32}
        @test sip_pro32.Cm === 0.0f0 && sip_pro32.Cprot === 0.001f0 && sip_pro32.Ccbc === 0.008f0

        sip_def = LeafSIPProperties()
        @test sip_def isa LeafSIPProperties{Float64}
        @test sip_def.Cm == 0.008 && sip_def.Cprot == 0.0 && sip_def.Ccbc == 0.0

        # Mutual exclusivity guard when both Cm > 0 and Cprot/Ccbc > 0 are supplied
        @test_throws ArgumentError LeafProspectProProperties(Cm = 0.005, Cprot = 0.001)
        @test_throws ArgumentError LeafMultiPigmentProperties(Cm = 0.005, Ccbc = 0.004)
        @test_throws ArgumentError LeafSIPProperties(Cm = 0.005, Cprot = 0.001)

        # ForwardDiff derivative wrt Cm right at Cm = 0.0 has zero Dual-vs-primal jump!
        d_cm0 = ForwardDiff.derivative(0.0) do cm
            lf = LeafProspectProProperties(Cm = cm)
            sum(prospect_budget(lf, opti64_full).absorptance)
        end
        d_cm_pos = ForwardDiff.derivative(1e-7) do cm
            lf = LeafProspectProProperties(Cm = cm)
            sum(prospect_budget(lf, opti64_full).absorptance)
        end
        @test isfinite(d_cm0) && d_cm0 > 0.0
        @test isapprox(d_cm0, d_cm_pos; rtol = 1e-4)

        # Short spectral grid (< 10 wavelength bins) no longer throws BoundsError in createLeafOpticalStruct
        opti_short = createLeafOpticalStruct((450.0:50.0:700.0) .* _u_nm)
        @test length(opti_short.nᵣ) == 5
    end

    @testset "7. Continuous N_tot >= 1 Dorsiventral Bridge & dorsiventral_from_traits" begin
        whole_pro = LeafProspectProProperties(N = 1.45, Ccab = 42.0, Ccar = 9.0, Cw = 0.012, Cm = 0.006)
        for nt in (1.0, 1.35, 1.999999, 2.0, 2.000001, 2.8)
            dv = dorsiventral_from_traits(
                whole_pro;
                N_tot = nt,
                mesophyll_fraction_palisade = 0.45,
                pigment_fraction_palisade = 0.75,
                surface_adaxial = RoughDielectricSurface(0.08, 1.03),
                surface_abaxial = RoughDielectricSurface(0.24, 0.98),
            )
            bd = dorsiventral_budget(dv, opti64_full; hemispherical = true)
            @test maximum(abs.(bd.adaxial.reflectance .+ bd.adaxial.transmission .+ bd.adaxial.absorptance .- 1.0)) <= 1e-13
            @test maximum(abs.(bd.adaxial.absorptance_palisade .+ bd.adaxial.absorptance_spongy .- bd.adaxial.absorptance)) <= 1e-13
            @test maximum(abs.(bd.adaxial.transmission .- bd.abaxial.transmission)) <= 1e-13
        end

        # Verify C^0 continuity across N_tot = 2.0
        bd_left = dorsiventral_budget(dorsiventral_from_traits(whole_pro; N_tot = 2.0 - 1e-7), opti64_full)
        bd_at2 = dorsiventral_budget(dorsiventral_from_traits(whole_pro; N_tot = 2.0), opti64_full)
        bd_right = dorsiventral_budget(dorsiventral_from_traits(whole_pro; N_tot = 2.0 + 1e-7), opti64_full)
        @test maximum(abs.(bd_left.adaxial.reflectance .- bd_at2.adaxial.reflectance)) <= 1e-6
        @test maximum(abs.(bd_right.adaxial.reflectance .- bd_at2.adaxial.reflectance)) <= 1e-6

        # ForwardDiff gradient wrt [N_tot, fm, fc] on sub-2 leaf (N_tot = 1.45)
        g_dv_sub2 = ForwardDiff.gradient([1.45, 0.45, 0.75]) do θ
            dv = dorsiventral_from_traits(
                whole_pro;
                N_tot = θ[1],
                mesophyll_fraction_palisade = θ[2],
                pigment_fraction_palisade = θ[3],
            )
            bd = dorsiventral_budget(dv, opti64_full)
            sum(bd.adaxial.reflectance) + sum(bd.adaxial.absorptance_palisade)
        end
        @test all(isfinite, g_dv_sub2) && all(!iszero, g_dv_sub2)

        # Multi-pigment whole leaf support + non-convex invalid (N_pal < 1 < N_spo) rejection
        whole_mp = LeafMultiPigmentProperties(N = 1.35f0, Cchla = 30.0f0, Cchlb = 10.0f0, Cm = 0.006f0)
        dv_mp32 = dorsiventral_from_traits(whole_mp)
        @test dv_mp32.palisade isa LeafMultiPigmentProperties{Float32}
        @test_throws ArgumentError DorsiventralLeafProperties(
            LeafProspectProProperties(N = 0.6, Cm = 0.003),
            LeafProspectProProperties(N = 1.2, Cm = 0.003),
        )
    end

    @testset "8. createMultiPigmentOpticalStruct Self-Contained Optical Bases" begin
        λ_bnds = (399.5:10.0:2499.5) .* _u_nm
        opt_pro = createLeafOpticalStruct(λ_bnds)
        for basis in (:calibrated_in_vivo, :prospect_mp_2017, :prospect_mp_plus_2022, :prospect_pro_pooled)
            mp64 = createMultiPigmentOpticalStruct(λ_bnds; basis = basis)
            mp32 = createMultiPigmentOpticalStruct(λ_bnds; basis = basis, FT = Float32)
            @test mp64 isa MultiPigmentOpticalProperties{Float64}
            @test mp32 isa MultiPigmentOpticalProperties{Float32}
            @test all(mp64.Kchla .>= 0.0) && all(mp64.Kchlb .>= 0.0)
            @test all(mp64.Kcar .>= 0.0) && all(mp64.Kant .>= 0.0)

            # Consistent [cm² µg⁻¹] anthocyanin units across all 4 bases (~0.07 at 535 nm)
            wl_nm = ustrip.(mp64.λ)
            idx_535 = argmin(abs.(wl_nm .- 535.0))
            @test 0.05 <= mp64.Kant[idx_535] <= 0.11

            # Clean splice to OPTI_2021 above 800 nm (zero Lorentzian wing leakage in NIR/SWIR)
            nir_mask = wl_nm .> 805.0
            @test isapprox(mp64.Kchla[nir_mask], opt_pro.Kcab[nir_mask]; rtol = 1e-12, atol = 1e-15)
            @test isapprox(mp64.Kchlb[nir_mask], opt_pro.Kcab[nir_mask]; rtol = 1e-12, atol = 1e-15)
            @test isapprox(mp64.Kcar[nir_mask], opt_pro.Kcar[nir_mask]; rtol = 1e-12, atol = 1e-15)
            @test isapprox(mp64.Kant[nir_mask], opt_pro.Kant[nir_mask]; rtol = 1e-12, atol = 1e-15)

            if basis !== :prospect_pro_pooled
                @test maximum(abs.(mp64.Kchla .- mp64.Kchlb)) > 0.01
            else
                @test mp64.Kchla == mp64.Kchlb
            end

            lf_mp = LeafMultiPigmentProperties(N = 1.6, Cchla = 32.0, Cchlb = 10.0, Ccar = 8.0, Canth = 2.0, Cw = 0.012, Cm = 0.006)
            bd_mp = prospect_budget(lf_mp, mp64)
            @test maximum(abs.(bd_mp.reflectance .+ bd_mp.transmission .+ bd_mp.absorptance .- 1.0)) <= 1e-13
        end
        @test_throws ArgumentError createMultiPigmentOpticalStruct(λ_bnds; basis = :unknown_basis)
    end

    @testset "9. Epicuticular Wax Dispersion & Surface-Coupled BSDF / Canopy Consumer Sites" begin
        λ_bnds = (399.5:100.0:2499.5) .* _u_nm
        opt_w = createLeafOpticalStruct(λ_bnds)
        n_wax = load_wax_refractive_index(λ_bnds)
        n_wax32 = load_wax_refractive_index(λ_bnds; FT = Float32)
        @test eltype(n_wax) === Float64 && eltype(n_wax32) === Float32
        @test length(n_wax) == length(opt_w.nᵣ)
        # Stuckens et al. (2009) epicuticular wax stays in [1.40, 1.50] across 400–2500 nm and well above SWIR n_r
        @test 1.40 <= n_wax[end] < n_wax[1] < 1.50
        @test n_wax[end] < n_wax[5] < 1.50
        @test n_wax[end] > opt_w.nᵣ[end] + 0.10

        # Conservation & SWIR reflectance increase across prospect_budget and leaf_sip_budget
        lf_w = LeafProspectProProperties(N = 1.6, Ccab = 40.0, Cw = 0.015, Cm = 0.006)
        bd_nowax = prospect_budget(lf_w, opt_w)
        bd_wax = prospect_budget(lf_w, opt_w; n_surface = n_wax)
        @test maximum(abs.(bd_wax.reflectance .+ bd_wax.transmission .+ bd_wax.absorptance .- 1.0)) <= 1e-13
        @test bd_wax.reflectance[end] > bd_nowax.reflectance[end]

        sip_w = LeafSIPProperties(Ccab = 40.0, Cw = 0.015, Cm = 0.006)
        bd_sip_wax = leaf_sip_budget(sip_w, opt_w; n_surface = n_wax)
        @test maximum(abs.(bd_sip_wax.reflectance .+ bd_sip_wax.transmission .+ bd_sip_wax.absorptance .- 1.0)) <= 1e-13

        # Asymmetric single-face wax (adaxial wax, abaxial mesophyll n_r) preserves exact hemispherical T_ad == T_ab!
        for N_val in (1.0, 1.45, 2.3)
            asym_wax_leaf = AsymmetricSurfaceLeafProperties(
                LeafProspectProProperties(N = N_val, Ccab = 40.0, Cw = 0.014, Cm = 0.006),
                RoughDielectricSurface(0.10, 1.0),
                RoughDielectricSurface(0.24, 1.0),
            )
            bd_aw = dorsiventral_budget(
                asym_wax_leaf,
                opt_w;
                hemispherical = true,
                n_surface_adaxial = n_wax,
                n_surface_abaxial = opt_w.nᵣ,
            )
            @test maximum(abs.(bd_aw.adaxial.transmission .- bd_aw.abaxial.transmission)) <= 1e-13
            @test maximum(abs.(bd_aw.adaxial.reflectance .+ bd_aw.adaxial.transmission .+ bd_aw.adaxial.absorptance .- 1.0)) <= 1e-13
        end

        # CoupledLeafBSDF consumer sites with asymmetric adaxial-only wax!
        # Single-zone CoupledLeafBSDF rejects distinct n_surface_adaxial vs n_surface_abaxial
        bsdf_single = CoupledLeafBSDF(lf_w, RoughDielectricSurface(0.15, 1.0); η_t = 0.35, γ_t = 4.5)
        @test_throws ArgumentError prepare_leaf_bsdf(
            bsdf_single, opt_w; n_surface_adaxial = n_wax, n_surface_abaxial = opt_w.nᵣ,
        )
        prep_single_wax = prepare_leaf_bsdf(bsdf_single, opt_w; n_surface = n_wax)
        prep_single_nowax = prepare_leaf_bsdf(bsdf_single, opt_w)
        sp_swir_wax = leaf_bsdf_spectrum(prep_single_wax, cosd(45.0), 0.0, cosd(45.0), π)
        sp_swir_nowax = leaf_bsdf_spectrum(prep_single_nowax, cosd(45.0), 0.0, cosd(45.0), π)
        @test sp_swir_wax.brdf[end] > sp_swir_nowax.brdf[end]

        # Dorsiventral CoupledLeafBSDF with adaxial wax (n_wax) and abaxial mesophyll (n_r):
        # verify directional conservation, cross-face BTDF reciprocity, 4x4 Hovenier reciprocity, and canopy Z reciprocity!
        dv_wax = dorsiventral_from_traits(
            lf_w;
            N_tot = 1.65,
            surface_adaxial = RoughDielectricSurface(0.10, 1.0),
            surface_abaxial = RoughDielectricSurface(0.24, 1.0),
        )
        prep_dv_wax = prepare_leaf_bsdf(
            CoupledLeafBSDF(dv_wax; η_t = 0.45, γ_t = 5.0),
            opt_w;
            n_surface_adaxial = n_wax,
            n_surface_abaxial = opt_w.nᵣ,
        )
        db_dv_wax = leaf_directional_budget(prep_dv_wax, 0.72; face = :adaxial)
        @test maximum(abs.(db_dv_wax.reflectance .+ db_dv_wax.transmission .+ db_dv_wax.absorptance .- 1.0)) <= 1e-13

        μ24, wμ24 = CanopyOptics.gauleg(24, 0.0, 1.0)
        ϕ32, wϕ32 = CanopyOptics.gauleg(32, 0.0, Float64(π))
        int_r_w = zeros(length(n_wax))
        int_t_w = zeros(length(n_wax))
        for km in eachindex(μ24), kp in eachindex(ϕ32)
            sp = leaf_bsdf_spectrum(prep_dv_wax, 0.72, 0.0, μ24[km], ϕ32[kp]; face = :adaxial, npol = 1)
            wt = wμ24[km] * 2.0 * wϕ32[kp] * μ24[km]
            int_r_w .+= wt .* sp.brdf
            int_t_w .+= wt .* sp.btdf
        end
        @test maximum(abs.(int_r_w .- db_dv_wax.reflectance)) <= 2.5e-3
        @test maximum(abs.(int_t_w .- db_dv_wax.transmission)) <= 2.5e-3

        s_w_ad = leaf_bsdf_spectrum(prep_dv_wax, 0.78, 0.3, 0.51, 2.1; face = :adaxial, npol = 4)
        s_w_ab = leaf_bsdf_spectrum(prep_dv_wax, 0.51, 2.1, 0.78, 0.3; face = :abaxial, npol = 4)
        @test maximum(abs.(s_w_ad.btdf .- s_w_ab.btdf)) <= 1e-13

        μ8, _ = CanopyOptics.gauleg(8, 0.0, 1.0)
        μ8 = collect(μ8)
        quad_w = CanopyQuadrature(n_leaf = 16, n_azimuth = 24)
        z_dv_wax = compute_Z_matrices_coupled_bsdf(prep_dv_wax, length(n_wax), μ8, planophile_leaves2(), 1; quadrature = quad_w)
        G8 = vec(CanopyOptics.G(μ8, planophile_leaves2(); quadrature = quad_w))
        for m in 0:1
            Wpp = z_dv_wax.Zpp[:, :, m + 1] .* reshape(G8, 1, :)
            Wmm = z_dv_wax.Zmm[:, :, m + 1] .* reshape(G8, 1, :)
            Wmp = z_dv_wax.Zmp[:, :, m + 1] .* reshape(G8, 1, :)
            Wpm = z_dv_wax.Zpm[:, :, m + 1] .* reshape(G8, 1, :)
            @test isapprox(Wmp, transpose(Wmp); rtol = 1e-12, atol = 1e-14)
            @test isapprox(Wpm, transpose(Wpm); rtol = 1e-12, atol = 1e-14)
            @test isapprox(Wpp, transpose(Wmm); rtol = 1e-12, atol = 1e-14)
        end
    end

    @testset "10. Dorsiventral 2-Zone Internal Fluorescence Escape (CPU, BigFloat & CUDA)" begin
        @test isbitstype(DorsiventralPlateEscape{Float64})
        @test isbitstype(DorsiventralPlateEscape{Float32})
        @test isbitstype(DorsiventralBirthProfile{Float64})
        @test isbitstype(DorsiventralBirthProfile{Float32})

        # Constructor guards
        @test_throws ArgumentError DorsiventralPlateEscape(-0.1, 0.2, 0.5, 0.5)
        @test_throws ArgumentError DorsiventralPlateEscape(0.1, 0.2, 0.0, 0.5)
        @test_throws ArgumentError DorsiventralPlateEscape(0.1, 0.2, 0.5, 0.5, 0.0)
        @test_throws ArgumentError DorsiventralBirthProfile(0.0, 0.0)
        @test_throws ArgumentError DorsiventralBirthProfile(1.0, 0.0, 0.5, 0.0, 0.0, 0.5, 1.0, 0.0)

        # Zero host allocations on scalar evaluation
        esc_alloc = DorsiventralPlateEscape(0.6, 0.4, 0.45, 0.42, 0.85)
        bp_alloc = DorsiventralBirthProfile(0.75, 0.25, 1.8, 1.0, 0.05, 0.6, 1.0, 0.12)
        leaf_emission_fates(esc_alloc, bp_alloc)
        @test (@allocated leaf_emission_fates(esc_alloc, bp_alloc)) == 0
        @test (@allocated leaf_emission_fates(esc_alloc, UniformBirthProfile())) == 0

        # Exact reduction to RegularPlateEscape(b, t) when splitting a homogeneous plate
        # at arbitrary fraction η ∈ (0, 1) with interface_transmission = 1!
        for η in (0.2, 0.45, 0.8), b in (0.0, 1e-6, 0.7, 5.5), t in (0.05, 0.48, 1.0), (a, f, g) in ((0.0, 1.0, 0.0), (0.35, 1.0, 0.2), (4.5, 1.0, 0.05))
            reg = RegularPlateEscape(b, t)
            dv_esc = DorsiventralPlateEscape(η * b, (1.0 - η) * b, t, t, 1.0)
            exp_bp = ExponentialBirthProfile(a, f, g)
            ref_fates = leaf_emission_fates(reg, exp_bp)

            L1 = CanopyOptics.emission_depth_integral(η * a)
            L2 = CanopyOptics.emission_depth_integral((1.0 - η) * a)
            f1_raw = f
            g1_raw = g * exp(-(1.0 - η) * a)
            f2_raw = f * exp(-η * a)
            g2_raw = g
            w1_raw = η * (f1_raw + g1_raw) * L1
            w2_raw = (1.0 - η) * (f2_raw + g2_raw) * L2
            dv_bp = DorsiventralBirthProfile(
                w1_raw, w2_raw,
                η * a, f1_raw, g1_raw,
                (1.0 - η) * a, f2_raw, g2_raw,
            )
            dv_fates = leaf_emission_fates(dv_esc, dv_bp)
            @test isapprox(dv_fates.front, ref_fates.front; rtol = 1e-12, atol = 1e-14)
            @test isapprox(dv_fates.back, ref_fates.back; rtol = 1e-12, atol = 1e-14)
            @test isapprox(dv_fates.reabsorbed, ref_fates.reabsorbed; rtol = 1e-12, atol = 1e-14)
        end

        # Exact anatomical reversal symmetry: swapping (adaxial <-> abaxial) swaps (front <-> back) and preserves reabsorbed!
        esc_fwd = DorsiventralPlateEscape(1.3, 0.4, 0.41, 0.49, 0.43)
        esc_rev = DorsiventralPlateEscape(0.4, 1.3, 0.49, 0.41, 0.43)
        bp_fwd = DorsiventralBirthProfile(0.78, 0.22, 2.4, 0.9, 0.1, 0.7, 0.85, 0.15)
        bp_rev = DorsiventralBirthProfile(0.22, 0.78, 0.7, 0.15, 0.85, 2.4, 0.1, 0.9)
        f_fwd = leaf_emission_fates(esc_fwd, bp_fwd)
        f_rev = leaf_emission_fates(esc_rev, bp_rev)
        @test f_fwd.front == f_rev.back
        @test f_fwd.back == f_rev.front
        @test isapprox(f_fwd.reabsorbed, f_rev.reabsorbed; rtol = 1e-14, atol = 1e-15)

        # Subnormal & extreme optical-depth stress tests against independent 256-bit BigFloat linear solver
        setprecision(BigFloat, 256) do
            for FT in (Float32, Float64)
                rtol = FT === Float32 ? FT(5e-5) : FT(1e-11)
                atol = FT === Float32 ? FT(5e-6) : FT(1e-13)
                sub_t = nextfloat(zero(FT))

                # Equal subnormal trapping limit (b1 = b2 = 0, t1 = t2 = t_int = nextfloat(0)):
                # exact random-walk escape from zone 1 is 2/3 front, 1/3 back!
                esc_sub = DorsiventralPlateEscape{FT}(zero(FT), zero(FT), sub_t, sub_t, sub_t)
                f_sub = leaf_emission_fates(esc_sub, DorsiventralBirthProfile{FT}(one(FT), zero(FT)))
                @test f_sub.reabsorbed === zero(FT)
                @test isapprox(f_sub.front, FT(2) / FT(3); rtol = rtol, atol = atol)
                @test isapprox(f_sub.back, FT(1) / FT(3); rtol = rtol, atol = atol)

                for (b1, b2, t1, t2, ti, w1, w2, a1, a2) in (
                    (zero(FT), zero(FT), FT(0.45), FT(0.45), FT(0.35), FT(0.8), FT(0.2), zero(FT), zero(FT)),
                    (FT(1e-5), FT(2e-5), sub_t, FT(0.4), FT(0.5), FT(0.7), FT(0.3), FT(1.5), FT(0.4)),
                    (FT(0.8), FT(0.3), sub_t, sub_t, FT(0.42), FT(0.75), FT(0.25), FT(3.0), FT(0.8)),
                    (FT(35.0), FT(0.2), FT(0.46), FT(0.44), FT(0.5), FT(0.6), FT(0.4), FT(25.0), FT(1.2)),
                )
                    esc_t = DorsiventralPlateEscape{FT}(b1, b2, t1, t2, ti)
                    bp_t = DorsiventralBirthProfile{FT}(w1, w2, a1, FT(0.85), FT(0.15), a2, FT(0.9), FT(0.1))
                    res_t = leaf_emission_fates(esc_t, bp_t)
                    @test res_t.front >= zero(FT) && res_t.back >= zero(FT) && res_t.reabsorbed >= zero(FT)
                    @test isapprox(res_t.front + res_t.back + res_t.reabsorbed, one(FT); rtol = rtol, atol = atol)

                    esc_bf = DorsiventralPlateEscape{BigFloat}(BigFloat(b1), BigFloat(b2), BigFloat(t1), BigFloat(t2), BigFloat(ti))
                    bp_bf = DorsiventralBirthProfile{BigFloat}(
                        BigFloat(w1), BigFloat(w2),
                        BigFloat(a1), BigFloat(FT(0.85)), BigFloat(FT(0.15)),
                        BigFloat(a2), BigFloat(FT(0.9)), BigFloat(FT(0.1)),
                    )
                    res_bf = leaf_emission_fates(esc_bf, bp_bf)
                    @test isapprox(res_t.front, FT(res_bf.front); rtol = rtol, atol = atol)
                    @test isapprox(res_t.back, FT(res_bf.back); rtol = rtol, atol = atol)
                    @test isapprox(res_t.reabsorbed, FT(res_bf.reabsorbed); rtol = rtol, atol = atol)
                end
            end
        end

        # Host helper dorsiventral_plate_optics + physical 685 nm vs 740 nm dorsiventral SIF escape asymmetry
        dpo_match = dorsiventral_plate_optics(
            0.65, 0.65, 1.41;
            surface_adaxial = ConeSurface(90.0, 1.0),
            surface_abaxial = ConeSurface(90.0, 1.0),
            hemispherical = true,
            internal_air_gap = false,
        )
        ref_dpo = diffuse_plate_optics(0.65, 1.41)
        @test isapprox(dpo_match.adaxial_effective_depth, ref_dpo.effective_depth; rtol = 1e-13)
        @test isapprox(dpo_match.adaxial_exit_transmission, ref_dpo.interface_exit_transmission; rtol = 1e-13)
        @test dpo_match.interface_transmission == 1.0

        dpo_gap = dorsiventral_plate_optics(0.65, 0.25, 1.41; internal_air_gap = true)
        @test 0.35 < dpo_gap.interface_transmission < 0.55

        # At 685 nm (strong palisade Chl reabsorption: k_pal = 1.8, k_spo = 0.6) vs 740 nm (weak reabsorption: k_pal = 0.04, k_spo = 0.015)
        # under adaxial blue excitation (w_pal = 0.88, w_spo = 0.12), front/back escape asymmetry is strongly amplified at 685 nm!
        dpo_685 = dorsiventral_plate_optics(1.8, 0.6, 1.41; internal_air_gap = true)
        dpo_740 = dorsiventral_plate_optics(0.04, 0.015, 1.41; internal_air_gap = true)
        bp_blue = DorsiventralBirthProfile(0.88, 0.12, 3.2, 1.0, 0.02, 0.9, 1.0, 0.05)
        f685 = leaf_emission_fates(dpo_685.escape, bp_blue)
        f740 = leaf_emission_fates(dpo_740.escape, bp_blue)
        @test (f685.front / f685.back) > 2.0 * (f740.front / f740.back)

        # CUDA device broadcast verification on NVIDIA L40S
        if CanopyOptics.CUDA.functional()
            CanopyOptics.CUDA.allowscalar(false)
            for FT in (Float32, Float64)
                escs_cpu = [
                    DorsiventralPlateEscape{FT}(FT(0.8), FT(0.3), FT(0.46), FT(0.44), one(FT)),
                    DorsiventralPlateEscape{FT}(FT(1.6), FT(0.5), FT(0.45), FT(0.42), FT(0.43)),
                    DorsiventralPlateEscape{FT}(zero(FT), zero(FT), nextfloat(zero(FT)), FT(0.45), FT(0.5)),
                ]
                bps_cpu = [
                    DorsiventralBirthProfile{FT}(FT(0.75), FT(0.25), FT(2.5), FT(0.9), FT(0.1), FT(0.8), FT(0.9), FT(0.1)),
                    DorsiventralBirthProfile{FT}(FT(0.85), FT(0.15), FT(3.5), one(FT), FT(0.02), FT(1.1), one(FT), FT(0.05)),
                    DorsiventralBirthProfile{FT}(FT(0.6), FT(0.4)),
                ]
                cpu_out = leaf_emission_fates.(escs_cpu, bps_cpu)
                gpu_out = Array(leaf_emission_fates.(CanopyOptics.CUDA.CuArray(escs_cpu), CanopyOptics.CUDA.CuArray(bps_cpu)))
                tol_gpu = FT === Float32 ? FT(5e-6) : FT(1e-13)
                for k in eachindex(cpu_out)
                    @test isapprox(gpu_out[k].front, cpu_out[k].front; rtol = tol_gpu, atol = tol_gpu)
                    @test isapprox(gpu_out[k].back, cpu_out[k].back; rtol = tol_gpu, atol = tol_gpu)
                    @test isapprox(gpu_out[k].reabsorbed, cpu_out[k].reabsorbed; rtol = tol_gpu, atol = tol_gpu)
                end
            end
        end
    end
end
