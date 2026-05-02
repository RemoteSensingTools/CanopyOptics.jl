using LinearAlgebra
using ForwardDiff

_relerr(A, B) = norm(A .- B) / max(norm(B), eps(Float64))

@testset "CanopyOptics" begin
    @testset "Public canopy API exports" begin
        LD = spherical_leaves()
        @test LD isa AbstractLeafDistribution
        @test G([0.5], LD)[1] > 0
        @test planophile_leaves2() isa LeafDistribution
    end

    @info "Testing spherical leaf G factor calculation ...";
    @testset "G function" begin
        for FT in [Float32, Float64]
            μ, w = CanopyOptics.gauleg(10, FT(0.0), FT(1.0));
            LD = CanopyOptics.spherical_leaves(FT)
            G = CanopyOptics.G(μ, LD)
            @test eltype(G) == FT
            @test CanopyOptics.G(-μ, LD) ≈ G
            # The default G quadrature (nLeg=40) is accurate to about
            # 1.3e-3 for this 10-point μ grid.
            tol = FT == Float32 ? FT(2e-3) : FT(1.5e-3)
            @test all(abs.(G .- FT(0.5)) .< tol)
        end
    end

    @testset "PROSPECT refl and gradient" begin
        function test_prospect(θ::AbstractVector{T}) where {T}
            opti = createLeafOpticalStruct((400.0:1.0:2500.0))
            leaf = LeafProspectProProperties{T}(Ccab = θ[1], Ccar = θ[2], Canth = θ[3])
            _, R = prospect(leaf, opti)
            return R
        end
        refl = test_prospect([40.0, 8.0, 4.0])
        @test all(refl .>= 0)
        @test all(refl .<= 1)
        jac = ForwardDiff.jacobian(test_prospect, [40.0, 8.0, 4.0])
        @test all(isfinite.(jac))
        # Increasing pigments decreases reflectance in visible.
        @test all(jac[1:100, :] .<= 0)
        # ...but has no effect in the SWIR.
        @test all(jac[(end - 100):end, :] .== 0)
    end

    @testset "BiLambertian Z-matrices" begin
        μ, w = CanopyOptics.gauleg(10, 0.0, 1.0)
        LD = CanopyOptics.spherical_leaves()
        mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.1, T = 0.05)
        Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(mod, μ, LD, 0)
        @test all(Z⁺⁺ .>= 0)
        @test all(Z⁻⁺ .>= 0)
        @test all(isfinite.(Z⁺⁺))
        @test all(isfinite.(Z⁻⁺))

        Z⁺⁺_m1, Z⁻⁺_m1 = CanopyOptics.compute_Z_matrices(mod, μ, LD, 1)
        Z⁺⁺_a1, Z⁻⁺_a1 = CanopyOptics.compute_Z_matrices_aniso(mod, μ, LD, 1)
        @test Z⁺⁺_m1 == Z⁺⁺_a1
        @test Z⁻⁺_m1 == Z⁻⁺_a1

        Z⁺⁺_stack, Z⁻⁺_stack = CanopyOptics.compute_Z_matrices(mod, μ, LD, 0:2)
        @test Z⁺⁺_stack[:, :, 2] == Z⁺⁺_m1
        @test Z⁻⁺_stack[:, :, 2] == Z⁻⁺_m1
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices(mod, μ, LD, -1:1)
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices(mod, μ, LD, 2:1)
    end

    @testset "Canopy quadrature controls" begin
        quad = CanopyOptics.CanopyQuadrature(n_leaf = 32, n_azimuth = 12)
        @test quad.n_leaf == 32
        @test quad.n_azimuth == 12
        @test CanopyOptics.CanopyQuadrature(nQuad = 7) == CanopyOptics.CanopyQuadrature(7, 7)
        @test_throws ArgumentError CanopyOptics.CanopyQuadrature(n_leaf = 0)
        @test_throws ArgumentError CanopyOptics.CanopyQuadrature(n_azimuth = 0)

        diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2, nQuad = 32)
        specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2, nQuad = 12)
        @test !hasproperty(diffuse, :nQuad)
        @test !hasproperty(specular, :nQuad)
    end

    @testset "Canopy clumping models" begin
        μ = [0.2, 0.6, 1.0]
        G = [0.45, 0.50, 0.55]

        no_clumping = CanopyOptics.NoClumping()
        @test CanopyOptics.clumping_index(no_clumping, μ) == ones(length(μ))
        @test CanopyOptics.effective_G(no_clumping, G, μ) == G

        constant = CanopyOptics.ConstantClumping(Ω = 0.7)
        @test CanopyOptics.clumping_index(constant, μ) == fill(0.7, length(μ))
        @test CanopyOptics.effective_G(constant, G, μ) == 0.7 .* G

        angular = CanopyOptics.EmpiricalDirectionalClumping(Ω₀ = 0.7, c = 2.0, e = 2.0)
        @test (CanopyOptics.ChenLeblancClumping(Ω₀ = 0.7, c = 2.0, e = 2.0) isa
               CanopyOptics.EmpiricalDirectionalClumping)
        Ω = CanopyOptics.clumping_index(angular, [1.0, 0.5, 0.0])
        @test Ω[1] ≈ 0.7
        @test Ω[1] < Ω[2] < Ω[3] <= 1
        @test CanopyOptics.clumping_index(angular, -0.5) ≈
              CanopyOptics.clumping_index(angular, 0.5)

        μ_grad = [0.2, 0.6, 0.9]
        G_grad = [0.45, 0.50, 0.55]
        clumping_sum(x) = begin
            model = CanopyOptics.ChenLeblancClumping(Ω₀ = x[1], c = x[2], e = x[3])
            sum(CanopyOptics.effective_G(model, G_grad, μ_grad))
        end
        @test all(isfinite, ForwardDiff.gradient(clumping_sum, [0.7, 2.0, 2.0]))
    end

    @testset "Canopy hotspot correction" begin
        FT = Float64
        k_s = FT(0.5) / FT(0.7)
        k_o = FT(0.5) / FT(0.7)
        L = FT(2)

        no_hotspot = CanopyOptics.NoHotSpot()
        base = exp(-(k_s + k_o) * L)
        @test CanopyOptics.joint_gap_probability(
            no_hotspot, k_s, k_o, 0.7, 0.7, 0.0, L) ≈ base

        disabled = CanopyOptics.KuuskHotSpot(h = 0.0)
        @test CanopyOptics.hotspot_correction(
            disabled, k_s, k_o, 0.7, 0.7, 0.0, L) ≈ 1

        hotspot = CanopyOptics.KuuskHotSpot(h = 0.1)
        exact_backscatter = CanopyOptics.joint_gap_probability(
            hotspot, k_s, k_o, 0.7, 0.7, 0.0, L)
        @test exact_backscatter ≈ exp(-k_s * L)

        near = CanopyOptics.hotspot_correction(
            hotspot, k_s, k_o, 0.7, 0.7, 0.0, L)
        far = CanopyOptics.hotspot_correction(
            hotspot, k_s, k_o, 0.7, 0.4, π, L)
        @test near > far > 1

        @test CanopyOptics.canopy_extinction(0.5f0, 0.25f0) === 2.0f0

        hs_sum(x) = begin
            model = CanopyOptics.KuuskHotSpot(h = x[1])
            CanopyOptics.joint_gap_probability(model, k_s, k_o, 0.7, 0.5, 0.2, L)
        end
        @test isfinite(ForwardDiff.derivative(h -> hs_sum([h]), 0.1))
    end

    @testset "Specular compute_reflection symmetry" begin
        mod = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2)
        LD = CanopyOptics.spherical_leaves()
        # Use same azimuth to avoid the known azimuth-sign edge case in getSpecularΩ.
        Ωᵢ = CanopyOptics.dirVector_μ(0.7, 0.0)
        Ωₒ = CanopyOptics.dirVector_μ(0.5, 0.0)
        fᵢₒ = CanopyOptics.compute_reflection(mod, Ωᵢ, Ωₒ, LD)
        fₒᵢ = CanopyOptics.compute_reflection(mod, Ωₒ, Ωᵢ, LD)
        @test fᵢₒ ≈ fₒᵢ
        @test fᵢₒ >= 0
    end

    @testset "Composite canopy scattering" begin
        μ, w = CanopyOptics.gauleg(5, 0.0, 1.0)
        LD = CanopyOptics.planophile_leaves2()
        diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2)
        diffuse2 = CanopyOptics.BiLambertianCanopyScattering(R = 0.1, T = 0.1)
        specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2)
        quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 32, n_azimuth = 12)
        composite = diffuse + diffuse2

        @test composite isa CanopyOptics.CompositeCanopyScattering
        @test composite.components == (diffuse, diffuse2)
        @test (diffuse + (diffuse2 + diffuse)).components == (diffuse, diffuse2, diffuse)

        Zpp, Zmp = CanopyOptics.compute_Z_matrices(composite, μ, LD, 0; quadrature)
        Zpp_d, Zmp_d = CanopyOptics.compute_Z_matrices(diffuse, μ, LD, 0; quadrature)
        Zpp_d2, Zmp_d2 = CanopyOptics.compute_Z_matrices(diffuse2, μ, LD, 0; quadrature)
        @test Zpp ≈ Zpp_d .+ Zpp_d2
        @test Zmp ≈ Zmp_d .+ Zmp_d2

        Zpp_a, Zmp_a = CanopyOptics.compute_Z_matrices_aniso(composite, μ, LD, 2; quadrature)
        Zpp_ad, Zmp_ad = CanopyOptics.compute_Z_matrices_aniso(diffuse, μ, LD, 2; quadrature)
        Zpp_ad2, Zmp_ad2 = CanopyOptics.compute_Z_matrices_aniso(diffuse2, μ, LD, 2; quadrature)
        @test Zpp_a ≈ Zpp_ad .+ Zpp_ad2
        @test Zmp_a ≈ Zmp_ad .+ Zmp_ad2

        Zpp_stack, Zmp_stack = CanopyOptics.compute_Z_matrices_aniso_analytic(
            composite, μ, LD, 3; quadrature)
        @test Zpp_stack[:, :, 3] ≈ Zpp_a
        @test Zmp_stack[:, :, 3] ≈ Zmp_a

        Zpp_public_stack, Zmp_public_stack = CanopyOptics.compute_Z_matrices(
            composite, μ, LD, 0:3; quadrature)
        @test Zpp_public_stack == Zpp_stack
        @test Zmp_public_stack == Zmp_stack

        mixed_convention = diffuse + specular
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices(
            mixed_convention, μ, LD, 0; quadrature)
    end

    @testset "Canopy Stokes expansion" begin
        μ = [0.35, 0.7]
        LD = CanopyOptics.spherical_leaves()
        quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 24, n_azimuth = 12)
        diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2)
        specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2)

        Zpp_I, Zmp_I = CanopyOptics.compute_Z_matrices(
            diffuse, μ, LD, 0:1; quadrature)
        Zpp_4, Zmp_4 = CanopyOptics.compute_Z_matrices(
            diffuse, μ, LD, 0:1; quadrature, npol = 4)

        @test size(Zpp_4) == (4 * length(μ), 4 * length(μ), 2)
        @test Zpp_4[1:4:end, 1:4:end, :] == Zpp_I
        @test Zmp_4[1:4:end, 1:4:end, :] == Zmp_I
        for si in 1:4, sj in 1:4
            (si, sj) == (1, 1) && continue
            @test all(iszero, Zpp_4[si:4:end, sj:4:end, :])
            @test all(iszero, Zmp_4[si:4:end, sj:4:end, :])
        end

        Ωin = CanopyOptics.dirVector_μ(0.7, 0.0)
        Ωout = CanopyOptics.dirVector_μ(-0.4, 0.5)
        M = CanopyOptics.compute_reflection_mueller(specular, Ωin, Ωout, LD, 4)
        @test M[1, 1] ≈ CanopyOptics.compute_reflection(specular, Ωin, Ωout, LD)
        @test abs(M[2, 1]) > 1e-12

        Zpp_s1, Zmp_s1 = CanopyOptics.compute_Z_matrices(
            specular, μ, LD, 0; quadrature)
        Zpp_s4, Zmp_s4 = CanopyOptics.compute_Z_matrices(
            specular, μ, LD, 0; quadrature, npol = 4)
        Zpp_sstack, Zmp_sstack = CanopyOptics.compute_Z_matrices_aniso_analytic(
            specular, μ, LD, 1; quadrature)
        @test Zpp_sstack[:, :, 1] ≈ Zpp_s1
        @test Zmp_sstack[:, :, 1] ≈ Zmp_s1
        @test Zpp_s4[1:4:end, 1:4:end] ≈ Zpp_s1
        @test Zmp_s4[1:4:end, 1:4:end] ≈ Zmp_s1
        @test any(abs.(Zpp_s4[2:4:end, 1:4:end]) .> 1e-12) ||
              any(abs.(Zmp_s4[2:4:end, 1:4:end]) .> 1e-12)

        composite = diffuse + specular
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices(
            composite, μ, LD, 0; quadrature, npol = 4)
    end

    @testset "Wood reflectance and Lambertian wood scattering" begin
        μ, _ = CanopyOptics.gauleg(5, 0.0, 1.0)
        LD = CanopyOptics.spherical_leaves()
        quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 24)

        constant = CanopyOptics.ConstantWoodReflectance(R = 0.23)
        @test CanopyOptics.wood_reflectance(constant) == 0.23
        @test CanopyOptics.wood_reflectance(constant, [450.0, 550.0]) == [0.23, 0.23]

        lut = CanopyOptics.LUTWoodReflectance([400.0, 500.0, 600.0],
                                             [0.10, 0.20, 0.40])
        @test CanopyOptics.wood_reflectance(lut, 450.0) ≈ 0.15
        @test CanopyOptics.wood_reflectance(lut, 1e7 / 500.0; grid_unit = :cm_inv) ≈ 0.20
        @test CanopyOptics.wood_reflectance(lut, 350.0) ≈ 0.10

        poly = CanopyOptics.PolynomialWoodReflectance(
            coeffs = [0.20, 0.05], x_offset = 400.0, x_scale = 100.0)
        @test CanopyOptics.wood_reflectance(poly, 500.0) ≈ 0.25

        poly_sum(c) = begin
            p = CanopyOptics.PolynomialWoodReflectance(
                coeffs = c, x_offset = 400.0, x_scale = 100.0)
            sum(CanopyOptics.wood_reflectance(p, [400.0, 500.0, 600.0]))
        end
        @test all(isfinite, ForwardDiff.gradient(poly_sum, [0.2, 0.05, 0.01]))

        wood = CanopyOptics.LambertianWoodCanopyScattering(constant)
        diffuse_ref = CanopyOptics.BiLambertianCanopyScattering(R = 0.23, T = 0.0)
        Zpp_w, Zmp_w = CanopyOptics.compute_Z_matrices(
            wood, μ, LD, 0:2; quadrature)
        Zpp_d, Zmp_d = CanopyOptics.compute_Z_matrices(
            diffuse_ref, μ, LD, 0:2; quadrature)
        @test Zpp_w == Zpp_d
        @test Zmp_w == Zmp_d

        spectral_wood = CanopyOptics.LambertianWoodCanopyScattering(reflectance = lut)
        @test_throws ArgumentError CanopyOptics.compute_Z_matrices(
            spectral_wood, μ, LD, 0; quadrature)
        Zpp_lut, Zmp_lut = CanopyOptics.compute_Z_matrices(
            spectral_wood, μ, LD, 0; quadrature, spectral_coordinate = 500.0)
        Zpp_ref, Zmp_ref = CanopyOptics.compute_Z_matrices(
            CanopyOptics.BiLambertianCanopyScattering(R = 0.20, T = 0.0),
            μ, LD, 0; quadrature)
        @test Zpp_lut == Zpp_ref
        @test Zmp_lut == Zmp_ref
    end

    @testset "Mixed canopy components" begin
        μ, _ = CanopyOptics.gauleg(4, 0.0, 1.0)
        μ = collect(μ)
        leaf_LD = CanopyOptics.spherical_leaves()
        wood_LD = CanopyOptics.erectophile_leaves()
        quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 24, n_azimuth = 8)

        leaf = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05)
        wood = CanopyOptics.LambertianWoodCanopyScattering(R = 0.20)

        leaf_component = CanopyOptics.CanopyComponent(
            scatterer = leaf, LAD = leaf_LD, area_index = 2.0)
        single = CanopyOptics.MixedCanopy(leaf_component)
        Zpp_single, Zmp_single = CanopyOptics.compute_Z_matrices(
            single, μ, 0:2; quadrature, npol = 4)
        Zpp_direct, Zmp_direct = CanopyOptics.compute_Z_matrices(
            leaf, μ, leaf_LD, 0:2; quadrature, npol = 4)
        @test Zpp_single ≈ Zpp_direct
        @test Zmp_single ≈ Zmp_direct

        wood_component = CanopyOptics.CanopyComponent(
            scatterer = wood, LAD = wood_LD, AI = 1.0,
            clumping = CanopyOptics.ConstantClumping(Ω = 0.5))
        mixed = CanopyOptics.MixedCanopy(leaf_component, wood_component)

        G_leaf = vec(CanopyOptics.G(μ, leaf_LD))
        G_wood = vec(CanopyOptics.G(μ, wood_LD))
        @test CanopyOptics.bulk_G(mixed, μ; clumped = false) ≈
              2.0 .* G_leaf .+ G_wood
        @test CanopyOptics.bulk_G(mixed, μ; clumped = true) ≈
              2.0 .* G_leaf .+ 0.5 .* G_wood

        Zpp_mix, Zmp_mix = CanopyOptics.compute_Z_matrices(
            mixed, μ, 0; quadrature)
        Zpp_leaf, Zmp_leaf = CanopyOptics.compute_Z_matrices(
            leaf, μ, leaf_LD, 0; quadrature)
        Zpp_wood, Zmp_wood = CanopyOptics.compute_Z_matrices(
            wood, μ, wood_LD, 0; quadrature)
        G_total = 2.0 .* G_leaf .+ G_wood

        Zpp_expected = zero.(Zpp_mix)
        Zmp_expected = zero.(Zmp_mix)
        for j in eachindex(μ)
            w_leaf = 2.0 * G_leaf[j] / G_total[j]
            w_wood = G_wood[j] / G_total[j]
            Zpp_expected[:, j] .= w_leaf .* Zpp_leaf[:, j] .+
                                  w_wood .* Zpp_wood[:, j]
            Zmp_expected[:, j] .= w_leaf .* Zmp_leaf[:, j] .+
                                  w_wood .* Zmp_wood[:, j]
        end

        @test Zpp_mix ≈ Zpp_expected
        @test Zmp_mix ≈ Zmp_expected

        unclumped_wood = CanopyOptics.CanopyComponent(
            scatterer = wood, LAD = wood_LD, AI = 1.0)
        mixed_unclumped = CanopyOptics.MixedCanopy(leaf_component, unclumped_wood)
        Zpp_unclumped, Zmp_unclumped = CanopyOptics.compute_Z_matrices(
            mixed_unclumped, μ, 0; quadrature)
        @test Zpp_mix ≈ Zpp_unclumped
        @test Zmp_mix ≈ Zmp_unclumped
    end

    @testset "Canopy Z supports ForwardDiff parameters" begin
        μ, _ = CanopyOptics.gauleg(4, 0.0, 1.0)
        LD = CanopyOptics.spherical_leaves()

        diffuse_sum(x) = begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = x[1], T = x[2])
            Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(
                mod, μ, LD, 0:2; quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 12))
            sum(Z⁺⁺) + sum(Z⁻⁺)
        end
        diffuse_grad = ForwardDiff.gradient(diffuse_sum, [0.4, 0.2])
        @test all(isfinite, diffuse_grad)

        specular_sum(x) = begin
            mod = CanopyOptics.SpecularCanopyScattering(nᵣ = x[1], κ = x[2])
            Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(
                mod, μ, LD, 0:1; quadrature = CanopyOptics.CanopyQuadrature(n_azimuth = 8))
            sum(Z⁺⁺) + sum(Z⁻⁺)
        end
        specular_grad = ForwardDiff.gradient(specular_sum, [1.5, 0.2])
        @test all(isfinite, specular_grad)
    end

    @testset "dielectric sanity" begin
        w = CanopyOptics.LiquidPureWater()
        ϵ_w = CanopyOptics.dielectric(w, 283.0, 10.0)
        @test real(ϵ_w) > 0
        @test imag(ϵ_w) > 0

        ice = CanopyOptics.PureIce()
        ϵ_i = CanopyOptics.dielectric(ice, 253.0, 10.0)
        @test real(ϵ_i) > 0
        @test imag(ϵ_i) > 0
        @test imag(ϵ_i) < imag(ϵ_w)
    end

    @testset "Analytic BiLambertian canopy Fourier moments" begin
        grid = collect(range(0.0, 1.0, length = 50))

        @testset "one-sided leaf projection moments" begin
            for μ in grid, μ_L in grid
                P, N = CanopyOptics._one_sided_projection_moments(μ, μ_L, 32)

                # P₀ is Shultis-Myneni's H function.
                @test isapprox(P[1], CanopyOptics.H(μ, μ_L); rtol = 1e-14, atol = 1e-14)

                a = μ * μ_L
                b = sqrt(max(0.0, 1 - μ^2)) * sqrt(max(0.0, 1 - μ_L^2))
                for m in 0:32
                    target = m == 0 ? a : (m == 1 ? b / 2 : 0.0)
                    @test isapprox(P[m + 1] - N[m + 1], target; rtol = 1e-14, atol = 1e-14)
                end
            end
        end

        μ, w = CanopyOptics.gauleg(8, 0.0, 1.0)
        LADs = (
            CanopyOptics.spherical_leaves(),
            CanopyOptics.planophile_leaves2(),
            CanopyOptics.erectophile_leaves(),
        )
        RTs = ((0.45, 0.05), (0.5, 0.5))

        @testset "m=0 matches Shultis-Myneni Eq. 45 assembly" begin
            for LD in LADs, (R, T) in RTs
                mod = CanopyOptics.BiLambertianCanopyScattering(R = R, T = T)
                quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                    mod, μ, LD, 0; quadrature)
                Zpp_ref, Zmp_ref = CanopyOptics.compute_Z_matrices(
                    mod, μ, LD, 0; quadrature)

                @test _relerr(Zpp[:, :, 1], Zpp_ref) ≤ 1e-12
                @test _relerr(Zmp[:, :, 1], Zmp_ref) ≤ 1e-12
            end
        end

        @testset "analytic path preserves ForwardDiff derivatives" begin
            LD = CanopyOptics.spherical_leaves()
            quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 24)

            function zsum_μ(x)
                mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                    mod, x, LD, 2; quadrature)
                return sum(Zpp) + sum(Zmp)
            end
            ∂μ = ForwardDiff.gradient(zsum_μ, [0.3, 0.7])
            @test all(isfinite, ∂μ)
            @test any(abs.(∂μ) .> 0)

            function zsum_leaf(ρτ)
                mod = CanopyOptics.BiLambertianCanopyScattering(R = ρτ[1], T = ρτ[2])
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                    mod, μ, LD, 2; quadrature)
                return sum(Zpp) + sum(Zmp)
            end
            ∂leaf = ForwardDiff.gradient(zsum_leaf, [0.45, 0.05])
            @test all(isfinite, ∂leaf)
        end

        @testset "legacy aniso signatures delegate to analytic closure" begin
            μb, _ = CanopyOptics.gauleg(4, 0.0, 1.0)
            for LD in LADs, (R, T) in RTs
                mod = CanopyOptics.BiLambertianCanopyScattering(R = R, T = T)
                quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                    mod, μb, LD, 4; quadrature)
                Zpp_public, Zmp_public = CanopyOptics.compute_Z_matrices(
                    mod, μb, LD, 0:4; quadrature)
                @test Zpp_public == Zpp
                @test Zmp_public == Zmp

                for m in 0:4
                    Zpp_single, Zmp_single = CanopyOptics.compute_Z_matrices(
                        mod, μb, LD, m; quadrature)
                    Zpp_aniso, Zmp_aniso = CanopyOptics.compute_Z_matrices_aniso(
                        mod, μb, LD, m; quadrature)
                    Zpp_compat, Zmp_compat = CanopyOptics.compute_Z_matrices_aniso(
                        mod, μb, LD, nothing, nothing, m; quadrature)

                    @test Zpp[:, :, m + 1] == Zpp_single == Zpp_aniso == Zpp_compat
                    @test Zmp[:, :, m + 1] == Zmp_single == Zmp_aniso == Zmp_compat
                end
            end
        end

        @testset "vSmartMOM normalization and reciprocity" begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.5, T = 0.5)
            LD = CanopyOptics.spherical_leaves()
            quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64)
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                mod, μ, LD, 32; quadrature)

            flux = vec(sum(w .* (Zpp[:, :, 1] .+ Zmp[:, :, 1]), dims = 1))
            @test all(abs.(flux .- 2) .≤ 2e-3)

            for LD in LADs
                mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                    mod, μ, LD, 32; quadrature)
                G = vec(CanopyOptics.G(Array(μ), LD))

                for m in 0:32
                    Wpp = Zpp[:, :, m + 1] .* reshape(G, 1, :)
                    Wmp = Zmp[:, :, m + 1] .* reshape(G, 1, :)
                    @test _relerr(Wpp, transpose(Wpp)) ≤ 1e-12
                    @test _relerr(Wmp, transpose(Wmp)) ≤ 1e-12
                end
            end
        end

        @testset "high-order smoke" begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05)
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                mod, μ, CanopyOptics.spherical_leaves(), 64;
                quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 48))

            @test all(isfinite, Zpp)
            @test all(isfinite, Zmp)
            @test maximum(abs.(Zpp)) < 10
            @test maximum(abs.(Zmp)) < 10
            @test all(diag(Zpp[:, :, 1]) .≥ 0)
            @test all(diag(Zmp[:, :, 1]) .≥ 0)
        end

        @testset "fully absorbing leaf limit" begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.0, T = 0.0)
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                mod, μ, CanopyOptics.spherical_leaves(), 4;
                quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 16))

            @test all(iszero, Zpp)
            @test all(iszero, Zmp)
        end
    end
end
