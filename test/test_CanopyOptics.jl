using LinearAlgebra

_relerr(A, B) = norm(A .- B) / max(norm(B), eps(Float64))

@testset "CanopyOptics" begin
    @info "Testing spherical leaf G factor calculation ...";
    @testset "G function" begin
        for FT in [Float32, Float64]
            μ, w = CanopyOptics.gauleg(10, FT(0.0), FT(1.0));
            LD = CanopyOptics.spherical_leaves(FT)
            G = CanopyOptics.G(μ, LD)
            @test eltype(G) == FT
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
        @test all(Z⁺⁺_m1 .== 0)
        @test all(Z⁻⁺_m1 .== 0)
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
        diffuse = CanopyOptics.BiLambertianCanopyScattering(R = 0.4, T = 0.2, nQuad = 32)
        specular = CanopyOptics.SpecularCanopyScattering(nᵣ = 1.5, κ = 0.2, nQuad = 12)
        composite = diffuse + specular

        @test composite isa CanopyOptics.CompositeCanopyScattering
        @test composite.components == (diffuse, specular)
        @test (diffuse + (specular + diffuse)).components == (diffuse, specular, diffuse)

        Zpp, Zmp = CanopyOptics.compute_Z_matrices(composite, μ, LD, 0)
        Zpp_d, Zmp_d = CanopyOptics.compute_Z_matrices(diffuse, μ, LD, 0)
        Zpp_s, Zmp_s = CanopyOptics.compute_Z_matrices(specular, μ, LD, 0)
        @test Zpp ≈ Zpp_d .+ Zpp_s
        @test Zmp ≈ Zmp_d .+ Zmp_s

        Zpp_a, Zmp_a = CanopyOptics.compute_Z_matrices_aniso(composite, μ, LD, 2)
        Zpp_ad, Zmp_ad = CanopyOptics.compute_Z_matrices_aniso(diffuse, μ, LD, 2)
        Zpp_as, Zmp_as = CanopyOptics.compute_Z_matrices_aniso(specular, μ, LD, 2)
        @test Zpp_a ≈ Zpp_ad .+ Zpp_as
        @test Zmp_a ≈ Zmp_ad .+ Zmp_as

        Zpp_stack, Zmp_stack = CanopyOptics.compute_Z_matrices_aniso_analytic(composite, μ, LD, 3)
        @test Zpp_stack[:, :, 3] ≈ Zpp_a
        @test Zmp_stack[:, :, 3] ≈ Zmp_a
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

        @testset "clipped projection moments" begin
            for μ in grid, μ_L in grid
                P, N = CanopyOptics._clipped_projection_moments(μ, μ_L, 32)

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
                mod = CanopyOptics.BiLambertianCanopyScattering(R = R, T = T, nQuad = 64)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(mod, μ, LD, 0)
                Zpp_ref, Zmp_ref = CanopyOptics.compute_Z_matrices(mod, μ, LD, 0)

                @test _relerr(Zpp[:, :, 1], Zpp_ref) ≤ 1e-12
                @test _relerr(Zmp[:, :, 1], Zmp_ref) ≤ 1e-12
            end
        end

        @testset "legacy aniso signatures delegate to analytic closure" begin
            μb, _ = CanopyOptics.gauleg(4, 0.0, 1.0)
            for LD in LADs, (R, T) in RTs
                mod = CanopyOptics.BiLambertianCanopyScattering(R = R, T = T, nQuad = 64)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(mod, μb, LD, 4)

                for m in 0:4
                    Zpp_single, Zmp_single = CanopyOptics.compute_Z_matrices_aniso(
                        mod, μb, LD, m)
                    Zpp_compat, Zmp_compat = CanopyOptics.compute_Z_matrices_aniso(
                        mod, μb, LD, nothing, nothing, m)

                    @test Zpp[:, :, m + 1] == Zpp_single == Zpp_compat
                    @test Zmp[:, :, m + 1] == Zmp_single == Zmp_compat
                end
            end
        end

        @testset "vSmartMOM normalization and reciprocity" begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.5, T = 0.5, nQuad = 64)
            LD = CanopyOptics.spherical_leaves()
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(mod, μ, LD, 32)

            flux = vec(sum(w .* (Zpp[:, :, 1] .+ Zmp[:, :, 1]), dims = 1))
            @test all(abs.(flux .- 2) .≤ 2e-3)

            for LD in LADs
                mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05, nQuad = 64)
                Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(mod, μ, LD, 32)
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
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.45, T = 0.05, nQuad = 48)
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                mod, μ, CanopyOptics.spherical_leaves(), 64)

            @test all(isfinite, Zpp)
            @test all(isfinite, Zmp)
            @test maximum(abs.(Zpp)) < 10
            @test maximum(abs.(Zmp)) < 10
            @test all(diag(Zpp[:, :, 1]) .≥ 0)
            @test all(diag(Zmp[:, :, 1]) .≥ 0)
        end

        @testset "fully absorbing leaf limit" begin
            mod = CanopyOptics.BiLambertianCanopyScattering(R = 0.0, T = 0.0, nQuad = 16)
            Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
                mod, μ, CanopyOptics.spherical_leaves(), 4)

            @test all(iszero, Zpp)
            @test all(iszero, Zmp)
        end
    end
end
