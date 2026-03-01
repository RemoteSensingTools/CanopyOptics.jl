@testset "CanopyOptics" begin
    @info "Testing spherical leaf G factor calculation ...";
    @testset "G function" begin
        for FT in [Float32, Float64]
            μ,w = CanopyOptics.gauleg(10,FT(0.0),FT(1.0));
            LD  = CanopyOptics.spherical_leaves(FT)
            G   = CanopyOptics.G(μ, LD)
            @test eltype(G) == FT;
            # Should be about 0.5 for spherical leaves:
            @test all(abs.(G .- 0.5) .< 0.001)

        end
    end

    @testset "PROSPECT refl and gradient" begin
        function test_prospect(θ::AbstractVector{T}) where {T}
            opti = createLeafOpticalStruct((400.0:1.0:2500.0))
            leaf = LeafProspectProProperties{T}(Ccab=θ[1], Ccar=θ[2], Canth=θ[3])
            _, R = prospect(leaf, opti)
            return R
        end
        refl = test_prospect([40.0, 8.0, 4.0])
        @test all(refl .>= 0)
        @test all(refl .<= 1)
        jac = ForwardDiff.jacobian(test_prospect, [40.0, 8.0, 4.0])
        @test all(isfinite.(jac))
        # Increasing pigments decreases reflectance in visible
        @test all(jac[1:100,:] .<= 0)
        # ...but has no effect in the SWIR
        @test all(jac[(end-100):end,:] .== 0)
    end

    @testset "BiLambertian Z-matrices" begin
        μ, w = CanopyOptics.gauleg(10, 0.0, 1.0)
        LD   = CanopyOptics.spherical_leaves()
        mod  = CanopyOptics.BiLambertianCanopyScattering(R=0.1, T=0.05)
        Z⁺⁺, Z⁻⁺ = CanopyOptics.compute_Z_matrices(mod, μ, LD, 0)
        # Phase matrix elements must be non-negative and finite
        @test all(Z⁺⁺ .>= 0)
        @test all(Z⁻⁺ .>= 0)
        @test all(isfinite.(Z⁺⁺))
        @test all(isfinite.(Z⁻⁺))
        # m>0 returns zero matrices (azimuthally uniform distribution)
        Z⁺⁺_m1, Z⁻⁺_m1 = CanopyOptics.compute_Z_matrices(mod, μ, LD, 1)
        @test all(Z⁺⁺_m1 .== 0)
        @test all(Z⁻⁺_m1 .== 0)
    end

    @testset "Specular compute_reflection symmetry" begin
        mod = CanopyOptics.SpecularCanopyScattering(nᵣ=1.5, κ=0.2)
        LD  = CanopyOptics.spherical_leaves()
        # Use same azimuth (ϕ=0) to avoid the known azimuth-sign edge case in getSpecularΩ
        Ωᵢ  = CanopyOptics.dirVector_μ(0.7, 0.0)
        Ωₒ  = CanopyOptics.dirVector_μ(0.5, 0.0)
        fᵢₒ = CanopyOptics.compute_reflection(mod, Ωᵢ, Ωₒ, LD)
        fₒᵢ = CanopyOptics.compute_reflection(mod, Ωₒ, Ωᵢ, LD)
        # Fresnel reflection is symmetric: f(Ωᵢ→Ωₒ) = f(Ωₒ→Ωᵢ)
        @test fᵢₒ ≈ fₒᵢ
        # Must be non-negative
        @test fᵢₒ >= 0
    end

    @testset "dielectric sanity" begin
        # Liquid water: real and imaginary parts both positive
        w = CanopyOptics.LiquidPureWater()
        ϵ_w = CanopyOptics.dielectric(w, 283.0, 10.0)
        @test real(ϵ_w) > 0
        @test imag(ϵ_w) > 0
        # Pure ice at -20°C: real part ≈ 3.17, very small imaginary part
        ice = CanopyOptics.PureIce()
        ϵ_i = CanopyOptics.dielectric(ice, 253.0, 10.0)
        @test real(ϵ_i) > 0
        @test imag(ϵ_i) > 0
        # Ice has far lower imaginary part than liquid water at 10 GHz
        @test imag(ϵ_i) < imag(ϵ_w)
    end
end
