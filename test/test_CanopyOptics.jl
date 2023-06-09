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
end
