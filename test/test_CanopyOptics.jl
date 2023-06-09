@info "Testing spherical leaf G factor calculation ...";
@testset "CanopyOptics --- G function" begin
    for FT in [Float32, Float64]
        μ,w = CanopyOptics.gauleg(10,FT(0.0),FT(1.0));
        LD  = CanopyOptics.spherical_leaves(FT)
        G   = CanopyOptics.G(μ, LD)
        @test eltype(G) == FT;
        # Should be about 0.5 for spherical leaves:
        @test all(abs.(G .- 0.5) .< 0.001)

    end
end
