using CanopyOptics, Test, LinearAlgebra

@testset "Shared leaf-area quadrature" begin
    for FT in (Float32, Float64)
        μ, wμ = CanopyOptics.gauleg(128, zero(FT), one(FT))
        # The outgoing-stream integration is separate from the leaf quadrature.
        tolerance = FT === Float32 ? FT(5e-5) : FT(5e-6)
        for (α, β) in ((FT(2), FT(2)), (FT(0.4), FT(0.7)))
            lad = beta_leaves(α, β)
            for n_leaf in (8, 24, 128)
                q = CanopyQuadrature(n_leaf=n_leaf)
                leaf = CanopyOptics._leaf_inclination_quadrature(lad,n_leaf,FT)
                @test eltype(leaf.w_measure) === FT
                @test sum(leaf.w_measure) ≈ one(FT)
                @test all(>=(zero(FT)),leaf.w_measure)
                @test CanopyOptics.G(μ,lad;quadrature=q) ≈ CanopyOptics.G(μ,lad;nLeg=n_leaf)
                # Conservatively scattering leaves must not create/absorb photons
                # because the inclination quadrature misses some PDF mass.
                scattering = BiLambertianCanopyScattering(R=FT(0.4),T=FT(0.6))
                Zpp,Zmp = CanopyOptics.compute_Z_matrices(scattering,μ,lad,0;quadrature=q)
                columns = vec(transpose(wμ)*(Zpp+Zmp))/FT(2)
                @test maximum(abs.(columns.-one(FT))) < tolerance
                @test eltype(Zpp) === FT
            end
        end
        # Default and explicit quadrature must describe the same leaf population.
        lad = beta_leaves(FT(2),FT(2))
        @test CanopyOptics.G(μ,lad) == CanopyOptics.G(μ,lad;quadrature=CanopyQuadrature())
        # One orientation PDF normalizes both the 4SAIL projection and analytic Z.
        q = CanopyQuadrature(n_leaf=24)
        geom = FourSAILGeometry(lad;sza_deg=FT(40),vza_deg=zero(FT),raa_deg=zero(FT),quadrature=q)
        @test geom.ks*geom.μs ≈ only(CanopyOptics.G(FT[geom.μs],lad;quadrature=q))
    end
end
