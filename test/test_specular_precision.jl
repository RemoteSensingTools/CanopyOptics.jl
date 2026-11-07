using Test, CanopyOptics, LinearAlgebra, TOML

@testset "Specular azimuth preserves floating-point type" begin
    settings = TOML.parsefile(joinpath(@__DIR__, "specular_precision.toml"))
    directions, scattering = settings["directions"], settings["scattering"]
    for FT in (Float32,Float64), mu_in in FT.(directions["incoming_cosines"]),
            mu_out in FT.(directions["outgoing_cosines"]),
            fraction in FT.(directions["outgoing_azimuth_pi"])
        phi = fraction*FT(pi)
        incoming = CanopyOptics.dirVector_μ(mu_in,zero(FT))
        outgoing = CanopyOptics.dirVector_μ(mu_out,phi)
        normal,alpha = @inferred CanopyOptics.getSpecularΩ(incoming,outgoing)
        @test normal isa CanopyOptics.dirVector_μ{FT}
        @test alpha isa FT
        # Independent Cartesian bisector, with incoming azimuth zero.
        direction(d) = FT[sqrt(one(FT)-d.μ^2)*cos(d.ϕ),
            sqrt(one(FT)-d.μ^2)*sin(d.ϕ),d.μ]
        expected = normalize(direction(incoming)+direction(outgoing))
        @test direction(normal) ≈ expected
        angles = @inferred CanopyOptics.getSpecularΩ(
            CanopyOptics.dirVector(acos(mu_in),zero(FT)),
            CanopyOptics.dirVector(acos(mu_out),phi))
        @test angles isa CanopyOptics.dirVector{FT}
    end
    for FT in (Float32,Float64)
        model=CanopyOptics.SpecularCanopyScattering(
            nᵣ=FT(scattering["refractive_index"]),κ=FT(scattering["roughness"]))
        mu=FT.(scattering["stream_cosines"])
        for m in scattering["fourier_orders"]
            pp,mp=CanopyOptics.compute_Z_matrices(model,mu,CanopyOptics.spherical_leaves(FT),m;
                quadrature=CanopyOptics.CanopyQuadrature(
                    n_leaf=scattering["leaf_nodes"],n_azimuth=scattering["azimuth_nodes"]))
            @test eltype(pp)===FT
            @test eltype(mp)===FT
            @test all(isfinite,pp)
            @test all(isfinite,mp)
        end
    end
end
