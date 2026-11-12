using Test, CanopyOptics, LinearAlgebra, TOML

@testset "Specular azimuth preserves floating-point type" begin
    settings = TOML.parsefile(joinpath(@__DIR__, "specular_precision.toml"))
    directions, scattering = settings["directions"], settings["scattering"]
    for FT in (Float32,Float64), mu_in in FT.(directions["incoming_cosines"]),
            mu_out in FT.(directions["outgoing_cosines"]),
            fraction in FT.(directions["outgoing_azimuth_pi"])
        phi = fraction*FT(pi)
        incoming = CanopyOptics.DirectionCosine(mu_in,zero(FT))
        outgoing = CanopyOptics.DirectionCosine(mu_out,phi)
        normal,alpha = @inferred CanopyOptics.getSpecularΩ(incoming,outgoing)
        @test normal isa CanopyOptics.DirectionCosine{FT}
        @test alpha isa FT
        # Independent Cartesian bisector, with incoming azimuth zero.
        direction(d) = FT[sqrt(one(FT)-d.μ^2)*cos(d.ϕ),
            sqrt(one(FT)-d.μ^2)*sin(d.ϕ),d.μ]
        expected = normalize(direction(incoming)+direction(outgoing))
        @test direction(normal) ≈ expected
        normal_ang, alpha_ang = @inferred CanopyOptics.getSpecularΩ(
            CanopyOptics.direction_from_angles(acos(mu_in),zero(FT)),
            CanopyOptics.direction_from_angles(acos(mu_out),phi))
        @test normal_ang isa CanopyOptics.DirectionCosine{FT}
        @test alpha_ang isa FT
        @test direction(normal_ang) ≈ expected
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
