using Test, CanopyOptics, TOML

@testset "Specular reflection is independent of direction representation" begin
    settings = TOML.parsefile(joinpath(@__DIR__,"specular_angle.toml"))
    directions,scattering,validation = settings["directions"],settings["scattering"],settings["validation"]
    for FT in (Float32,Float64)
        model=SpecularCanopyScattering(nᵣ=FT(scattering["refractive_index"]),κ=FT(scattering["roughness"]))
        distributions=(spherical_leaves(FT),beta_leaves(FT(scattering["beta_alpha"]),FT(scattering["beta_beta"])))
        rtol=FT(validation[FT===Float32 ? "float32_rtol" : "float64_rtol"])
        for lad in distributions,mu_in in FT.(directions["incoming_cosines"]),
                mu_out in FT.(directions["outgoing_cosines"]),fraction in FT.(directions["outgoing_azimuth_pi"])
            phi=fraction*FT(pi)
            incoming=CanopyOptics.dirVector_μ(mu_in,zero(FT))
            outgoing=CanopyOptics.dirVector_μ(mu_out,phi)
            expected=compute_reflection(model,incoming,outgoing,lad)
            actual=compute_reflection(model,
                CanopyOptics.dirVector(acos(mu_in),zero(FT)),
                CanopyOptics.dirVector(acos(mu_out),phi),lad)
            @test actual isa FT
            @test isfinite(actual)
            @test actual>=zero(FT)
            @test isapprox(actual,expected;rtol)
        end
    end
end
