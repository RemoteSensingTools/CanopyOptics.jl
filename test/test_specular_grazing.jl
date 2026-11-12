using Test, TOML, CanopyOptics, ForwardDiff

@testset "Specular roughness grazing limit" begin
    settings = TOML.parsefile(joinpath(@__DIR__,"specular_grazing.toml"))
    r,s,v = settings["roughness"],settings["scattering"],settings["validation"]
    for FT in (Float32,Float64)
        halfpi = FT(pi)/FT(2)
        atol = FT(v[FT===Float32 ? "float32_atol" : "float64_atol"])
        rtol = FT(v[FT===Float32 ? "float32_rtol" : "float64_rtol"])
        for roughness in FT.(r["values"])
            expected = iszero(roughness) ? one(FT) : zero(FT)
            @test (@inferred CanopyOptics.K(roughness,halfpi)) === expected
            @test CanopyOptics.K(roughness,-halfpi) === expected
            near = CanopyOptics.K(roughness,prevfloat(halfpi))
            @test isfinite(near)
            @test zero(FT)<=near<=one(FT)
            for fraction in FT.(r["interior_angles_pi"])
                angle = fraction*FT(pi)
                # Independent high-precision trigonometric ratio away from grazing.
                reference = FT(exp(-BigFloat(roughness)*sin(BigFloat(angle))/cos(BigFloat(angle))))
                @test isapprox(CanopyOptics.K(roughness,angle),reference;rtol,atol)
            end
            # The fix must retain ordinary AD behavior away from the endpoint.
            angle = FT(r["derivative_angle_pi"])*FT(pi)
            value = CanopyOptics.K(roughness,angle)
            @test isapprox(ForwardDiff.derivative(k->CanopyOptics.K(k,angle),roughness),
                -tan(angle)*value;rtol,atol)
            @test isapprox(ForwardDiff.derivative(a->CanopyOptics.K(roughness,a),angle),
                -roughness*value/cos(angle)^2;rtol,atol)
        end
        model = CanopyOptics.SpecularCanopyScattering(nᵣ=FT(s["refractive_index"]),κ=FT(s["roughness"]))
        mu = FT.(s["stream_cosines"])
        lad = CanopyOptics.beta_leaves(FT(s["leaf_beta_alpha"]),FT(s["leaf_beta_beta"]))
        quadrature = CanopyOptics.CanopyQuadrature(n_leaf=s["leaf_nodes"],n_azimuth=s["azimuth_nodes"])
        for npol in s["stokes_components"],m in s["fourier_orders"]
            pp,mp = CanopyOptics.compute_Z_matrices(model,mu,lad,m;quadrature,npol)
            @test eltype(pp)===FT
            @test eltype(mp)===FT
            @test all(isfinite,pp)
            @test all(isfinite,mp)
            # For vertical streams (last_block, μ = 1):
            # - Same-sign transmission `pp` requires a vertical leaf (θ* = π/2), so K(κ, π/2) = 0 removes it.
            # - Opposite-sign reflection `mp` requires a horizontal leaf (θ* = 0), where Beta(2, 2) has pdf(0) = 0.
            last_block = (length(mu)-1)*npol+1:length(mu)*npol
            @test all(iszero,pp[last_block,last_block])
            @test all(iszero,mp[last_block,last_block])
        end
        lad_plano = CanopyOptics.planophile_leaves2(FT)
        for npol in s["stokes_components"]
            pp_p,mp_p = CanopyOptics.compute_Z_matrices(model,mu,lad_plano,0;quadrature,npol)
            last_block = (length(mu)-1)*npol+1:length(mu)*npol
            @test all(iszero,pp_p[last_block,last_block])
            @test mp_p[first(last_block),first(last_block)] > zero(FT)
        end
    end
end

if get(ENV,"CANOPY_SPECULAR_GPU","false")=="true"
    using CUDA
    CUDA.functional() || error("CUDA is required for the requested roughness test")
    CUDA.allowscalar(false)
    @testset "Specular roughness on CUDA" begin
        settings = TOML.parsefile(joinpath(@__DIR__,"specular_grazing.toml"))
        r,v = settings["roughness"],settings["validation"]
        for FT in (Float32,Float64)
            halfpi = FT(pi)/FT(2)
            angles = [FT.(r["interior_angles_pi"]).*FT(pi);prevfloat(halfpi);halfpi;-halfpi]
            roughness = reshape(FT.(r["values"]),:,1)
            angles = reshape(angles,1,:)
            reference = CanopyOptics.K.(roughness,angles)
            device_result = CanopyOptics.K.(CuArray(roughness),CuArray(angles))
            CUDA.synchronize()
            actual = Array(device_result)
            atol = FT(v[FT===Float32 ? "float32_atol" : "float64_atol"])
            rtol = FT(v[FT===Float32 ? "float32_rtol" : "float64_rtol"])
            @test eltype(actual)===FT
            @test all(isfinite,actual)
            @test all(x->zero(FT)<=x<=one(FT),actual)
            @test isapprox(actual,reference;rtol,atol)
        end
        println("Specular roughness verified on ",CUDA.name(CUDA.device()))
    end
end
