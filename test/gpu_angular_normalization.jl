# Explicit CUDA acceptance check; the normalizer and projection geometry are host setup.
using CanopyOptics,CUDA,Test,TOML,LinearAlgebra
CUDA.functional(true) || error("CUDA is required for this acceptance check")
CUDA.allowscalar(false)
settings = TOML.parsefile(joinpath(@__DIR__,"angular_normalization.toml"))

@testset "Discrete leaf-source geometry on CPU/CUDA" begin
    for FT in (Float32,Float64)
        prefix = FT===Float32 ? "float32" : "float64"
        atol,rtol = (FT(settings["validation"][prefix*"_"*name]) for name in ("atol","rtol"))
        μ,w = CanopyOptics.gauleg(first(settings["numerics"]["streams"]),zero(FT),one(FT))
        append!(μ,FT.(settings["numerics"]["view_cosines"]))
        append!(w,zeros(FT,length(settings["numerics"]["view_cosines"])))
        normalization = DiscreteAngularNormalization(μ,w;atol)
        cosines = FT.(settings["numerics"]["leaf_cosines"])
        rays = vcat(μ,-μ)
        nd,nl,nm = length(rays),length(cosines),settings["numerics"]["fourier_order"]+1
        coefficients = zeros(FT,nd*nm,2nl)
        for i in eachindex(rays),l in eachindex(cosines)
            P,N = leaf_projection_moments(rays[i],cosines[l],nm-1;normalization)
            for m in 1:nm
                coefficients[i+nd*(m-1),l] = FT(2)*P[m]
                coefficients[i+nd*(m-1),nl+l] = FT(2)*N[m]
            end
        end
        front,back = FT.(settings["source"]["front"]),FT.(settings["source"]["back"])
        scales = FT.(settings["source"]["sample_scales"])
        sources = vcat(front,back)*transpose(scales)
        expected = coefficients*sources
        output = CUDA.zeros(FT,size(expected))
        CUDA.@sync mul!(output,CuArray(coefficients),CuArray(sources))
        @test eltype(output)===FT
        @test isapprox(Array(output),expected;atol,rtol)
        @test all(iszero,Array(output[:,iszero.(scales)]))
        @test any(<(zero(FT)),Array(output))
        flux = Array(sum(CuArray(vcat(w,w)).*view(output,1:nd,:);dims=1))
        @test isapprox(vec(flux),sum(front+back).*scales;atol,rtol)
    end
    println((device=CUDA.name(CUDA.device()),cuda_version=Base.pkgversion(CUDA),julia_version=VERSION))
end
