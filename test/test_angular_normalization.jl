using CanopyOptics,Test,LinearAlgebra,TOML,ForwardDiff

@testset "Opt-in discrete angular normalization" begin
    settings = TOML.parsefile(joinpath(@__DIR__,"angular_normalization.toml"))
    for FT in (Float32,Float64)
        prefix = FT===Float32 ? "float32" : "float64"
        atol,rtol = (FT(settings["validation"][prefix*"_"*name]) for name in ("atol","rtol"))
        continuous = ContinuousAngularNormalization()
        for streams in settings["numerics"]["streams"]
            μ,w = CanopyOptics.gauleg(streams,zero(FT),one(FT))
            append!(μ,FT.(settings["numerics"]["view_cosines"]))
            append!(w,zeros(FT,length(settings["numerics"]["view_cosines"])))
            normalization = DiscreteAngularNormalization(μ,w;atol)
            # Closed-form horizontal/vertical leaf limits are independent of
            # the production one-sided projection routine.
            @test (@inferred leaf_angular_scale(continuous,one(FT)))===one(FT)
            @test isapprox(@inferred(leaf_angular_scale(normalization,one(FT))),one(FT);atol,rtol)
            vertical = inv(FT(4)/FT(pi)*sum(w.*sqrt.(one(FT).-μ.^2)))
            @test isapprox(leaf_angular_scale(normalization,zero(FT)),vertical;atol,rtol)
            for leaf_cosine in FT.(settings["numerics"]["leaf_cosines"])
                sum_front,sum_back = zero(FT),zero(FT)
                for i in eachindex(μ),sign in (-one(FT),one(FT))
                    P,N = leaf_projection_moments(sign*μ[i],leaf_cosine,settings["numerics"]["fourier_order"];normalization)
                    @test eltype(P)===eltype(N)===FT
                    @test all(isfinite,P)&&all(isfinite,N)
                    sum_front += FT(2)*w[i]*P[1]
                    sum_back += FT(2)*w[i]*N[1]
                end
                @test isapprox(sum_front,one(FT);atol,rtol)
                @test isapprox(sum_back,one(FT);atol,rtol)
            end
            for (a,b) in values(settings["distributions"]),n_leaf in settings["numerics"]["leaf_nodes"]
                lad = beta_leaves(FT(a),FT(b))
                quadrature = CanopyQuadrature(n_leaf=n_leaf)
                G = vec(CanopyOptics.G(μ,lad;quadrature))
                for (r,t) in values(settings["optics"])
                    model = BiLambertianCanopyScattering(R=FT(r),T=FT(t))
                    modes = 0:settings["numerics"]["fourier_order"]
                    default = compute_Z_matrices(model,μ,lad,modes;quadrature)
                    reference = compute_Z_matrices(model,μ,lad,modes;quadrature,angular_normalization=continuous)
                    @test default==reference
                    Zpp,Zmp = compute_Z_matrices(model,μ,lad,modes;quadrature,angular_normalization=normalization)
                    @test eltype(Zpp)===eltype(Zmp)===FT
                    if iszero(r+t)
                        @test all(iszero,Zpp)&&all(iszero,Zmp)
                    else
                        columns = vec(transpose(w)*(Zpp[:,:,1]+Zmp[:,:,1]))/FT(2)
                        @test isapprox(columns,ones(FT,length(μ));atol,rtol)
                        @test minimum(Zpp[:,:,1])>=zero(FT)&&minimum(Zmp[:,:,1])>=zero(FT)
                        for Z in (Zpp,Zmp),m in axes(Z,3)
                            K = Z[:,:,m].*transpose(G)
                            @test isapprox(K,transpose(K);atol,rtol)
                        end
                    end
                    single = compute_Z_matrices(model,μ,lad,last(modes);quadrature,angular_normalization=normalization)
                    legacy = CanopyOptics.compute_Z_matrices_aniso(model,μ,lad,last(modes);quadrature,angular_normalization=normalization)
                    @test single==(Zpp[:,:,end],Zmp[:,:,end])==legacy
                end
            end
            # Stokes expansion is still scalar diffuse I→I, not new polarization physics.
            lad = beta_leaves(FT.(settings["distributions"]["central"])...)
            R,T = FT.(settings["optics"]["red"])
            model = BiLambertianCanopyScattering(R,T)
            scalar = compute_Z_matrices(model,μ,lad,0;angular_normalization=normalization)
            legacy_buffers = CanopyOptics.compute_Z_matrices_aniso(model,μ,lad,
                nothing,nothing,0;angular_normalization=normalization)
            @test legacy_buffers==scalar
            for npol in (3,4)
                blocks = compute_Z_matrices(model,μ,lad,0;npol,angular_normalization=normalization)
                for (block,plain) in zip(blocks,scalar)
                    @test block[1:npol:end,1:npol:end]==plain
                    block[1:npol:end,1:npol:end] .= zero(FT)
                    @test all(iszero,block)
                end
            end
            @test_throws ArgumentError compute_Z_matrices(model,reverse(μ),lad,0;angular_normalization=normalization)
            @test_throws ArgumentError compute_Z_matrices(model,(FT===Float32 ? Float64 : Float32).(μ),lad,0;angular_normalization=normalization)
        end
        μ,w = FT[one(FT)],FT[one(FT)]
        owned = DiscreteAngularNormalization(view(μ,:),view(w,:);atol)
        μ[1]=zero(FT);w[1]=zero(FT)
        @test owned.μ==owned.weights==FT[one(FT)]
        @test_throws ArgumentError DiscreteAngularNormalization(FT[zero(FT)],FT[one(FT)];atol)
        @test_throws ArgumentError DiscreteAngularNormalization(FT[one(FT)],FT[-one(FT)];atol)
        @test_throws ArgumentError DiscreteAngularNormalization(FT[one(FT)],FT[zero(FT)];atol)
        @test_throws DimensionMismatch DiscreteAngularNormalization(FT[one(FT)],FT[];atol)
        @test_throws ArgumentError leaf_projection_moments(FT(Inf),one(FT),0)
        @test_throws ArgumentError leaf_projection_moments(one(FT),-one(FT),0)
        @test_throws ArgumentError leaf_projection_moments(one(FT),one(FT),-1)
    end
    # Existing optical-parameter AD remains available with Float64 angular geometry.
    μ,w = CanopyOptics.gauleg(first(settings["numerics"]["streams"]),0.0,1.0)
    norm = DiscreteAngularNormalization(μ,w;atol=settings["validation"]["float64_atol"])
    lad = beta_leaves(settings["distributions"]["central"]...)
    R,T = settings["optics"]["red"]
    f(r) = sum(first(compute_Z_matrices(BiLambertianCanopyScattering(r,T),μ,lad,0;angular_normalization=norm)))
    step = settings["validation"]["finite_difference_step"]
    @test isapprox(ForwardDiff.derivative(f,R),(f(R+step)-f(R-step))/(2step);rtol=settings["validation"]["derivative_rtol"])
end
