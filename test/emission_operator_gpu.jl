using CanopyOptics, CUDA, Test
include("emission_operator_batches.jl")
include("piecewise_birth_batches.jl")
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
@testset "Excitation allocation and emission operators CUDA application" begin
    check_emission_operator_batches(CUDA.CuArray, CUDA.synchronize)
end
@testset "Piecewise birth profile CUDA application" begin
    check_piecewise_birth_batches(CUDA.CuArray, CUDA.synchronize)
end
