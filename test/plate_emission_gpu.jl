using CanopyOptics, CUDA, Test
include("plate_emission_reference.jl")
include("plate_emission_batches.jl")
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
CUDA.versioninfo()
@testset "Regular plate emission CUDA independent references" begin
    check_plate_emission_batches(CUDA.CuArray, CUDA.synchronize)
end
