using CanopyOptics, Test
const CUDA = CanopyOptics.CUDA
include("emission_spectrum_batches.jl")
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
CUDA.versioninfo()
@testset "Conditional spectrum CUDA application" begin
    check_emission_spectrum_batches(CUDA.CuArray, CUDA.synchronize)
end
