using CanopyOptics, CUDA, Test
include("pigment_batches.jl")
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
@testset "Pigment allocation CUDA batches" begin
    check_pigment_batches(CUDA.CuArray, CUDA.synchronize)
end
include("pigment_spectral_fixture.jl")
@testset "Actual PROSPECT pigment spectra CUDA" begin
    check_pigment_spectra(CUDA.CuArray, CUDA.synchronize)
end
