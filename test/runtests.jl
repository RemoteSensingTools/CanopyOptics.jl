using CanopyOptics
using ForwardDiff
using Test

include("test_CanopyOptics.jl")

include("test_leaf_quadrature.jl")
include("test_specular_precision.jl")
include("test_specular_grazing.jl")
include("test_specular_angle.jl")
include("test_prospect_budget.jl")
include("test_plate_emission.jl")
include("test_plate_preparation.jl")
include("test_emission_spectrum.jl")
include("emission_spectrum_batches.jl")
@testset "Conditional spectrum CPU application" begin
    check_emission_spectrum_batches(identity, () -> nothing)
end
include("test_excitation_allocation.jl")
include("test_emission_operators.jl")
include("emission_operator_batches.jl")
@testset "Excitation allocation and emission operators CPU application" begin
    check_emission_operator_batches(identity, () -> nothing)
end
include("test_piecewise_birth.jl")
@testset "Piecewise birth profile CPU batches" begin
    check_piecewise_birth_batches(identity, () -> nothing)
end
include("plate_emission_batches.jl")
@testset "Regular plate emission CPU batches" begin
    check_plate_emission_batches(identity,()->nothing)
end

include("test_pigment_absorption.jl")
include("pigment_batches.jl")
@testset "Pigment allocation CPU batches" begin
    check_pigment_batches(identity,()->nothing)
end
include("pigment_spectral_fixture.jl")
@testset "Actual PROSPECT pigment spectra CPU" begin
    check_pigment_spectra(identity,()->nothing)
end

include("test_angular_normalization.jl")
include("test_leaf_model_families.jl")
include("test_qa_and_indexing.jl")
