module CanopyOptics

###### Julia packages to import  ############
using Distributions            # Distributions for angular distributions of canopy elements
using FastGaussQuadrature      # Quadrature points for numerical integration
using Unitful                  # Units attached to some variable
using UnitfulEquivalences      # Spectral conversions
using SpecialFunctions:expint  # expint in Prospect
using DelimitedFiles           # File IO (Prospect csv File)
using DocStringExtensions      # Documentation
using LazyArtifacts            # Artifacts
using LinearAlgebra            # Well, guess...
using Polynomials              # Polynomials for some empirical functions
using YAML                     # YAML input files 
using QuadGK                   # Numerical Integration
using CUDA 
using ForwardDiff

import SpecialFunctions.expint 
#"Definition of Stokes vector types:"
#using vSmartMOM.Scattering: Stokes_I, Stokes_IQU, Stokes_IQUV

# Filename for ProspectPro optical properties
const OPTI_2021 = artifact"Prospect" * "/dataSpec_PRO.csv";

# Needs to be more modular later:
FT = Float64 
###### Own files to include #################
include("initialization/constants.jl")

include("types/canopy_types.jl")
include("types/angle_types.jl")
include("types/material_types.jl")

include("utils/quadrature.jl")
include("utils/canopy_angles.jl")
include("utils/fresnel.jl")

include("core/prospect.jl")
include("core/leafAngleRoutines.jl")

include("initialization/loadProspect.jl")
include("initialization/default_constructors.jl")

include("forest_prototyping/types.jl")
include("forest_prototyping/parameters_from_yaml.jl")
include("forest_prototyping/probabilities.jl")
include("forest_prototyping/subroutines.jl")
#include("forest_prototyping/output_check.jl")

include("utils/dielectric.jl")


export prospect
export createLeafOpticalStruct, LeafProspectProProperties, LeafOpticalProperties, dirVector, dirVector_μ
export AbstractCanopyScatteringType, BiLambertianCanopyScattering, SpecularCanopyScattering
export PureIce, LiquidPureWater, LiquidSaltWater
export LeafProspectProProperties, LeafOpticalProperties, dielectric
# Functions:
export compute_Z_matrices, prospect, compute_specular_reflection
# MW stuff
export wood_forward, wood_backward, afsal, asal, abs_components 

end # module
