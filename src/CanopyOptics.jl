module CanopyOptics

###### Julia packages to import  ############
using Distributions            # Distributions for angular distributions of canopy elements
using FastGaussQuadrature      # Quadrature points for numerical integration
using Unitful                  # Units attached to some variable
using UnitfulEquivalences      # Spectral conversions
using SpecialFunctions:expint  # expint in Prospect
using DelimitedFiles           # File IO (Prospect csv File)
using DocStringExtensions      # Documentation
using UnPack                   # @unpack and stuff
using LazyArtifacts            # Artifacts
using LinearAlgebra

# Filename for ProspectPro optical properties
const OPTI_2021 = artifact"Prospect" * "/dataSpec_PRO.csv";

###### Own files to include #################
include("types/canopy_types.jl")
include("types/angle_types.jl")
include("utils/quadrature.jl")
include("utils/canopy_angles.jl")
include("utils/fresnel.jl")


include("core/prospect.jl")
include("core/leafAngleRoutines.jl")

include("initialization/loadProspect.jl")
include("initialization/default_constructors.jl")

export prospect
export createLeafOpticalStruct, LeafProspectProProperties, LeafOpticalProperties, dirVector

end # module
