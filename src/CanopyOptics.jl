"""
    CanopyOptics

Canopy radiative-transfer utilities.

Canopy scattering matrices use `Z[i_out, j_in]`: rows are outgoing streams and
columns are incoming streams; `Z⁺⁺` is same-sign transmission and `Z⁻⁺` is
sign-change reflection. Bi-Lambertian kernels follow vSmartMOM's atmospheric
phase-matrix convention: the scalar single-scattering albedo is not folded into
`Z`, so conservative bi-Lambertian leaves satisfy

```math
\\sum_i w_i\\,(Z^{++}_{ij} + Z^{-+}_{ij}) \\simeq 2 .
```

Specular kernels carry their Fresnel/roughness strength in the returned
contribution. Layer solvers that multiply by a separate `ϖ` need a matching
effective albedo when using specular components.
"""
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

###### Own files to include #################
include("initialization/constants.jl")

include("types/canopy_types.jl")
include("types/angle_types.jl")
include("types/material_types.jl")

include("utils/quadrature.jl")
include("utils/canopy_angles.jl")
include("utils/fresnel.jl")

include("leaf_optics/prospect.jl")
include("wood_optics/reflectance.jl")
include("canopy_structure/clumping.jl")
include("canopy_structure/hotspot.jl")

include("canopy_scattering/projection_geometry.jl")
include("canopy_structure/components.jl")
include("canopy_scattering/quadrature_controls.jl")
include("canopy_scattering/stokes.jl")
include("canopy_scattering/specular.jl")
include("canopy_scattering/bilambertian_fourier.jl")
include("canopy_scattering/z_matrices.jl")

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
export AbstractCanopyScatteringType, CanopyQuadrature, BiLambertianCanopyScattering,
       CompositeCanopyScattering, SpecularCanopyScattering,
       LambertianWoodCanopyScattering, LambertianWood
export AbstractLeafDistribution, LeafDistribution,
       planophile_leaves, planophile_leaves2, uniform_leaves,
       plagiophile_leaves, erectophile_leaves, spherical_leaves,
       flat_leaves, beta_leaves, βparameters
export AbstractWoodReflectance, AbstractLUTWoodReflectance,
       ConstantWoodReflectance, LUTWoodReflectance, PolynomialWoodReflectance,
       wood_reflectance
export AbstractClumping, NoClumping, ConstantClumping,
       EmpiricalDirectionalClumping, ChenLeblancClumping,
       clumping_index, effective_G
export AbstractHotSpot, NoHotSpot, KuuskHotSpot, canopy_extinction,
       hotspot_separation, hotspot_correction, joint_gap_probability
export CanopyComponent, MixedCanopy, component_G, bulk_G
export PureIce, LiquidPureWater, LiquidSaltWater
export LeafProspectProProperties, LeafOpticalProperties, dielectric
# Functions:
export compute_Z_matrices, compute_Z_matrices_aniso_analytic,
       prospect, compute_reflection, compute_reflection_mueller,
       G, G2, bfG
# MW stuff
export wood_forward, wood_backward, afsal, asal, abs_components 

end # module
