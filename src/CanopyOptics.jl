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
using DelimitedFiles: readdlm
using Distributions: Beta, Uniform, UnivariateDistribution, params, pdf, support
using FastGaussQuadrature: gausslegendre
using ForwardDiff: ForwardDiff
using KernelAbstractions: KernelAbstractions, @index, @kernel
using LazyArtifacts: LazyArtifacts, @artifact_str
using LinearAlgebra: ⋅
using SpecialFunctions: erfcx, expint
using Statistics: mean
using TOML: TOML
using Unitful: @u_str, NoUnits, uconvert, unit, ustrip
using UnitfulEquivalences: Spectral

import LinearAlgebra: dot

@static if VERSION >= v"1.11"
    eval(Expr(:public,
        :gauleg,
        :A,
        :Asm,
        :H,
        :compute_Ψ,
        :compute_lambertian_Γ,
        :getSpecularΩ,
        :dot_product,
        :Fᵣ,
        :fresnel_components,
        :K,
        :calctav,
        :prospect_primal,
        :compute_Γ,
        :compute_Γ_isotropic,
        :compute_Z_matrices_aniso,
    ))
end

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

include("leaf_optics/prospect_layers.jl")
include("leaf_optics/prospect.jl")
include("leaf_optics/pigment_absorption.jl")
include("leaf_optics/excitation_selection.jl")
include("leaf_optics/surface_optics.jl")
include("leaf_optics/pigment_bases.jl")
include("leaf_optics/dorsiventral.jl")
include("leaf_optics/leaf_bsdf.jl")
include("leaf_optics/plate_emission.jl")
include("leaf_optics/plate_preparation.jl")
include("leaf_optics/dorsiventral_emission.jl")
include("leaf_optics/emission_spectrum.jl")
include("wood_optics/reflectance.jl")
include("canopy_structure/clumping.jl")
include("canopy_structure/hotspot.jl")

include("canopy_scattering/projection_geometry.jl")
include("canopy_structure/components.jl")
include("canopy_scattering/quadrature_controls.jl")
include("canopy_scattering/stokes.jl")
include("canopy_scattering/specular.jl")
include("canopy_scattering/bilambertian_fourier.jl")
include("canopy_scattering/angular_normalization.jl")
include("canopy_scattering/z_matrices.jl")
include("canopy_scattering/dorsiventral_bilambertian.jl")
include("foursail/foursail.jl")

include("initialization/loadProspect.jl")
include("initialization/default_constructors.jl")

include("utils/dielectric.jl")


export prospect, prospect_budget, prospect_constituent_coefficients,
       constituent_optical_depths, LeafMultiPigmentProperties,
       MultiPigmentOpticalProperties, MultiChlorophyllAntenna,
       createMultiPigmentOpticalStruct, load_wax_refractive_index,
       LeafSIPProperties, leaf_sip_budget,
       AbstractLeafSurface, ConeSurface, RoughDielectricSurface,
       calctav_smooth, beckmann_D, smith_beckmann_G1,
       microfacet_directional_reflectance, microfacet_hemispherical_reflectance,
       surface_interface_optics,
       DorsiventralLeafProperties, AsymmetricSurfaceLeafProperties,
       dorsiventral_from_traits, single_slab_dorsiventral_budget, dorsiventral_budget,
       CoupledLeafBSDF, PreparedLeafBSDFSpectrum, prepare_leaf_bsdf,
       surface_brdf_mueller, leaf_directional_budget, leaf_bsdf_spectrum,
       DorsiventralBiLambertianCanopyScattering, compute_Γ_directional,
       compute_Z_matrices_directional, compute_Z_matrices_coupled_bsdf
export TabulatedPhotonSpectrum, AbstractEmissionSupport, StrictStokes,
       emission_allowed, emission_probabilities
export AbstractLeafBirthProfile, UniformBirthProfile, ExponentialBirthProfile,
       RegularPlateEscape, DorsiventralPlateEscape, DorsiventralBirthProfile,
       leaf_emission_fates, diffuse_plate_optics, dorsiventral_plate_optics
export AbstractAbsorptionPartition, HomogeneousAbsorption, AbsorptionPartitionResult,
       AbsorptionPartitionStatus, AbsorptionPartitionValid, AbsorptionNoMaterial,
       AbsorptionInvalidInput, AbsorptionIncompatibleInput, partition_absorption,
       absorption_valid, AbstractAntennaSelection, ChlorophyllOnly, CarotenoidAntenna,
       antenna_excitation, select_absorbed_photons, photosystem_excitation
export createLeafOpticalStruct, LeafProspectProProperties, LeafOpticalProperties,
       DirectionCosine, direction_from_angles, dirVector, dirVector_μ
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
export LeafComponent, StemComponent, BranchComponent, TreeCanopy
export FourSAILGeometry, FourSAILGeometrySet, FourSAILResult, foursail, foursail!
export AbstractMaterial, AbstractWater, AbstractSoil, AbstractVegetation
export PureIce, LiquidPureWater, LiquidSaltWater, SoilMW, LeafUlabyElRayes1987
export dielectric
# Functions:
export compute_Z_matrices, compute_Z_matrices_aniso_analytic,
       compute_reflection, compute_reflection_mueller,
       G, G2, bfG
export AbstractAngularNormalization, ContinuousAngularNormalization,
       DiscreteAngularNormalization, leaf_angular_scale, leaf_projection_moments

end # module
