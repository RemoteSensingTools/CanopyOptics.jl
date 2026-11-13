# # 4SAIL Polar BRDF and the Angular NDVI Effect

# Where [the principal-plane tutorial](foursail_hotspot.md) sweeps a
# single line through view-angle space, this one runs 4SAIL over the
# full `(VZA, VAZ)` half-sphere in a single batched call and uses the
# result to look at the *angular* dependence of NDVI.
#
# Two ideas are worth highlighting:
#
# * [`FourSAILGeometrySet`](@ref) lets you evaluate 4SAIL for an entire
#   directional grid in one go — no Julia-level loop over view angles.
# * NDVI is not a single number for a canopy: even a vegetation-only
#   scene with a fixed leaf spectrum shows several percent of angular
#   spread, with a dip on the hotspot side.
using CanopyOptics
include(joinpath(@__DIR__, "..", "..", "..", "d3_plots.jl")) # hide

# ## Scene
const sza_deg     = 30.0
const lai         = 4.0
const lambda_nm   = [685.0, 800.0]
const soil_albedo = 0.0
const hotspot_h   = 0.05

# Dense (VZA, VAZ) grid. `FourSAILGeometrySet` consumes flattened
# (VZA, VAZ) vectors; we reshape the results back to a 2D grid after.
vzas_deg = collect(0.0:5.0:80.0)
vazs_deg = collect(0.0:5.0:355.0)
n_vza, n_vaz = length(vzas_deg), length(vazs_deg)

vza_all = repeat(vzas_deg; inner = n_vaz)
raa_all = repeat(vazs_deg; outer = n_vza)   # 4SAIL's raa = azimuth difference

# ## Leaf optics from PROSPECT-PRO
leaf = CanopyOptics.LeafProspectProProperties{Float64}(
    N = 1.5, Ccab = 40.0, Ccar = 8.0,
    Canth = 0.0, Cbrown = 0.0,
    Cw = 0.012, Cm = 0.009,
    Cprot = 0.0, Ccbc = 0.0,
)

opti = CanopyOptics.LeafOpticalProperties(400.0:1.0:2500.0)
T_grid, R_grid = CanopyOptics.prospect(leaf, opti)
λ_grid   = [Float64(v.val) for v in opti.λ]
grid_idx = [argmin(abs.(λ_grid .- λ)) for λ in lambda_nm]
R_leaf   = R_grid[grid_idx]
T_leaf   = T_grid[grid_idx]

# ## One batched 4SAIL call for the entire half-sphere
LD    = CanopyOptics.planophile_leaves2(Float64)
geoms = FourSAILGeometrySet(LD;
    sza_deg = sza_deg,
    vza_deg = vza_all,
    raa_deg = raa_all,
    quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64, n_azimuth = 16),
)

result = foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                  hotspot = hotspot_h, mode = :full)

# `result.rsot` is (n_wavelength × n_directions). Reshape each
# wavelength row back to (n_vza × n_vaz):
brf = Array{Float64}(undef, length(lambda_nm), n_vza, n_vaz)
for iλ in eachindex(lambda_nm)
    brf[iλ, :, :] .= permutedims(reshape(result.rsot[iλ, :], n_vaz, n_vza))
end

R_red = @view brf[1, :, :]
R_nir = @view brf[2, :, :]
NDVI  = (R_nir .- R_red) ./ (R_nir .+ R_red)

println("685 nm BRF extrema: ", extrema(R_red))
println("800 nm BRF extrema: ", extrema(R_nir))
println("NDVI extrema:        ", extrema(NDVI))

# Anchor the angular NDVI effect to nadir so we can see *how much* it
# varies across the view hemisphere relative to the most common
# space-based reference geometry.
nadir_ndvi = NDVI[1, 1]
ΔNDVI = NDVI .- nadir_ndvi
println("NDVI at nadir: ", nadir_ndvi,
        "   ΔNDVI extrema (relative to nadir): ", extrema(ΔNDVI))

# ## Polar BRDF plots
#
# Interactive polar hemisphere plots over `(VZA, VAZ)` with `0°` azimuth at top
# (same azimuth as the sun in 4SAIL's convention). The red star marks the hotspot
# direction at `(vza = sza, vaz = 0°)`.
Δlim = maximum(abs, ΔNDVI)
hs_marker = [(; vza = sza_deg, vaz = 0.0, color = "#ef4444", label = "hotspot")]

D3Figure(
    title = "4SAIL polar BRDF and NDVI angular effect",
    subtitle = "SZA = $(Int(sza_deg))°, LAI = $(lai), planophile LAD, soil albedo = $(soil_albedo), h = $(hotspot_h)",
    cols = 2,
    rowHeight = 310,
    panels = [
        d3_polar(vzas_deg, vazs_deg, brf[1, :, :];
                 title = "685 nm BRF (rsot)", cbarLabel = "BRF",
                 colormap = "viridis", markers = hs_marker),
        d3_polar(vzas_deg, vazs_deg, brf[2, :, :];
                 title = "800 nm BRF (rsot)", cbarLabel = "BRF",
                 colormap = "viridis", markers = hs_marker),
        d3_polar(vzas_deg, vazs_deg, NDVI;
                 title = "NDVI", cbarLabel = "NDVI",
                 colormap = "viridis", markers = hs_marker),
        d3_polar(vzas_deg, vazs_deg, ΔNDVI;
                 title = "NDVI − NDVI(nadir)", cbarLabel = "ΔNDVI",
                 colormap = "balance", clim = [-Δlim, Δlim], markers = hs_marker),
    ],
)

# ## Reading the figure
#
# * The 685 nm BRDF is dim and dominated by the hotspot peak: most of
#   the red signal *is* the directly back-scattered branch, because
#   chlorophyll absorbs almost everything that scatters into the diffuse
#   field.
# * The 800 nm BRDF is much brighter and noticeably flatter, with a
#   broad bowl shape — multiple scattering smears the angular contrast
#   out across the upper hemisphere.
# * The two effects combine in the NDVI panel: NDVI dips on the hotspot
#   side (where the red BRF is locally lifted) and rises on the
#   forward-scatter side. Even with a fixed leaf spectrum and a
#   horizontally homogeneous canopy this is a several-percent angular
#   variation — large enough to matter for cross-instrument comparison
#   if view geometries differ.
