# # 4SAIL Hotspot in the Solar Principal Plane

# This tutorial walks through a red/NIR principal-plane BRF diagnostic
# with [`foursail`](@ref), using a realistic PROSPECT-PRO leaf and the
# [`FourSAILGeometrySet`](@ref) batched entry point. It shows three
# things in one figure:
#
# * the full BRF `rsot` with hotspot off (`h = 0`)
# * the full BRF `rsot` with a Kuusk-style hotspot turned on
# * the direct/direct branch `rsost` from the same hotspot solve
#
# The shape `rsot = rsost + rsodt` is what makes 4SAIL convenient as a
# diagnostic tool — the embedded `rsost` is the path that carries the
# hotspot, so you can read off the hotspot increment without rerunning
# the model in a reduced mode.
using CanopyOptics
using CairoMakie

# ## Geometry — signed VZA in the principal plane
#
# The "principal plane" is the vertical plane containing both the sun
# and the observer. We sweep view zenith from −80° to +80° and split
# the sign by relative azimuth:
#
# * `signed_vza < 0` → `raa = 0°`: backscatter side, same azimuth as
#   the sun, the side that carries the hotspot.
# * `signed_vza > 0` → `raa = 180°`: forward-scatter side.
#
# (4SAIL's `raa_deg` argument is the *azimuth difference* between sun
# and view directions, so RAA = 0° puts the observer in the same
# half-plane as the sun.)
const sza_deg = 30.0
const lai     = 4.0
const hotspot_h = 0.05
const lambda_nm = [685.0, 800.0]

signed_vza = collect(-80.0:2.0:80.0)
vza_abs    = abs.(signed_vza)
raa_deg    = ifelse.(signed_vza .< 0, 0.0, 180.0)

LD    = CanopyOptics.planophile_leaves2(Float64)
geoms = FourSAILGeometrySet(LD;
    sza_deg = sza_deg,
    vza_deg = vza_abs,
    raa_deg = raa_deg,
    quadrature = CanopyOptics.CanopyQuadrature(n_leaf = 64, n_azimuth = 16),
)

# ## Leaf optics from PROSPECT-PRO
#
# Evaluate PROSPECT once on a 1 nm grid in the 400–2500 nm range, then
# pull the two requested wavelengths out by nearest-neighbor index.
leaf = CanopyOptics.LeafProspectProProperties{Float64}(
    N = 1.5, Ccab = 40.0, Ccar = 8.0,
    Canth = 0.0, Cbrown = 0.0,
    Cw = 0.012, Cm = 0.009,
    Cprot = 0.0, Ccbc = 0.0,
)

opti = CanopyOptics.createLeafOpticalStruct(400.0:1.0:2500.0)
T_grid, R_grid = CanopyOptics.prospect(leaf, opti)
λ_grid   = [Float64(v.val) for v in opti.λ]
grid_idx = [argmin(abs.(λ_grid .- λ)) for λ in lambda_nm]
R_leaf   = R_grid[grid_idx]
T_leaf   = T_grid[grid_idx]

println("PROSPECT leaf optics at the two probe wavelengths:")
for (λ, ρ, τ) in zip(lambda_nm, R_leaf, T_leaf)
    println("  λ = ", λ, " nm   ρ = ", round(ρ, digits = 4),
            "   τ = ", round(τ, digits = 4))
end

# ## Two 4SAIL solves — hotspot off, hotspot on
#
# `foursail` returns a [`FourSAILResult`](@ref). For the batched
# `FourSAILGeometrySet` path each field is a (n_wavelength × n_angle)
# matrix.
soil_albedo = 0.0

sail_nohot = foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                      hotspot = 0.0, mode = :full)
sail_hot   = foursail(R_leaf, T_leaf, soil_albedo, geoms, lai;
                      hotspot = hotspot_h, mode = :full)

size(sail_hot.rsot)

# Quick numeric readout at three reference angles: the hotspot
# (`signed_vza = -sza`), nadir, and the symmetric forward-scatter angle.
hotspot_idx = findfirst(==(-sza_deg), signed_vza)
nadir_idx   = findfirst(==(0.0), signed_vza)
forward_idx = findfirst(==(sza_deg), signed_vza)

for (iλ, λ) in pairs(lambda_nm)
    println("λ = ", λ, " nm")
    for (label, idx) in (("hotspot (-sza)", hotspot_idx),
                         ("nadir",          nadir_idx),
                         ("forward (+sza)", forward_idx))
        idx === nothing && continue
        println("  ", rpad(label, 16),
                "  rsot(h=0) = ",      round(sail_nohot.rsot[iλ, idx], digits = 4),
                "  rsot(h=$(hotspot_h)) = ", round(sail_hot.rsot[iλ, idx], digits = 4),
                "  rsost = ",          round(sail_hot.rsost[iλ, idx], digits = 4))
    end
end

# ## Plot — full BRF and the hotspot increment
#
# Top row: BRF curves themselves (with and without hotspot, plus the
# direct/direct branch `rsost`). Bottom row: the hotspot increment
# `rsot(h) − rsot(0)` is concentrated near `vza = −sza`, exactly where
# the camera looks back along the sun direction.
fig = Figure(size = (1100, 700))
colors = (sail_nohot = :dodgerblue3,
          sail_hot   = :seagreen4,
          sail_ss    = :darkorange3,
          increment  = :purple4)

for (iλ, λ) in pairs(lambda_nm)
    ax = Axis(fig[iλ, 1];
        xlabel = iλ == length(lambda_nm) ?
                 "signed VZA in solar principal plane (deg)" : "",
        ylabel = "BRF",
        title  = "$(Int(round(λ))) nm — principal-plane BRF")
    lines!(ax, signed_vza, sail_nohot.rsot[iλ, :];
           color = colors.sail_nohot, linewidth = 2.2, linestyle = :dash,
           label = "rsot, h = 0")
    lines!(ax, signed_vza, sail_hot.rsot[iλ, :];
           color = colors.sail_hot,   linewidth = 2.5,
           label = "rsot, h = $(hotspot_h)")
    lines!(ax, signed_vza, sail_hot.rsost[iλ, :];
           color = colors.sail_ss,    linewidth = 2.0, linestyle = :dot,
           label = "rsost (direct/direct), h = $(hotspot_h)")
    vlines!(ax, [-sza_deg]; color = :gray30, linestyle = :dashdot, linewidth = 1.2,
            label = "hotspot direction")
    vlines!(ax, [0.0];      color = :gray60, linestyle = :dot,    linewidth = 1.0)
    axislegend(ax; position = :lt, framevisible = false)

    axd = Axis(fig[iλ, 2];
        xlabel = iλ == length(lambda_nm) ?
                 "signed VZA in solar principal plane (deg)" : "",
        ylabel = "Δ BRF",
        title  = "$(Int(round(λ))) nm — hotspot increment")
    lines!(axd, signed_vza, sail_hot.rsot[iλ, :] .- sail_nohot.rsot[iλ, :];
           color = colors.increment, linewidth = 2.4,
           label = "rsot(h) − rsot(0)")
    hlines!(axd, [0.0]; color = :gray70, linewidth = 1.0)
    vlines!(axd, [-sza_deg]; color = :gray30, linestyle = :dashdot, linewidth = 1.2)
    axislegend(axd; position = :rt, framevisible = false)
end

Label(fig[0, 1:2],
      "4SAIL principal-plane BRF and hotspot increment " *
      "(SZA = $(Int(sza_deg))°, LAI = $(lai), planophile LAD, soil albedo = $(soil_albedo))";
      fontsize = 16)
fig

# ## Reading the figure
#
# * Without the hotspot (`h = 0`) the BRF is smooth across the principal
#   plane. The slight back–forward asymmetry is just the LAD-weighted
#   single-scattering shape.
# * Turning the hotspot on lifts the curve sharply near `vza = -sza`.
#   The increment is essentially zero on the forward-scatter side, so
#   the multiple-scattering remainder `rsodt` is hardly perturbed — the
#   hotspot lives in the `rsost` branch.
# * `rsost` itself sits well below the full BRF in the NIR because most
#   of the 800 nm signal arrives through diffuse multiple scattering,
#   not single scattering. In the red band leaves absorb most of what
#   they scatter once, so `rsost` and the full `rsot` are much closer.
