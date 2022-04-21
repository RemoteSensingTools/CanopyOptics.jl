# calculate scattering coefficient for forest based on distorted born approximation

using YAML
using SpecialFunctions
using QuadGK
using Parameters
using Test
using BenchmarkTools
# using LinearAlgebra Import later if needed

include("types.jl")
include("parameters_from_yaml.jl")
include("probabilities.jl")
include("subroutines.jl")
include("output_check.jl")

## 
## Load input parameters
## 

const INPUT_FILE = "deccheckin.yaml"
const input_params = parameters_from_yaml(INPUT_FILE)

@unpack bfrghz, leaf, branch_1, branch_2, trunk, 
d_c, d_t, ϵ_g, l, sig = input_params

# Conversion to standard metric

const c           = 3e8             # Speed of light m/s
const bfr         = bfrghz*1e9      # GHz to Hz
const lm          = l*1e-2          # cm to m
const s           = sig*1e-2        # cm to m (surface rms height)
const k₀          = 2π*bfr/c        # Free space wave number (2π/λ)
const n_ϕ         = 41              # No. of ϕ
const n_θ         = 37              # No. of θ
const ej          = 0.0 + 1.0im     # Complex unit vector

# Loop over angle of incidence θ_i (not currently looped)
const ip = 1

const θ_iᵈ  = 40
const θ_iʳ  = deg2rad(θ_iᵈ)
const ϕ_iᵈ  = 0
const ϕ_iʳ  = deg2rad(ϕ_iᵈ) 
const δθ   = π/n_θ
const δϕ   = 2π/n_ϕ
const Δϕ   = 1/n_ϕ 
const Δθ   = 1/n_θ

# Roughness factor for ground scattering in direct-reflected term
# r_g defined right above 3.1.9
const r_g = exp(-4*(k₀*s*cos(θ_iʳ))^2)

# dim1 = polarization [v, h]
# dim2 = layer [crown, trunk]
at = zeros(2, 2, 20)

# Compute scattering amplitudes from leaves (forward and backward)
afhhl, afvvl = afsal(θ_iʳ, leaf)
σ_l = asal(θ_iʳ, leaf)

# Forward scattering of all wood types
# Each output is a 2x2 matrix (vv vh ; hv hh)
afb1, afb2, aft = wood_forward.([branch_1, branch_2, trunk])

# Backward scattering of all wood types
# Each output is a WoodBackscatter type, with 3 fields: d, dr, vh
# Each of these fields is an array 
# (d: [vv, vh, hh], dr: [vv, hh], vh: [vh1, vh3])
σ_b1, σ_b2, σ_t = wood_backward.([branch_1, branch_2, trunk])

############################

# CALCULATION OF PROPAGATION CONSTANT IN LAYER 1(TOP)

# 3.1.8??? 
K_vc = k₀*cos(θ_iʳ)+(2π*(leaf.ρ*afvvl + branch_1.ρ*abs_components(afb1[1,1]) + branch_2.ρ*abs_components(afb2[1,1])))/(k₀*cos(θ_iʳ))
K_hc = k₀*cos(θ_iʳ)+(2π*(leaf.ρ*afhhl + branch_1.ρ*abs_components(afb1[2,2]) + branch_2.ρ*abs_components(afb2[2,2])))/(k₀*cos(θ_iʳ))

ath1, atv1 = abs.(imag.([K_hc, K_vc]))
K_hc=complex(real(K_hc),abs(imag(K_hc)))
K_vc=complex(real(K_vc),abs(imag(K_vc)))

atv1 = (atv1 <= 1.0E-20 ? 0.0001 : atv1)

skdh1= 1/ath1
skdv1= 1/atv1

at[1, 1, ip]=atv1
at[2, 1, ip]=ath1

############################
# CALCULATION OF PROPAGATION CONSTANT IN LAYER 2 (BOTTOM) 

K_ht = k₀*cos(θ_iʳ)+(2π*trunk.ρ*aft[2,2])/(k₀*cos(θ_iʳ))
K_vt = k₀*cos(θ_iʳ)+(2π*trunk.ρ*aft[1,1])/(k₀*cos(θ_iʳ))

ath2, atv2 = abs.(imag.([K_ht, K_vt]))
K_ht = complex(real(K_ht),abs(imag(K_ht)))
K_vt = complex(real(K_vt),abs(imag(K_vt)))

at[1, 2, ip]=abs(imag(K_vc+K_vt))
at[2, 2, ip]=abs(imag(K_hc+K_ht))

############################
# Calculation of reflection coefficient from ground 

rgh = (cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))
rgv = (ϵ_g*cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(ϵ_g*cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))

K_hcⁱ, K_vcⁱ, K_htⁱ, K_vtⁱ = imag.((K_hc, K_vc, K_ht, K_vt))

reflhh = rgh*exp(ej*(K_hc+K_hc)*d_c+ej*(K_ht+K_ht)*d_t)
reflhc = conj(reflhh)
reflvv = rgv*exp(ej*(K_vc+K_vc)*d_c+ej*(K_vt+K_vt)*d_t)
reflha = abs(reflhh)
reflva = abs(reflvv)

dattenh1 = exp(-2*(K_hcⁱ+K_hcⁱ)*d_c)
dattenh2 = exp(-2*(K_htⁱ+K_htⁱ)*d_t)
dattenv1 = exp(-2*(K_vcⁱ+K_vcⁱ)*d_c)
dattenv2 = exp(-2*(K_vtⁱ+K_vtⁱ)*d_t)
dattenvh1 = exp(-2*(K_vcⁱ+K_hcⁱ)*d_c)
dattenvh2 = exp(-2*(K_vtⁱ+K_htⁱ)*d_t)
a1     = (ej)*(K_vc-K_hc)
e1     = (ej)*(conj(K_vc)-conj(K_hc))
a2     = (ej)*(K_vt-K_ht)
e2     = (ej)*(conj(K_vt)-conj(K_ht))

############################

# Row is polarization [vv, vh, hh]
# Column is layer [branch+leaf (d1), trunk (d2), leaf (d3)]
σ_d = zeros(3,3)

# Perform 3.1.2 over all polarizations, using computed subexpressions (term_c, term_t)
term_c = map((x,y) -> funcm(2*x, 2*y, d_c), [K_vcⁱ, K_vcⁱ, K_hcⁱ], [K_vcⁱ, K_hcⁱ, K_hcⁱ])
term_t = map((x,y) -> funcm(2*x, 2*y, d_t), [K_vtⁱ, K_vtⁱ, K_htⁱ], [K_vtⁱ, K_htⁱ, K_htⁱ])

σ_d[:,1] = (branch_1.ρ*σ_b1.d + branch_2.ρ*σ_b2.d + leaf.ρ*σ_l.d) .* term_c
σ_d[:,2] = trunk.ρ*σ_t.d .* term_t .* [dattenv1, dattenvh1, dattenh1]
σ_d[:,3] = leaf.ρ*σ_l.d .* term_c

############################

# Row is polarization [vv, hh]
# Column is layer [branch+leaf (d1), trunk (d2), leaf (d3)]
σ_dr = zeros(2,3)
σ_dr[:,1] = 4*d_c*(branch_1.ρ*σ_b1.dr + branch_2.ρ*σ_b2.dr + leaf.ρ*σ_l.dr)*r_g
σ_dr[:,2] = 4*d_t*trunk.ρ*σ_t.dr*r_g
σ_dr[:,3] = 4*d_c*leaf.ρ*σ_l.dr*r_g

σ_dr[1,:] *= reflva^2
σ_dr[2,:] *= reflha^2

############################

# TODO: Perform this vh polarization cleanup (cmd-F sgvh)

# Row is 
# Columns is 
σ_vh = zeros(3,3)

############################
# Backscat. cross section, hh pol.
        
sghhd = σ_d[3,1]+σ_d[3,2]    
sghhr = 0
sghhdr = σ_dr[2,1]+σ_dr[2,2]
sghdri = (σ_dr[2,1]+σ_dr[2,2])/2
sghh  = sghhd+sghhr+sghhdr
sghhi=sghhd+sghhr+sghdri

############################
# Backscat. cross section, vv pol.

sgvvd = σ_d[1,1]+σ_d[1,2]
sgvvr = 0
sgvvdr = σ_dr[1,1]+σ_dr[1,2]
sgvdri = (σ_dr[1,1]+σ_dr[1,2])/2
sgvv  = sgvvd+sgvvr+sgvvdr
sgvvi=sgvvd+sgvvr+sgvdri

############################
# Backscat. cross section, vh pol.

sgvhd = σ_d[2,1]+σ_d[2,2]

factvh1 = exp(-2*(K_hcⁱ-K_vcⁱ)*d_c)
factvh2 = exp(-2*(K_vcⁱ-K_hcⁱ)*d_c)
factvh3 = exp((-a1-e1)*d_c)

# vh1 
σ_vh[1,1] = (branch_1.ρ*σ_b1.vh[1]+branch_2.ρ*σ_b2.vh[1]+leaf.ρ*σ_l.vh[1])*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
σ_vh[2,1] = (branch_1.ρ*σ_b1.vh[1]+branch_2.ρ*σ_b2.vh[1]+leaf.ρ*σ_l.vh[1])*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
σ_vh[3,1] = abs(2*real((branch_1.ρ*σ_b1.vh[2] + branch_2.ρ*σ_b2.vh[2] +leaf.ρ*σ_l.vh[2])*reflvv*reflhc*cfun(a1,e1,d_c)*r_g))
# vh2
σ_vh[1,3] = (leaf.ρ*σ_l.vh[1])*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
σ_vh[2,3] = (leaf.ρ*σ_l.vh[1])*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
σ_vh[3,3] = abs(2*real((leaf.ρ*σ_l.vh[2])*reflvv*reflhc*cfun(a1,e1,d_c)*r_g))
# vh3
σ_vh[1,2] = factvh1*(trunk.ρ*σ_t.vh[1])*(reflha)^2*funcm(2*K_vtⁱ,-2*K_htⁱ,d_t)*r_g
σ_vh[2,2] = factvh2*(trunk.ρ*σ_t.vh[1])*(reflva)^2*funcp(2*K_vtⁱ,-2*K_htⁱ,d_t)*r_g
σ_vh[3,2] = abs(2*real(factvh3*(trunk.ρ*σ_t.vh[2])*reflvv*reflhc*cfun(a2,e2,d_t)*r_g))

sgvhdr1=σ_vh[1,1]+σ_vh[2,1]+σ_vh[3,1]
sgvh1  =   σ_d[2,1]+σ_vh[1,1]+σ_vh[2,1]+σ_vh[3,1]
sgvhi1 =   σ_d[2,1]+σ_vh[1,1]+σ_vh[2,1]
sgvhidr1=σ_vh[1,1]+σ_vh[2,1]

sgvhdr3=σ_vh[1,3]+σ_vh[2,3]+σ_vh[3,3]
sgvh3  =   σ_d[2,3]+σ_vh[1,3]+σ_vh[2,3]+σ_vh[3,3]

sgvhi3=σ_d[2,3]+σ_vh[1,3]+σ_vh[2,3]
sgvhidr3=σ_vh[1,3]+σ_vh[2,3]

sgvhdr2=σ_vh[1,2]+σ_vh[2,2]+σ_vh[3,2]

sgvh2  =   σ_d[2,2]+σ_vh[1,2]+σ_vh[2,2]+σ_vh[3,2]
                
sgvhi2=σ_d[2,2]+σ_vh[1,2]+σ_vh[2,2]
sgvhidr2=σ_vh[1,2]+σ_vh[2,2]

sgvh = sgvh1+sgvh2
sgvhi = sgvhi1+sgvhi2
sgvhidr=sgvhidr1+sgvhidr2

### 

σ = [sgvv sgvh sghh]
σ_i = [sgvvi sgvhi sghhi]

############################
# Calculation of backscat cross sections in db

svhi = zeros(20)
svhi1 = zeros(20)
svhi2 = zeros(20)
svhi3 = zeros(20)

sgvhr = 0

σ_d_db = 10 * log10.(σ_d)
σ_dr_db = 10 * log10.(σ_dr)

svhi[ip]  = 10.0*log10(sgvhidr)
svhi1[ip] = 10.0*log10(sgvhidr1)
svhi2[ip] = 10.0*log10(sgvhidr2)
svhi3[ip] = 10.0*log10(sgvhidr3)

sghho  = 10.0*log10(sghh)
sghhoi = 10.0*log10(sghhi)
sgvvo  = 10.0*log10(sgvv)
sgvvoi = 10.0*log10(sgvvi)
sgvho  = 10.0*log10(abs(sgvh))
sgvhoi = 10.0*log10(abs(sgvhi))
shhdd  = 10.0*log10(sghhd)

shhdrd	= 10.0*log10(sghhdr)

shhdri	= 10.0*log10(sghdri)
shhrd	= 10.0*log10(sghhr)

svhdd	= 10.0*log10(abs(sgvhd))
sgvhdr	= sgvhdr1+sgvhdr2+sgvhdr3

svhdrd	= 10.0*log10(abs(sgvhdr))
svhdrd1	= 10.0*log10(abs(sgvhdr1))
svhdrd2	= 10.0*log10(abs(sgvhdr2))
svhdrd3	= 10.0*log10(abs(sgvhdr3))            
sgvho	= 10.0*log10(sgvhd+sgvhdr+sgvhr)

svhrd	= 10.0*log10(abs(sgvhr))
svvdd	= 10.0*log10(sgvvd)

svvdrd	= 10.0*log10(sgvvdr)

svvdri	= 10.0*log10(sgvdri)
svvrd	= 10.0*log10(sgvvr)

################################
# Add the effect of rough ground

ksig=k₀*s

σ_g, gd = grdoh(ksig)

σ_t = σ + σ_g
σ_t_i = σ_i + σ_g

###############################

output = Forest_Scattering_Output(θ_iᵈ, σ_d_db[3,1], σ_d_db[3,2], σ_d_db[3,3], 
                                  σ_dr_db[2,1], σ_dr_db[2,2], σ_dr_db[2,3],
                                  10.0*log10(σ_t[3]), shhdd, shhdrd,
                                  
                                  σ_d_db[1,1], σ_d_db[1,2], σ_d_db[1,3],
                                  σ_dr_db[1,1], σ_dr_db[1,2], σ_dr_db[1,3],
                                  10.0*log10(σ_t[1]), svvdd, svvdrd,
                                  
                                  σ_d_db[2,1], σ_d_db[2,2], σ_d_db[2,3],
                                  svhdrd1, svhdrd2, svhdrd3,
                                  10.0*log10(σ_t[2]), svhdd, svhdrd,
                                  
                                  10.0*log10(σ_g[3]), 10.0*log10(σ_g[2]), 10.0*log10(σ_g[1]),
                                  
                                  sghhoi, sgvhoi, sgvvoi,
                                  
                                  10.0*log10(sgvhidr1), 10.0*log10(sgvhidr2), 10.0*log10(sgvhidr3),
                                  
                                  σ_t[3]/σ_t_i[3], σ_t[2]/σ_t_i[2], σ_t[1]/σ_t_i[1], 
                                  ath1, atv1, imag(K_hc+K_ht), imag(K_vc+K_vt))

check_output_matches_fortran(output)