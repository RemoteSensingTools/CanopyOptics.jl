# calculate scattering coefficient for forest based on distorted born approximation

using YAML
using SpecialFunctions
using QuadGK
using Parameters
using Test
using BenchmarkTools

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

const c           = 3e8             # Speed of light 
const bfr         = bfrghz*1e9      # GHz to Hz
const lm          = l*1e-2          # cm to m
const s           = sig*1e-2        # cm to m (surface rms height)
const k₀          = 2π*bfr/c        # Free space wave number (2π/λ)
const zk          = k₀              # 
const n_ϕ         = 41              # No. of ϕ
const n_θ         = 37              # No. of θ
const ej          = 0.0 + 1.0im     # Complex unit vector

## 
## Calculation of Parameters
## 

# Calculate integral of p(θ)*sin(θ)^2 from 0 to pi
# That is, the average integral over inclinations
const ail = sum(x -> prob(x, leaf.pdf_num, leaf.pdf_param) * sin(x)^2 * π/n_θ, collect(1:n_θ) * π/n_θ)

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

# Calculation of skin depth skdh skdv and the bistatic cross sections sghh,sghv,sgvv

afhhl, afvvl = afsal(θ_iʳ, ail, leaf)
σ_d_l, σ_dr_l, σ_vh1_l, σ_vh3_l = asal(θ_iʳ, leaf)

# Compute scattering amplitudes from primary branches

afb1 = woodf(branch_1)

afhhb1 = complex(abs(real(afb1[2,2])),abs(imag(afb1[2,2])))
afvvb1 = complex(abs(real(afb1[1,1])),abs(imag(afb1[1,1])))

σ_d_b1, σ_dr_b1, σ_vh1_b1, σ_vh3_b1 = woodb(branch_1)

# compute scattering amplitudes from secondary branches

afb2 = woodf(branch_2)

afhhb2 = complex(abs(real(afb2[2,2])),abs(imag(afb2[2,2])))
afvvb2 = complex(abs(real(afb2[1,1])),abs(imag(afb2[1,1])))

σ_d_b2, σ_dr_b2, σ_vh1_b2, σ_vh3_b2 = woodb(branch_2)

# Compute scattering amplitudes from trunks

aft = woodf(trunk)
σ_d_t, σ_dr_t, σ_vh1_t, σ_vh3_t = woodb(trunk)

# Using reciprocity and scatterer symmetry to calculate rho*sigma

# TODO: Remove these variables

# Branch / trunk / leaf (VH1)
bsgvh11 = branch_1.ρ*σ_vh1_b1 + branch_2.ρ*σ_vh1_b2 +leaf.ρ*σ_vh1_l
bsgvh12 = trunk.ρ*σ_vh1_t
bsgvh13 = leaf.ρ*σ_vh1_l
        
# Branch / trunk / leaf (VH3)
bsgvh31 = branch_1.ρ*σ_vh3_b1 + branch_2.ρ*σ_vh3_b2 +leaf.ρ*σ_vh3_l
bsgvh32 = trunk.ρ*σ_vh3_t
bsgvh33 = leaf.ρ*σ_vh3_l 

bsgvh21=bsgvh11
bsgvh22=bsgvh12
bsgvh23=bsgvh13

afhh1 = leaf.ρ*afhhl + branch_1.ρ*afhhb1 + branch_2.ρ*afhhb2
afhh2 = trunk.ρ*aft[2,2]

afvv1 = leaf.ρ*afvvl + branch_1.ρ*afvvb1 + branch_2.ρ*afvvb2
afvv2 = trunk.ρ*aft[1,1]

############################

# CALCULATION OF PROPAGATION CONSTANT IN LAYER 1(TOP)

# 3.1.8??? 
K_vc = k₀*cos(θ_iʳ)+(2π*afvv1)/(k₀*cos(θ_iʳ))
K_hc = k₀*cos(θ_iʳ)+(2π*afhh1)/(k₀*cos(θ_iʳ))

ath1= abs(imag(K_hc))
atv1= abs(imag(K_vc))
K_hc=complex(real(K_hc),abs(imag(K_hc)))
K_vc=complex(real(K_vc),abs(imag(K_vc)))

atv1 = (atv1 <= 1.0E-20 ? 0.0001 : atv1)

skdh1= 1/ath1
skdv1= 1/atv1

at[1, 1, ip]=atv1
at[2, 1, ip]=ath1

############################
# CALCULATION OF PROPAGATION CONSTANT IN LAYER 2 (BOTTOM) 

K_ht = k₀*cos(θ_iʳ)+(2π*afhh2)/(k₀*cos(θ_iʳ))
K_vt = k₀*cos(θ_iʳ)+(2π*afvv2)/(k₀*cos(θ_iʳ))
K_ht = complex(real(K_ht),abs(imag(K_ht)))
K_vt = complex(real(K_vt),abs(imag(K_vt)))

ath2= abs(imag(K_ht))
atv2= abs(imag(K_vt))
temp0h=abs(imag(K_hc+K_ht))
temp0v=abs(imag(K_vc+K_vt))
skdh2= 1/ath2
skdv2= 1/atv2

at[1, 2, ip]=temp0v
at[2, 2, ip]=temp0h


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

factvh1 = exp(-2*(K_hcⁱ-K_vcⁱ)*d_c)
factvh2 = exp(-2*(K_vcⁱ-K_hcⁱ)*d_c)
factvh3 = exp((-a1-e1)*d_c)

############################

# Row is polarization [vv, vh, hh]
# Column is layer [branch+leaf (d1), trunk (d2), leaf (d3)]
σ_d = zeros(3,3)

# Perform 3.1.2 over all polarizations, using computed subexpressions (term_c, term_t)
term_c = map((x,y) -> funcm(2*x, 2*y, d_c), [K_vcⁱ, K_vcⁱ, K_hcⁱ], [K_vcⁱ, K_hcⁱ, K_hcⁱ])
term_t = map((x,y) -> funcm(2*x, 2*y, d_t), [K_vtⁱ, K_vtⁱ, K_htⁱ], [K_vtⁱ, K_htⁱ, K_htⁱ])

σ_d[:,1] = (branch_1.ρ*σ_d_b1 + branch_2.ρ*σ_d_b2 + leaf.ρ*σ_d_l) .* term_c
σ_d[:,2] = trunk.ρ*σ_d_t .* term_t .* [dattenv1, dattenvh1, dattenh1]
σ_d[:,3] = leaf.ρ*σ_d_l .* term_c

############################

# Row is polarization [vv, hh]
# Column is layer [branch+leaf (d1), trunk (d2), leaf (d3)]
σ_dr = zeros(2,3)
σ_dr[:,1] = 4*d_c*(branch_1.ρ*σ_dr_b1 + branch_2.ρ*σ_dr_b2 + leaf.ρ*σ_dr_l)*r_g
σ_dr[:,2] = 4*d_t*trunk.ρ*σ_dr_t*r_g
σ_dr[:,3] = 4*d_c*leaf.ρ*σ_dr_l*r_g

σ_dr[1,:] *= reflva^2
σ_dr[2,:] *= reflha^2

############################

# TODO: Perform this vh polarization cleanup (cmd-F sgvh)

# Row is 
# Columns is 
σ_vh = zeros(3,3)

############################
# Backscat. cross section, hh pol.

sghdri1=σ_dr[2,1]/2
sghdri2=σ_dr[2,2]/2
sghdri3=σ_dr[2,3]/2

sghhr1 = 0
sghhr2 = 0
        
sghhd = σ_d[3,1]+σ_d[3,2]    
sghhr = sghhr1+sghhr2
sghhdr = σ_dr[2,1]+σ_dr[2,2]
sghdri = sghdri1+sghdri2
sghh  = sghhd+sghhr+sghhdr
sghhi=sghhd+sghhr+sghdri

############################
# Backscat. cross section, vv pol.
sgvdri1=σ_dr[1,1]/2.0
sgvdri2=σ_dr[1,2]/2.0
sgvdri3=σ_dr[1,3]/2.0

sgvvr1 = 0
sgvvr2 = 0

sgvvd = σ_d[1,1]+σ_d[1,2]
sgvvr = sgvvr1+sgvvr2
sgvvdr = σ_dr[1,1]+σ_dr[1,2]
sgvdri = sgvdri1+sgvdri2
sgvv  = sgvvd+sgvvr+sgvvdr
sgvvi=sgvvd+sgvvr+sgvdri

############################
# Backscat. cross section, vh pol.

sgvhd = σ_d[2,1]+σ_d[2,2]

sgvh11 = bsgvh11*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
sgvh21 = bsgvh21*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g

avh31  =   bsgvh31*reflvv*reflhc*cfun(a1,e1,d_c)*r_g
sgvh31 =  abs(2.0*real(avh31))

sgvhr1 = 0
sgvhr2 = 0
sgvhr3 = 0

sgvhdr1=sgvh11+sgvh21+sgvh31
sgvh1  =   σ_d[2,1]+sgvhr1+sgvh11+sgvh21+sgvh31
sgvhi1=σ_d[2,1]+sgvhr1+sgvh11+sgvh21
sgvhidr1=sgvh11+sgvh21

sgvh13 = bsgvh13*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
sgvh23 = bsgvh23*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g

avh33  =   bsgvh33*reflvv*reflhc*cfun(a1,e1,d_c)*r_g
sgvh33 =  abs(2.0*real(avh33))

sgvhdr3=sgvh13+sgvh23+sgvh33
sgvh3  =   σ_d[2,3]+sgvhr3+sgvh13+sgvh23+sgvh33

sgvhi3=σ_d[2,3]+sgvhr3+sgvh13+sgvh23
sgvhidr3=sgvh13+sgvh23

sgvh12 = factvh1*bsgvh12*(reflha)^2*funcm(2*K_vtⁱ,-2*K_htⁱ,d_t)
sgvh12=sgvh12*r_g

sgvh22 = factvh2*bsgvh22*(reflva)^2*funcp(2*K_vtⁱ,-2*K_htⁱ,d_t)
sgvh22 = sgvh22*r_g

avh32  =   factvh3*bsgvh32*reflvv*reflhc*cfun(a2,e2,d_t)*r_g
sgvh32 =  abs(2*real(avh32))
sgvhdr2=sgvh12+sgvh22+sgvh32

sgvh2  =   σ_d[2,2]+sgvhr2+sgvh12+sgvh22+sgvh32
                
sgvhi2=σ_d[2,2]+sgvhr2+sgvh12+sgvh22
sgvhidr2=sgvh12+sgvh22

sgvh = sgvh1+sgvh2
sgvhi = sgvhi1+sgvhi2
sgvhidr=sgvhidr1+sgvhidr2

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

shhg,svvg,svhg,gd = grdoh(ksig)

shht = sghh + shhg
svvt = sgvv + svvg
svht = sgvh + svhg
shhti=sghhi+shhg
svvti=sgvvi+svvg
svhti=sgvhi+svhg

###############################

output = Forest_Scattering_Output(θ_iᵈ, σ_d_db[3,1], σ_d_db[3,2], σ_d_db[3,3], 
                                  σ_dr_db[2,1], σ_dr_db[2,2], σ_dr_db[2,3],
                                  10.0*log10(shht), shhdd, shhdrd,
                                  
                                  σ_d_db[1,1], σ_d_db[1,2], σ_d_db[1,3],
                                  σ_dr_db[1,1], σ_dr_db[1,2], σ_dr_db[1,3],
                                  10.0*log10(svvt), svvdd, svvdrd,
                                  
                                  σ_d_db[2,1], σ_d_db[2,2], σ_d_db[2,3],
                                  svhdrd1, svhdrd2, svhdrd3,
                                  10.0*log10(svht), svhdd, svhdrd,
                                  
                                  10.0*log10(shhg), 10.0*log10(svhg), 10.0*log10(svvg),
                                  
                                  sghhoi, sgvhoi, sgvvoi,
                                  
                                  10.0*log10(sgvhidr1), 10.0*log10(sgvhidr2), 10.0*log10(sgvhidr3),
                                  
                                  shht/shhti, svht/svhti, svvt/svvti, 
                                  ath1, atv1, temp0h, temp0v)

check_output_matches_fortran(output)