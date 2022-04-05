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

@unpack bfrghz, amajcm, bmincm, tmm, ρ_l, ϵ_l, 
ntypel, parml, r_b1, l_b1, ρ_b1, ϵ_b1, ntypeb1, 
parmb1, r_b2, l_b2, ρ_b2, ϵ_b2, ntypeb2, parmb2, 
r_t, l_t, ρ_t, ϵ_t, ntypet, parmt, 
d_c, d_t, ϵ_g, l, sig = input_params

# Conversion to standard metric

const c           = 3e8             # Speed of light 
const bfr         = bfrghz*1e9      # GHz to Hz
const amaj_temp   = amajcm*1e-2     # cm to m
const bmin_temp   = bmincm*1e-2     # cm to m
const t_temp      = tmm*1e-3        # mm to m
const radb1m      = 0.5*r_b1*1e-2  # diameter to radius, cm to m
const radb2m      = 0.5*r_b2*1e-2  # diameter to radius, cm to m
const radtm_temp  = 0.5*r_t*1e-2   # diameter to radius, cm to m
const lm          = l*1e-2          # cm to m
const s           = sig*1e-2        # cm to m (surface rms height)
const k₀          = 2π*bfr/c        # Free space wave number (2π/λ)
const zk          = k₀              # 
const n_ϕ         = 41              # No. of ϕ
const n_θ         = 37              # No. of θ
const ej          = 0.0 + 1.0im 

## 
## Calculation of Parameters
## 

# Calculate integral of p(θ)*sin(θ)^2 from 0 to pi
# That is, the average integral over inclinations
const ail = sum(x -> prob(x, ntypel, parml) * sin(x)^2 * π/n_θ, collect(1:n_θ) * π/n_θ)

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

athc = zeros(20)
atvc = zeros(20)
atht = zeros(20)
atvt = zeros(20)

## Common Block functionality 

branch1 = wood(ϵ_b1, radb1m, l_b1)
branch2 = wood(ϵ_b2, radb2m, l_b2)
leaf_common = leaf(ϵ_l, amaj_temp, bmin_temp, t_temp)
trunk_common = wood(ϵ_t, radtm_temp, l_t)
parm2_common = parm2(0.0, 0.0, 0.0, 0.0)

# Calculation of skin depth skdh skdv and the bistatic cross sections sghh,sghv,sgvv

afhhl, afvvl = afsal(θ_iʳ, ail, leaf_common)
sbdl, sbdrl, sbvh1l, sbvh3l = asal(θ_iʳ, ntypel, parml, leaf_common)

# Compute scattering amplitudes from primary branches

afb1 = woodf(ntypeb1, parmb1, branch1, parm2_common)

afhhb1 = complex(abs(real(afb1[2,2])),abs(imag(afb1[2,2])))
afvvb1 = complex(abs(real(afb1[1,1])),abs(imag(afb1[1,1])))

sbd_b1, sbdr_b1, sbvh1b1,sbvh3b1 = woodb(ntypeb1, parmb1, branch1, parm2_common)

# compute scattering amplitudes from secondary branches

afb2 = woodf(ntypeb2,parmb2, branch2, parm2_common)

afhhb2 = complex(abs(real(afb2[2,2])),abs(imag(afb2[2,2])))
afvvb2 = complex(abs(real(afb2[1,1])),abs(imag(afb2[1,1])))

sbd_b2, sbdr_b2, sbvh1b2,sbvh3b2 = woodb(ntypeb2, parmb2, branch2, parm2_common)

# Compute scattering amplitudes from trunks

aft = woodf(ntypet,parmt, trunk_common, parm2_common)
sbd_t, sbdr_t, sbvh1t,sbvh3t = woodb(ntypet,parmt, trunk_common, parm2_common)

# Using reciprocity and scatterer symmetry to calculate rho*sigma

bsgd1 = ρ_b1*sbd_b1 + ρ_b2*sbd_b2 + ρ_l*sbdl
bsgd2 = ρ_t*sbd_t
bsgd3 = ρ_l*sbdl

bsgdr1 = ρ_b1*sbdr_b1 + ρ_b2*sbdr_b2 + ρ_l*sbdrl
bsgdr2 = ρ_t*sbdr_t
bsgdr3 = ρ_l*sbdrl

bsgvh11 = ρ_b1*sbvh1b1 + ρ_b2*sbvh1b2 +ρ_l*sbvh1l
bsgvh12 = ρ_t*sbvh1t
bsgvh13 = ρ_l*sbvh1l
        
bsgvh31 = ρ_b1*sbvh3b1 + ρ_b2*sbvh3b2 +ρ_l*sbvh3l
bsgvh32 = ρ_t*sbvh3t
bsgvh33 = ρ_l*sbvh3l 

bsgvh21=bsgvh11
bsgvh22=bsgvh12
bsgvh23=bsgvh13

afhh1 = ρ_l*afhhl + ρ_b1*afhhb1 + ρ_b2*afhhb2
afhh2 = ρ_t*aft[2,2]

afvv1 = ρ_l*afvvl + ρ_b1*afvvb1 + ρ_b2*afvvb2
afvv2 = ρ_t*aft[1,1]

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
athc[ip]=ath1
atvc[ip]=atv1

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
atht[ip]=temp0h
atvt[ip]=temp0v

############################
# Calculation of reflection coefficient from ground 

rgh = (cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))
rgv = (ϵ_g*cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(ϵ_g*cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))

K_hcⁱ, K_vcⁱ, K_htⁱ, K_vtⁱ = imag.((K_hc, K_vc, K_ht, K_vt))

kpvc1 = conj(K_vc)
kmhc1 = conj(K_hc)

kpvc2 = conj(K_vt)
kmhc2 = conj(K_ht)

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
e1     = (ej)*(kpvc1-kmhc1)
a2     = (ej)*(K_vt-K_ht)
e2     = (ej)*(kpvc2-kmhc2)

factvh1 = exp(-2*(K_hcⁱ-K_vcⁱ)*d_c)
factvh2 = exp(-2*(K_vcⁱ-K_hcⁱ)*d_c)
factvh3 = exp((-a1-e1)*d_c)

############################
# Backscat. cross section, hh pol.

sghhd1 = bsgd1[3]*funcm(2*K_hcⁱ,2*K_hcⁱ,d_c)
sghhd2 = bsgd2[3]*dattenh1*funcm(2*K_htⁱ,2*K_htⁱ,d_t)
sghhd3 = bsgd3[3]*funcm(2*K_hcⁱ,2*K_hcⁱ,d_c)


sghhdr1=4*d_c*bsgdr1[2]*reflha^2*r_g
sghhdr2=4*d_t*bsgdr2[2]*reflha^2*r_g
sghhdr3=4*d_c*bsgdr3[2]*reflha^2*r_g

sghdri1=sghhdr1/2
sghdri2=sghhdr2/2
sghdri3=sghhdr3/2

sghhr1 = 0
sghhr2 = 0
        
sghhd = sghhd1+sghhd2    
sghhr = sghhr1+sghhr2
sghhdr = sghhdr1+sghhdr2
sghdri = sghdri1+sghdri2
sghh  = sghhd+sghhr+sghhdr
sghhi=sghhd+sghhr+sghdri

############################
# Backscat. cross section, vv pol.

sgvvd1 = bsgd1[1]*funcm(2*K_vcⁱ, 2*K_vcⁱ, d_c)
sgvvd2 = bsgd2[1]*dattenv1*funcm(2*K_vtⁱ, 2*K_vtⁱ,d_t)
sgvvd3 = bsgd3[1]*funcm(2*K_vcⁱ,2*K_vcⁱ,d_c)

sgvvdr1=4*d_c*bsgdr1[1]*reflva^2*r_g
sgvvdr2=4*d_t*bsgdr2[1]*reflva^2*r_g
sgvvdr3=4*d_c*bsgdr3[1]*reflva^2*r_g
sgvdri1=sgvvdr1/2.0
sgvdri2=sgvvdr2/2.0
sgvdri3=sgvvdr3/2.0

sgvvr1 = 0
sgvvr2 = 0

sgvvd = sgvvd1+sgvvd2
sgvvr = sgvvr1+sgvvr2
sgvvdr = sgvvdr1+sgvvdr2
sgvdri = sgvdri1+sgvdri2
sgvv  = sgvvd+sgvvr+sgvvdr
sgvvi=sgvvd+sgvvr+sgvdri

############################
# Backscat. cross section, vh pol.

sgvhd1 = bsgd1[2]*funcm(2*K_vcⁱ,2*K_hcⁱ,d_c)
sgvhd3 = bsgd3[2]*funcm(2*K_vcⁱ,2*K_hcⁱ,d_c)
sgvhd2 = bsgd2[2]*dattenvh1*funcm(2*K_vtⁱ,2*K_htⁱ,d_t)
sgvhd = sgvhd1+sgvhd2

sgvh11 = bsgvh11*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
sgvh21 = bsgvh21*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g

avh31  =   bsgvh31*reflvv*reflhc*cfun(a1,e1,d_c)*r_g
sgvh31 =  abs(2.0*real(avh31))

sgvhr1 = 0
sgvhr2 = 0
sgvhr3 = 0

sgvhdr1=sgvh11+sgvh21+sgvh31
sgvh1  =   sgvhd1+sgvhr1+sgvh11+sgvh21+sgvh31
sgvhi1=sgvhd1+sgvhr1+sgvh11+sgvh21
sgvhidr1=sgvh11+sgvh21

sgvh13 = bsgvh13*(reflha)^2*funcm(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g
sgvh23 = bsgvh23*(reflva)^2*funcp(2*K_vcⁱ,-2*K_hcⁱ,d_c)*r_g

avh33  =   bsgvh33*reflvv*reflhc*cfun(a1,e1,d_c)*r_g
sgvh33 =  abs(2.0*real(avh33))

sgvhdr3=sgvh13+sgvh23+sgvh33
sgvh3  =   sgvhd3+sgvhr3+sgvh13+sgvh23+sgvh33

sgvhi3=sgvhd3+sgvhr3+sgvh13+sgvh23
sgvhidr3=sgvh13+sgvh23

sgvh12 = factvh1*bsgvh12*(reflha)^2*funcm(2*K_vtⁱ,-2*K_htⁱ,d_t)
sgvh12=sgvh12*r_g

sgvh22 = factvh2*bsgvh22*(reflva)^2*funcp(2*K_vtⁱ,-2*K_htⁱ,d_t)
sgvh22 = sgvh22*r_g

avh32  =   factvh3*bsgvh32*reflvv*reflhc*cfun(a2,e2,d_t)*r_g
sgvh32 =  abs(2.0*real(avh32))
sgvhdr2=sgvh12+sgvh22+sgvh32

sgvh2  =   sgvhd2+sgvhr2+sgvh12+sgvh22+sgvh32
                
sgvhi2=sgvhd2+sgvhr2+sgvh12+sgvh22
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
shhdd1 = 10.0*log10(sghhd1)
shhdd2 = 10.0*log10(sghhd2)
shhdd3 = 10.0*log10(sghhd3)

shhdrd	= 10.0*log10(sghhdr)
shhdrd1	= 10.0*log10(sghhdr1)
shhdrd2	= 10.0*log10(sghhdr2)
shhdrd3	= 10.0*log10(sghhdr3)

shhdri	= 10.0*log10(sghdri)
shhrd	= 10.0*log10(sghhr)

svhdd	= 10.0*log10(abs(sgvhd))
svhdd1	= 10.0*log10(abs(sgvhd1))
svhdd2	= 10.0*log10(abs(sgvhd2))
svhdd3	= 10.0*log10(abs(sgvhd3))
sgvhdr	= sgvhdr1+sgvhdr2+sgvhdr3

svhdrd	= 10.0*log10(abs(sgvhdr))
svhdrd1	= 10.0*log10(abs(sgvhdr1))
svhdrd2	= 10.0*log10(abs(sgvhdr2))
svhdrd3	= 10.0*log10(abs(sgvhdr3))            
sgvho	= 10.0*log10(sgvhd+sgvhdr+sgvhr)

svhrd	= 10.0*log10(abs(sgvhr))
svvdd	= 10.0*log10(sgvvd)
svvdd1	= 10.0*log10(sgvvd1)
svvdd2	= 10.0*log10(sgvvd2)
svvdd3	= 10.0*log10(sgvvd3)

svvdrd	= 10.0*log10(sgvvdr)
svvdrd1	= 10.0*log10(sgvvdr1)
svvdrd2	= 10.0*log10(sgvvdr2)
svvdrd3	= 10.0*log10(sgvvdr3)

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

output = Forest_Scattering_Output(θ_iᵈ, shhdd1, shhdd2, shhdd3, 
                                  shhdrd1, shhdrd2, shhdrd3,
                                  10.0*log10(shht), shhdd, shhdrd,
                                  
                                  svvdd1, svvdd2, svvdd3, 
                                  svvdrd1, svvdrd2, svvdrd3,
                                  10.0*log10(svvt), svvdd, svvdrd,
                                  
                                  svhdd1, svhdd2, svhdd3, 
                                  svhdrd1, svhdrd2, svhdrd3,
                                  10.0*log10(svht), svhdd, svhdrd,
                                  
                                  10.0*log10(shhg), 10.0*log10(svhg), 10.0*log10(svvg),
                                  
                                  sghhoi, sgvhoi, sgvvoi,
                                  
                                  10.0*log10(sgvhidr1), 10.0*log10(sgvhidr2), 10.0*log10(sgvhidr3),
                                  
                                  shht/shhti, svht/svhti, svvt/svvti, 
                                  ath1, atv1, temp0h, temp0v)

check_output_matches_fortran(output)