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

INPUT_FILE = "deccheckin.yaml"
input_params = parameters_from_yaml(INPUT_FILE)

@unpack bfrghz, amajcm, bmincm, tmm, rhol, epsl_temp, 
ntypel, parml, radb1, lb1, rhob1, epsb1, ntypeb1, 
parmb1, radb2, lb2, rhob2, epsb2, ntypeb2, parmb2, 
radt, lt_temp, rhot, epst_temp, ntypet, parmt, 
d_c, d_t, epsg, l, sig = input_params

# Conversion to standard metric

const c      = 3.0e+08            # Speed of light 
const bfr    = bfrghz*1.0e09      # GHz to Hz
const amaj_temp   = amajcm*1.0e-02     # cm to m
const bmin_temp   = bmincm*1.0e-02     # cm to m
t_temp      = tmm*1.0e-3         # mm to m
radb1m = 0.5*radb1*1.0e-02  # diameter to radius, cm to m
radb2m = 0.5*radb2*1.0e-02  # diameter to radius, cm to m
radtm_temp  = 0.5*radt*1.0e-02   # diameter to radius, cm to m
lm     = l*1.0e-02          # cm to m
sigm   = sig*1.0e-02        # cm to m
    
ak0_temp    = 2.0*π*bfr/c       # Free space wave number (2π/λ)
zk     = ak0_temp                # 
ksig   = ak0_temp*sigm           # 

const n_ϕ    = 41   #  
const n_θ    = 37   #  

## 
## Calculation of Parameters
## 

# Calculate integral of p(θ)*sin(θ)^2 from 0 to pi
# That is, the average integral over inclinations
ail = sum(x -> prob(x, ntypel, parml) * sin(x)^2 * π/n_θ, collect(1:n_θ) * π/n_θ)

# Loop over angle of incidence-θ_iᵈ 
ip = 1

θ_iᵈ  = 40
θ_iᵈ  = (θ_iᵈ < 0.001 ? θ_iᵈ=0.1 : θ_iᵈ) 
θ_iʳ  = deg2rad(θ_iᵈ)
ϕ_iᵈ  = 0.0 
ϕ_iʳ  = deg2rad(0.0) 

# Incidence geometry in radians
geom_i = IncidentGeometry(θ_iʳ, ϕ_iʳ)

# Roughness factor for ground scattering in direct-reflected term
grough = exp(-4.0*(ak0_temp*sigm*cos(θ_iʳ))^2)

athc = zeros(20)
atvc = zeros(20)
atht = zeros(20)
atvt = zeros(20)

## Common Block functionality 

index=1
epsb_temp=epsb1
lb_temp=lb1
radbm_temp=radb1m

data_common = data(ak0_temp, epsb_temp, epst_temp, radbm_temp, radtm_temp, lb_temp, lt_temp)
integ_common = integ(n_ϕ, n_θ)
leaf_common = leaf(epsl_temp, amaj_temp, bmin_temp, t_temp)
parm1_common = parm1(0.0, 0.0, 0.0, 0.0)
parm2_common = parm2(0.0, 0.0, 0.0, 0.0)
a_common = a(0.0, 0.0, 0.0, 0.0, bfr, epsg)
b_common = b(zk, sig, 0.0, 0.0)

# Calculation of skin depth skdh skdv and the bistatic cross sections sghh,sghv,sgvv

afhhl, afvvl = afsal(θ_iʳ, ail, leaf_common, data_common)

sbdl, sbdrl, sbvh1l, sbvh3l = asal(θ_iʳ, ntypel, parml, leaf_common, data_common, integ_common)

# Compute scattering amplitudes from primary branches

ntypeb1 = 11

afb1 = woodf(index,geom_i,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

afhhb1 = complex(abs(real(afb1[2,2])),abs(imag(afb1[2,2])))
afvvb1 = complex(abs(real(afb1[1,1])),abs(imag(afb1[1,1])))


sbd_b1, sbdr_b1, sbvh1b1,sbvh3b1 = woodb(index,geom_i,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

# @show sbd_b1, sbdr_b1, sbvh1b1,sbvh3b1 

# exit()

# compute scattering amplitudes from secondary branches

index=1
data_common.epsb=epsb2
data_common.lb=lb2
data_common.radbm=radb2m

afb2 = woodf(index,geom_i,ntypeb2,parmb2, data_common, integ_common, parm1_common, parm2_common)

afhhb2 = complex(abs(real(afb2[2,2])),abs(imag(afb2[2,2])))
afvvb2 = complex(abs(real(afb2[1,1])),abs(imag(afb2[1,1])))

sbd_b2, sbdr_b2, sbvh1b2,sbvh3b2 = woodb(index,geom_i,ntypeb2,parmb2, data_common, integ_common, parm1_common, parm2_common)

# Compute scattering amplitudes from trunks

index = 2

aft = woodf(index,geom_i,ntypet,parmt, data_common, integ_common, parm1_common, parm2_common)

sbd_t, sbdr_t, sbvh1t,sbvh3t = woodb(index,geom_i,ntypet,parmt, data_common, integ_common, parm1_common, parm2_common)

# Using reciprocity and scatterer symmetry to calculate rho*sigma

bsgd1 = rhob1*sbd_b1 + rhob2*sbd_b2 + rhol*sbdl
bsgd2 = rhot*sbd_t
bsgd3 = rhol*sbdl

bsgdr1 = rhob1*sbdr_b1 + rhob2*sbdr_b2 + rhol*sbdrl
bsgdr2 = rhot*sbdr_t
bsgdr3 = rhol*sbdrl

bsgvh11 = rhob1*sbvh1b1 + rhob2*sbvh1b2 +rhol*sbvh1l
bsgvh12 = rhot*sbvh1t
bsgvh13 = rhol*sbvh1l
        
bsgvh31 = rhob1*sbvh3b1 + rhob2*sbvh3b2 +rhol*sbvh3l
bsgvh32 = rhot*sbvh3t
bsgvh33 = rhol*sbvh3l 

bsgvh21=bsgvh11
bsgvh22=bsgvh12
bsgvh23=bsgvh13

afhh1 = rhol*afhhl + rhob1*afhhb1 + rhob2*afhhb2
afhh2 = rhot*aft[2,2]

afvv1 = rhol*afvvl + rhob1*afvvb1 + rhob2*afvvb2
afvv2 = rhot*aft[1,1]

############################

# 3.1.8??? 
kv1 = data_common.ak0*cos(θ_iʳ)+(2π*afvv1)/(data_common.ak0*cos(θ_iʳ))
kh1 = data_common.ak0*cos(θ_iʳ)+(2π*afhh1)/(data_common.ak0*cos(θ_iʳ))

ath1= abs(imag(kh1))
atv1= abs(imag(kv1))
kh1=complex(real(kh1),abs(imag(kh1)))
kv1=complex(real(kv1),abs(imag(kv1)))

atv1 = (atv1 <= 1.0E-20 ? 0.0001 : atv1)

skdh1= 1/ath1
skdv1= 1/atv1
athc[ip]=ath1
atvc[ip]=atv1

############################

kh2 = data_common.ak0*cos(θ_iʳ)+(2π*afhh2)/(data_common.ak0*cos(θ_iʳ))
kv2 = data_common.ak0*cos(θ_iʳ)+(2π*afvv2)/(data_common.ak0*cos(θ_iʳ))
kh2=complex(real(kh2),abs(imag(kh2)))
kv2=complex(real(kv2),abs(imag(kv2)))

ath2= abs(imag(kh2))
atv2= abs(imag(kv2))
temp0h=abs(imag(kh1+kh2))
temp0v=abs(imag(kv1+kv2))
skdh2= 1/ath2
skdv2= 1/atv2
atht[ip]=temp0h
atvt[ip]=temp0v

############################

rgh = (cos(θ_iʳ)-sqrt(epsg-sin(θ_iʳ)^2))/(cos(θ_iʳ)+sqrt(epsg-sin(θ_iʳ)^2))
rgv = (epsg*cos(θ_iʳ)-sqrt(epsg-sin(θ_iʳ)^2))/(epsg*cos(θ_iʳ)+sqrt(epsg-sin(θ_iʳ)^2))
a_common.kvim1= imag(kv1)
a_common.khim1= imag(kh1)

K_hc = kh1
K_vc = kv1
K_ht = kh2
K_vt = kv2

kpvc1 = conj(kv1)
kmhc1 = conj(kh1)

a_common.kvim2= imag(kv2)
a_common.khim2= imag(kh2)
kpvc2 = conj(kv2)
kmhc2 = conj(kh2)

reflhh = rgh*exp(parm1_common.ej*(kh1+kh1)*d_c+parm1_common.ej*(kh2+kh2)*d_t)
reflhc = conj(reflhh)
reflvv = rgv*exp(parm1_common.ej*(kv1+kv1)*d_c+parm1_common.ej*(kv2+kv2)*d_t)
reflha = abs(reflhh)
reflva = abs(reflvv)

x1     = a_common.khim1
y1     = a_common.kvim1
x2     = a_common.khim2
y2     = a_common.kvim2
dattenh1 = exp(-2*(x1+x1)*d_c)
dattenh2 = exp(-2*(x2+x2)*d_t)
dattenv1 = exp(-2*(y1+y1)*d_c)
dattenv2 = exp(-2*(y2+y2)*d_t)
dattenvh1 = exp(-2*(y1+x1)*d_c)
dattenvh2 = exp(-2*(y2+x2)*d_t)
a1     = (parm1_common.ej)*(kv1-kh1)
e1     = (parm1_common.ej)*(kpvc1-kmhc1)
a2     = (parm1_common.ej)*(kv2-kh2)
e2     = (parm1_common.ej)*(kpvc2-kmhc2)

factvh1 = exp(-2*(x1-y1)*d_c)
factvh2 = exp(-2*(y1-x1)*d_c)
factvh3 = exp((-a1-e1)*d_c)

sghhd1 = bsgd1[3]*funcm(2*x1,2*x1,d_c)
sghhd2 = bsgd2[3]*dattenh1*funcm(2*x2,2*x2,d_t)
sghhd3 = bsgd3[3]*funcm(2*x1,2*x1,d_c)

sghhdr1=4.0*d_c*bsgdr1[2]*reflha^2*grough
sghhdr2=4.0*d_t*bsgdr2[2]*reflha^2*grough
sghhdr3=4.0*d_c*bsgdr3[2]*reflha^2*grough
sghdri1=sghhdr1/2.0
sghdri2=sghhdr2/2.0
sghdri3=sghhdr3/2.0

sghhr1 = 0
sghhr2 = 0
        
sghhd = sghhd1+sghhd2    
sghhr = sghhr1+sghhr2
sghhdr = sghhdr1+sghhdr2
sghdri = sghdri1+sghdri2
sghh  = sghhd+sghhr+sghhdr
sghhi=sghhd+sghhr+sghdri

############################

sgvvd1 = bsgd1[1]*funcm(2.0*y1,2.0*y1,d_c)
sgvvd2 = bsgd2[1]*dattenv1*funcm(2.0*y2,2.0*y2,d_t)
sgvvd3 = bsgd3[1]*funcm(2.0*y1,2.0*y1,d_c)

sgvvdr1=4.0*d_c*bsgdr1[1]*reflva^2*grough
sgvvdr2=4.0*d_t*bsgdr2[1]*reflva^2*grough
sgvvdr3=4.0*d_c*bsgdr3[1]*reflva^2*grough
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

sgvhd1 = bsgd1[2]*funcm(2.0*y1,2.0*x1,d_c)
sgvhd3 = bsgd3[2]*funcm(2.0*y1,2.0*x1,d_c)
sgvhd2 = bsgd2[2]*dattenvh1*funcm(2.0*y2,2.0*x2,d_t)
sgvhd = sgvhd1+sgvhd2

sgvh11 = bsgvh11*(reflha)^2*funcm(2.0*y1,-2.0*x1,d_c)*grough
sgvh21 = bsgvh21*(reflva)^2*funcp(2.0*y1,-2.0*x1,d_c)*grough

avh31  =   bsgvh31*reflvv*reflhc*cfun(a1,e1,d_c)*grough
sgvh31 =  abs(2.0*real(avh31))

sgvhr1 = 0
sgvhr2 = 0
sgvhr3 = 0

sgvhdr1=sgvh11+sgvh21+sgvh31
sgvh1  =   sgvhd1+sgvhr1+sgvh11+sgvh21+sgvh31
sgvhi1=sgvhd1+sgvhr1+sgvh11+sgvh21
sgvhidr1=sgvh11+sgvh21

sgvh13 = bsgvh13*(reflha)^2*funcm(2.0*y1,-2.0*x1,d_c)*grough
sgvh23 = bsgvh23*(reflva)^2*funcp(2.0*y1,-2.0*x1,d_c)*grough

avh33  =   bsgvh33*reflvv*reflhc*cfun(a1,e1,d_c)*grough
sgvh33 =  abs(2.0*real(avh33))

sgvhdr3=sgvh13+sgvh23+sgvh33
sgvh3  =   sgvhd3+sgvhr3+sgvh13+sgvh23+sgvh33

sgvhi3=sgvhd3+sgvhr3+sgvh13+sgvh23
sgvhidr3=sgvh13+sgvh23

sgvh12 = factvh1*bsgvh12*(reflha)^2*funcm(2.0*y2,-2.0*x2,d_t)
sgvh12=sgvh12*grough

sgvh22 = factvh2*bsgvh22*(reflva)^2*funcp(2.0*y2,-2.0*x2,d_t)
sgvh22 = sgvh22*grough

avh32  =   factvh3*bsgvh32*reflvv*reflhc*cfun(a2,e2,d_t)*grough
sgvh32 =  abs(2.0*real(avh32))
sgvhdr2=sgvh12+sgvh22+sgvh32

sgvh2  =   sgvhd2+sgvhr2+sgvh12+sgvh22+sgvh32
                
sgvhi2=sgvhd2+sgvhr2+sgvh12+sgvh22
sgvhidr2=sgvh12+sgvh22

sgvh = sgvh1+sgvh2
sgvhi = sgvhi1+sgvhi2
sgvhidr=sgvhidr1+sgvhidr2


############################

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

# klx=data_common.ak0*lm
# kly=klx
ksig=data_common.ak0*sigm

shhg,svvg,svhg,gd = grdoh(ksig, θ_iᵈ, a_common, b_common)

shht = sghh + shhg
svvt = sgvv + svvg
svht = sgvh + svhg
shhti=sghhi+shhg
svvti=sgvvi+svvg
svhti=sgvhi+svhg

###############################

println("--------------------------------------------------------------------------------")
theti = θ_iʳ
println("\nbackscat. cross section of forest in db     hh polarization coherent\n")
println("theti \t sighhd1 \t\t sighhd2 \t\t sighhd3")
println("$theti \t $shhdd1 \t $shhdd2 \t $shhdd3 \n")
println("theti \t sighhdr1 \t\t sighhdr2 \t\t sighhdr3")
println("$theti \t $shhdrd1 \t $shhdrd2 \t $shhdrd3 \n")
println("theti \t sighht \t\t sighhd \t\t sighhdr")
println("$theti \t $(10.0*log10(shht)) \t $shhdd \t $shhdrd \n")

println("--------------------------------------------------------------------------------")

println("\nbackscat. cross section of forest in db     vv polarization coherent\n")
println("theti \t sigvvd1 \t\t sigvvd2 \t\t sigvvd3")
println("$theti \t $svvdd1 \t $svvdd2 \t $svvdd3 \n")
println("theti \t sigvvdr1 \t\t sigvvdr2 \t\t sigvvdr3")
println("$theti \t $svvdrd1 \t $svvdrd2 \t $svvdrd3 \n")
println("theti \t sigvvt \t\t sigvvd \t\t sigvvdr")
println("$theti \t $(10.0*log10(svvt)) \t $svvdd \t $svvdrd \n")

println("--------------------------------------------------------------------------------")

println("\nbackscat. cross section of forest in db     vh polarization coherent\n")
println("theti \t sigvhd1 \t\t sigvhd2 \t\t sigvhd3")
println("$theti \t $svhdd1 \t $svhdd2 \t $svhdd3 \n")
println("theti \t sigvhdr1 \t\t sigvhdr2 \t\t sigvhdr3")
println("$theti \t $svhdrd1 \t $svhdrd2 \t $svhdrd3 \n")
println("theti \t sigvht \t\t sigvhd \t\t sigvhdr")
println("$theti \t $(10.0*log10(svht)) \t $svhdd \t $svhdrd \n")

println("--------------------------------------------------------------------------------")

println("\nbackscat. cross section of forest in db     total coherent\n")
println("theti \t sighht \t\t sigvht \t\t sigvvt")
println("$theti \t $(10.0*log10(shht)) \t $(10.0*log10(svht)) \t $(10.0*log10(svvt)) \n")
println("theti \t ghhd \t\t\t gvhd \t\t\t gvvd")
println("$theti \t $(10.0*log10(shhg)) \t $(10.0*log10(svhg)) \t $(10.0*log10(svvg)) \n")

println("--------------------------------------------------------------------------------")

println("\nbackscat. cross section of forest in db     total incoherent\n")
println("theti \t sighhi \t\t sigvhi \t\t sigvvi")
println("$theti \t $sghhoi \t $sgvhoi \t $sgvvoi \n")

println("--------------------------------------------------------------------------------")

println("\nVH incoherent terms due to double bounce from     both layers in dB\n")
println("theti \t sigvhi1 \t\t sigvhi2 \t\t sigvhi3")
println("$theti \t $(10.0*log10(sgvhidr1)) \t $(10.0*log10(sgvhidr2)) \t $(10.0*log10(sgvhidr3)) \n")

println("--------------------------------------------------------------------------------")

println("\nEnhancement factors for total canopy backscatter\n")
println("theti \t cofhh \t\t cofvh \t\t cofvv")
println("$theti \t $(shht/shhti) \t $(svht/svhti) \t $(svvt/svvti) \n")

println("theti \t athc \t\t atvc \t\t atht \t\t atvt")
println("$theti \t $ath1 \t $atv1 \t $temp0h \t $temp0v \n")

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