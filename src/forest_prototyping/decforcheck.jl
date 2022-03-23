# calculate scattering coefficient for forest based on distorted born approximation

using YAML
using SpecialFunctions
using QuadGK
using Parameters
using Test

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
d1, d2, epsg, l, sig = input_params

println(input_params)

# Convert complex numbers 
epsl_temp = complex(epsl_temp[1], epsl_temp[2])
epsb1 = complex(epsb1[1], epsb1[2])
epsb2 = complex(epsb2[1], epsb2[2])
epst_temp = complex(epst_temp[1], epst_temp[2])
epsg = complex(epsg[1], epsg[2])

# Conversion to standard metric

c      = 3.0e+08            # Speed of light 

bfr    = bfrghz*1.0e09      # GHz to Hz
amaj_temp   = amajcm*1.0e-02     # cm to m
bmin_temp   = bmincm*1.0e-02     # cm to m
t_temp      = tmm*1.0e-3         # mm to m
radb1m = 0.5*radb1*1.0e-02  # diameter to radius, cm to m
radb2m = 0.5*radb2*1.0e-02  # diameter to radius, cm to m
radtm_temp  = 0.5*radt*1.0e-02   # diameter to radius, cm to m
lm     = l*1.0e-02          # cm to m
sigm   = sig*1.0e-02        # cm to m
    
ak0_temp    = 2.0*π*bfr/c       # Free space wave number
zk     = ak0_temp                # 
ksig   = ak0_temp*sigm           # 
nph_temp    = 41                 # 
nth_temp    = 37                 # 

## 
## Calculation of Parameters
## 

# Calculate integral of p(th)*sin(th)^2 from 0 to pi
# That is, the average integral over inclinations
ail = sum(x -> prob(x, ntypel, parml) * sin(x)^2 * π/nth_temp, collect(1:nth_temp) * π/nth_temp)

# Loop over angle of incidence-theti 
ip = 0
phiir=0.0
ip=ip+1
θ_i  = 40
θ_i = (θ_i < 0.001 ? θ_i=0.1 : θ_i) 
thetir  = deg2rad(θ_i)
si = sin(thetir)
ci = cos(thetir)
si2 = si^2 

# Incidence geometry in radians
geom_i = IncidentGeometry(thetir, phiir)

# Roughness factor for ground scattering in direct-reflected term
grough = exp(-4.0*(ak0_temp*sigm*ci)^2)

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
integ_common = integ(nph_temp, nth_temp)
leaf_common = leaf(epsl_temp, amaj_temp, bmin_temp, t_temp)
parm1_common = parm1(0.0, 0.0, 0.0, 0.0)
parm2_common = parm2(0.0, 0.0, 0.0, 0.0)
a_common = a(0.0, 0.0, 0.0, 0.0, bfr, epsg)
b_common = b(zk, sig, 0.0, 0.0)

# Calculation of skin depth skdh skdv and the bistatic cross sections sghh,sghv,sgvv

afhhl, afvvl = afsal(thetir, ail, leaf_common, data_common)

sbhhdl, sbvhdl, sbvvdl, sbhdrl, sbvdrl, sbvh1l, sbvh3l = asal(thetir, ntypel, parml, leaf_common, data_common, integ_common)

# Compute scattering amplitudes from primary branches

ntypeb1 = 11

afb1 = woodf(index,geom_i,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

afhhb1 = complex(abs(real(afb1[2,2])),abs(imag(afb1[2,2])))
afvvb1 = complex(abs(real(afb1[1,1])),abs(imag(afb1[1,1])))

sbvvdb1,sbvhdb1,sbhhdb1,sbhdrb1,sbvdrb1,sbvh1b1,sbvh3b1 = woodb(index,geom_i,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

# compute scattering amplitudes from secondary branches

index=1
data_common.epsb=epsb2
data_common.lb=lb2
data_common.radbm=radb2m

afb2 = woodf(index,geom_i,ntypeb2,parmb2, data_common, integ_common, parm1_common, parm2_common)

afhhb2 = complex(abs(real(afb2[2,2])),abs(imag(afb2[2,2])))
afvvb2 = complex(abs(real(afb2[1,1])),abs(imag(afb2[1,1])))

sbvvdb2,sbvhdb2,sbhhdb2,sbhdrb2,sbvdrb2,sbvh1b2,sbvh3b2 = woodb(index,geom_i,ntypeb2,parmb2, data_common, integ_common, parm1_common, parm2_common)

# Compute scattering amplitudes from trunks

index = 2

aft = woodf(index,geom_i,ntypet,parmt, data_common, integ_common, parm1_common, parm2_common)

sbvvdt,sbvhdt,sbhhdt,sbhdrt,sbvdrt,sbvh1t,sbvh3t = woodb(index,geom_i,ntypet,parmt, data_common, integ_common, parm1_common, parm2_common)

# Using reciprocity and scatterer symmetry to calculate rho*sigma

bsghhd1 = rhob1*sbhhdb1 + rhob2*sbhhdb2 +rhol*sbhhdl
bsghhd2 = rhot*sbhhdt
bsghhd3 = rhol*sbhhdl

bsgvhd1 = rhob1*sbvhdb1 +rhob2*sbvhdb2 + rhol*sbvhdl
bsgvhd2 = rhot*sbvhdt
bsgvhd3 = rhol*sbvhdl

bsgvvd1 = rhob1*sbvvdb1 + rhob2*sbvvdb2 + rhol*sbvvdl
bsgvvd2 = rhot*sbvvdt
bsgvvd3 = rhol*sbvvdl

bsghdr1 = rhob1*sbhdrb1 + rhob2*sbhdrb2 + rhol*sbhdrl
bsghdr2 = rhot*sbhdrt
bsghdr3 = rhol*sbhdrl

bsgvdr1 = rhob1*sbvdrb1 + rhob2*sbvdrb2 + rhol*sbvdrl
bsgvdr2 = rhot*sbvdrt
bsgvdr3 = rhol*sbvdrl

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

kv1 = data_common.ak0*ci+(2.0*pi*afvv1)/(data_common.ak0*ci)
kh1 = data_common.ak0*ci+(2.0*pi*afhh1)/(data_common.ak0*ci)

ath1= abs(imag(kh1))
atv1= abs(imag(kv1))
kh1=complex(real(kh1),abs(imag(kh1)))
kv1=complex(real(kv1),abs(imag(kv1)))

atv1 = (atv1 <= 1.0E-20 ? 0.0001 : atv1)

skdh1= 1.0/ath1
skdv1= 1.0/atv1
athc[ip]=ath1
atvc[ip]=atv1

############################

kh2 = data_common.ak0*ci+(2.0*pi*afhh2)/(data_common.ak0*ci)
kv2 = data_common.ak0*ci+(2.0*pi*afvv2)/(data_common.ak0*ci)
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

rgh = (ci-sqrt(epsg-si2))/(ci+sqrt(epsg-si2))
rgv = (epsg*ci-sqrt(epsg-si2))/(epsg*ci+sqrt(epsg-si2))
a_common.kvim1= imag(kv1)
a_common.khim1= imag(kh1)
kpvc1 = conj(kv1)
kmhc1 = conj(kh1)

a_common.kvim2= imag(kv2)
a_common.khim2= imag(kh2)
kpvc2 = conj(kv2)
kmhc2 = conj(kh2)

reflhh = rgh*exp(parm1_common.ej*(kh1+kh1)*d1+parm1_common.ej*(kh2+kh2)*d2)
reflhc = conj(reflhh)
reflvv = rgv*exp(parm1_common.ej*(kv1+kv1)*d1+parm1_common.ej*(kv2+kv2)*d2)
reflha = abs(reflhh)
reflva = abs(reflvv)

x1     = a_common.khim1
y1     = a_common.kvim1
x2     = a_common.khim2
y2     = a_common.kvim2
dattenh1 = exp(-2.0*(x1+x1)*d1)
dattenh2 = exp(-2.0*(x2+x2)*d2)
dattenv1 = exp(-2.0*(y1+y1)*d1)
dattenv2 = exp(-2.0*(y2+y2)*d2)
dattenvh1 = exp(-2.0*(y1+x1)*d1)
dattenvh2 = exp(-2.0*(y2+x2)*d2)
a1     = (parm1_common.ej)*(kv1-kh1)
e1     = (parm1_common.ej)*(kpvc1-kmhc1)
a2     = (parm1_common.ej)*(kv2-kh2)
e2     = (parm1_common.ej)*(kpvc2-kmhc2)

factvh1 = exp(-2.0*(x1-y1)*d1)
factvh2 = exp(-2.0*(y1-x1)*d1)
factvh3 = exp((-a1-e1)*d1)

sghhd1 = bsghhd1*funcm(2.0*x1,2.0*x1,d1)
sghhd2 = bsghhd2*dattenh1*funcm(2.0*x2,2.0*x2,d2)
sghhd3 = bsghhd3*funcm(2.0*x1,2.0*x1,d1)

sghhdr1=4.0*d1*bsghdr1*reflha^2*grough
sghhdr2=4.0*d2*bsghdr2*reflha^2*grough
sghhdr3=4.0*d1*bsghdr3*reflha^2*grough
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

sgvvd1 = bsgvvd1*funcm(2.0*y1,2.0*y1,d1)
sgvvd2 = bsgvvd2*dattenv1*funcm(2.0*y2,2.0*y2,d2)
sgvvd3 = bsgvvd3*funcm(2.0*y1,2.0*y1,d1)

sgvvdr1=4.0*d1*bsgvdr1*reflva^2*grough
sgvvdr2=4.0*d2*bsgvdr2*reflva^2*grough
sgvvdr3=4.0*d1*bsgvdr3*reflva^2*grough
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

sgvhd1 = bsgvhd1*funcm(2.0*y1,2.0*x1,d1)
sgvhd3 = bsgvhd3*funcm(2.0*y1,2.0*x1,d1)
sgvhd2 = bsgvhd2*dattenvh1*funcm(2.0*y2,2.0*x2,d2)
sgvhd = sgvhd1+sgvhd2

sgvh11 = bsgvh11*(reflha)^2*funcm(2.0*y1,-2.0*x1,d1)*grough
sgvh21 = bsgvh21*(reflva)^2*funcp(2.0*y1,-2.0*x1,d1)*grough

avh31  =   bsgvh31*reflvv*reflhc*cfun(a1,e1,d1)*grough
sgvh31 =  abs(2.0*real(avh31))

sgvhr1 = 0
sgvhr2 = 0
sgvhr3 = 0

sgvhdr1=sgvh11+sgvh21+sgvh31
sgvh1  =   sgvhd1+sgvhr1+sgvh11+sgvh21+sgvh31
sgvhi1=sgvhd1+sgvhr1+sgvh11+sgvh21
sgvhidr1=sgvh11+sgvh21

sgvh13 = bsgvh13*(reflha)^2*funcm(2.0*y1,-2.0*x1,d1)*grough
sgvh23 = bsgvh23*(reflva)^2*funcp(2.0*y1,-2.0*x1,d1)*grough

avh33  =   bsgvh33*reflvv*reflhc*cfun(a1,e1,d1)*grough
sgvh33 =  abs(2.0*real(avh33))

sgvhdr3=sgvh13+sgvh23+sgvh33
sgvh3  =   sgvhd3+sgvhr3+sgvh13+sgvh23+sgvh33

sgvhi3=sgvhd3+sgvhr3+sgvh13+sgvh23
sgvhidr3=sgvh13+sgvh23

sgvh12 = factvh1*bsgvh12*(reflha)^2*funcm(2.0*y2,-2.0*x2,d2)
sgvh12=sgvh12*grough

sgvh22 = factvh2*bsgvh22*(reflva)^2*funcp(2.0*y2,-2.0*x2,d2)
sgvh22 = sgvh22*grough

avh32  =   factvh3*bsgvh32*reflvv*reflhc*cfun(a2,e2,d2)*grough
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

klx=data_common.ak0*lm
kly=klx
ksig=data_common.ak0*sigm

shhg,svvg,svhg,ghhd,gvvd,gvhd = grdoh(ksig, θ_i, a_common, b_common)

shht = sghh + shhg
svvt = sgvv + svvg
svht = sgvh + svhg
shhti=sghhi+shhg
svvti=sgvvi+svvg
svhti=sgvhi+svhg

###############################

output = Forest_Scattering_Output(θ_i, shhdd1, shhdd2, shhdd3, 
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

# ak0  - free space wave number.
# d - slab thickness
# bfr - frequency in hz
# eps1-relative dielctric constant of scatterer.sth
# eps2-relative dielectric constant of water droplets
# epsg- relative dielectic constant of ground
# rho- density of scatterers- no./m^3.
# kh - horizontal propagation constant
# kmv - vertical propagation constant
# theti- angle of incidence-degrees.
# atmh - attenuation of mean horizontal wave
# atmv - attenuation of mean vertical   wave
# skdh,skdv     - horizontal and vertical skin depth in m.
# sghh ,sgvh ,sgvv - average backscattering coefficients.
# sghho ,sghvo ,sgvvo - average backscattering in lb.