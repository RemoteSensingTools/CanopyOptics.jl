# calculate scattering coefficient for forest based on distorted born approximation

using YAML
using SpecialFunctions
using QuadGK

include("types.jl")
include("probabilities.jl")
include("subroutines.jl")

## 
## Load input parameters
## 

INPUT_FILE = "deccheckin.yaml"
input_params = YAML.load_file(INPUT_FILE)

# Frequency
bfrghz  = input_params["frequency"]

# Leaf Parameters 
amajcm  = input_params["leaf_parameters"]["major_axis"]
bmincm  = input_params["leaf_parameters"]["minor_axis"]
tmm     = input_params["leaf_parameters"]["thickness"]
rhol    = input_params["leaf_parameters"]["density"]
epsl_temp    = input_params["leaf_parameters"]["dielectric_constant"]
ntypel  = input_params["leaf_parameters"]["leaf_inclination_pdf_type"]
parml   = input_params["leaf_parameters"]["pdf_parameter"]

# Primary Branch Parameters
radb1   = input_params["primary_branch_parameters"]["diameter"]
lb1     = input_params["primary_branch_parameters"]["length"]
rhob1   = input_params["primary_branch_parameters"]["density"]
epsb1   = input_params["primary_branch_parameters"]["dielectric_constant"]
ntypeb1 = input_params["primary_branch_parameters"]["branch_inclination_pdf_type"]
parmb1  = input_params["primary_branch_parameters"]["pdf_parameter"]

# Secondary Branch Parameters
radb2   = input_params["secondary_branch_parameters"]["diameter"]
lb2     = input_params["secondary_branch_parameters"]["length"]
rhob2   = input_params["secondary_branch_parameters"]["density"]
epsb2   = input_params["secondary_branch_parameters"]["dielectric_constant"]
ntypeb2 = input_params["secondary_branch_parameters"]["branch_inclination_pdf_type"]
parmb2  = input_params["secondary_branch_parameters"]["pdf_parameter"]

# Trunk Parameters
radt    = input_params["trunk_parameters"]["diameter"]
lt_temp      = input_params["trunk_parameters"]["length"]
rhot    = input_params["trunk_parameters"]["density"]
epst_temp    = input_params["trunk_parameters"]["dielectric_constant"]
ntypet  = input_params["trunk_parameters"]["branch_inclination_pdf_type"]
parmt   = input_params["trunk_parameters"]["pdf_parameter"]
d1      = input_params["trunk_parameters"]["crown_height"]
d2      = input_params["trunk_parameters"]["trunk_height"]
epsg    = input_params["trunk_parameters"]["soil_dielectric"] 
l       = input_params["trunk_parameters"]["corr_length"]
sig     = input_params["trunk_parameters"]["rms_height"]

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
theti  = 40
theti = (theti < 0.001 ? theti=0.1 : theti) 
thetir  = theti*π/180.0
si = sin(thetir)
ci = cos(thetir)
si2 = si^2 
grough = exp(-4.0*(ak0_temp*sigm*ci)^2)

# print(epsl)

## Common Block functionality 

# data_struct = data(0, 0, 0, 0, 0, 0, 0)

index=1
epsb_temp=epsb1
lb_temp=lb1
radbm_temp=radb1m

## 

# a_common = a()
# b_common = b()
# ds_common = ds(nph_temp, nth_temp)
data_common = data(ak0_temp, epsb_temp, epst_temp, radbm_temp, radtm_temp, lb_temp, lt_temp)
integ_common = integ(nph_temp, nth_temp)
leaf_common = leaf(epsl_temp, amaj_temp, bmin_temp, t_temp)
parm1_common = parm1(0.0, 0.0, 0.0, 0.0)
parm2_common = parm2(0.0, 0.0, 0.0, 0.0)


afhhl, afvvl = afsal(thetir, ail, leaf_common, data_common)

println("afsal")
println(afhhl)
println(afvvl)

sgbhhd, sgbvhd, sgbvvd, sgbhdr, sgbvdr, sgbvh1, sgbvh3 = asal(thetir, ntypel, parml, leaf_common, data_common, integ_common)

println("asal")
println(sgbhhd)
println(sgbvhd)
println(sgbvvd)
println(sgbhdr)
println(sgbvdr)
println(sgbvh1)
println(sgbvh3)


println(data_common)

ntypeb1 = 11

afhhb1, afvhb1, afhvb1, afvvb1 = woodf(index,thetir,phiir,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

println("woodf")
println(afhhb1)
println(afvhb1)
println(afhvb1)
println(afvvb1)

afhhb1 = complex(abs(real(afhhb1)),abs(imag(afhhb1)))
afvvb1 = complex(abs(real(afvvb1)),abs(imag(afvvb1)))

println("two complex statements")
println(afhhb1)
println(afvvb1)

sbhhdb1,sbvhdb1,sbvvdb1,sbhdrb1,sbvdrb1,sbvh1b1,sbvh3b1 = woodb(index,thetir,phiir,ntypeb1,parmb1, data_common, integ_common, parm1_common, parm2_common)

println("woodb")
println(sbhhdb1)
println(sbvhdb1)
println(sbvvdb1)
println(sbhdrb1)
println(sbvdrb1)
println(sbvh1b1)
println(sbvh3b1)

# ak0  - free space wave number.
# d - slab thickness
# bfr - frequency in hz
# eps1-relative dielctric constant of scatterer.sth
# eps2-relative dielectric constant of water droplets
# epsg- relative dielectic constant of ground
# rho- density of scatterers- no./m**3.
# kh - horizontal propagation constant
# kmv - vertical propagation constant
# theti- angle of incidence-degrees.
# atmh - attenuation of mean horizontal wave
# atmv - attenuation of mean vertical   wave
# skdh,skdv     - horizontal and vertical skin depth in m.
# sghh ,sgvh ,sgvv - average backscattering coefficients.
# sghho ,sghvo ,sgvvo - average backscattering in lb.