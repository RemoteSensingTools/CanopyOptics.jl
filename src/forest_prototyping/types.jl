mutable struct Forest_Scattering_Parameters 

    # Frequency
    bfrghz

    # Leaf Parameters 
    amajcm
    bmincm
    tmm
    rhol
    epsl_temp
    ntypel
    parml

    # Primary Branch Parameters
    radb1   
    lb1     
    rhob1   
    epsb1   
    ntypeb1 
    parmb1  

    # Secondary Branch Parameters
    radb2   
    lb2     
    rhob2   
    epsb2   
    ntypeb2 
    parmb2  

    # Trunk Parameters
    radt    
    lt_temp 
    rhot    
    epst_temp
    ntypet  
    parmt   
    d1      
    d2      
    epsg    
    l       
    sig     

end
mutable struct a 

    khim1
    khim2
    kvim1
    kvim2
    bfr
    epsg

end

mutable struct b

    zk
    sig
    zlx
    zly

end

mutable struct ds

    d1
    d2

end

mutable struct data

    ak0
    epsb
    epst
    radbm
    radtm
    lb
    lt

end

mutable struct integ 

    nph
    nth

end

mutable struct leaf 

    epsl
    amaj
    bmin
    t

end


mutable struct parm1

    r0
    h
    ej
    epsi

end

mutable struct parm2

    thetai
    phyi
    thetas
    phys

end

mutable struct Forest_Scattering_Output

    theti

    sighhd1
    sighhd2
    sighhd3

    sighhdr1
    sighhdr2
    sighhdr3

    sighht
    sighhd
    sighhdr

    sigvvd1
    sigvvd2
    sigvvd3

    sigvvdr1
    sigvvdr2
    sigvvdr3

    sigvvt
    sigvvd
    sigvvdr

    sigvhd1
    sigvhd2
    sigvhd3

    sigvhdr1
    sigvhdr2
    sigvhdr3

    sigvht
    sigvhd
    sigvhdr

    ghhd
    gvhd
    gvvd

    sighhi
    sigvhi
    sigvvi

    sigvhi1
    sigvhi2
    sigvhi3

    cofhh
    cofvh
    cofvv

    athc
    atvc
    atht
    atvt

end