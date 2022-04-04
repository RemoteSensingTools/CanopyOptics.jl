mutable struct Forest_Scattering_Parameters 

    # Frequency
    bfrghz

    # Leaf Parameters 
    amajcm
    bmincm
    tmm
    ρ_l
    ϵ_l
    ntypel
    parml

    # Primary Branch Parameters
    radb1   
    l_b1     
    ρ_b1
    ϵ_b1   
    ntypeb1 
    parmb1  

    # Secondary Branch Parameters
    radb2   
    l_b2     
    ρ_b2   
    ϵ_b2   
    ntypeb2 
    parmb2  

    # Trunk Parameters
    radt    
    l_t
    ρ_t    
    ϵ_t
    ntypet  
    parmt   

    d_c
    d_t  
    ϵ_g    
    l       
    sig     

end
mutable struct a 

    khim1
    khim2
    kvim1
    kvim2
    
    bfr
    ϵ_g

end

mutable struct b

    zk
    sig
    zlx
    zly

end

mutable struct ds

    d_t
    d_c

end

mutable struct leaf 

    ϵ_l
    amaj
    bmin
    t

end

mutable struct branch 

    ϵ_b
    radbm
    l_b

end

mutable struct trunk 

    ϵ_t
    radtm
    l_t

end

mutable struct parm1

    r0
    h
    ej
    ϵ_i

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