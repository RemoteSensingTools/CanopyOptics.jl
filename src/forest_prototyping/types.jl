"""
All leaf-related parameters
"""
mutable struct Leaf 

    ϵ           # Dielectric constant
    a_maj       # Major axis (m)
    b_min       # Minor axis (m)
    t           # Thickness (m)
    ρ           # Density (m3)

    pdf_num     # PDF number
    pdf_param   # PDF parameter

end

"""
All wood-related parameters (branch or trunk)
"""
mutable struct Wood

    ϵ           # Dielectric constant
    r           # Radius (m)
    l           # Length (m)
    ρ           # Density (m3)

    pdf_num     # PDF number
    pdf_param   # PDF parameter

end

mutable struct parm2

    thetai
    phyi
    thetas
    phys

end

mutable struct Forest_Scattering_Parameters 

    # Frequency
    bfrghz

    # Leaf object
    leaf::Leaf

    # Two branch objects
    branch_1::Wood
    branch_2::Wood

    # Trunk object
    trunk::Wood

    d_c
    d_t  
    ϵ_g    
    l       
    sig     

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