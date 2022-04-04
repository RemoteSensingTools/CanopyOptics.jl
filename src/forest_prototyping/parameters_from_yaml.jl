
function parameters_from_yaml(filepath)

    input_params = YAML.load_file(filepath)

    # Frequency
    bfrghz  = input_params["frequency"]

    # Leaf Parameters 
    amajcm  = input_params["leaf_parameters"]["major_axis"]
    bmincm  = input_params["leaf_parameters"]["minor_axis"]
    tmm     = input_params["leaf_parameters"]["thickness"]
    ρ_l     = input_params["leaf_parameters"]["density"]
    ϵ_l     = input_params["leaf_parameters"]["dielectric_constant"]
    ntypel  = input_params["leaf_parameters"]["leaf_inclination_pdf_type"]
    parml   = input_params["leaf_parameters"]["pdf_parameter"]

    # Primary Branch Parameters
    radb1   = input_params["primary_branch_parameters"]["diameter"]
    lb1     = input_params["primary_branch_parameters"]["length"]
    ρ_b1    = input_params["primary_branch_parameters"]["density"]
    ϵ_b1    = input_params["primary_branch_parameters"]["dielectric_constant"]
    ntypeb1 = input_params["primary_branch_parameters"]["branch_inclination_pdf_type"]
    parmb1  = input_params["primary_branch_parameters"]["pdf_parameter"]

    # Secondary Branch Parameters
    radb2   = input_params["secondary_branch_parameters"]["diameter"]
    lb2     = input_params["secondary_branch_parameters"]["length"]
    ρ_b2    = input_params["secondary_branch_parameters"]["density"]
    ϵ_b2    = input_params["secondary_branch_parameters"]["dielectric_constant"]
    ntypeb2 = input_params["secondary_branch_parameters"]["branch_inclination_pdf_type"]
    parmb2  = input_params["secondary_branch_parameters"]["pdf_parameter"]

    # Trunk Parameters
    radt    = input_params["trunk_parameters"]["diameter"]
    lt_temp      = input_params["trunk_parameters"]["length"]
    ρ_t     = input_params["trunk_parameters"]["density"]
    ϵ_t     = input_params["trunk_parameters"]["dielectric_constant"]
    ntypet  = input_params["trunk_parameters"]["branch_inclination_pdf_type"]
    parmt   = input_params["trunk_parameters"]["pdf_parameter"]
    d_c      = input_params["trunk_parameters"]["crown_height"]
    d_t      = input_params["trunk_parameters"]["trunk_height"]
    ϵ_g     = input_params["trunk_parameters"]["soil_dielectric"] 
    l       = input_params["trunk_parameters"]["corr_length"]
    sig     = input_params["trunk_parameters"]["rms_height"]

    # Convert complex numbers 
    ϵ_l  = eval(Meta.parse(ϵ_l))
    ϵ_b1 = eval(Meta.parse(ϵ_b1))
    ϵ_b2 = eval(Meta.parse(ϵ_b2))
    ϵ_t  = eval(Meta.parse(ϵ_t))
    ϵ_g  = eval(Meta.parse(ϵ_g))

    return Forest_Scattering_Parameters(bfrghz, amajcm, bmincm, tmm, ρ_l, ϵ_l, 
                                        ntypel, parml, radb1, lb1, ρ_b1, ϵ_b1, ntypeb1, 
                                        parmb1, radb2, lb2, ρ_b2, ϵ_b2, ntypeb2, parmb2, 
                                        radt, lt_temp, ρ_t, ϵ_t, ntypet, parmt, 
                                        d_c, d_t, ϵ_g, l, sig)

end

function Base.show(io::IO, x::Forest_Scattering_Parameters)


    println(io, "Frequency: $(x.bfrghz)")

    println(io, "------------------------------")
    println(io, "Leaf Parameters")
    println(io, "------------------------------")
    println(io, "amajcm: $(x.amajcm)")
    println(io, "bmincm: $(x.bmincm)")
    println(io, "tmm: $(x.tmm)")
    println(io, "ρ_l: $(x.ρ_l)")
    println(io, "ϵ_l: $(x.ϵ_l)")
    println(io, "ntypel: $(x.ntypel)")
    println(io, "parml: $(x.parml)")

    println(io, "------------------------------")
    println(io, "Primary Branch Parameters")
    println(io, "------------------------------")
    println(io, "radb1: $(x.radb1)")
    println(io, "lb1: $(x.lb1)")
    println(io, "ρ_b1: $(x.ρ_b1)")
    println(io, "ϵ_b1: $(x.ϵ_b1)")
    println(io, "ntypeb1: $(x.ntypeb1)")
    println(io, "parmb1: $(x.parmb1)")

    println(io, "------------------------------")
    println(io, "Secondary Branch Parameters")
    println(io, "------------------------------")
    println(io, "radb2: $(x.radb2)")
    println(io, "lb2: $(x.lb2)")
    println(io, "ρ_b2: $(x.ρ_b2)")
    println(io, "ϵ_b2: $(x.ϵ_b2)")
    println(io, "ntypeb2: $(x.ntypeb2)")
    println(io, "parmb2: $(x.parmb2)")

    println(io, "------------------------------")
    println(io, "Trunk Parameters")
    println(io, "------------------------------")
    println(io, "radt: $(x.radt)")
    println(io, "lt_temp: $(x.lt_temp)")
    println(io, "ρ_t: $(x.ρ_t)")
    println(io, "ϵ_t: $(x.ϵ_t)")
    println(io, "ntypet: $(x.ntypet)")
    println(io, "parmt: $(x.parmt)")

    println(io, "------------------------------")
    println(io, "Other")
    println(io, "------------------------------")
    println(io, "d_c: $(x.d_c)")
    println(io, "d_t: $(x.d_t)")
    println(io, "ϵ_g: $(x.ϵ_g)")
    println(io, "l: $(x.l)")
    println(io, "sig: $(x.sig)")

end