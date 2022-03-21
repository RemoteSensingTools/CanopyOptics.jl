
function parameters_from_yaml(filepath)

    input_params = YAML.load_file(filepath)

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

    return Forest_Scattering_Parameters(bfrghz, amajcm, bmincm, tmm, rhol, epsl_temp, 
                                        ntypel, parml, radb1, lb1, rhob1, epsb1, ntypeb1, 
                                        parmb1, radb2, lb2, rhob2, epsb2, ntypeb2, parmb2, 
                                        radt, lt_temp, rhot, epst_temp, ntypet, parmt, 
                                        d1, d2, epsg, l, sig)

end

function Base.show(io::IO, x::Forest_Scattering_Parameters)


    println(io, "Frequency: $(x.bfrghz)")

    println(io, "------------------------------")
    println(io, "Leaf Parameters")
    println(io, "------------------------------")
    println(io, "amajcm: $(x.amajcm)")
    println(io, "bmincm: $(x.bmincm)")
    println(io, "tmm: $(x.tmm)")
    println(io, "rhol: $(x.rhol)")
    println(io, "epsl_temp: $(x.epsl_temp)")
    println(io, "ntypel: $(x.ntypel)")
    println(io, "parml: $(x.parml)")

    println(io, "------------------------------")
    println(io, "Primary Branch Parameters")
    println(io, "------------------------------")
    println(io, "radb1: $(x.radb1)")
    println(io, "lb1: $(x.lb1)")
    println(io, "rhob1: $(x.rhob1)")
    println(io, "epsb1: $(x.epsb1)")
    println(io, "ntypeb1: $(x.ntypeb1)")
    println(io, "parmb1: $(x.parmb1)")

    println(io, "------------------------------")
    println(io, "Secondary Branch Parameters")
    println(io, "------------------------------")
    println(io, "radb2: $(x.radb2)")
    println(io, "lb2: $(x.lb2)")
    println(io, "rhob2: $(x.rhob2)")
    println(io, "epsb2: $(x.epsb2)")
    println(io, "ntypeb2: $(x.ntypeb2)")
    println(io, "parmb2: $(x.parmb2)")

    println(io, "------------------------------")
    println(io, "Trunk Parameters")
    println(io, "------------------------------")
    println(io, "radt: $(x.radt)")
    println(io, "lt_temp: $(x.lt_temp)")
    println(io, "rhot: $(x.rhot)")
    println(io, "epst_temp: $(x.epst_temp)")
    println(io, "ntypet: $(x.ntypet)")
    println(io, "parmt: $(x.parmt)")

    println(io, "------------------------------")
    println(io, "Other")
    println(io, "------------------------------")
    println(io, "d1: $(x.d1)")
    println(io, "d2: $(x.d2)")
    println(io, "epsg: $(x.epsg)")
    println(io, "l: $(x.l)")
    println(io, "sig: $(x.sig)")

end