

function parameters_from_yaml(filepath)

    input_params  = YAML.load_file(filepath)

    # Frequency
    bfrghz        = input_params["frequency"]

    # Leaf Parameters and Conversions to Metric
    a_maj_cm      = input_params["leaf_parameters"]["major_axis"]
    b_min_cm      = input_params["leaf_parameters"]["minor_axis"]
    tmm           = input_params["leaf_parameters"]["thickness"]
    ρ_l           = input_params["leaf_parameters"]["density"]
    ϵ_l           = eval(Meta.parse(input_params["leaf_parameters"]["dielectric_constant"]))
    pdf_num_l     = input_params["leaf_parameters"]["leaf_inclination_pdf_type"]
    pdf_param_l   = input_params["leaf_parameters"]["pdf_parameter"]
    a_maj         = a_maj_cm * 1e-2    # cm to m
    b_min         = b_min_cm * 1e-2    # cm to m
    t             = tmm * 1e-3         # mm to m

    leaf          = Leaf(ϵ_l, a_maj, b_min, t, ρ_l, pdf_num_l, pdf_param_l)

    # Primary Branch Parameters
    r_b1_cm       = input_params["primary_branch_parameters"]["diameter"]
    l_b1          = input_params["primary_branch_parameters"]["length"]
    ρ_b1          = input_params["primary_branch_parameters"]["density"]
    ϵ_b1          = eval(Meta.parse(input_params["primary_branch_parameters"]["dielectric_constant"]))
    pdf_num_b1    = input_params["primary_branch_parameters"]["branch_inclination_pdf_type"]
    pdf_param_b1  = input_params["primary_branch_parameters"]["pdf_parameter"]
    r_b1          = 0.5*r_b1_cm*1e-2  # diameter to radius, cm to m

    branch_1      = Wood(ϵ_b1, r_b1, l_b1, ρ_b1, pdf_num_b1, pdf_param_b1)

    # Secondary Branch Parameters
    r_b2_cm       = input_params["secondary_branch_parameters"]["diameter"]
    l_b2          = input_params["secondary_branch_parameters"]["length"]
    ρ_b2          = input_params["secondary_branch_parameters"]["density"]
    ϵ_b2          = eval(Meta.parse(input_params["secondary_branch_parameters"]["dielectric_constant"]))
    pdf_num_b2    = input_params["secondary_branch_parameters"]["branch_inclination_pdf_type"]
    pdf_param_b2  = input_params["secondary_branch_parameters"]["pdf_parameter"]
    r_b2          = 0.5*r_b2_cm*1e-2  # diameter to radius, cm to m

    branch_2      = Wood(ϵ_b2, r_b2, l_b2, ρ_b2, pdf_num_b2, pdf_param_b2)

    # Trunk Parameters
    r_t_cm        = input_params["trunk_parameters"]["diameter"]
    l_t           = input_params["trunk_parameters"]["length"]
    ρ_t           = input_params["trunk_parameters"]["density"]
    ϵ_t           = eval(Meta.parse(input_params["trunk_parameters"]["dielectric_constant"]))
    pdf_num_t     = input_params["trunk_parameters"]["branch_inclination_pdf_type"]
    pdf_param_t   = input_params["trunk_parameters"]["pdf_parameter"]
    r_t           = 0.5*r_t_cm*1e-2  # diameter to radius, cm to m

    trunk         = Wood(ϵ_t, r_t, l_t, ρ_t, pdf_num_t, pdf_param_t)

    d_c      = input_params["trunk_parameters"]["crown_height"]
    d_t      = input_params["trunk_parameters"]["trunk_height"]
    ϵ_g     = eval(Meta.parse(input_params["trunk_parameters"]["soil_dielectric"] ))
    l       = input_params["trunk_parameters"]["corr_length"]
    sig     = input_params["trunk_parameters"]["rms_height"]

    return Forest_Scattering_Parameters(bfrghz, leaf, branch_1, branch_2, trunk,
                                        d_c, d_t, ϵ_g, l, sig)

end

function Base.show(io::IO, x::Forest_Scattering_Parameters)


    println(io, "Frequency: $(x.bfrghz)")

    println(io, "------------------------------")
    println(io, "Leaf Parameters")
    println(io, "------------------------------")
    println(io, "a_maj: $(x.leaf.a_maj * 1e2) cm")
    println(io, "b_min: $(x.leaf.b_min * 1e2) cm")
    println(io, "t: $(x.leaf.t * 1e3) mm")
    println(io, "ρ: $(x.leaf.ρ) m3")
    println(io, "ϵ: $(x.leaf.ϵ)")
    println(io, "pdf_num: $(x.leaf.pdf_num)")
    println(io, "pdf_param: $(x.leaf.pdf_param)")

    println(io, "------------------------------")
    println(io, "Primary Branch Parameters")
    println(io, "------------------------------")
    println(io, "r: $(x.branch_1.r  * 1e2) cm")
    println(io, "l: $(x.branch_1.l) m")
    println(io, "ρ: $(x.branch_1.ρ) m3")
    println(io, "ϵ: $(x.branch_1.ϵ)")
    println(io, "pdf_num: $(x.branch_1.pdf_num)")
    println(io, "pdf_param: $(x.branch_1.pdf_param)")

    println(io, "------------------------------")
    println(io, "Secondary Branch Parameters")
    println(io, "------------------------------")
    println(io, "r_b2: $(x.branch_2.r * 1e2) cm")
    println(io, "l_b2: $(x.branch_2.l) m")
    println(io, "ρ_b2: $(x.branch_2.ρ) m3")
    println(io, "ϵ_b2: $(x.branch_2.ϵ)")
    println(io, "pdf_num: $(x.branch_2.pdf_num)")
    println(io, "pdf_param: $(x.branch_2.pdf_param)")

    println(io, "------------------------------")
    println(io, "Trunk Parameters")
    println(io, "------------------------------")
    println(io, "r_t: $(x.trunk.r * 1e2) cm")
    println(io, "l_t: $(x.trunk.l) m")
    println(io, "ρ_t: $(x.trunk.ρ) m3")
    println(io, "ϵ_t: $(x.trunk.ϵ)")
    println(io, "pdf_num: $(x.trunk.pdf_num)")
    println(io, "pdf_param: $(x.trunk.pdf_param)")

    println(io, "------------------------------")
    println(io, "Other")
    println(io, "------------------------------")
    println(io, "d_c: $(x.d_c)")
    println(io, "d_t: $(x.d_t)")
    println(io, "ϵ_g: $(x.ϵ_g)")
    println(io, "l: $(x.l)")
    println(io, "sig: $(x.sig)")

end