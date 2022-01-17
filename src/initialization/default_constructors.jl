"""
    $(FUNCTIONNAME)(FT=Float64)

Creates a `LeafDistribution` following a planophile (mostly horizontal) distribution using a Beta distribution

Standard values for θₗ (26.76) and s (18.51) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1
"""
function planophile_leaves(FT=Float64)
    θ = FT(26.76)
    s = FT(18.51)
    α, β =  βparameters(deg2rad(θ), deg2rad(s));
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2/π));
end

function planophile_leaves2(FT=Float64)
    return LeafDistribution(Beta(FT(1), FT(2.768)), FT(2/π));
end

"""
    $(FUNCTIONNAME)(FT=Float64)

Creates a `LeafDistribution` following a uniform distribution (all angles equally likely)
"""
function uniform_leaves(FT=Float64)
    return LeafDistribution(Uniform(FT(0), FT(1)),FT(2/π));
end

"""
    $(FUNCTIONNAME)(FT=Float64)

Creates a `LeafDistribution` following a plagiophile (mostly 45degrees) distribution using a Beta distribution

Standard values for θₗ (45.0) and s (16.27) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1
"""
function plagiophile_leaves(FT=Float64)
    θ = FT(45)
    s = FT(16.27)
    α, β =  βparameters(deg2rad(θ), deg2rad(s));
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2/π));
end

"""
    $(FUNCTIONNAME)(FT=Float64)

Creates a `LeafDistribution` following a erectophile (mostly vertical) distribution using a Beta distribution

Standard values for θₗ (63.24) and s (18.5) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1
"""
function erectophile_leaves(FT=Float64)
    θ = FT(63.24)
    s = FT(18.5)
    α, β =  βparameters(deg2rad(θ), deg2rad(s));
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2/π));
end

"""
    $(FUNCTIONNAME)(FT=Float64)

Creates a `LeafDistribution` following a spherical (distributed as facets on a sphere) distribution using a Beta distribution

Standard values for θₗ (57.3) and s (21.55) from Bonan https://doi.org/10.1017/9781107339217.003, Table 2.1
"""
function spherical_leaves(FT=Float64)
    θ = FT(57.3)
    s = FT(21.55)
    α, β =  βparameters(deg2rad(θ), deg2rad(s));
    return LeafDistribution(Beta(FT(α), FT(β)), FT(2/π));
end
