"""
    $(FUNCTIONNAME)(λ_bnds)

    Loads in the PROSPECT-PRO database of pigments (and other) absorption cross section in leaves, returns a [`LeafOpticalProperties`](@ref) type struct with spectral units attached.
# Arguments
    - `λ_bnds` an array (with or without a spectral grid unit) that defines the upper and lower limits over which to average the absorption cross sections

# Examples
```julia-repl
julia> using Unitful                                                               # Need to include Unitful package 
julia> opti = createLeafOpticalStruct((400.0:5:2400)*u"nm");                       # in nm
julia> opti = createLeafOpticalStruct((0.4:0.1:2.4)*u"μm");                        # in μm
julia> opti = CanopyOptics.createLeafOpticalStruct((10000.0:100:25000.0)u"1/cm");  # in wavenumber (cm⁻¹)
```
"""
function createLeafOpticalStruct(λ_bnds)
    FT = eltype(ustrip(λ_bnds[1]))
    # Reference input grid converted to nm:
    λ_ref = unit2nm(λ_bnds)
    
    KS = readdlm(OPTI_2021, '\t',FT)
    N  = length(λ_bnds)-1
    λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc = [zeros(FT,N) for i = 1:N];
    vars = (λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)
    λ     = KS[:,1]*u"nm";
    @inbounds for i=1:N
        start = min(λ_ref[i],λ_ref[i+1])
        stop  = max(λ_ref[i],λ_ref[i+1])
        ind_all = findall((λ .≥ start) .& (λ .< stop) )
        isempty(ind_all) ? (@warn "Some λ_bnds bins are empty, better coarsen the input grid $(λ_bnds[i]) -> $(λ_bnds[i+1])" ) : nothing
        for j in eachindex(vars)
            vars[j][i] =  mean(view(KS,ind_all,j)); 
        end
    end
    # Get output unit:
    ou = unit(λ_bnds[1])
    (typeof(ou) <: Unitful.FreeUnits{(), NoDims, nothing}) ? out_unit = u"nm"  : out_unit = ou
    return LeafOpticalProperties(uconvSpectral(λ_out*u"nm",out_unit), 
                                nᵣ, 
                                Kcab,
                                Kcar,
                                Kant,
                                Kb,
                                Kw,
                                Km,
                                Kp,
                                Kcbc)
end

"""
    $(FUNCTIONNAME)()
    As in $(FUNCTIONNAME)(λ_bnds) but reads in the in Prospect-PRO database at original resolution (400-2500nm in 1nm steps)
"""
function createLeafOpticalStruct()
    λ_bnds = (399.5:2500.5)u"nm"
    createLeafOpticalStruct(λ_bnds)
end

"Convert an input array with Spectral units (or none) to a grid in nm"
function unit2nm(λ_bnds)
    try
        return uconvSpectral(λ_bnds, u"nm");
    catch e # Assume without unit is in "nm" already
        @warn "Assuming unit of [nm] here for optical struct"
        return λ_bnds*u"nm"
    end 
end

"Avoiding annoying issues with uconvert for Spectra() (can't do nm->nm!)"
# From unit of `in` to `out_unit` (just attaching out_unit if unit is missing)
function uconvSpectral(in,out_unit)
    unit(in[1]) == out_unit ? in : uconvert.(out_unit, in, Spectral())
end

