"""
    $(FUNCTIONNAME)(leaf::LeafProspectProProperties{FT},
                    optis) where {FT<:AbstractFloat}

Computes leaf optical properties (reflectance and transittance) based on PROSPECT-PRO
# Arguments
- `leaf`  [`LeafProspectProProperties`](@ref) type struct which provides leaf composition
- `optis` [`LeafOpticalProperties`](@ref) type struct, which provides absorption cross sections and spectral grid
# Examples
```julia-repl
julia> opti = createLeafOpticalStruct((400.0:5:2400)*u"nm");
julia> leaf = LeafProspectProProperties{Float64}(Ccab=30.0);
julia> T,R = prospect(leaf,opti);
```
"""
function prospect(
            leaf::LeafProspectProProperties{FT},
            optis
) where {FT}
    # ***********************************************************************
    # Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
    # properties spectra; Remote Sens. Environ.; 34:75-91.
    # Reference:
    # Féret, Gitelson, Noble & Jacquemoud [2017]. PROSPECT-D: Towards modeling
    # leaf optical properties through a complete lifecycle
    # Remote Sensing of Environment; 193:204215
    # DOI: http://doi.org/10.1016/j.rse.2017.03.004
    # The specific absorption coefficient corresponding to brown pigment is()
    # provided by Frederic Baret [EMMAH, INRA Avignon, baret@avignon.inra.fr]
    # & used with his autorization.
    # ***********************************************************************

    (;N, Ccab, Ccar, Cbrown, Canth, Cw, Cm, Cprot, Ccbc) = leaf;
    (;Kcab, Kant, Kb, Kcar, Km, Kw, nᵣ, Kp, Kcbc)        = optis;

    # This can go into a separate multiple dispatch function as the rest remains constant across versions!
    Kall=(Ccab*Kcab + Ccar*Kcar + Canth*Kant + Cbrown*Kb + Cw*Kw + Cm*Km + Cprot*Kp +Ccbc * Kcbc) / N;
   
    # Adding eps() here to keep it stable and NOT set to 1 manually when Kall=0 (ForwardDiff won't work otherwise)
    tau = (FT(1) .-Kall).*exp.(-Kall) .+ Kall.^2 .*real.(expint.(Kall.+eps(FT)))

    # ***********************************************************************
    # reflectance & transmittance of one layer
    # ***********************************************************************
    # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969)
    # Interaction of isotropic ligth with a compact plant leaf; J. Opt.
    # Soc. Am., 59[10]:1376-1379.
    # ***********************************************************************
    # reflectivity & transmissivity at the interface
    #-------------------------------------------------
    # From Prospect-D, uses 40 here instead of 59 from CVT)

    #talf    = calctav.(59.,nr)
    talf    = calctav.(FT(40),nᵣ)
    ralf    = FT(1) .-talf

    t12     = calctav.(FT(90), nᵣ)
    r12     = FT(1) .-t12
    t21     = t12./(nᵣ.^2)
    r21     = FT(1) .-t21

    # top surface side
    denom   = FT(1) .-r21.*r21.*tau.^2
    Ta      = talf.*tau.*t21./denom
    Ra      = ralf.+r21.*tau.*Ta
    # bottom surface side
    t       = t12.*tau.*t21./denom
    r       = r12+r21.*tau.*t

    # ***********************************************************************
    # reflectance & transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers [N real]
    # Normal case()
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # | transmitted through a pile of plates; Proc. Roy. Soc. Lond.
    # 11:545-556.
    # ***********************************************************************
    D       = sqrt.(((FT(1) .+r.+t).*(FT(1) .+r.-t).*(FT(1) .-r.+t).*(FT(1) .-r.-t)).+5eps(FT))
    #println(typeof(D), typeof(r), typeof(t))
    rq      = r.^2
    tq      = t.^2
    a       = (FT(1) .+rq.-tq.+D)./(2r)
    b       = (FT(1) .-rq.+tq.+D)./(2t)

    bNm1    = b.^(N-1);                  #
    bN2     = bNm1.^2
    a2      = a.^2
    denom   = a2.*bN2.-1
    Rsub    = a.*(bN2.-1)./denom
    Tsub    = bNm1.*(a2.-1)./denom

    # Case of zero absorption
    j       = findall(r.+t .>= 1)
    Tsub[j] = t[j]./(t[j]+(1 .-t[j])*(leaf.N-1))
    Rsub[j] = 1 .-Tsub[j]

    # Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    denom   = FT(1) .-Rsub.*r
    # lambertian Tranmsission
    T    = Ta.*Tsub./denom
    # lambertian Reflectance
    R    = Ra.+Ta.*Rsub.*t./denom
    
    return T,R
end



"""
    calctav(α::FT, nr::FT) where {FT<:AbstractFloat}

Computes transmission of isotropic radiation across a dielectric surface 
(Stern F., 1964; Allen W.A.,Appl. Opt., 3(1):111-113 1973)). 
From calctav.m in PROSPECT-D
# Arguments
- `α` angle of incidence [degrees]
- `nr` Index of refraction
"""
function calctav(α::FT, nr::FT2) where {FT,FT2}
    a   = ((nr+1) ^ 2) / 2;
    a3  = a  ^ 3;
    n2  = nr ^ 2;
    n4  = nr ^ 4;
    n6  = nr ^ 6;
    np  = n2 + 1;
    np2 = np ^ 2;
    np3 = np ^ 3;
    nm2 = (n2 - 1) ^2;
    k   = ((n2-1) ^ 2) / -4;
    k2  = k  ^ 2;
    sa2 = sind(α) ^ 2;

    _b1 = (α==90 ? 0 : sqrt( (sa2 - np/2)^2 + k ));
    _b2 = sa2 - np/2;
    b   = _b1 - _b2;
    b3  = b ^ 3;
    ts  = ( k2 / (6*b3) + k/b - b/2 ) - ( k2 / (6*a3) + k/a - a/2 );
    tp1 = -2 * n2 * (b-a) / (np2);
    tp2 = -2 * n2 * np * log(b/a) / (nm2);
    tp3 = n2 * (1/b - 1/a) / 2;
    tp4 = 16 * n4 * (n4+1) * log((2*np*b - nm2) / (2*np*a - nm2)) / (np3*nm2);
    tp5 = 16 * n6 * (1 / (2*np*b-nm2) - 1 / (2*np*a-nm2)) / (np3);
    tp  = tp1 + tp2 + tp3 + tp4 + tp5;
    tav = (ts + tp) / (2 * sa2);

    return tav
end
