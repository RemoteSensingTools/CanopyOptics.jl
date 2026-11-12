"""
    prospect(leaf, optics[, surface]; n_surface = nothing, μ_in = 1) -> (T, R)

Compute leaf spectral transmittance `T` and reflectance `R` from a
[`LeafProspectProProperties`](@ref) or [`LeafMultiPigmentProperties`](@ref) state
and an [`AbstractProspectProperties`](@ref) optical database.

Optionally supply an [`AbstractLeafSurface`](@ref) (`ConeSurface` or
`RoughDielectricSurface`) and cuticular refractive-index vector `n_surface`.
See [`prospect_budget`](@ref) for the full conserved `(; transmission, reflectance, absorptance)` budget.

# Examples
```jldoctest
julia> using Unitful

julia> opti = LeafOpticalProperties((400.0:500.0:2400.0)u"nm");

julia> leaf = LeafProspectProProperties(Ccab = 30.0);

julia> T, R = prospect(leaf, opti);

julia> length(T) == length(R) == 4 && all(0 .< T .+ R .<= 1)
true
```
"""
function prospect(leaf::LeafProspectProProperties, optis)
    budget = prospect_budget(leaf, optis)
    return budget.transmission, budget.reflectance
end

"""
    prospect_budget(leaf, optis)

PROSPECT-PRO transmission, reflection and directly integrated material absorption.
Returns a named tuple of spectral arrays. Unlike `1-R-T`, the local absorptance
remains meaningful in the nonabsorbing/weak-absorption limit. No flux is clipped
or rescaled. `prospect` retains its `(T, R)` return interface.
Absorptance is the monochromatic fraction of incident energy (equivalently photon
count) absorbed by all material, not chlorophyll-only absorption or energy retained
after fluorescence emission.

Requires finite `N >= 1`, finite nonnegative constituent contents/absorption coefficients
and refractive indices `> 1`. Both `prospect` and `prospect_budget` enforce `N >= 1`. The
PROSPECT first-surface 40-degree and internal diffuse 90-degree conventions and
Fresnel routine are used. Its accuracy is tested across the optical database's
index range; the nearly index-matched limit is ill-conditioned and is not covered
by that accuracy claim. Invalid computed surface transmissions are rejected.
CPU and first-order ForwardDiff paths are
supported; this function does not claim CUDA execution. The physical absorption
Hessian is singular at exactly zero material optical depth.
"""
function prospect_budget(
    leaf::LeafProspectProProperties{FT},
    optis;
    n_surface = nothing,
) where {FT}
    if n_surface !== nothing
        BFT = float(promote_type(_leaf_base_float(FT), _leaf_base_float(eltype(optis.nᵣ))))
        b = prospect_budget(leaf, optis, ConeSurface{BFT}(BFT(40), one(BFT)); n_surface = n_surface)
        return (;
            transmission = b.transmission,
            reflectance = b.reflectance,
            absorptance = b.absorptance,
        )
    end
    (; N, Ccab, Ccar, Cbrown, Canth, Cw, Cm, Cprot, Ccbc)=leaf
    (; Kcab, Kant, Kb, Kcar, Km, Kw, nᵣ, Kp, Kcbc)=optis
    n=prospect_primal(N)
    isfinite(n) && n>=one(n) || throw(ArgumentError("PROSPECT requires finite N>=1"))
    nonnegative(x) =
        isfinite(prospect_primal(x)) && prospect_primal(x)>=zero(prospect_primal(x))
    all(nonnegative, (Ccab, Ccar, Canth, Cbrown, Cw, Cm, Cprot, Ccbc)) ||
        throw(ArgumentError("PROSPECT constituent contents must be finite and nonnegative"))
    all(
        coefficients->all(nonnegative, coefficients),
        (Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc),
    ) || throw(
        ArgumentError("PROSPECT absorption coefficients must be finite and nonnegative"),
    )
    all(
        x->isfinite(prospect_primal(x)) && prospect_primal(x)>one(prospect_primal(x)),
        nᵣ,
    ) || throw(ArgumentError("PROSPECT requires finite refractive indices >1"))
    products=prospect_constituent_coefficients(leaf,optis)
    # Preserve the original operation order while sharing exactly the same
    # concentration × coefficient products with pigment attribution.
    Kall=(products.chlorophyll+products.carotenoids+products.anthocyanins+
        products.brown_material+products.water+products.dry_matter+
        products.protein+products.carbon_constituents)/N
    talf=calctav.(FT(40), nᵣ)
    t12=calctav.(FT(90), nᵣ)
    surface_transmission(x) =
        isfinite(prospect_primal(x)) &&
        zero(prospect_primal(x))<prospect_primal(x)<=one(prospect_primal(x))
    all(surface_transmission, talf) && all(surface_transmission, t12) || throw(
        ArgumentError(
            "Legacy PROSPECT Fresnel calculation produced an invalid surface transmission; nearly matched refractive indices are unsupported",
        ),
    )
    values=prospect_spectral_budget.(Kall, N, talf, t12, nᵣ)
    return (;
        transmission = getproperty.(values, :transmission),
        reflectance = getproperty.(values, :reflectance),
        absorptance = getproperty.(values, :absorptance),
    )
end

"""
    prospect_constituent_coefficients(leaf, optics)

Host-side PROSPECT-PRO composition products for each spectral bin, named by
physical constituent. These are the exact concentration × specific-absorption
products used by `prospect_budget`, before dividing their sum by the common plate
count N. That common scale cancels in homogeneous constituent attribution.

No spectrum, plant profile or units are selected here. Contents and coefficients
must use the matching units documented by `LeafProspectProProperties` and
`LeafOpticalProperties`. Bulk dry matter overlaps protein/carbon constituents;
set Cm=0 when those PRO constituents already describe the same mass. Validation
occurs in `prospect_budget` and the scalar `partition_absorption` kernel. This
array assembly is separate from the CPU/GPU scalar partition/selection rules.
"""
function prospect_constituent_coefficients(leaf::LeafProspectProProperties,optics)
    return (;chlorophyll=leaf.Ccab*optics.Kcab,
        carotenoids=leaf.Ccar*optics.Kcar,
        anthocyanins=leaf.Canth*optics.Kant,
        brown_material=leaf.Cbrown*optics.Kb,
        water=leaf.Cw*optics.Kw,
        dry_matter=leaf.Cm*optics.Km,
        protein=leaf.Cprot*optics.Kp,
        carbon_constituents=leaf.Ccbc*optics.Kcbc)
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
