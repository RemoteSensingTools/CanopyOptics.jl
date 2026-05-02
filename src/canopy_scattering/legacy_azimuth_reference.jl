"""
Legacy and reference azimuth-integration paths for bi-Lambertian canopy
scattering.

The production bi-Lambertian Z assembly is
[`compute_Z_matrices_aniso_analytic`](@ref). The methods in this file are kept
for compatibility with older callers and for cross-checking analytic Fourier
moments against direct azimuth integration.
"""

"""
    compute_Γ_isotropic(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT},
                        Ωᵒᵘᵗ::dirVector_μ{FT}) where FT

Analytic area scattering transfer function assuming an isotropic (spherical) leaf
angle distribution, following Shultis & Myneni (1988) Eq. (40):

``Γ_{\\mathrm{iso}}(β) = \\frac{ω}{3π}(\\sin β - β \\cos β) + \\frac{T}{3}\\cos β``

where `β = arccos(Ωⁱⁿ ⋅ Ωᵒᵘᵗ)` is the scattering angle and `ω = R + T`.

Used for validation against the general anisotropic [`compute_Γ`](@ref).
"""
function compute_Γ_isotropic(mod::BiLambertianCanopyScattering, Ωⁱⁿ::dirVector_μ{FT}, Ωᵒᵘᵗ::dirVector_μ{FT}) where FT
    (;R,T) = mod
    β = acos( Ωᵒᵘᵗ ⋅ Ωⁱⁿ)
    ω = R + T
    
    # Eq 40 in Shultis and Myneni
    Γ = (ω/3π) * (sin(β) - β * cos(β)) + T/3*cos(β)
    return Γ
end

"""
    compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering,
                             μ, LD, Zup, Zdown, m; quadrature)

Compatibility overload for callers that still pass precomputed azimuth grids.

`Zup` and `Zdown` are ignored because the current implementation delegates to
the analytic Fourier path through [`compute_Z_matrices_aniso`](@ref). Keeping
this method avoids breaking older vSmartMOM-side call sites while preventing
the stale precomputed convention from affecting results.
"""
function compute_Z_matrices_aniso(mod::BiLambertianCanopyScattering,
                                  μ::AbstractArray{FT,1},
                                  LD::AbstractLeafDistribution,
                                  Zup, Zdown, m::Int;
                                  quadrature::CanopyQuadrature = CanopyQuadrature(),
                                  nQuad = nothing) where FT
    # Zup/Zdown are retained only for API compatibility with callers that
    # previously used the brute-force precomputed-azimuth path.
    q = _resolve_quadrature(quadrature, nQuad)
    return compute_Z_matrices_aniso(mod, μ, LD, m; quadrature = q)
end


"""
    precompute_Zazi(mod::BiLambertianCanopyScattering, μ, LD; quadrature)

Build direct-azimuth samples of the bi-Lambertian transfer function.

This is a reference/compatibility path for the old
`compute_Z_matrices_aniso(mod, μ, LD, Zup, Zdown, m)` API. New production code
should call [`compute_Z_matrices`](@ref) or
[`compute_Z_matrices_aniso_analytic`](@ref) instead.
"""
function precompute_Zazi(mod::BiLambertianCanopyScattering,
                         μ::AbstractArray{FT,1},
                         LD::AbstractLeafDistribution;
                         quadrature::CanopyQuadrature = CanopyQuadrature(),
                         nQuad = nothing) where FT
    q = _resolve_quadrature(quadrature, nQuad)
    n_azimuth = q.n_azimuth
    # Quadrature points in the azimuth:
    ϕ, w_azi = gauleg(n_azimuth,FT(0),FT(π));
    # Fourier weights (cosine decomposition)
    
    # Transmission (same direction)
    Zup = zeros(length(μ),length(μ), n_azimuth)
    # Reflection (change direction)
    Zdown = zeros(length(μ),length(μ), n_azimuth)
    
    # Create outgoing vectors in θ and ϕ
    dirOutꜛ = [dirVector_μ(a,b) for a in -μ, b in ϕ];
    dirOutꜜ = [dirVector_μ(a,b) for a in μ, b in ϕ];

    Threads.@threads for i in eachindex(μ)
        # Incoming beam at ϕ = 0
        Ωⁱⁿ = dirVector_μ(μ[i], FT(0));
        # Compute over μ and μ_azi:
        Zup[i,:,:]   = [compute_Γ(mod, Ωⁱⁿ, dirOutꜛ[j, k], LD; quadrature = q)
                        for j in axes(dirOutꜛ, 1), k in axes(dirOutꜛ, 2)]
        Zdown[i,:,:] = [compute_Γ(mod, Ωⁱⁿ, dirOutꜜ[j, k], LD; quadrature = q)
                        for j in axes(dirOutꜜ, 1), k in axes(dirOutꜜ, 2)]
        #Zup[:,:,i]   = compute_Γ_isotropic.((mod,),(Ωⁱⁿ,),dirOutꜛ);
        #Zdown[:,:,i] = compute_Γ_isotropic.((mod,),(Ωⁱⁿ,),dirOutꜜ);
    end
    return Zup, Zdown
end

"""
    precompute_Zazi_(mod::BiLambertianCanopyScattering, μ, LD; quadrature)

Vectorized experimental version of [`precompute_Zazi`](@ref).

This helper is retained only as an implementation reference for the direct
azimuth path. It is not used by the canonical Z API.
"""
function precompute_Zazi_(mod::BiLambertianCanopyScattering,
                          μ::AbstractArray{FT,1},
                          LD::AbstractLeafDistribution;
                          quadrature::CanopyQuadrature = CanopyQuadrature(),
                          nQuad = nothing) where FT
    (;R,T) = mod
    q = _resolve_quadrature(quadrature, nQuad)
    n_leaf = q.n_leaf
    n_azimuth = q.n_azimuth
    # Quadrature points in μ
    n_μ  = length(μ);
    arr_type = typeof(μ);

    #μ,w         = CanopyOptics.gauleg(n_μ,   FT(0),  FT(1.0));
    dϕ,  w_azi  = CanopyOptics.gauleg(n_azimuth, FT(0),  FT(π));
    dϕᴸ, w_aziᴸ = CanopyOptics.gauleg(n_azimuth+1, FT(0),  FT(2π));
    θᴸ,wᴸ       = CanopyOptics.gauleg(n_leaf, FT(0),FT(π/2));
    
    μᴸ = cos.(θᴸ)
    Fᵢ = pdf.(LD.LD,2θᴸ/π) * LD.scaling
    # Reshape stuff:
    μⁱⁿ   = reshape(arr_type(μ),  n_μ,  1,     1,     1,     1   );
    μᵒᵘᵗ  = reshape(arr_type(deepcopy(μ)), 1,   n_μ,    1,     1,     1   );
    _dϕ   = reshape(arr_type(dϕ), 1,    1,   n_azimuth,   1,     1   );
    _μᴸ   = reshape(arr_type(μᴸ), 1,    1,     1,   n_leaf,   1   );
    _dϕᴸ  = reshape(arr_type(dϕᴸ),1,    1,     1,     1,   n_azimuth+1 );

    # Quadrature points
    wᴸ  = wᴸ .*  Fᵢ
    #_w_azi  = reshape(arr_type(w_azi),  1,  1,   n_azimuth,   1,     1   );
    _wᴸ     = reshape(arr_type(wᴸ),     1,  1,    1,    n_leaf,   1   );
    _w_aziᴸ = reshape(arr_type(w_aziᴸ), 1,  1,    1,      1,    n_azimuth+1);

    integrand  =  CanopyOptics.leaf_dot_products.(μⁱⁿ, -μᵒᵘᵗ, _dϕ,_μᴸ, _dϕᴸ);
    iPos       = (integrand+abs.(integrand))./2;
    iNeg       = (integrand-abs.(integrand))./2;
    Γ⁻         = -1/2π * sum(sum(iNeg.*_w_aziᴸ, dims=5).*_wᴸ,dims=4); 
    Γ⁺         =  1/2π * sum(sum(iPos.*_w_aziᴸ, dims=5).*_wᴸ,dims=4);

    Γ = R .* Γ⁻ .+ T .* Γ⁺; 
    Γup = reshape(Γ, n_μ,n_μ, n_azimuth);
    
    integrand  =  CanopyOptics.leaf_dot_products.(μⁱⁿ, μᵒᵘᵗ, _dϕ,_μᴸ, _dϕᴸ);
    iPos       = (integrand+abs.(integrand))./2;
    iNeg       = (integrand-abs.(integrand))./2;
    Γ⁻         = -1/2π * sum(sum(iNeg.*_w_aziᴸ, dims=5).*_wᴸ,dims=4); 
    Γ⁺         =  1/2π * sum(sum(iPos.*_w_aziᴸ, dims=5).*_wᴸ,dims=4);
    
    Γ = R .* Γ⁻ .+ T .* Γ⁺; 
    Γdown = reshape(Γ, n_μ,n_μ, n_azimuth);
    return Γup,Γdown
end

"""
    leaf_dot_products(μⁱⁿ, μᵒᵘᵗ, dϕᵒᵘᵗ, μᴸ, dϕᴸ)

Return `(Ω_in ⋅ Ω_leaf) * (Ω_out ⋅ Ω_leaf)` from cosine/azimuth variables.

Used by the vectorized legacy azimuth precompute path.
"""
function leaf_dot_products(μⁱⁿ::FT, μᵒᵘᵗ::FT, dϕᵒᵘᵗ::FT, μᴸ::FT, dϕᴸ::FT) where FT
    (μⁱⁿ  * μᴸ + sqrt(1-μⁱⁿ^2)   * sqrt(1-μᴸ^2) * cos(dϕᴸ)) *
    (μᵒᵘᵗ * μᴸ + sqrt(1-μᵒᵘᵗ^2)  * sqrt(1-μᴸ^2) * cos(dϕᴸ - dϕᵒᵘᵗ))
end
