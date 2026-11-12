"""
    DorsiventralBiLambertianCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}

Dorsiventral bi-Lambertian canopy leaf scattering model with distinct adaxial
(upper, sky-facing: `R_ad`, `T_ad`) and abaxial (lower, soil-facing: `R_ab`, `T_ab`)
hemispherical reflectance and transmittance.

Note: For true elastic leaf radiative transfer (e.g. [`dorsiventral_budget`](@ref)
with `hemispherical = true`), bi-hemispherical Helmholtz reciprocity guarantees
`T_ad == T_ab` even when `R_ad != R_ab`; unequal `T_ad != T_ab` is permitted here
only for effective/upscaled asymmetric canopy elements (such as clumped shoots).

In `CanopyOptics`'s canopy coordinate convention (`+μ` = downward stream into the
canopy, `-μ` = upward stream toward the sky), a leaf with inclination `θ_L ∈ [0, π/2]`
(`μ_L = cos θ_L >= 0`) has downward abaxial normal `Ω_L = dirVector_μ(+μ_L, ϕ_L)`
and upward adaxial normal `-Ω_L`. Because downward and upward streams illuminate
opposite faces on average (`Z⁺⁺ != Z⁻⁻` and `Z⁻⁺ != Z⁺⁻` whenever `R_ad != R_ab`
or `T_ad != T_ab`), use [`compute_Z_matrices_directional`](@ref) to assemble all
four scattering blocks `(; Zpp, Zmp, Zpm, Zmm, ω_down, ω_up)`.
"""
struct DorsiventralBiLambertianCanopyScattering{FT<:Real} <: AbstractCanopyScatteringType{FT}
    R_ad::FT
    T_ad::FT
    R_ab::FT
    T_ab::FT

    function DorsiventralBiLambertianCanopyScattering{FT}(
        R_ad::FT, T_ad::FT, R_ab::FT, T_ab::FT,
    ) where {FT<:Real}
        BFT = _leaf_base_float(FT)
        for x in (R_ad, T_ad, R_ab, T_ab)
            v = prospect_primal(x)
            isfinite(v) && v >= zero(BFT) ||
                throw(ArgumentError("Dorsiventral leaf R and T must be finite and nonnegative"))
        end
        tol = BFT(16) * eps(BFT)
        (prospect_primal(R_ad + T_ad) <= one(BFT) + tol &&
         prospect_primal(R_ab + T_ab) <= one(BFT) + tol) ||
            throw(ArgumentError("Single-face R + T cannot exceed 1"))
        return new{FT}(R_ad, T_ad, R_ab, T_ab)
    end
end

function DorsiventralBiLambertianCanopyScattering(
    R_ad::Real, T_ad::Real, R_ab::Real, T_ab::Real,
)
    FT = _canopy_parameter_type(R_ad, T_ad, R_ab, T_ab)
    return DorsiventralBiLambertianCanopyScattering{FT}(FT(R_ad), FT(T_ad), FT(R_ab), FT(T_ab))
end

function DorsiventralBiLambertianCanopyScattering(;
    R_ad::Real = 0.30,
    T_ad::Real = 0.10,
    R_ab::Real = 0.35,
    T_ab::Real = 0.10,
)
    return DorsiventralBiLambertianCanopyScattering(R_ad, T_ad, R_ab, T_ab)
end

function DorsiventralBiLambertianCanopyScattering{FT}(;
    R_ad = FT(0.30),
    T_ad = FT(0.10),
    R_ab = FT(0.35),
    T_ab = FT(0.10),
) where {FT<:Real}
    return DorsiventralBiLambertianCanopyScattering{FT}(FT(R_ad), FT(T_ad), FT(R_ab), FT(T_ab))
end

@inline function _is_face_symmetric(mod::DorsiventralBiLambertianCanopyScattering)
    return _same_val_and_partials(mod.R_ad, mod.R_ab) &&
           _same_val_and_partials(mod.T_ad, mod.T_ab)
end

"""
    compute_Γ_directional(mod::DorsiventralBiLambertianCanopyScattering,
                          Ωⁱⁿ::dirVector_μ, Ωᵒᵘᵗ::dirVector_μ, LD; quadrature)

Direct 2D `(θ_L, ϕ_L)` quadrature reference for the dorsiventral area scattering
transfer function `Γ(Ωⁱⁿ → Ωᵒᵘᵗ)`. Reduces identically to `compute_Γ`
when `R_ad == R_ab` and `T_ad == T_ab`.
"""
function compute_Γ_directional(
    mod::DorsiventralBiLambertianCanopyScattering,
    Ωⁱⁿ::dirVector_μ{FT},
    Ωᵒᵘᵗ::dirVector_μ{FT},
    LD::AbstractLeafDistribution;
    quadrature::CanopyQuadrature = CanopyQuadrature(),
) where {FT<:Real}
    (; R_ad, T_ad, R_ab, T_ab) = mod
    ZFT = promote_type(FT, typeof(R_ad))
    leaf_quad = _leaf_inclination_quadrature(LD, quadrature.n_leaf, FT)
    θₗ = leaf_quad.θ
    ϕ, w_azi = gauleg(quadrature.n_azimuth, FT(0), FT(2π))

    acc = zero(ZFT)
    inv_2pi = inv(ZFT(2) * ZFT(π))

    for l in eachindex(θₗ)
        μ_l = FT(cos(θₗ[l]))
        wl = ZFT(leaf_quad.w_measure[l])
        az_sum = zero(ZFT)
        for ia in eachindex(ϕ)
            Ω_l = dirVector_μ(μ_l, ϕ[ia])
            d_in = ZFT(Ωⁱⁿ ⋅ Ω_l)
            d_out = ZFT(Ωᵒᵘᵗ ⋅ Ω_l)

            # Incoming hits adaxial when d_in > 0, abaxial when d_in < 0
            in_ad = d_in > zero(ZFT) ? d_in : zero(ZFT)
            in_ab = d_in < zero(ZFT) ? -d_in : zero(ZFT)
            # Outgoing exits adaxial when d_out < 0, abaxial when d_out > 0
            out_ad = d_out < zero(ZFT) ? -d_out : zero(ZFT)
            out_ab = d_out > zero(ZFT) ? d_out : zero(ZFT)

            kern = ZFT(R_ad) * in_ad * out_ad +
                   ZFT(R_ab) * in_ab * out_ab +
                   ZFT(T_ad) * in_ad * out_ab +
                   ZFT(T_ab) * in_ab * out_ad
            az_sum += ZFT(w_azi[ia]) * kern
        end
        acc += wl * az_sum * inv_2pi
    end
    return acc
end

@inline function _psi_dorsiventral(
    R_ad::FT, T_ad::FT, R_ab::FT, T_ab::FT,
    P_in::FT, N_in::FT, P_out::FT, N_out::FT,
) where {FT<:Real}
    two = one(FT) + one(FT)
    return two * (R_ad * P_in * N_out +
                  R_ab * N_in * P_out +
                  T_ad * P_in * P_out +
                  T_ab * N_in * N_out)
end

"""
    compute_Z_matrices_directional(mod::DorsiventralBiLambertianCanopyScattering,
                                   μ::AbstractVector, LD::AbstractLeafDistribution,
                                   m_max::Int; quadrature = CanopyQuadrature(),
                                   npol = 1, normalize_by_ω = false)

Analytic closed-form cosine Fourier moments `m = 0:m_max` for all four directional
canopy scattering blocks `(; Zpp, Zmp, Zpm, Zmm, ω_down, ω_up)` with zero
azimuthal quadrature:
- `Zpp`: downward → downward (`+μ_j → +μ_i`)
- `Zmp`: downward → upward (`+μ_j → -μ_i`, canopy top reflection)
- `Zpm`: upward → downward (`-μ_j → +μ_i`, canopy bottom reflection)
- `Zmm`: upward → upward (`-μ_j → -μ_i`)
- `ω_down`, `ω_up`: exact stream-dependent single-scattering albedos for `+μ` and `-μ`.

By default (`normalize_by_ω = false`), columns are normalized by `G(μ_in)` only
(so at `m = 0`, `0.5 * sum(w .* (Zpp[:, j, 1] .+ Zmp[:, j, 1])) ≈ ω_down[j]`),
preventing `ϖ`-mixing bugs when combined with specular or coupled BSDF kernels.
Setting `normalize_by_ω = true` divides out `ω_down[j]` and `ω_up[j]`, matching
legacy `BiLambertianCanopyScattering` identically when `R_ad == R_ab` and `T_ad == T_ab`.
"""
function compute_Z_matrices_directional(
    mod::DorsiventralBiLambertianCanopyScattering,
    μ::AbstractVector{FT},
    LD::AbstractLeafDistribution,
    m_max::Int;
    quadrature::CanopyQuadrature = CanopyQuadrature(),
    nQuad = nothing,
    npol::Integer = 1,
    normalize_by_ω::Bool = false,
) where {FT<:Real}
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))
    n = _validate_npol(npol)
    q = _resolve_quadrature(quadrature, nQuad)

    ZFT = float(promote_type(FT, typeof(mod.R_ad), typeof(mod.T_ad), typeof(mod.R_ab), typeof(mod.T_ab)))
    BFT = _leaf_base_float(ZFT)
    !isempty(μ) && all(x -> (v = prospect_primal(x); isfinite(v) && zero(BFT) < v <= one(BFT)), μ) ||
        throw(ArgumentError("compute_Z_matrices_directional requires non-empty μ in (0, 1]"))
    R_ad = ZFT(mod.R_ad); T_ad = ZFT(mod.T_ad)
    R_ab = ZFT(mod.R_ab); T_ab = ZFT(mod.T_ab)
    ω_ad = R_ad + T_ad
    ω_ab = R_ab + T_ab

    nμ = length(μ)
    nm = m_max + 1
    Zpp = zeros(ZFT, nμ, nμ, nm)
    Zmp = zeros(ZFT, nμ, nμ, nm)
    Zpm = zeros(ZFT, nμ, nμ, nm)
    Zmm = zeros(ZFT, nμ, nμ, nm)
    ω_down = zeros(ZFT, nμ)
    ω_up = zeros(ZFT, nμ)

    if prospect_primal(ω_ad + ω_ab) <= zero(BFT)
        return (;
            Zpp = _expand_intensity_Z(Zpp, n),
            Zmp = _expand_intensity_Z(Zmp, n),
            Zpm = _expand_intensity_Z(Zpm, n),
            Zmm = _expand_intensity_Z(Zmm, n),
            ω_down,
            ω_up,
        )
    end

    μ_vec = collect(μ)
    μ_work = ZFT.(μ_vec)
    leaf_quad = _leaf_inclination_quadrature(LD, q.n_leaf, FT)
    θₗ = leaf_quad.θ
    w_measure = leaf_quad.w_measure
    G_vec = ZFT.(vec(CanopyOptics.G(μ_vec, LD; quadrature = q)))

    Pꜜ = zeros(ZFT, nμ, nm)
    Nꜜ = zeros(ZFT, nμ, nm)
    Pꜛ = zeros(ZFT, nμ, nm)
    Nꜛ = zeros(ZFT, nμ, nm)
    P = zeros(ZFT, nm)
    N = zeros(ZFT, nm)
    Pm = zeros(ZFT, nm)
    Nm = zeros(ZFT, nm)

    for l in eachindex(θₗ)
        μ_L = ZFT(cos(θₗ[l]))
        lw = ZFT(w_measure[l])

        for i in 1:nμ
            μ_i = μ_work[i]
            _one_sided_projection_moments!(P, N, μ_i, μ_L, m_max)
            _one_sided_projection_moments!(Pm, Nm, -μ_i, μ_L, m_max)
            @inbounds for k in 1:nm
                Pꜜ[i, k] = P[k]
                Nꜜ[i, k] = N[k]
                Pꜛ[i, k] = Pm[k]
                Nꜛ[i, k] = Nm[k]
            end
            @inbounds ω_down[i] += lw * (ω_ad * P[1] + ω_ab * N[1])
            @inbounds ω_up[i]   += lw * (ω_ad * Pm[1] + ω_ab * Nm[1])
        end

        for k in 1:nm, j in 1:nμ
            @inbounds begin
                Pd_j = Pꜜ[j, k]; Nd_j = Nꜜ[j, k]
                Pu_j = Pꜛ[j, k]; Nu_j = Nꜛ[j, k]
                for i in 1:nμ
                    Zpp[i, j, k] += lw * _psi_dorsiventral(R_ad, T_ad, R_ab, T_ab, Pd_j, Nd_j, Pꜜ[i, k], Nꜜ[i, k])
                    Zmp[i, j, k] += lw * _psi_dorsiventral(R_ad, T_ad, R_ab, T_ab, Pd_j, Nd_j, Pꜛ[i, k], Nꜛ[i, k])
                    Zpm[i, j, k] += lw * _psi_dorsiventral(R_ad, T_ad, R_ab, T_ab, Pu_j, Nu_j, Pꜜ[i, k], Nꜜ[i, k])
                    Zmm[i, j, k] += lw * _psi_dorsiventral(R_ad, T_ad, R_ab, T_ab, Pu_j, Nu_j, Pꜛ[i, k], Nꜛ[i, k])
                end
            end
        end
    end

    @inbounds for j in 1:nμ
        gj = G_vec[j]
        if prospect_primal(gj) > zero(BFT)
            ω_down[j] /= gj
            ω_up[j]   /= gj
        end
    end

    two = one(ZFT) + one(ZFT)
    four = two + two
    for k in 1:nm
        ff = k == 1 ? two : four
        for j in 1:nμ
            gj = G_vec[j]
            if prospect_primal(gj) <= zero(BFT)
                Zpp[:, j, k] .= zero(ZFT); Zmp[:, j, k] .= zero(ZFT)
                Zpm[:, j, k] .= zero(ZFT); Zmm[:, j, k] .= zero(ZFT)
            else
                scale_dn = if normalize_by_ω
                    prospect_primal(ω_down[j]) > zero(BFT) ? ff / (ω_down[j] * gj) : zero(ZFT)
                else
                    ff / gj
                end
                scale_up = if normalize_by_ω
                    prospect_primal(ω_up[j]) > zero(BFT) ? ff / (ω_up[j] * gj) : zero(ZFT)
                else
                    ff / gj
                end
                @inbounds for i in 1:nμ
                    Zpp[i, j, k] *= scale_dn
                    Zmp[i, j, k] *= scale_dn
                    Zpm[i, j, k] *= scale_up
                    Zmm[i, j, k] *= scale_up
                end
            end
        end
    end

    return (;
        Zpp = _expand_intensity_Z(Zpp, n),
        Zmp = _expand_intensity_Z(Zmp, n),
        Zpm = _expand_intensity_Z(Zpm, n),
        Zmm = _expand_intensity_Z(Zmm, n),
        ω_down,
        ω_up,
    )
end

function compute_Z_matrices(
    mod::DorsiventralBiLambertianCanopyScattering,
    μ::AbstractVector{FT},
    LD::AbstractLeafDistribution,
    m::Union{Integer,AbstractUnitRange{<:Integer}};
    quadrature::CanopyQuadrature = CanopyQuadrature(),
    nQuad = nothing,
    npol::Integer = 1,
) where {FT<:Real}
    _is_face_symmetric(mod) || throw(ArgumentError(
        "Asymmetric DorsiventralBiLambertianCanopyScattering breaks up/down symmetry (Z⁺⁺ != Z⁻⁻, Z⁻⁺ != Z⁺⁻); call compute_Z_matrices_directional instead",
    ))
    sym = BiLambertianCanopyScattering(mod.R_ad, mod.T_ad)
    return compute_Z_matrices(sym, μ, LD, m; quadrature, nQuad, npol)
end

function compute_Z_matrices_aniso(
    mod::DorsiventralBiLambertianCanopyScattering,
    μ::AbstractArray{FT,1},
    LD::AbstractLeafDistribution,
    m::Int;
    quadrature::CanopyQuadrature = CanopyQuadrature(),
    nQuad = nothing,
    npol::Integer = 1,
) where {FT<:Real}
    return compute_Z_matrices(mod, μ, LD, m; quadrature, nQuad, npol)
end

function compute_Z_matrices_aniso_analytic(
    mod::DorsiventralBiLambertianCanopyScattering,
    μ::AbstractVector{FT},
    LD::AbstractLeafDistribution,
    m_max::Int;
    quadrature::CanopyQuadrature = CanopyQuadrature(),
    nQuad = nothing,
    npol::Integer = 1,
) where {FT<:Real}
    _is_face_symmetric(mod) || throw(ArgumentError(
        "Asymmetric DorsiventralBiLambertianCanopyScattering breaks up/down symmetry; call compute_Z_matrices_directional instead",
    ))
    sym = BiLambertianCanopyScattering(mod.R_ad, mod.T_ad)
    return compute_Z_matrices_aniso_analytic(sym, μ, LD, m_max; quadrature, nQuad, npol)
end

@inline function _interp_lut_01(tab::AbstractVector{T}, n_lut::Int, μ_val::T) where {T<:Real}
    u = clamp(μ_val, zero(T), one(T)) * T(n_lut)
    idx = clamp(floor(Int, prospect_primal(u)) + 1, 1, n_lut)
    frac = u - T(idx - 1)
    @inbounds return (one(T) - frac) * tab[idx] + frac * tab[idx + 1]
end

"""
    compute_Z_matrices_coupled_bsdf(prep::PreparedLeafBSDFSpectrum, idx_λ::Integer,
                                    μ::AbstractVector, LD::AbstractLeafDistribution,
                                    m_max::Int; quadrature = CanopyQuadrature(n_leaf=16, n_azimuth=24))

Project a coupled directional leaf BSDF ([`CoupledLeafBSDF`](@ref), prepared via
[`prepare_leaf_bsdf`](@ref) at wavelength index `idx_λ`) onto the 4-block canopy
Fourier phase matrices `(; Zpp, Zmp, Zpm, Zmm, ω_down, ω_up)` (`G_in`-normalized,
`normalize_by_ω = false`).

Uses rotationally invariant 3D inner products for both surface microfacet
reflection (`H = k_out - k_in`) and internal volume/forward-lobe scattering,
preserving exact leaf-level energy conservation `R(μ_i) + T(μ_i) + A(μ_i) == 1`
inside canopy radiative transfer. For very glossy leaves (`σ ≲ 0.05`, specular
lobe width `~ σ` rad), increase `quadrature` (e.g. `n_leaf >= 32, n_azimuth >= 48`)
to resolve the narrow specular peak without angular aliasing.
"""
function compute_Z_matrices_coupled_bsdf(
    prep::PreparedLeafBSDFSpectrum{FTP},
    idx_λ::Integer,
    μ::AbstractVector{FT},
    LD::AbstractLeafDistribution,
    m_max::Int;
    quadrature::CanopyQuadrature = CanopyQuadrature(n_leaf = 16, n_azimuth = 24),
) where {FTP<:Real,FT<:Real}
    @boundscheck checkbounds(prep.nᵣ, idx_λ)
    m_max < 0 && throw(ArgumentError("m_max must be non-negative"))

    ZFT = float(promote_type(FTP, FT))
    BFT = _leaf_base_float(ZFT)
    !isempty(μ) && all(x -> (v = prospect_primal(x); isfinite(v) && zero(BFT) < v <= one(BFT)), μ) ||
        throw(ArgumentError("compute_Z_matrices_coupled_bsdf requires non-empty μ in (0, 1]"))
    nμ = length(μ)
    nm = m_max + 1

    nsurf_ad = ZFT(prep.n_surf_ad[idx_λ])
    nsurf_ab = ZFT(prep.n_surf_ab[idx_λ])
    σ_ad = ZFT(prep.σ_ad); ns_ad = ZFT(prep.n_scale_ad); n_eff_ad = ns_ad * nsurf_ad
    σ_ab = ZFT(prep.σ_ab); ns_ab = ZFT(prep.n_scale_ab); n_eff_ab = ns_ab * nsurf_ab
    th_ad = ZFT(prep.t_hem_ext_ad[idx_λ])
    th_ab = ZFT(prep.t_hem_ext_ab[idx_λ])
    rv_ad = ZFT(prep.r_vol_hem_ad[idx_λ])
    rv_ab = ZFT(prep.r_vol_hem_ab[idx_λ])
    t_hem = ZFT(prep.t_hem[idx_λ])
    a_ad = ZFT(prep.a_hem_ad[idx_λ])
    a_ab = ZFT(prep.a_hem_ab[idx_λ])
    η_t = ZFT(prep.η_t)
    γ_t = ZFT(prep.γ_t)
    cn = ZFT(prep.c_norm[idx_λ])
    d_bar = ZFT(prep.d_bar[idx_λ])
    has_lobe = prospect_primal(η_t) > zero(BFT)
    surf_ad_obj = RoughDielectricSurface{ZFT}(σ_ad, ns_ad)
    surf_ab_obj = RoughDielectricSurface{ZFT}(σ_ab, ns_ab)

    μ_vec = collect(μ)
    G_vec = ZFT.(vec(CanopyOptics.G(μ_vec, LD; quadrature)))
    leaf_quad = _leaf_inclination_quadrature(LD, quadrature.n_leaf, FT)
    # Periodic midpoint quadrature on [0, 2π) with even N_ϕ: spectrally accurate on S^1
    # AND invariant under ϕ_L -> π - ϕ_L, preserving discrete 4-block reciprocity to ~1e-14!
    nL_ϕ = 2 * max(1, div(quadrature.n_azimuth + 1, 2))
    dϕ_out, w_dϕ = gauleg(quadrature.n_azimuth, FT(0), FT(π)) # even symmetry in Δϕ!

    # Helper evaluating scalar leaf BSDF f(k_in -> k_out; n_ad) * |k_in ⋅ n_ad| * |k_out ⋅ n_ad|
    two = one(ZFT) + one(ZFT)
    four = two + two
    inv_pi = inv(ZFT(π))
    inv_2pi = inv(two * ZFT(π))

    # Precompute 257-node 1D lookup tables on μ_loc ∈ [0, 1] for r_ext(μ) and L_s'(μ)
    # so the 5D (j, i, Δϕ, ϕ_L, θ_L) canopy integration never repeats 2D quadrature in the inner loop.
    n_lut = 256
    r_ext_ad_tab = Vector{ZFT}(undef, n_lut + 1)
    r_ext_ab_tab = Vector{ZFT}(undef, n_lut + 1)
    L_ad_tab     = Vector{ZFT}(undef, n_lut + 1) # L_{surf_ad}(μ)
    L_ab_tab     = Vector{ZFT}(undef, n_lut + 1) # L_{surf_ab}(μ)
    inv_nlut     = inv(ZFT(n_lut))
    for kl in 0:n_lut
        μ_q = max(ZFT(BFT(1e-6)), ZFT(kl) * inv_nlut)
        r_ext_ad_tab[kl + 1] = microfacet_directional_reflectance(μ_q, n_eff_ad, σ_ad)
        r_ext_ab_tab[kl + 1] = microfacet_directional_reflectance(μ_q, n_eff_ab, σ_ab)
        if has_lobe
            s_q = sqrt(max(zero(ZFT), one(ZFT) - μ_q * μ_q))
            L_ad_tab[kl + 1] = _btdf_lobe_escape_integral(μ_q, s_q, surf_ad_obj, nsurf_ad, th_ad, γ_t, cn)
            L_ab_tab[kl + 1] = _btdf_lobe_escape_integral(μ_q, s_q, surf_ab_obj, nsurf_ab, th_ab, γ_t, cn)
        else
            L_ad_tab[kl + 1] = zero(ZFT)
            L_ab_tab[kl + 1] = zero(ZFT)
        end
    end

    @inline function eval_leaf_kernel(kin, kout, om_l)
        # om_l is downward abaxial normal; n_ad = -om_l
        d_in = _dot3(kin, om_l)   # > 0 => hits adaxial; < 0 => hits abaxial
        d_out = _dot3(kout, om_l) # < 0 => exits adaxial; > 0 => exits abaxial
        (iszero(prospect_primal(d_in)) || iszero(prospect_primal(d_out))) && return zero(ZFT)

        hit_ad = prospect_primal(d_in) > zero(BFT)
        μi_loc = min(one(ZFT), hit_ad ? d_in : -d_in)

        is_refl = (prospect_primal(d_in) * prospect_primal(d_out)) < zero(BFT)
        μo_loc = min(one(ZFT), abs(d_out))

        σ_i = hit_ad ? σ_ad : σ_ab
        ne_i = hit_ad ? n_eff_ad : n_eff_ab
        th_i = hit_ad ? th_ad : th_ab
        r_ext_i = _interp_lut_01(hit_ad ? r_ext_ad_tab : r_ext_ab_tab, n_lut, μi_loc)
        E_i = (one(ZFT) - r_ext_i) / th_i

        if is_refl
            # H = ω_i + ω_o = k_out - k_in
            Hx = kout[1] - kin[1]; Hy = kout[2] - kin[2]; Hz = kout[3] - kin[3]
            H2 = Hx * Hx + Hy * Hy + Hz * Hz
            Hnorm = sqrt(H2)
            Hz_loc = μi_loc + μo_loc
            Hz_loc2 = Hz_loc * Hz_loc
            tan2_th = max(zero(ZFT), H2 - Hz_loc2) / Hz_loc2
            cos2_th = Hz_loc2 / H2
            D = beckmann_D(tan2_th, cos2_th, σ_i)
            G2 = smith_beckmann_G1(μi_loc, σ_i) * smith_beckmann_G1(μo_loc, σ_i)
            c_loc = Hnorm / two
            s_loc2 = max(zero(ZFT), one(ZFT) - c_loc * c_loc)
            F = fresnel_unpolarized_cos(c_loc, s_loc2, ne_i)
            # Note: μi_loc * μo_loc * f_surf == (D * G2 * F) / 4 analytically (no near-horizon 0/0)!
            w_surf = (D * G2 * F) / four

            r_ext_o = _interp_lut_01(hit_ad ? r_ext_ad_tab : r_ext_ab_tab, n_lut, μo_loc)
            E_o = (one(ZFT) - r_ext_o) / th_i
            rv_h = hit_ad ? rv_ad : rv_ab
            f_vol = rv_h * inv_pi * E_i * E_o
            return w_surf + μi_loc * μo_loc * f_vol
        else
            th_t = hit_ad ? th_ab : th_ad
            r_ext_o = _interp_lut_01(hit_ad ? r_ext_ab_tab : r_ext_ad_tab, n_lut, μo_loc)
            E_o = (one(ZFT) - r_ext_o) / th_t
            if !has_lobe
                return μi_loc * μo_loc * (t_hem * inv_pi * E_i * E_o)
            else
                c_dot = _dot3(kin, kout)
                p_tilde = exp(γ_t * (c_dot - one(ZFT))) / cn
                L_tr_i = _interp_lut_01(hit_ad ? L_ab_tab : L_ad_tab, n_lut, μi_loc)
                L_in_o = _interp_lut_01(hit_ad ? L_ad_tab : L_ab_tab, n_lut, μo_loc)
                closure = (one(ZFT) - L_tr_i) * (one(ZFT) - L_in_o) * inv_pi / d_bar
                f_t = t_hem * E_i * E_o * ((one(ZFT) - η_t) * inv_pi + η_t * (p_tilde + closure))
                return μi_loc * μo_loc * f_t
            end
        end
    end

    Zpp = zeros(ZFT, nμ, nμ, nm)
    Zmp = zeros(ZFT, nμ, nμ, nm)
    Zpm = zeros(ZFT, nμ, nμ, nm)
    Zmm = zeros(ZFT, nμ, nμ, nm)
    ω_down = zeros(ZFT, nμ)
    ω_up = zeros(ZFT, nμ)

    # Precompute leaf normals (θ_L, ϕ_L) and weights on uniform periodic midpoint grid
    nL_θ = length(leaf_quad.θ)
    om_L_grid = Array{NTuple{3,ZFT}}(undef, nL_θ, nL_ϕ)
    w_L_grid = Matrix{ZFT}(undef, nL_θ, nL_ϕ)
    w_az_unif = inv(ZFT(nL_ϕ))
    two_pi = two * ZFT(π)
    for l in 1:nL_θ
        μL = ZFT(cos(leaf_quad.θ[l]))
        sL = ZFT(sin(leaf_quad.θ[l]))
        wl = ZFT(leaf_quad.w_measure[l]) * w_az_unif
        for ip in 1:nL_ϕ
            p = two_pi * (ZFT(ip) - inv(two)) * w_az_unif
            om_L_grid[l, ip] = (sL * cos(p), sL * sin(p), μL)
            w_L_grid[l, ip] = wl
        end
    end

    # Compute exact stream-dependent albedos ω_down[j], ω_up[j] from leaf directional albedos
    for j in 1:nμ
        μj = ZFT(μ_vec[j])
        sj = sqrt(max(zero(ZFT), one(ZFT) - μj * μj))
        kin_dn = (sj, zero(ZFT), μj)
        kin_up = (sj, zero(ZFT), -μj)
        acc_dn = zero(ZFT)
        acc_up = zero(ZFT)
        acc_G_dn = zero(ZFT)
        acc_G_up = zero(ZFT)
        for ip in 1:nL_ϕ, l in 1:nL_θ
            om_l = om_L_grid[l, ip]
            wL = w_L_grid[l, ip]
            for (kin, is_dn) in ((kin_dn, true), (kin_up, false))
                d_in = _dot3(kin, om_l)
                if !iszero(prospect_primal(d_in))
                    hit_ad = prospect_primal(d_in) > zero(BFT)
                    μi_loc = min(one(ZFT), abs(d_in))
                    r_s = _interp_lut_01(hit_ad ? r_ext_ad_tab : r_ext_ab_tab, n_lut, μi_loc)
                    E_i = (one(ZFT) - r_s) / (hit_ad ? th_ad : th_ab)
                    a_loc = E_i * (hit_ad ? a_ad : a_ab)
                    ω_loc = one(ZFT) - a_loc
                    if is_dn
                        acc_dn += wL * μi_loc * ω_loc
                        acc_G_dn += wL * μi_loc
                    else
                        acc_up += wL * μi_loc * ω_loc
                        acc_G_up += wL * μi_loc
                    end
                end
            end
        end
        if prospect_primal(acc_G_dn) > zero(BFT)
            ω_down[j] = acc_dn / acc_G_dn
        end
        if prospect_primal(acc_G_up) > zero(BFT)
            ω_up[j] = acc_up / acc_G_up
        end
    end

    # Assemble Z blocks using symmetric half-azimuths (-Δϕ/2 -> +Δϕ/2) so discrete time-reversal
    # maps ϕ_L -> π - ϕ_L and preserves 4-block Helmholtz reciprocity to machine precision!
    half = inv(two)
    for j in 1:nμ
        μj = ZFT(μ_vec[j])
        sj = sqrt(max(zero(ZFT), one(ZFT) - μj * μj))

        for i in 1:nμ
            μi = ZFT(μ_vec[i])
            si = sqrt(max(zero(ZFT), one(ZFT) - μi * μi))
            for idp in eachindex(dϕ_out)
                dϕ = ZFT(dϕ_out[idp])
                wd = two * ZFT(w_dϕ[idp]) # factor 2 maps [0, π] -> [0, 2π]
                ch, sh = cos(dϕ * half), sin(dϕ * half)
                kin_dn  = (sj * ch, -sj * sh, μj)
                kin_up  = (sj * ch, -sj * sh, -μj)
                kout_dn = (si * ch,  si * sh, μi)
                kout_up = (si * ch,  si * sh, -μi)

                f_pp = zero(ZFT); f_mp = zero(ZFT)
                f_pm = zero(ZFT); f_mm = zero(ZFT)
                @inbounds for ip in 1:nL_ϕ, l in 1:nL_θ
                    om_l = om_L_grid[l, ip]
                    wL = w_L_grid[l, ip]
                    f_pp += wL * eval_leaf_kernel(kin_dn, kout_dn, om_l)
                    f_mp += wL * eval_leaf_kernel(kin_dn, kout_up, om_l)
                    f_pm += wL * eval_leaf_kernel(kin_up, kout_dn, om_l)
                    f_mm += wL * eval_leaf_kernel(kin_up, kout_up, om_l)
                end

                @inbounds for k in 1:nm
                    cm = cos(ZFT(k - 1) * dϕ) * wd
                    Zpp[i, j, k] += cm * f_pp
                    Zmp[i, j, k] += cm * f_mp
                    Zpm[i, j, k] += cm * f_pm
                    Zmm[i, j, k] += cm * f_mm
                end
            end
        end
    end

    for k in 1:nm
        ff = k == 1 ? two : four
        for j in 1:nμ
            gj = G_vec[j]
            scale = prospect_primal(gj) > zero(BFT) ? ff / gj : zero(ZFT)
            @inbounds for i in 1:nμ
                Zpp[i, j, k] *= scale
                Zmp[i, j, k] *= scale
                Zpm[i, j, k] *= scale
                Zmm[i, j, k] *= scale
            end
        end
    end

    return (; Zpp, Zmp, Zpm, Zmm, ω_down, ω_up)
end
