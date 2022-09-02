"""
Given a complex number, return the complex number with absolute value of both parts
"""
abs_components(x::Complex) = complex(abs(real(x)),abs(imag(x)))

"""
Given a complex number, return the complex number with same real part and absolute value of imaginary part
"""
abs_imag_only(x::Complex) = complex(real(x),abs(imag(x)))

"""
Calculate forward scattering amplitudes
"""
function afsal(Θᵢʳ, leaf::Leaf, k₀)

    # Calculate integral of p(θ)*sin(θ)^2 from 0 to pi
    # That is, the average integral over inclinations
    # prob needs to be replaced with Distributions.jl (which can integrate better anyhow)
    f(x) = prob(x, leaf.pdf_num, leaf.pdf_param) * sin(x)^2
    ai1  = quadgk(f, 0, π)[1]
    # Volume of leaf 
    vol  = π * leaf.a_maj * leaf.b_min * leaf.t / 4
    # ν² * vol / 4π      
    β    = k₀^2 * vol / (4π) 

    xi   = leaf.ϵ - 1
    xii  = xi/(2 * (1 + xi))
    afhhl = β * xi * (1 - xii * ai1)
    afvvl = β * xi * (1 - xii * (2 * sin(Θᵢʳ)^2 + (cos(Θᵢʳ)^2 - 2 * sin(Θᵢʳ)^2) * ai1))
    return (afhhl, afvvl)
end

"""
Calculate scattering coefficients for a thin disk
"""
function asal(Θᵢʳ, leaf::Leaf, k₀)

    # Equations 14 and 29, leading term
    a2 = abs(k₀^2 * (leaf.ϵ - 1) / sqrt(4π))^2 

    # 
    eb = (leaf.ϵ - 1) / leaf.ϵ
    b2 = abs(eb)^2

    sm = [0.0, 0.0, 0.0]
    sidr = [0.0, 0.0]
    sivh = [0.0, 0.0]

    cnst3 = 2π*leaf.a_maj*leaf.b_min/4

    cnti=1.0
	cntj=1.0

    # Double integral trapezoidal rule? 

    # Loop over θ
    for i = 1:n_θ+1

        θ_curr = (i - 1) * δθ
        
        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,leaf.pdf_num,leaf.pdf_param)
        pdf < 1.0e-03 && break

        sumx = [0.0, 0.0, 0.0]
        sjdr = [0.0, 0.0]
        sjvh = [0.0, 0.0]

        for j = 1:n_ϕ+1

            ϕ_curr=(j-1)*δϕ

            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj

            # calculate swig--beam pattern
            alphb = (cos(θ_curr)*sin(Θᵢʳ))*sin(ϕ_curr)
            alphd = alphb - (sin(θ_curr)*cos(Θᵢʳ))
            betad = sin(Θᵢʳ)*cos(ϕ_curr)

            nud   = sqrt(alphd^2*(leaf.a_maj/2)^2+betad^2*(leaf.b_min/2)^2)* k₀/π
            nub   = sqrt(alphb^2*(leaf.a_maj/2)^2+betad^2*(leaf.b_min/2)^2)* k₀/π

            zetad = 2π*nud
            zetab = 2π*nub
            swigd = cnst3 * (besselj1(zetad) / zetad)
            swigb = cnst3 * (besselj1(zetab) / zetab)
            
            # calculate integrands
            sthcp2 = (sin(θ_curr)*cos(ϕ_curr))^2
            ss2    = (sin(θ_curr)*cos(Θᵢʳ)*sin(ϕ_curr)-(cos(θ_curr)*sin(Θᵢʳ)))^2
            scsph2 = (sin(θ_curr)*cos(Θᵢʳ))^2*sin(ϕ_curr)^2
            scs2   = (cos(θ_curr)*sin(Θᵢʳ))^2 - scsph2
            scs1   = scsph2+(cos(θ_curr)*sin(Θᵢʳ))^2
            sscci  = cos(Θᵢʳ)^2-sin(Θᵢʳ)^2
            sccs   = (cos(θ_curr)*sin(Θᵢʳ))^2-scsph2
            
            cnt    = cnti*cntj

            sumx[1] += (abs(1.0-eb*sthcp2))^2*cnt*swigd^2
            sumx[2] += sthcp2*ss2*cnt*swigd^2
            sumx[3] += (abs(1.0-eb*ss2))^2*cnt*swigd^2

            sjdr[1] +=  (abs(1.0-eb*sthcp2))^2*cnt*swigb^2 # hh
            sjdr[2] +=  (abs(sscci+eb*sccs))^2*cnt*swigb^2 # vv

            sjvh[1] += sthcp2*scs1*cnt*swigb^2
            sjvh[2] += sthcp2*scs2*cnt*swigb^2

            cntj=1.0
        end

        sm   += sumx * pdf
        sidr += sjdr * pdf 

		sivh[1]+=sjvh[1]*pdf
		sivh[2]+=sjvh[2]*pdf
		cnti=1.0
    end

    sgbd = a2*leaf.t^2 * sm * Δθ * Δϕ
    sgbd[2] *= b2

    sgbdr  =  a2*leaf.t^2 * sidr*Δθ*Δϕ
    sgbvh  =  a2*leaf.t^2 * b2 * abs.(sivh)*Δθ*Δϕ

    return BackscatterFields(sgbd, sgbdr, sgbvh)
end

"""
Calculation of average forward scattering amplitudes of a finite length cylinder, 
exact series solution
"""
function wood_forward(wood::Wood)

    # Maximum n for scattering loop 
    nmax= min(20, Integer(floor(k₀*wood.r+4*(k₀*wood.r)^(1/3)+2)))

    # To store total averaged scattering matrix
    S = zeros(Float64, 2, 2)

    # Loop over θs
    for i = 1:n_θ+1

        # Current θ
        θ_curr = (i-1)*δθ
        cnti = (i == 1 || i == n_θ+1) ? 0.5 : 1.0

        # Probability of θ 
        pdf = prob(θ_curr,wood.pdf_num,wood.pdf_param)
        (pdf < 1.0e-03) && break

        # To store averaged scattering matrix for current θ
        S_θ = zeros(Float64, 2, 2)

        # Loop over ϕs
        for j = 1:n_ϕ+1

            # Current ϕ
            ϕ_curr = (j-1)*δϕ
            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : 1.0
            cnt=cnti*cntj

            # Get scattering at current angle
            dsi, S_curr = scattering(θ_curr, ϕ_curr, nmax, 1, x -> -x, x->x, x -> π-x, x->x, wood)
            S_θ += (cnt*S_curr/dsi) * Δθ * Δϕ
        end

        # Weight scattering sum by θ pdf
        S += S_θ * pdf

    end

    # Return total S
    return S
end

"""
Calculation of average backscatter cross-section of a finite length cylinder, exact series solution
(KARAM, FUNG, AND ANTAR, 1988)
"""
function wood_backward(wood::Wood)

    # Maximum n for scattering loop
    nmax=min(20, Integer(floor(k₀*wood.r+4.0*(k₀*wood.r)^(1/3)+2.0)))

    # To store total averaged scattering matrix
    S_d  = zeros(3)  # hh, vh, vv
    S_dr = zeros(2)  # hh, vv
    S_vh = zeros(2)  # smvh[1] is smvh1, smvh[2] is **smvh3**

    # Loop over θs
    for i = 1:n_θ+1

        # Current θ
        θ_curr = (i-1) * δθ
        cnti = (i == 1 || i == n_θ+1) ? 0.5 : 1.0

        # Probability of θ 
        pdf = prob(θ_curr,wood.pdf_num,wood.pdf_param)
        (pdf < 1.0e-03) && break

        # To store averaged scattering matrix for current θ
        S_d_θ  = zeros(3)  # hh, vh, vv
        S_dr_θ = zeros(2)  # hh, vv
        S_vh_θ = zeros(2)  # smvh[1] is smvh1, smvh[2] is **smvh3**

        # Loop over ϕs
        for j = 1:n_ϕ+1

            # Current ϕ
            ϕ_curr= (j-1) * δϕ
            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : 1.0
            cnt = cnti*cntj

            # Get direct scattering at current angle
            dsi, S = scattering(θ_curr, ϕ_curr, nmax, 1, x -> -x, x->-x, x -> x, x->x+π, wood)
            S_d_θ += abs.([S[1,1] ; S[1,2] ; S[2,2]] / dsi).^2*cnt

            # Get direct-reflected scattering at current angle
            dsi, S = scattering(θ_curr, ϕ_curr, nmax, 1, x -> x, x->-x, x -> π-x, x->x+π, wood)
            S_dr_θ += abs.([S[1,1] ; S[2,2]] / dsi).^2*cnt
            S_vh_θ[1] += abs(S[1,2]/(dsi))^2*cnt
            fvhc=conj(S[1,2])

            # Get vh scattering at current angle
            dsi, S = scattering(θ_curr, ϕ_curr, nmax, -1, x -> -x, x->-x, x -> π-x, x->x+π, wood)
            S_vh_θ[2] += abs(S[1,2]*fvhc/(dsi*dsi))*cnt
        end

        # Weight scattering sum by θ pdf
        S_d  += S_d_θ  * pdf
        S_dr += S_dr_θ * pdf
        S_vh += S_vh_θ * pdf

    end

    # Why are we multiplying by 4π here? 
    S_d  = 4π * Δϕ * Δθ * S_d  
    S_dr = 4π * Δϕ * Δθ * S_dr 
    S_vh = 4π * Δϕ * Δθ * S_vh 

    return BackscatterFields(S_d, S_dr, S_vh)

end


"""
Calculate scattering cross-section of a finite length cylinder with specified orientation 
"""
function scattering(θ_curr, ϕ_curr, nmax, sign_i, new_tvs, new_ths, new_thetas, new_phys, wood::Wood)

    # B.44 from Zyl and Kim book Synthetic Aperture Radar Polarimetry
    tvi = -sin(θ_curr)*cos(Θᵢʳ)*cos(ϕ_curr-ϕ_iʳ)-sign_i*cos(θ_curr)*sin(Θᵢʳ)
    # B.45 from Zyl and Kim book Synthetic Aperture Radar Polarimetry
    thi = sign_i*sin(θ_curr)*sin(ϕ_curr-ϕ_iʳ)
    ti  = sqrt(tvi^2+thi^2)

    ths = new_ths(thi)
    tvs = new_tvs(tvi)
    ts  = sqrt(ths^2+tvs^2)

    cthi = sin(θ_curr)*sin(Θᵢʳ)*cos(ϕ_curr-ϕ_iʳ)-sign_i*cos(θ_curr)*cos(Θᵢʳ)
    cphi = (sin(Θᵢʳ)*cos(θ_curr)*cos(ϕ_curr-ϕ_iʳ)+sign_i*sin(θ_curr)*cos(Θᵢʳ))/sqrt(1-cthi^2)

    cthi,cphi = max.(min.([cthi, cphi], 1.0), -1.0)

    thetai, phyi = acos.([cthi, cphi])

    thetas = new_thetas(thetai)
    phys   = new_phys(phyi)

    S = scattering_amplitude(nmax, wood, thetai, thetas, phyi, phys)

    dsi=ti*ts  

    # B.35 from Zyl and Kim book Synthetic Aperture Radar Polarimetry
    Ts = [tvs  ths ; -ths tvs]
    Ti = [tvi -thi ;  thi tvi]
    S_new = Ts * S * Ti

    return dsi, S_new

end

"""
Calculate bistatic scattering amplitude in the cylinder coordinate

Equations 26, 27 in 
Electromagnetic Wave Scattering from Some Vegetation Samples
(KARAM, FUNG, AND ANTAR, 1988)
http://www2.geog.ucl.ac.uk/~plewis/kuusk/karam.pdf
"""
function scattering_amplitude(nmax, wood_c::Wood, θᵢ, Θₛ, ϕᵢ, ϕₛ)

    μᵢ = cos(θᵢ)
    μₛ = cos(Θₛ)
    sin_θₛ = sqrt(1-μₛ^2)

    # [hh hv ; vh vv]
    s = zeros(Complex, 2,2)

    # Equation 26 (Karam, Fung, Antar)
    argum = k₀*(wood_c.l/2)*( μᵢ+μₛ )
    μ_si = (argum == 0.0) ? 1.0 : sin(argum)/argum

    # Calculate scattering amplitude for n = 0
    z₀, A₀, B₀, e_0v, e_0h, ηh_0v, ηh_0h = cylinder(0, wood_c, θᵢ, Θₛ)
    s0 = [B₀*ηh_0h 0 ; 0 e_0v*(B₀ * μₛ * μᵢ - sin_θₛ * z₀)]

    # Loop over all needed n values
    for n = 1:nmax

        # Cylindrical scattering coefficients at current n 
        zₙ,Aₙ,Bₙ,e_nv,e_nh,ηh_nv,ηh_nh = cylinder(n, wood_c, θᵢ, Θₛ)

        sphn = 0.0
        π_δ  = 0.0001
        cphn = (π - π_δ < ϕₛ - ϕᵢ < π + π_δ) ? (-1)^n : 1.0
        
        # Inside summation term in each line of Equation 27
        s[1,1] += 2*(ηh_nh*Bₙ + ej*e_nh*Aₙ*μᵢ) * cphn
        s[1,2] +=   (ηh_nv*Bₙ + ej*e_nv*Aₙ*μᵢ) * sphn
        s[2,1] +=   ((e_nh* μᵢ *Bₙ - ej*ηh_nh*Aₙ)* μₛ - sin_θₛ * e_nh*zₙ) * sphn
        s[2,2] += 2*((e_nv* μᵢ *Bₙ - ej*ηh_nv*Aₙ)* μₛ - sin_θₛ * e_nv*zₙ) * cphn
        
    end

    # Added coefficients to finish Equation 27
    h = (wood_c.l/2)
    F = k₀^2 * h * (wood_c.ϵ-1.0) * μ_si * (s0 + s)
    F[1,2] *= 2 * ej
    F[2,1] *= 2 * ej

    return F

end

"""
Compute the coefficients of the series expansion of cylinder scattering amplitude.

Equation 24 in 
Electromagnetic Wave Scattering from Some Vegetation Samples
(KARAM, FUNG, AND ANTAR, 1988)
http://www2.geog.ucl.ac.uk/~plewis/kuusk/karam.pdf
"""
function cylinder(n::Integer, wood_c::Wood, Θᵢ::Real, Θₛ::Real)

    # Constants
    a  = wood_c.r
    ϵᵣ = wood_c.ϵ
    μᵢ  = cos(Θᵢ)
    μₛ  = cos(Θₛ)
    coef1   = sqrt(ϵᵣ-μᵢ^2)	
    sin_Θᵢ  = sqrt( 1-μᵢ^2)
    
    # Last equations in 24
    u  = k₀ * a * coef1
    vᵢ = k₀ * a * sqrt(1-μᵢ^2)
    vₛ = k₀ * a * sqrt(1-μₛ^2)

    # Bessel functions and derivatives at u, v
    Jₙu=besselj(n,u)
    Hₙv=besselj(n,vᵢ) - ej*bessely(n,vᵢ)
    d_Jₙu=dbessj(n,u)
    d_Hₙv=dbessj(n,vᵢ) - ej*dbessy(n,vᵢ)
    
    # Equation 26 (Karam, Fung, Antar)
    zₙ, znp1, znm1 = zeta(n, a, u, vₛ)

    # Equation 28 (Karam, Fung, Antar)
    Aₙ = (znm1-znp1)/(2*coef1)
    Bₙ = (znm1+znp1)/(2*coef1)

    # Compute Rₙ with sub-expressions
    term1   = (d_Hₙv/(vᵢ*Hₙv)-d_Jₙu/(u*Jₙu))
    term2   = (d_Hₙv/(vᵢ*Hₙv)-ϵᵣ*d_Jₙu/(u*Jₙu))
    term3   = (1/vᵢ^2-1/u^2)*n*μᵢ
    Rₙ      = (π*vᵢ^2*Hₙv/2) * (term1*term2-term3^2)

    # Calculate remaining terms, using above sub-expressions
    e_nv   = ej*sin_Θᵢ*term1/(Rₙ*Jₙu)
    e_nh   = -sin_Θᵢ*term3/(Rₙ*Jₙu)	
    ηh_nv = -e_nh
    ηh_nh = ej*sin_Θᵢ*term2/(Rₙ*Jₙu)

    return (zₙ, Aₙ, Bₙ, e_nv, e_nh, ηh_nv, ηh_nh)
end

"""
Calculate zeta function using bessel functions

Equation 26 in 
Electromagnetic Wave Scattering from Some Vegetation Samples
(KARAM, FUNG, AND ANTAR, 1988)
http://www2.geog.ucl.ac.uk/~plewis/kuusk/karam.pdf
"""
function zeta(n::Integer, a, u, vs)

    # Pre-compute bessel function outputs at n-1, n, n+1, n+2
    bjnu=besselj(n,u)
    bjnv=besselj(n,vs)

    bjnp1u=besselj(n+1,u)
    bjnp1v=besselj(n+1,vs)

    bjnp2v=besselj(n+2,vs)
    bjnp2u=besselj(n+2,u)

    bjnm1u=besselj(n-1,u)
	bjnm1v=besselj(n-1,vs)

    # Compute z_n values
    znp1 = (a^2 / (u^2-vs^2)) * (u * bjnp1v * bjnp2u - vs * bjnp1u * bjnp2v) 
    zn   = (a^2 / (u^2-vs^2)) * (u * bjnv   * bjnp1u - vs * bjnu   * bjnp1v) 
    znm1 = (a^2 / (u^2-vs^2)) * (u * bjnm1v * bjnu   - vs * bjnm1u * bjnv)   

    return zn, znp1, znm1
end

# Bessel function derivatives 
dbessj(n, x) = -n/x*besselj(n,x)+besselj(n-1,x)
dbessy(n, x) = -n/x*bessely(n,x)+bessely(n-1,x)

"""
3.1.2 in SAATCHI, MCDONALD
"""
function funcm(x, y, d)

    dag = x+y
    dagd = dag * d
    z = abs(dagd)

    # Edge cases 
    z < 1.0e-03 && return d
    z > 1.6e02 && return 1.0/dag

    return (1.0  - exp(-(dagd)))/dag

end

function funcp(x, y, d)

    dag=x+y
    dagd=dag*d
    z=abs(dagd)

    # Edge cases 
    z < 1e-3 && return d
    z > 1.6e2 && return 1.0/dag

    return (exp(dagd)-1.0)/dag

end

function cfun(a, e, d)

    dag=a+e
    dagd=dag*d
    a3=abs(dagd)

    # Edge cases 
    a3 < 1.0e-03  && return d
    a3 > 1.6e02 && return 1.0/dag

    return (1.0 - exp(-(dagd)))/dag

end

"""
Oh et al, 1992 semi-emprical model for rough surface scattering

Equation 3.1.5 in 
Coherent Effects in Microwave Backscattering Models for Forest Canopies
(SAATCHI, MCDONALD)
https://www.researchgate.net/publication/3202927_Semi-empirical_model_of_the_
ensemble-averaged_differential_Mueller_matrix_for_microwave_backscattering_from_bare_soil_surfaces
"""
function grdoh(ksig)

    r0 = abs((1-sqrt(ϵ_g))/(1+sqrt(ϵ_g)))^2
    g=0.7*(1-exp(-0.65*ksig^1.8))
    q=0.23*(1-exp(-ksig))*sqrt(r0)
    sp1=(2*Θᵢʳ/π)^(1/(3*r0))
    sp=1-sp1*exp(-ksig)

    rgh  = (cos(Θᵢʳ)-sqrt(ϵ_g-sin(Θᵢʳ)^2))/(cos(Θᵢʳ)+sqrt(ϵ_g-sin(Θᵢʳ)^2))
    rgv  = (ϵ_g*cos(Θᵢʳ)-sqrt(ϵ_g-sin(Θᵢʳ)^2))/(ϵ_g*cos(Θᵢʳ)+sqrt(ϵ_g-sin(Θᵢʳ)^2))

    rh0, rv0 = abs.([rgh rgv]).^2

    sighh=g*sp*cos(Θᵢʳ)^3*(rh0+rv0)
    sigvv=g*cos(Θᵢʳ)^3*(rh0+rv0)/sp
    sigvh=q*sigvv	

    # 3.1.5
    shhg=sighh*exp(-4*imag(K_hc*d_c + K_ht*d_t)) # Factor of 2 pulled out in exp
    svvg=sigvv*exp(-4*imag(K_vc*d_c + K_vt*d_t)) # Factor of 2 pulled out in exp
    svhg=sigvh*exp(-2*imag((K_hc+K_vc)*d_c + (K_ht+K_vt)*d_t))

    sg = [shhg svhg svvg]
    gd = 10 * log10.([sighh sigvh sigvv])

    return sg, gd 

end