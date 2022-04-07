"""
Calculate forward scattering amplitudes
"""
function afsal(θ_iʳ, ai1, leaf::Leaf)
    vol = π * leaf.a_maj * leaf.b_min * leaf.t / 4      # Volume of leaf 
    β = k₀^2 * vol / (4π)                   # ν² * vol / 4π
    xi = leaf.ϵ - 1
    xii = xi/(2 * (1 + xi))
    afhhl = β * xi * (1 - xii * ai1)
    afvvl = β * xi * (1 - xii * (2 * sin(θ_iʳ)^2 + (cos(θ_iʳ)^2 - 2 * sin(θ_iʳ)^2) * ai1))
    return (afhhl, afvvl)
end

"""
Calculate scattering coefficients for a thin disk
"""
function asal(θ_iʳ, leaf::Leaf)

    # Needed constants
    exi = leaf.ϵ - 1.0
    ea = k₀^2 * exi / sqrt(4π)
    a2 = abs(ea)^2
    eb = exi / (1 + exi)
    b2 = abs(eb)^2

    sm = [0.0, 0.0, 0.0]
    sidr = [0.0, 0.0]

    sivh1 = 0.0
    sivh3 = 0.0
    cnst3 = 2π*leaf.a_maj*leaf.b_min/4.0

    cnti=1.0
	cntj=1.0

    # Loop over θ
    for i = 1:n_θ+1

        θ_curr = (i - 1) * δθ
        
        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,leaf.pdf_num,leaf.pdf_param)
        pdf < 1.0e-03 && break

        sumx = [0.0, 0.0, 0.0]
        sjdr = [0.0, 0.0]

        sjvh1=0.0
        sjvh3=0.0

        for j = 1:n_ϕ+1

            ϕ_curr=(j-1)*δϕ

            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj

            # calculate swig--beam pattern
            alphb=(cos(θ_curr)*sin(θ_iʳ))*sin(ϕ_curr)
            alphd=alphb - (sin(θ_curr)*cos(θ_iʳ))
            betad=sin(θ_iʳ)*cos(ϕ_curr)

            nud=sqrt(alphd^2*(leaf.a_maj/2)^2+betad^2*(leaf.b_min/2)^2)* k₀/π
            nub=sqrt(alphb^2*(leaf.a_maj/2)^2+betad^2*(leaf.b_min/2)^2)* k₀/π

            zetad=2π*nud
            zetab=2π*nub
            swigd=cnst3 * (besselj1(zetad) / zetad)
            swigb=cnst3 * (besselj1(zetab) / zetab)
            
            # calculate integrands
            sthcp2=(sin(θ_curr)*cos(ϕ_curr))^2
            ss2=((sin(θ_curr)*cos(θ_iʳ))*sin(ϕ_curr)-(cos(θ_curr)*sin(θ_iʳ)))^2
            scsph2=(sin(θ_curr)*cos(θ_iʳ))^2*sin(ϕ_curr)^2
            scs2 = (cos(θ_curr)*sin(θ_iʳ))^2 - scsph2
            scs1 = scsph2+(cos(θ_curr)*sin(θ_iʳ))^2
            sscci = cos(θ_iʳ)^2-sin(θ_iʳ)^2
            sccs = (cos(θ_curr)*sin(θ_iʳ))^2-scsph2
            
            cnt=cnti*cntj

            sumx[3] += (abs(1.0-eb*sthcp2))^2*cnt*swigd^2
            sumx[2] += sthcp2*ss2*cnt*swigd^2
            sumx[1] += (abs(1.0-eb*ss2))^2*cnt*swigd^2

            sjdr[1] +=  (abs(sscci+eb*sccs))^2*cnt*swigb^2 # vv
            sjdr[2] +=  (abs(1.0-eb*sthcp2))^2*cnt*swigb^2 # hh

            sjvh1=sthcp2*scs1*cnt*swigb^2+sjvh1
            sjvh3=sthcp2*scs2*cnt*swigb^2+sjvh3

            cntj=1.0
        end

        sm   += sumx * pdf
        sidr += sjdr * pdf 

		sivh1+=sjvh1*pdf
		sivh3+=sjvh3*pdf
		cnti=1.0
    end

    sgbd = a2*leaf.t^2 * sm * Δθ * Δϕ
    sgbd[2] *= b2

    sgbdr  =  a2*leaf.t^2 * sidr*Δθ*Δϕ
    sgbvh1 =  a2*leaf.t^2 * b2  *abs(sivh1)*Δθ*Δϕ
    sgbvh3 =  a2*leaf.t^2 * b2  *abs(sivh3)*Δθ*Δϕ

    return (sgbd, sgbdr, sgbvh1, sgbvh3)
end

"""
Calculation of average forward scattering amplitudes of a
finite length cylinder, exact series solution
"""
function woodf(wood::Wood)

    nmax= Integer(floor(k₀*wood.r+4.0*(k₀*wood.r)^(1/3)+2.0))

    nmax = min(nmax, 20)

    fsm = [0.0 0.0 ; 0.0 0.0]

    cnti=1.0
	cntj=1.0

    for i = 1:n_θ+1

        θ_curr = (i-1)*δθ

        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,wood.pdf_num,wood.pdf_param)

        (pdf < 1.0e-03) && break

        fsum = [0.0 0.0 ; 0.0 0.0]

        for j = 1:n_ϕ+1

            ϕ_curr = (j-1)*δϕ
            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj

            cnt=cnti*cntj

            cthi =  sin(θ_curr) * sin(θ_iʳ) * cos(ϕ_curr-ϕ_iʳ) - cos(θ_curr) * cos(θ_iʳ)
            tvi  = -sin(θ_curr) * cos(θ_iʳ) * cos(ϕ_curr-ϕ_iʳ) - cos(θ_curr) * sin(θ_iʳ)

            thi = sin(θ_curr) * sin(ϕ_curr-ϕ_iʳ)
            ti=sqrt(tvi^2+thi^2)
            cphi = (sin(θ_iʳ)*cos(θ_curr)*cos(ϕ_curr-ϕ_iʳ)+sin(θ_curr)*cos(θ_iʳ))/sqrt(1.0-cthi^2)

            tvs=-tvi
            ths=thi
            ts=sqrt(tvs^2+ths^2)

            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)
            
            thetai=acos(cthi)
            thetas=π-thetai
            phyi=acos(cphi)
            phys=phyi

            fp = scat(nmax, wood, thetai, thetas, phyi, phys)

            # SCATTERING OUTPUTS MATCH # 

            dsi=ti*ts  

            fvv = tvs*fp[1,1]*tvi-tvs*fp[1,2]*thi-ths*fp[2,1]*tvi+ths*fp[2,2]*thi
            fvh = tvs*fp[1,1]*thi+tvs*fp[1,2]*tvi-ths*fp[2,1]*thi-ths*fp[2,2]*tvi
            fhv = ths*fp[1,1]*tvi+tvs*fp[2,1]*tvi-ths*fp[1,2]*thi-tvs*fp[2,2]*thi
            fhh = ths*fp[1,1]*thi+tvs*fp[2,1]*thi+ths*fp[1,2]*tvi+tvs*fp[2,2]*tvi

            fsum = cnt*[fvv fvh ; fhv fhh]/dsi + fsum

            cntj=1.0

        end

        fsm += fsum * pdf
        cnti=1.0

    end

    return fsm * Δϕ * Δθ

end

function backscattering(θ_curr, ϕ_curr, nmax, sign_i, new_tvs, new_thetas, wood::Wood)

    cthi=sin(θ_curr)*sin(θ_iʳ)*cos(ϕ_curr-ϕ_iʳ)-sign_i*cos(θ_curr)*cos(θ_iʳ)
    tvi=-sin(θ_curr)*cos(θ_iʳ)*cos(ϕ_curr-ϕ_iʳ)-sign_i*cos(θ_curr)*sin(θ_iʳ)
    thi=sign_i*sin(θ_curr)*sin(ϕ_curr-ϕ_iʳ)
    ti=sqrt(tvi^2+thi^2)
    cphi=(sin(θ_iʳ)*cos(θ_curr)*cos(ϕ_curr-ϕ_iʳ)+sign_i*sin(θ_curr)*cos(θ_iʳ))/sqrt(1.0-cthi^2)

    ths=-thi
    tvs= new_tvs(tvi)
    ts=sqrt(ths^2+tvs^2)

    cthi = max(min(cthi, 1.0), -1.0)
    cphi = max(min(cphi, 1.0), -1.0)

    thetai=acos(cthi)
    thetas=new_thetas(thetai)
    phyi=acos(cphi)
    phys=phyi+π

    fp = scat(nmax, wood, thetai, thetas, phyi, phys)

    dsi=ti*ts  
    fvv = tvs*fp[1,1]*tvi-tvs*fp[1,2]*thi-ths*fp[2,1]*tvi+ths*fp[2,2]*thi
    fvh = tvs*fp[1,1]*thi+tvs*fp[1,2]*tvi-ths*fp[2,1]*thi-ths*fp[2,2]*tvi
    fhv = ths*fp[1,1]*tvi+tvs*fp[1,2]*tvi-ths*fp[2,1]*thi-tvs*fp[2,2]*thi
    fhh = ths*fp[1,1]*thi+tvs*fp[1,2]*thi+ths*fp[2,1]*tvi+tvs*fp[2,2]*tvi

    return dsi, fvv, fvh, fhv, fhh

end

"""
Calculation of average backscatter cross-section of a finite length cylinder, exact series solution
(KARAM, FUNG, AND ANTAR, 1988)
"""
function woodb(wood::Wood)

    nmax=min(20, Integer(floor(k₀*wood.r+4.0*(k₀*wood.r)^(1/3)+2.0)))

    smd = [0.0, 0.0, 0.0]
    smdr = [0.0, 0.0]

    smvh1=0.0
    smvh3=0.0

    cnti=1.0
	cntj=1.0

    for i = 1:n_θ+1

        θ_curr = (i-1) * δθ

        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,wood.pdf_num,wood.pdf_param)

        (pdf < 1.0e-03) && break

        sumd = [0.0, 0.0, 0.0]
        sumdr = [0.0, 0.0]

        sumvh1=0.0
        sumvh3=0.0

        for j = 1:n_ϕ+1

            ϕ_curr= (j-1) * δϕ

            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj
            cnt = cnti*cntj

            dsi, fvv, fvh, fhv, fhh = backscattering(θ_curr, ϕ_curr, nmax, 1, x -> -x, x -> x, wood)

            sumd += abs.([fvv ; fvh ; fhh] / dsi).^2*cnt

            ############### 

            dsi, fvv, fvh, fhv, fhh = backscattering(θ_curr, ϕ_curr, nmax, 1, x -> x, x -> π-x, wood)

            sumdr += abs.([fvv ; fhh] / dsi).^2*cnt
            sumvh1 = abs(fvh/(dsi))^2*cnt + sumvh1
            fvhc=conj(fvh)
            
            ################## 

            dsi, fvv, fvh, fhv, fhh = backscattering(θ_curr, ϕ_curr, nmax, -1, x -> -x, x -> π-x, wood)

            sumvh3 = abs(fvh*fvhc/(dsi*dsi))*cnt + sumvh3
            cntj = 1.0
        end

        smd += sumd * pdf
        smdr += sumdr * pdf
        smvh1 = sumvh1*pdf + smvh1
        smvh3 = sumvh3*pdf + smvh3
        cnti = 1.0

    end

    sbd   = 4π * smd   * Δϕ * Δθ
    sbdr  = 4π * smdr  * Δϕ * Δθ
    sbvh1 = 4π * smvh1 * Δϕ * Δθ
    sbvh3 = 4π * smvh3 * Δϕ * Δθ

    return sbd, sbdr, sbvh1, sbvh3

end

function scat(nmax, wood_c::Wood, θ_i, θ_s, phyi, phys)

    k02=k₀^2
    cos_θ_i=cos(θ_i)
    cos_θ_s=cos(θ_s)
    sths=sqrt(1.0-cos_θ_s^2)

    # Equation 26 (Karam, Fung, Antar)
    argum = k₀*(wood_c.l/2)*(cos_θ_i+cos_θ_s)
    μ_si = (argum == 0.0) ? 1.0 : sin(argum)/argum

    z0, a0, b0, e0h, e0v, eta0h, eta0v = cylinder(0, wood_c, θ_i, θ_s)

    s0 = [e0v*(b0*cos_θ_s*cos_θ_i-sths*z0) 0 ; 0 b0*eta0h]
    s = zeros(Complex, 2,2)

    for n = 1:nmax

        zn,an,bn,enh,env,etanh,etanv = cylinder(n, wood_c, θ_i, θ_s)

        sphn=0.0
        π_δ = 0.0001
        cphn = (π - π_δ < phys-phyi < π + π_δ) ? (-1)^n : 1.0
        
        s[1,1] += +2.0*((env*cos_θ_i*bn-ej*etanv*an)*cos_θ_s-sths*env*zn)*cphn
        s[1,2] +=     +((enh*cos_θ_i*bn-ej*etanh*an)*cos_θ_s-sths*enh*zn)*sphn
        s[2,1] +=     +(etanv*bn+ej*env*an*cos_θ_i)*sphn
        s[2,2] += +2.0*(etanh*bn+ej*enh*an*cos_θ_i)*cphn	
        
    end

    f = (s0 + s) * k02 * μ_si * (wood_c.l/2) * (wood_c.ϵ-1.0)
    f[1,2] *= 2 * μ_si* ej
    f[2,1] *= 2 * μ_si* ej

    return f

end

"""
Compute the coefficients of the series expansion of cylinder scattering amplitude.

Equation 24 in 
Electromagnetic Wave Scattering from Some Vegetation Samples
(KARAM, FUNG, AND ANTAR, 1988)
http://www2.geog.ucl.ac.uk/~plewis/kuusk/karam.pdf
"""
function cylinder(n::Integer, wood_c::Wood, θ_i::Real, θ_s::Real)

    # Constants
    a = wood_c.r
    ϵᵣ = wood_c.ϵ
    cos_θ_i  = cos(θ_i)
    cos_θ_s  = cos(θ_s)
    coef1    = sqrt(ϵᵣ-cos_θ_i^2)	
    sin_θ_i  = sqrt( 1-cos_θ_i^2)
    
    # Last equations in 24
    u  = k₀ * a * coef1
    vᵢ = k₀ * a * sqrt(1-cos_θ_i^2)
    vₛ = k₀ * a * sqrt(1-cos_θ_s^2)

    # Bessel functions and derivatives at u, v
    Jₙu=besselj(n,u)
    Hₙv=besselj(n,vᵢ) - ej*bessely(n,vᵢ)
    d_Jₙu=dbessj(n,u)
    d_Hₙv=dbessj(n,vᵢ) - ej*dbessy(n,vᵢ)
    
    # Equation 26 (Karam, Fung, Antar)
    zn, znp1, znm1 = zeta(n, a, u, vₛ)

    # Equation 28 (Karam, Fung, Antar)
    aₙ = (znm1-znp1)/(2*coef1)
    bₙ = (znm1+znp1)/(2*coef1)

    # Compute Rₙ with sub-expressions
    term1   = (d_Hₙv/(vᵢ*Hₙv)-d_Jₙu/(u*Jₙu))
    term2   = (d_Hₙv/(vᵢ*Hₙv)-ϵᵣ*d_Jₙu/(u*Jₙu))
    term3   = (1.0/(vᵢ^2)-1.0/(u^2))*n*cos_θ_i
    Rₙ      = (π*vᵢ^2*Hₙv/2) * (term1*term2-term3^2)

    # Calculate remaining terms, using above sub-expressions
    e_nv   = ej*sin_θ_i*term1/(Rₙ*Jₙu)
    e_nh   = -sin_θ_i*term3/(Rₙ*Jₙu)	
    e_tanv = -e_nh
    e_tanh = ej*sin_θ_i*term2/(Rₙ*Jₙu)

    return (zn, aₙ, bₙ, e_nh, e_nv, e_tanh, e_tanv)
end

"""
Zeta function described in Equation 26 of 

Electromagnetic Wave Scattering from Some Vegetation Samples
(KARAM, FUNG, AND ANTAR, 1988)
http://www2.geog.ucl.ac.uk/~plewis/kuusk/karam.pdf
"""
function zeta(n::Integer, a, u, vs)

    bjnu=besselj(n,u)
    bjnv=besselj(n,vs)

    bjnp1u=besselj(n+1,u)
    bjnp1v=besselj(n+1,vs)

    bjnp2v=besselj(n+2,vs)
    bjnp2u=besselj(n+2,u)

    bjnm1u=besselj(n-1,u)
	bjnm1v=besselj(n-1,vs)

    znp1 = (a^2 / (u^2-vs^2)) * (u * bjnp1v * bjnp2u - vs * bjnp1u * bjnp2v) 
    zn   = (a^2 / (u^2-vs^2)) * (u * bjnv   * bjnp1u - vs * bjnu   * bjnp1v) 
    znm1 = (a^2 / (u^2-vs^2)) * (u * bjnm1v * bjnu   - vs * bjnm1u * bjnv)   

    return zn, znp1, znm1
end

dbessj(n, x) = -n/x*besselj(n,x)+besselj(n-1,x)
dbessy(n, x) = -n/x*bessely(n,x)+bessely(n-1,x)

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

URL=https://www.researchgate.net/publication/3202927_Semi-empirical_model_of_the_
ensemble-averaged_differential_Mueller_matrix_for_microwave_backscattering_from_bare_soil_surfaces
"""
function grdoh(ksig)

    r0 = abs((1-sqrt(ϵ_g))/(1+sqrt(ϵ_g)))^2
    g=0.7*(1-exp(-0.65*ksig^1.8))
    q=0.23*(1-exp(-ksig))*sqrt(r0)
    sp1=(2*θ_iʳ/π)^(1/(3*r0))
    sp=1-sp1*exp(-ksig)

    rgh  = (cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))
    rgv  = (ϵ_g*cos(θ_iʳ)-sqrt(ϵ_g-sin(θ_iʳ)^2))/(ϵ_g*cos(θ_iʳ)+sqrt(ϵ_g-sin(θ_iʳ)^2))

    rh0=abs(rgh)^2
    rv0=abs(rgv)^2

    sighh=g*sp*cos(θ_iʳ)^3*(rh0+rv0)
    sigvv=g*cos(θ_iʳ)^3*(rh0+rv0)/sp
    sigvh=q*sigvv	

    # 3.1.5
    shhg=sighh*exp(-4*imag(K_hc*d_c + K_ht*d_t)) # Factor of 2 pulled out in exp
    svvg=sigvv*exp(-4*imag(K_vc*d_c + K_vt*d_t)) # Factor of 2 pulled out in exp
    svhg=sigvh*exp(-2*imag((K_hc+K_vc)*d_c + (K_ht+K_vt)*d_t))

    gd = 10 * log10.([sigvv sigvh sighh])

    return shhg, svvg, svhg, gd 

end