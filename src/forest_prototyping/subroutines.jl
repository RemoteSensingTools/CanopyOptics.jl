"""
Calculate forward scattering amplitudes
"""
function afsal(θ_iʳ, ai1, l_c::leaf)
    vol = π * l_c.amaj * l_c.bmin * l_c.t / 4      # Volume of leaf 
    β = k₀^2 * vol / (4π)                   # ν² * vol / 4π
    xi = l_c.ϵ_l-1
    xii = xi/(2 * (1 + xi))
    afhhl = β * xi * (1 - xii * ai1)
    afvvl = β * xi * (1 - xii * (2 * sin(θ_iʳ)^2 + (cos(θ_iʳ)^2 - 2 * sin(θ_iʳ)^2) * ai1))
    return (afhhl, afvvl)
end

"""
Calculate scattering coefficients for a thin disk
"""
function asal(θ_iʳ, ntype, parm, l_c::leaf)

    # Needed constants
    exi = l_c.ϵ_l - 1.0
    ea = k₀^2 * exi / sqrt(4π)
    a2 = abs(ea)^2
    eb = exi / (1 + exi)
    b2 = abs(eb)^2

    δθ =  π/n_θ
    δϕ = 2π/n_ϕ

    sm = [0.0, 0.0, 0.0]
    sidr = [0.0, 0.0]

    sivh1 = 0.0
    sivh3 = 0.0
    cnst1 = k₀/π
    cnst2 = 2π
    cnst3 = cnst2*l_c.amaj*l_c.bmin/4.0

    cnti=1.0
	cntj=1.0

    # Loop over θ
    for i = 1:n_θ+1

        θ_curr = (i - 1) * δθ
        
        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,ntype,parm)
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

            nud=sqrt(alphd^2*(l_c.amaj/2)^2+betad^2*(l_c.bmin/2)^2)*cnst1
            nub=sqrt(alphb^2*(l_c.amaj/2)^2+betad^2*(l_c.bmin/2)^2)*cnst1

            zetad=cnst2*nud
            zetab=cnst2*nub
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

    Δϕ = δϕ/(2*π)
    Δθ = δθ/π

    sgbd = a2*l_c.t^2 * sm * Δθ * Δϕ
    sgbd[2] *= b2

    sgbdr  =  a2*l_c.t^2 * sidr*Δθ*Δϕ
    sgbvh1 =  a2*l_c.t^2 * b2  *abs(sivh1)*Δθ*Δϕ
    sgbvh3 =  a2*l_c.t^2 * b2  *abs(sivh3)*Δθ*Δϕ

    return (sgbd, sgbdr, sgbvh1, sgbvh3)
end

function woodf(index, ntype, parm, branch_c::branch, trunk_c::trunk, p1_c::parm1, p2_c::parm2)

    p1_c.ej = 0.0 + 1.0im

    if (index == 1) 
        p1_c.r0 = branch_c.radbm
        p1_c.h = 0.5 * branch_c.l_b
        p1_c.ϵ_i = branch_c.ϵ_b
    end

    if (index == 2)
        p1_c.r0 = trunk_c.radtm
        p1_c.h = 0.5 * trunk_c.l_t
        p1_c.ϵ_i = trunk_c.ϵ_t
    end

    nmax= Integer(floor(k₀*p1_c.r0+4.0*(k₀*p1_c.r0)^(1/3)+2.0))

    nmax = min(nmax, 20)

    ci=cos(θ_iʳ)
	si=sin(θ_iʳ)
    
    ϕ_i = 0.0
    δθ = π/n_θ
    δϕ = 2.0*π/n_ϕ
    Δϕ = δϕ/(2π)
    Δθ = δθ/π

    fsm = [0.0 0.0 ; 0.0 0.0]

    cnti=1.0
	cntj=1.0

    for i = 1:n_θ+1

        θ_curr = (i-1)*δθ

        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,ntype,parm)

        (pdf < 1.0e-03) && break

        cth=cos(θ_curr)
        sth=sin(θ_curr)

        fsum = [0.0 0.0 ; 0.0 0.0]

        for j = 1:n_ϕ+1

            ϕ_curr = (j-1)*δϕ
            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj
            sph=sin(ϕ_curr-ϕ_i)
            cph=cos(ϕ_curr-ϕ_i)
            cnt=cnti*cntj

            cthi = +sth * si * cph - cth * ci
            tvi  = -sth * ci * cph - cth * si
            thi=sth*sph
            ti=sqrt(tvi^2+thi^2)
            cphi=(+si*cth*cph+sth*ci)/sqrt(1.0-cthi^2)

            tvs=-tvi
            ths=thi
            ts=sqrt(tvs^2+ths^2)

            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)
            
            p2_c.thetai=acos(cthi)
            p2_c.thetas=π-p2_c.thetai
            p2_c.phyi=acos(cphi)
            p2_c.phys=p2_c.phyi

            fpvv,fpvh,fphv,fphh = scat(nmax, p1_c, p2_c, i=i)

            # SCATTERING OUTPUTS MATCH # 

            dsi=ti*ts  
            fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
            fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
            fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
            fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

            fsum = cnt*[fvv fvh ; fhv fhh]/dsi + fsum

            cntj=1.0

        end

        fsm += fsum * pdf
        cnti=1.0

    end

    return fsm * Δϕ * Δθ

end

function backscattering(ϕ_curr, sth, cth, nmax, sign_i, new_tvs, new_thetas, p1_c::parm1, p2_c::parm2)

    si = sin(θ_iʳ)
    ci = cos(θ_iʳ)

    ϕ_i = 0.0

    sph=sin(ϕ_curr-ϕ_i)
    cph=cos(ϕ_curr-ϕ_i)
    cthi=sth*si*cph-sign_i*cth*ci
    tvi=-sth*ci*cph-sign_i*cth*si
    thi=sign_i*sth*sph
    ti=sqrt(tvi^2+thi^2)
    cphi=(si*cth*cph+sign_i*sth*ci)/sqrt(1.0-cthi^2)

    ths=-thi
    tvs= new_tvs(tvi)
    ts=sqrt(ths^2+tvs^2)

    cthi = max(min(cthi, 1.0), -1.0)
    cphi = max(min(cphi, 1.0), -1.0)

    p2_c.thetai=acos(cthi)
    p2_c.thetas=new_thetas(p2_c.thetai)
    p2_c.phyi=acos(cphi)
    p2_c.phys=p2_c.phyi+π

    fpvv,fpvh,fphv,fphh = scat(nmax, p1_c, p2_c)

    dsi=ti*ts  
    fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
    fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
    fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
    fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

    return dsi, fvv, fvh, fhv, fhh

end

"""
Calculation of average backscatter cross-section of a finite length cylinder, exact series solution
(KARAM, FUNG, AND ANTAR, 1988)
"""
function woodb(index,ntype,parm, branch_c, trunk_c, p1_c::parm1, p2_c::parm2)

    p1_c.ej=complex(0.0,1.0)

    if (index == 1)
        p1_c.r0 = branch_c.radbm
        p1_c.h=0.5*branch_c.l_b
        p1_c.ϵ_i=branch_c.ϵ_b
    elseif index == 2
        p1_c.r0 = trunk_c.radtm
        p1_c.ϵ_i = trunk_c.ϵ_t
    end

    nmax=min(20, Integer(floor(k₀*p1_c.r0+4.0*(k₀*p1_c.r0)^(1/3)+2.0)))

    δθ    = π/n_θ
    δϕ   = 2π/n_ϕ
    Δϕ = 1/n_ϕ # δϕ/(2π)
    Δθ = 1/n_θ

    smd = [0.0, 0.0, 0.0]
    smdr = [0.0, 0.0]

    smvh1=0.0
    smvh3=0.0

    cnti=1.0
	cntj=1.0

    for i = 1:n_θ+1

        θ_curr = (i-1) * δθ

        cnti = (i == 1 || i == n_θ+1) ? 0.5 : cnti

        pdf = prob(θ_curr,ntype,parm)

        (pdf < 1.0e-03) && break

        cth=cos(θ_curr)
        sth=sin(θ_curr)

        sumd = [0.0, 0.0, 0.0]
        sumdr = [0.0, 0.0]

        sumvh1=0.0
        sumvh3=0.0

        for j = 1:n_ϕ+1

            ϕ_curr= (j-1) * δϕ

            cntj = (j == 1 || j == n_ϕ+1) ? 0.5 : cntj
            cnt = cnti*cntj

            dsi, fvv, fvh, fhv, fhh = backscattering(ϕ_curr, sth, cth, nmax, 1, x -> -x, x -> x, p1_c, p2_c)

            sumd += abs.([fvv ; fvh ; fhh] / dsi).^2*cnt

            ############### 

            dsi, fvv, fvh, fhv, fhh = backscattering(ϕ_curr, sth, cth, nmax, 1, x -> x, x -> π-x, p1_c, p2_c)

            sumdr += abs.([fvv ; fhh] / dsi).^2*cnt
            sumvh1 = abs(fvh/(dsi))^2*cnt + sumvh1
            fvhc=conj(fvh)
            
            ################## 

            dsi, fvv, fvh, fhv, fhh = backscattering(ϕ_curr, sth, cth, nmax, -1, x -> -x, x -> π-x, p1_c, p2_c)

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

function scat(nmax, p1_c::parm1, p2_c::parm2; i=-1)

    k02=k₀^2
    cthi=cos(p2_c.thetai)
    cths=cos(p2_c.thetas)

    sths=sqrt(1.0-cths^2)
    argum=k₀*p1_c.h*(cthi+cths)

    q = (argum == 0.0) ? 1.0 : sin(argum)/argum

    z0, a0, b0, e0h, e0v, eta0h, eta0v = cylinder(0, p1_c, p2_c)

    shh0=b0*eta0h
    shv0=0
    svh0=0
    svv0=e0v*(b0*cths*cthi-sths*z0)
    shh=0.0
    shv=0.0
    svv=0.0
    svh=0.0

    for n = 1:nmax

        zn,an,bn,enh,env,etanh,etanv = cylinder(n, p1_c, p2_c)

        sphn=0.0
        π_δ = 0.0001
        cphn = (π - π_δ < p2_c.phys-p2_c.phyi < π + π_δ) ? (-1.0)^n : 1.0

        shh += +2.0*(etanh*bn+p1_c.ej*enh*an*cthi)*cphn	
        shv += +(etanv*bn+p1_c.ej*env*an*cthi)*sphn
        svh += +((enh*bn*cthi-p1_c.ej*etanh*an)*cths-sths*enh*zn)*sphn
        svv += +2.0*((env*cthi*bn-p1_c.ej*etanv*an)*cths-sths*env*zn)*cphn
        
    end

    fhh=(shh0+shh)*k02*q*p1_c.h*(p1_c.ϵ_i-1.0)
    fhv=2.0*(shv0+shv)*k02*p1_c.h*q*p1_c.ej*(p1_c.ϵ_i-1.0)
    fvh=2.0*(svh0+svh)*k02*p1_c.h*q*p1_c.ej*(p1_c.ϵ_i-1.0)
    fvv=(svv0+svv)*k02*p1_c.h*q*(p1_c.ϵ_i-1.0)

    return fvv,fvh,fhv,fhh

end

function cylinder(n, p1_c::parm1, p2_c::parm2)

    cthi=cos(p2_c.thetai)
    cths=cos(p2_c.thetas)

    cthi2=cthi^2
    coef1= sqrt(p1_c.ϵ_i-cthi2)	
    u = k₀*p1_c.r0*coef1
    vi = k₀*p1_c.r0*sqrt(1.0-cthi2)
    vs = k₀*p1_c.r0*sqrt(1.0-(cths^2))
    sthi=sqrt(1.0-cthi2)

    bjnu=besselj(n,u)
    dbjnu=dbessj(n,u)
    hnv=besselj(n,vi)-p1_c.ej*bessely(n,vi)
    dhnv=dbessj(n,vi)-p1_c.ej*dbessy(n,vi)

    zn,znp1,znm1 = zeta(n, u, vs)

    zn=zn*p1_c.r0^2
    znp1=znp1*p1_c.r0^2
    znm1=znm1*p1_c.r0^2
    an = (znm1-znp1)/(2.0*coef1)
    bn = (znm1+znp1)/(2.0*coef1)

    rn1=0.5*π*(vi^2)*hnv
    coefenv=(dhnv/(vi*hnv)-dbjnu/(u*bjnu))
    coefetanh = (dhnv/(vi*hnv)-p1_c.ϵ_i*dbjnu/(u*bjnu))
    coefenh=(1.0/(vi^2)-1.0/(u^2))*n*cthi

    rn=rn1*(coefenv*coefetanh-coefenh^2)

    env=p1_c.ej*sthi*coefenv/(rn*bjnu)
    enh=-sthi*coefenh/(rn*bjnu)	

    etanv=-enh
    etanh=p1_c.ej*sthi*coefetanh/(rn*bjnu)

    return (zn,an,bn,enh,env,etanh,etanv)
end

function zeta(n, u, vs)

    bjnu=besselj(n,u)
    bjnv=besselj(n,vs)

    bjnp1u=besselj(n+1,u)
    bjnp1v=besselj(n+1,vs)

    bjnp2v=besselj(n+2,vs)
    bjnp2u=besselj(n+2,u)

    bjnm1u=besselj(n-1,u)
	bjnm1v=besselj(n-1,vs)

    zn   = (u * bjnv   * bjnp1u - vs * bjnu   * bjnp1v) / (u^2-vs^2)
    znp1 = (u * bjnp1v * bjnp2u - vs * bjnp1u * bjnp2v) / (u^2-vs^2)
    znm1 = (u * bjnm1v * bjnu   - vs * bjnm1u * bjnv)   / (u^2-vs^2)

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
    z < 1.0e-03 && return d
    z > 1.6e02 && return 1.0/dag

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
function grdoh(ksig, a_c::a, ::b)

    er=sqrt(a_c.ϵ_g)
    r0 = abs((1.0-er)/(1.0+er))^2
    g=0.7*(1.0-exp(-0.65*ksig^1.8))
    q=0.23*(1.0-exp(-ksig))*sqrt(r0)
    sp1=(2.0*θ_iʳ/π)^(1.0/(3.0*r0))
    sp=1.0-sp1*exp(-ksig)

    rgh  = (cos(θ_iʳ)-sqrt(a_c.ϵ_g-sin(θ_iʳ)^2))/(cos(θ_iʳ)+sqrt(a_c.ϵ_g-sin(θ_iʳ)^2))
    rgv  = (a_c.ϵ_g*cos(θ_iʳ)-sqrt(a_c.ϵ_g-sin(θ_iʳ)^2))/(a_c.ϵ_g*cos(θ_iʳ)+sqrt(a_c.ϵ_g-sin(θ_iʳ)^2))

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