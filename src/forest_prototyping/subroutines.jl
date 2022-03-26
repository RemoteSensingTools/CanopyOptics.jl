"""
Calculate forward scattering amplitudes
"""
function afsal(θ_i_rad, ai1, l_c::leaf, dat_c::data)
    vol = π * (l_c.amaj * l_c.bmin) * l_c.t / 4.0
    beta = dat_c.ak0^2 * vol / (4.0 * π)
    xi = l_c.epsl-1.0
    xii = xi/(2.0 * (1.0 + xi))
    afhhl = beta * xi * (1.0 - xii * ai1)
    afvvl = beta * xi * (1.0 - xii * (2.0 * sin(θ_i_rad)^2 + (cos(θ_i_rad) ^ 2 - 2.0 * sin(θ_i_rad)^2) * ai1))
    return (afhhl, afvvl)
end

"""
Calculate scattering coefficients for a thin disk
"""
function asal(thetir, ntype, parm, l_c::leaf, dat_c::data, i_c::integ)

    # Needed constants
    exi = l_c.epsl - 1.0
    ea = dat_c.ak0^2 * exi / (sqrt(4.0 * π))
    a2 = abs(ea)^2
    eb = exi / (1.0 + exi)
    b2 = abs(eb) ^ 2

    cthi  = cos(thetir)
    sthi  = sin(thetir)
    cthi2 = cthi^2
    sthi2 = sthi^2

    δθ      = π/i_c.nth
    δϕ      = 2.0*π/i_c.nph
    θ_curr  = 0.0

    sm = [0.0, 0.0, 0.0]
    sidr = [0.0, 0.0]

    sivh1= 0.0
    sivh3 = 0.0
    cnst1 = dat_c.ak0/π
    cnst2 = 2.0*π
    cnst3 = cnst2*l_c.amaj*l_c.bmin/4.0

    cnti=1.0
	cntj=1.0

    # Loop over θ
    for i = 1:i_c.nth+1

        θ_curr = (i - 1) * δθ
        
        cnti = (i == 1 || i == i_c.nth+1) ? 0.5 : cnti

        pdf = prob(θ_curr,ntype,parm)
        pdf < 1.0e-03 && break

        cth = cos(θ_curr)
        sth = sin(θ_curr)
        cs=cth*sthi
        sc=sth*cthi
        cs2=cs^2
        sc2=sc^2
        ph=0.0

        sumx = [0.0, 0.0, 0.0]
        sjdr = [0.0, 0.0]

        sjvh1=0.0
        sjvh3=0.0

        for j = 1:i_c.nph+1

            ph=(j-1)*δϕ

            cntj = (j == 1 || j == i_c.nph+1) ? 0.5 : cntj

            sph=sin(ph)
            cph=cos(ph)
            sph2=sph^2
            # calculate swig--beam pattern
            alphb=cs*sph
            alphd=alphb - sc
            betad=sthi*cph
            alphd2=alphd^2
            alphb2=alphb^2
            betad2=betad^2
            amaj2=(l_c.amaj/2)^2
            bmin2=(l_c.bmin/2)^2
            nud=sqrt(alphd2*amaj2+betad2*bmin2)*cnst1
            nub=sqrt(alphb2*amaj2+betad2*bmin2)*cnst1

            zetad=cnst2*nud
            zetab=cnst2*nub
            swigd=cnst3* (besselj1(zetad) / zetad)
            swigb=cnst3* (besselj1(zetab) / zetab)
            swigd2=swigd^2
            swigb2=swigb^2
            
            # calculate integrands
            sthcp2=(sth*cph)^2
            ss2=(sc*sph-cs)^2
            scsph2=sc2*sph2
            scs2 = cs2 - scsph2
            scs1 = scsph2+cs2
            sscci = cthi2-sthi2
            sccs=cs2-scsph2
            cnt=cnti*cntj

            sumx[3] += (abs(1.0-eb*sthcp2))^2*cnt*swigd2
            sumx[2] += sthcp2*ss2*cnt*swigd2
            sumx[1] += (abs(1.0-eb*ss2))^2*cnt*swigd2

            sjdr[1] +=  (abs(sscci+eb*sccs))^2*cnt*swigb2 # vv
            sjdr[2] +=  (abs(1.0-eb*sthcp2))^2*cnt*swigb2 # hh

            sjvh1=sthcp2*scs1*cnt*swigb2+sjvh1
            sjvh3=sthcp2*scs2*cnt*swigb2+sjvh3

            cntj=1.0
        end

        sm += sumx .* pdf
        sidr += sjdr * pdf 

		sivh1=sjvh1*pdf+sivh1
		sivh3=sjvh3*pdf+sivh3
		cnti=1.0
    end

    t2=l_c.t^2

    delph=δϕ/(2*π)
    delth=δθ/π

    sgbd = a2*t2 * sm * delth * delph
    sgbd[2] *= b2

    sgbdr = a2*t2*sidr *delth*delph
    sgbvh1=a2*t2*b2*abs(sivh1)*delth*delph
    sgbvh3=a2*t2*b2*abs(sivh3)*delth*delph

    return (sgbd, sgbdr, sgbvh1, sgbvh3)
end

function woodf(index, geom_i_rad::IncidentGeometry, ntype, parm, dat_c::data, i_c::integ, p1_c::parm1, p2_c::parm2)

    @unpack θ_i, ϕ_i = geom_i_rad

    p1_c.ej = 0.0 + 1.0im

    if (index == 1) 
        p1_c.r0 = dat_c.radbm
        p1_c.h = 0.5 * dat_c.lb
        p1_c.epsi = dat_c.epsb
    end

    if (index == 2)
        p1_c.r0 = dat_c.radtm
        p1_c.h = 0.5 * dat_c.lt
        p1_c.epsi = dat_c.epst
    end

    nmax= Integer(floor(dat_c.ak0*p1_c.r0+4.0*(dat_c.ak0*p1_c.r0)^(1/3)+2.0))

    nmax = (nmax > 20) ? 20 : nmax
    ci=cos(θ_i)
	si=sin(θ_i)
    phi = 0.0
    δθ = π/i_c.nth
    δϕ = 2.0*π/i_c.nph
    delph = δϕ/(2.0*π)
    delth = δθ/π
    th    = 0.0

    fsm = [0.0 0.0 ; 0.0 0.0]

    cnti=1.0
	cntj=1.0

    for i = 1:i_c.nth+1

        th = (i-1)*δθ

        cnti = (i == 1 || i == i_c.nth+1) ? 0.5 : cnti

        pdf = prob(th,ntype,parm)

        (pdf < 1.0e-03) && break

        cth=cos(th)
        sth=sin(th)
        ph=0.0

        fsum = [0.0 0.0 ; 0.0 0.0]

        for j = 1:i_c.nph+1
            ph=(j-1)*δϕ
            cntj = (j == 1 || j == i_c.nph+1) ? 0.5 : cntj
            sph=sin(ph-phi)
            cph=cos(ph-phi)
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

            fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c, i=i)

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

    return fsm * delph * delth

end

function backscattering(ph, sth, cth, nmax, sign_i, new_tvs, new_thetas, dat_c::data, p1_c::parm1, p2_c::parm2)

    phi=0.0
    sph=sin(ph-phi)
    cph=cos(ph-phi)
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

    fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c)

    dsi=ti*ts  
    fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
    fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
    fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
    fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

    return dsi, fvv, fvh, fhv, fhh

end

function woodb(index,geom_i::IncidentGeometry,ntype,parm, dat_c::data, i_c::integ, p1_c::parm1, p2_c::parm2)

    @unpack θ_i, ϕ_i = geom_i

    p1_c.ej=complex(0.0,1.0)

    if (index == 1)
        p1_c.r0 = dat_c.radbm
        p1_c.h=0.5*dat_c.lb
        p1_c.epsi=dat_c.epsb
    elseif index == 2
        p1_c.r0 = dat_c.radtm
        p1_c.epsi = dat_c.epst
    end

    k02 = dat_c.ak0^2

    nmax=min(20, Integer(floor(dat_c.ak0*p1_c.r0+4.0*(dat_c.ak0*p1_c.r0)^(1/3)+2.0)))

    δθ    = π/i_c.nth
    δϕ   = 2.0*π/i_c.nph
    delph = δϕ/(2.0*π)
    delth = δθ/π
    th    = 0.0

    smd = [0.0, 0.0, 0.0]
    smdr = [0.0, 0.0]

    smvh1=0.0
    smvh3=0.0

    cnti=1.0
	cntj=1.0

    for i = 1:i_c.nth+1

        θ_curr = (i-1) * δθ

        cnti = (i == 1 || i == i_c.nth+1) ? 0.5 : cnti

        pdf = prob(θ_curr,ntype,parm)

        (pdf < 1.0e-03) && break

        cth=cos(θ_curr)
        sth=sin(θ_curr)
        ph=0.0

        sumd = [0.0, 0.0, 0.0]
        sumdr = [0.0, 0.0]

        sumvh1=0.0
        sumvh3=0.0

        for j = 1:i_c.nph+1

            ph=(j-1)*δϕ

            cntj = (j == 1 || j == i_c.nph+1) ? 0.5 : cntj
            cnt = cnti*cntj

            dsi, fvv, fvh, fhv, fhh = backscattering(ph, sth, cth, nmax, 1, x -> -x, x -> x, dat_c, p1_c, p2_c)

            sumd += abs.([fvv ; fvh ; fhh] / dsi).^2*cnt

            ############### 

            dsi, fvv, fvh, fhv, fhh = backscattering(ph, sth, cth, nmax, 1, x -> x, x -> π-x, dat_c, p1_c, p2_c)

            sumdr += abs.([fvv ; fhh] / dsi).^2*cnt
            sumvh1 = abs(fvh/(dsi))^2*cnt + sumvh1
            fvhc=conj(fvh)
            
            ################## 

            dsi, fvv, fvh, fhv, fhh = backscattering(ph, sth, cth, nmax, -1, x -> -x, x -> π-x, dat_c, p1_c, p2_c)

            sumvh3 = abs(fvh*fvhc/(dsi*dsi))*cnt + sumvh3
            cntj=1.0
        end

        smd += sumd * pdf
        smdr += sumdr * pdf
        smvh1 = sumvh1*pdf + smvh1
        smvh3 = sumvh3*pdf + smvh3
        cnti=1.0

    end

    sbd   = 4 * π * smd * delph * delth
    sbdr  = 4 * π * smdr * delph * delth
    sbvh1 = 4 * π * smvh1 * delph * delth
    sbvh3 = 4 * π * smvh3 * delph * delth

    return sbd, sbdr, sbvh1, sbvh3

end

function scat(nmax, dat_c::data, p1_c::parm1, p2_c::parm2; i=-1)

    k02=dat_c.ak0^2
    cthi=cos(p2_c.thetai)
    cths=cos(p2_c.thetas)
    cths2=cths^2

    sths=sqrt(1.0-cths2)
    argum=dat_c.ak0*p1_c.h*(cthi+cths)

    q = (argum == 0.0) ? 1.0 : sin(argum)/argum

    z0, a0, b0, e0h, e0v, eta0h, eta0v = cylinder(0, dat_c, p1_c, p2_c)

    shh0=b0*eta0h
    shv0=0
    svh0=0
    svv0=e0v*(b0*cths*cthi-sths*z0)
    shh=0.0
    shv=0.0
    svv=0.0
    svh=0.0

    for n = 1:nmax

        zn,an,bn,enh,env,etanh,etanv = cylinder(n, dat_c, p1_c, p2_c)

        sphn=0.0
        π_δ = 0.0001
        cphn = (π - π_δ < p2_c.phys-p2_c.phyi < π + π_δ) ? (-1.0)^n : 1.0

        shh=shh+2.0*(etanh*bn+p1_c.ej*enh*an*cthi)*cphn	
        shv=shv+(etanv*bn+p1_c.ej*env*an*cthi)*sphn
        svh=svh+((enh*bn*cthi-p1_c.ej*etanh*an)*cths-sths*enh*zn)*sphn
        svv=svv+2.0*((env*cthi*bn-p1_c.ej*etanv*an)*cths-sths*env*zn)*cphn
        
    end

    fhh=(shh0+shh)*k02*q*p1_c.h*(p1_c.epsi-1.0)
    fhv=2.0*(shv0+shv)*k02*p1_c.h*q*p1_c.ej*(p1_c.epsi-1.0)
    fvh=2.0*(svh0+svh)*k02*p1_c.h*q*p1_c.ej*(p1_c.epsi-1.0)
    fvv=(svv0+svv)*k02*p1_c.h*q*(p1_c.epsi-1.0)

    return fvv,fvh,fhv,fhh

end

function cylinder(n, dat_c::data, p1_c::parm1, p2_c::parm2)

    cthi=cos(p2_c.thetai)
    cths=cos(p2_c.thetas)

    cthi2=cthi^2
    coef1= sqrt(p1_c.epsi-cthi2)	
    u = dat_c.ak0*p1_c.r0*coef1
    vi = dat_c.ak0*p1_c.r0*sqrt(1.0-cthi2)
    vs = dat_c.ak0*p1_c.r0*sqrt(1.0-(cths^2))
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
    coefetanh = (dhnv/(vi*hnv)-p1_c.epsi*dbjnu/(u*bjnu))
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

function grdoh(ksig, θ_i, a_c::a, b_c::b)

    θ_i_rad= deg2rad(θ_i)

    er=sqrt(a_c.epsg)
    r0 = abs((1.0-er)/(1.0+er))^2
    g=0.7*(1.0-exp(-0.65*ksig^1.8))
    q=0.23*(1.0-exp(-ksig))*sqrt(r0)
    sp1=(2.0*θ_i_rad/π)^(1.0/(3.0*r0))
    sp=1.0-sp1*exp(-ksig)

    rgh  = (cos(θ_i_rad)-sqrt(a_c.epsg-sin(θ_i_rad)^2))/(cos(θ_i_rad)+sqrt(a_c.epsg-sin(θ_i_rad)^2))
    rgv  = (a_c.epsg*ci-sqrt(a_c.epsg-sin(θ_i_rad)^2))/(a_c.epsg*ci+sqrt(a_c.epsg-sin(θ_i_rad)^2))

    rh0=abs(rgh)^2
    rv0=abs(rgv)^2

    sighh=g*sp*ci^3*(rh0+rv0)
    sigvv=g*ci^3*(rh0+rv0)/sp
    sigvh=q*sigvv	
    shhg=sighh*exp(-4.0*(a_c.khim1*d1+a_c.khim2*d2))
    svvg=sigvv*exp(-4.0*(a_c.kvim1*d1+a_c.kvim2*d2))
    svhg=sigvh*exp(-2.0*((a_c.khim1+a_c.kvim1)*d1+(a_c.khim2+a_c.kvim2)*d2))

    ghhd=10*log10(sighh)
    gvvd=10*log10(sigvv)
    gvhd=10*log10(sigvh)

    return shhg, svvg, svhg, ghhd, gvvd, gvhd

end