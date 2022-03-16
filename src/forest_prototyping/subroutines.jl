function afsal(thetir, ai1, l_c::leaf, dat_c::data)
    vol = π * (l_c.amaj * l_c.bmin) * l_c.t / 4.0
    beta = dat_c.ak0^2 * vol / (4.0 * π)
    xi = l_c.epsl-1.0
    xii = xi/(2.0 * (1.0 + xi))
    afhhl = beta * xi * (1.0 - xii * ai1)
    si1 = sin(thetir) ^ 2
    ci2 = cos(thetir) ^ 2
    afvvl = beta * xi * (1.0 - xii * (2.0 * si2 + (ci2 - 2.0 * si2) * ai1))
    return (afhhl, afvvl)
end

function asal(thetir, ntype, parm, l_c::leaf, dat_c::data, i_c::integ)
    exi = l_c.epsl - 1.0
    ea = dat_c.ak0^2 * exi / (sqrt(4.0 * π))
    a2 = abs(ea)^2
    eb = exi / (1.0 + exi)
    b2 = abs(eb) ^ 2

    cthi  = cos(thetir)
    sthi  = sin(thetir)
    cthi2 = cthi^2
    sthi2 = sthi^2
    nth1 = i_c.nth+1
    nph1 = i_c.nph+1
    dth   = pi/i_c.nth
    dph   = 2.0*pi/i_c.nph
    th    = 0.0
    smhh  = 0.0
    smvh  = 0.0
    smvv  = 0.0
    sihhd = 0.0
    sivhd = 0.0
    sivvd = 0.0
    sihhdr= 0.0
    sivvdr= 0.0
    sivh1= 0.0
    sivh3 = 0.0
    cnst1 = dat_c.ak0/pi
    cnst2 = 2.0*pi
    cnst3 = cnst2*l_c.amaj*l_c.bmin/4.0

    cnti=1.0
	cntj=1.0

    for i = 1:nth1
        th = (i - 1) * dth
        if (i == 1 || i == nth1) 
            cnti = .5
        end
        pdf = prob(th,ntype,parm)
        if (pdf < 1.0e-03) 
            break
        end 
        cth=cos(th)
        sth=sin(th)
        cs=cth*sthi
        sc=sth*cthi
        cs2=cs^2
        sc2=sc^2
        ph=0.0
        sumhh=0.0
        sumvh=0.0
        sumvv=0.0
        sjhhd=0.0
        sjvhd=0.0
        sjvvd=0.0
        sjhhdr=0.0
        sjvvdr=0.0
        sjvh1=0.0
        sjvh3=0.0

        for j = 1:nph1

            ph=(j-1)*dph
            if (j == 1 || j == nph1) 
                cntj=.5
            end
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
            sumhh=(abs(1.0-eb*sthcp2))^2*cnt*swigd2+sumhh
            sumvh=sthcp2*ss2*cnt*swigd2+sumvh
            sumvv=(abs(1.0-eb*ss2))^2*cnt*swigd2+sumvv
            sjhhdr=(abs(1.0-eb*sthcp2))^2*cnt*swigb2+sjhhdr
            sjvvdr=(abs(sscci+eb*sccs))^2*cnt*swigb2+sjvvdr
            sjvh1=sthcp2*scs1*cnt*swigb2+sjvh1
            sjvh3=sthcp2*scs2*cnt*swigb2+sjvh3

            cntj=1.0
        end

        smhh=sumhh*pdf+smhh
		smvh=sumvh*pdf+smvh
		smvv=sumvv*pdf+smvv
		sihhdr=sjhhdr*pdf+sihhdr
		sivvdr=sjvvdr*pdf+sivvdr
		sivh1=sjvh1*pdf+sivh1
		sivh3=sjvh3*pdf+sivh3
		cnti=1.0
    end

    t2=l_c.t^2
    delph=dph/(2*pi)
    delth=dth/(pi)

    sgbhhd=a2*t2*smhh*delth*delph    
    sgbvhd=a2*b2*t2*smvh*delth*delph
    sgbvvd=a2*t2*smvv*delth*delph
    sgbhdr=a2*t2*sihhdr*delth*delph
    sgbvdr=a2*t2*sivvdr*delth*delph
    sgbvh1=a2*t2*b2*abs(sivh1)*delth*delph
    sgbvh3=a2*t2*b2*abs(sivh3)*delth*delph

    return (sgbhhd, sgbvhd, sgbvvd, sgbhdr, sgbvdr, sgbvh1, sgbvh3)
end

function woodf(index, theti, phii, ntype, parm, dat_c::data, i_c::integ, p1_c::parm1, p2_c::parm2)

    shh = 0.0
    svh = 0.0
    shv = 0.0
    svv = 0.0

    p1_c.ej = 0.0 + 1.0im

    if (index == 1) 
        p1_c.r0 = dat_c.radbm
        p1_c.h=0.5*dat_c.lb
        p1_c.epsi=dat_c.epsb
    end
    if (index == 2)
        p1_c.r0=dat_c.radtm
        p1_c.h=0.5*dat_c.lt
        p1_c.epsi=dat_c.epst
    end

    nmax= Integer(floor(dat_c.ak0*p1_c.r0+4.0*(dat_c.ak0*p1_c.r0)^0.33333+2.000001))

    nmax = (nmax > 20) ? 20 : nmax
    ci=cos(theti)
	si=sin(theti)

    thets = π - theti
    phi = 0.0
    phs = π
    cs = cos(thets)
    ss = sin(thets)
    nth1 = i_c.nth+1
    nph1 = i_c.nph+1
    dth = π/i_c.nth
    dph = 2.0*π/i_c.nph
    delph = dph/(2.0*π)
    delth = dth/π
    th    = 0.0
    fsmhh  = 0.0
    fsmvh  = 0.0
    fsmvv = 0.0
    fsmhv = 0.0
    cnti=1.0
	cntj=1.0

    for i = 1:nth1
        th = (i-1)*dth

        if (i == 1 || i == nth1)
            cnti = 0.5
        end
        pdf = prob(th,ntype,parm)
        if (pdf < 1.0e-03) 
            break
        end 

        cth=cos(th)
        sth=sin(th)
        ph=0.0
        fsumhh=0.0
        fsumvh=0.0
        fsumvv=0.0
        fsumhv=0.0

        for j = 1:nph1
            ph=(j-1)*dph
            if (j == 1 || j == nph1) 
                cntj=.5
            end
            sph=sin(ph-phi)
            cph=cos(ph-phi)
            cnt=cnti*cntj

            cthi=+sth*si*cph-cth*ci
            tvi=-sth*ci*cph-cth*si
            thi=sth*sph
            ti=sqrt(tvi^2+thi^2)
            cphi=(+si*cth*cph+sth*ci)/sqrt(1.0-cthi^2)

            tvs=-tvi
            ths=thi
            ts=sqrt(tvs^2+ths^2)

            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)
            
            p2_c.thetai=acos(cthi)
            p2_c.thetas=pi-p2_c.thetai
            p2_c.phyi=acos(cphi)
            p2_c.phys=p2_c.phyi

            fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c)

            # SCATTERING OUTPUTS MATCH # 

            dsi=ti*ts  
            fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
            fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
            fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
            fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi


            fsumvv=cnt*fvv/dsi + fsumvv
            fsumvh=cnt*fvh/dsi + fsumvh
            fsumhv=cnt*fhv/dsi + fsumhv
            fsumhh=cnt*fhh/dsi + fsumhh

            cntj=1.0

        end

        fsmhh = fsumhh*pdf + fsmhh
        fsmvh = fsumvh*pdf + fsmvh
        fsmhv = fsumhv*pdf + fsmhv
        fsmvv = fsumvv*pdf + fsmvv
        cnti=1.0		
    end

    shh= fsmhh*delph*delth
    svh= fsmvh*delph*delth
    shv= fsmhv*delph*delth
    svv= fsmvv*delph*delth

    return (shh, svh, shv, svv)

end

function woodb(index,theti,phii,ntype,parm, dat_c::data, i_c::integ, p1_c::parm1, p2_c::parm2)

    ci=cos(theti)
    si=sin(theti)
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

    nmax=min(20, Integer(floor(dat_c.ak0*p1_c.r0+4.0*(dat_c.ak0*p1_c.r0)^0.33333+2.000001)))

    nth1=i_c.nth+1
    nph1=i_c.nph+1
    dth   = pi/i_c.nth
    dph   = 2.0*pi/i_c.nph
    delph=dph/(2.0*pi)
    delth=dth/pi
    th    = 0.0
    smhhd = 0.0
    smvhd = 0.0
    smvvd=0.0
    smhhdr=0.0
    smvvdr=0.0
    smhhdr2=0.0
    smvvdr2=0.0
    smvh1=0.0
    smvh3=0.0

    cnti=1.0
	cntj=1.0

    for i = 1:nth1

        th = (i-1)*dth

        if (i == 1 || i == nth1)
            cnti = 0.5
        end
        pdf = prob(th,ntype,parm)
        if (pdf < 1.0e-03) 
            break
        end 

        cth=cos(th)
        sth=sin(th)
        ph=0.0

        sumhhd = 0.0
        sumvhd = 0.0
        sumvvd=0.0
        sumhhdr=0.0
        sumvvdr=0.0
        sumhhdr2=0.0
        sumvvdr2=0.0
        sumvh1=0.0
        sumvh3=0.0

        for j = 1:nph1

            ph=(j-1)*dph
            if (j == 1 || j == nph1) 
                cntj=.5
            end

            cnt = cnti*cntj

            phi=0.0
            sph=sin(ph-phi)
            cph=cos(ph-phi)
            cthi=sth*si*cph-cth*ci
            tvi=-sth*ci*cph-cth*si
            thi=sth*sph
            ti=sqrt(tvi^2+thi^2)
            cphi=(si*cth*cph+sth*ci)/sqrt(1.0-cthi^2)

            ths=-thi
            tvs=-tvi
            ts=sqrt(ths^2+tvs^2)

            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)

            p2_c.thetai=acos(cthi)
            p2_c.thetas=p2_c.thetai
            p2_c.phyi=acos(cphi)
            p2_c.phys=p2_c.phyi+pi


            fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c)

            dsi=ti*ts  
            fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
            fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
            fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
            fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

            sumhhd=abs(fhh/(dsi))^2*cnt + sumhhd
            sumvvd=abs(fvv/(dsi))^2*cnt + sumvvd
            sumvhd=abs(fvh/(dsi))^2*cnt + sumvhd

            ############### 

            phi=0.0
            sph=sin(ph-phi)
            cph=cos(ph-phi)		
            cthi=sth*si*cph-cth*ci
            tvi=-sth*ci*cph-cth*si
            thi=sth*sph
            ti=sqrt(tvi^2+thi^2)
            cphi=(si*cth*cph+sth*ci)/sqrt(1.0-cthi^2)
            tvs=tvi
            ths=-thi
            ts=sqrt(tvs^2+ths^2)
            
            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)            

            p2_c.thetai=acos(cthi)
            p2_c.thetas=pi-p2_c.thetai
            p2_c.phyi=acos(cphi)
            p2_c.phys=p2_c.phyi+pi
            
            fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c)
            
            dsi=ti*ts  
            fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
            fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
            fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
            fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

            sumhhdr= abs(fhh/(dsi))^2*cnt + sumhhdr
            sumvvdr= abs(fvv/(dsi))^2*cnt + sumvvdr
            sumvh1= abs(fvh/(dsi))^2*cnt + sumvh1
            fvhc=conj(fvh)

            
            ################## 

            phi=0.0
            sph=sin(ph-phi)
            cph=cos(ph-phi)
            cthi=sth*si*cph+cth*ci
            tvi=-sth*ci*cph+cth*si
            thi=-sth*sph
            ti=sqrt(tvi^2+thi^2)
            cphi=(si*cth*cph-sth*ci)/sqrt(1.0-cthi^2)
            tvs=-tvi
            ths=-thi
            ts=sqrt(tvs^2+ths^2)

            cthi = max(min(cthi, 1.0), -1.0)
            cphi = max(min(cphi, 1.0), -1.0)    

            p2_c.thetai=acos(cthi)
            p2_c.thetas=pi-p2_c.thetai
            p2_c.phyi=acos(cphi)
            p2_c.phys=p2_c.phyi+pi

            fpvv,fpvh,fphv,fphh = scat(nmax, dat_c, p1_c, p2_c, i = i)

            dsi=ti*ts 
            fvv = tvs*fpvv*tvi-tvs*fpvh*thi-ths*fphv*tvi+ths*fphh*thi
            fvh = tvs*fpvv*thi+tvs*fpvh*tvi-ths*fphv*thi-ths*fphh*tvi
            fhv = ths*fpvv*tvi+tvs*fphv*tvi-ths*fpvh*thi-tvs*fphh*thi
            fhh = ths*fpvv*thi+tvs*fphv*thi+ths*fpvh*tvi+tvs*fphh*tvi

            sumvh3 = abs(fvh*fvhc/(dsi*dsi))*cnt + sumvh3
            sumvvdr2 = abs(fvv/(dsi))^2*cnt + sumvvdr2
            sumhhdr2 = abs(fhh/(dsi))^2*cnt + sumhhdr2
            cntj=1.0
        end

        smhhd = sumhhd*pdf + smhhd
        smvhd = sumvhd*pdf + smvhd
        smvvd = sumvvd*pdf + smvvd
        smhhdr= sumhhdr*pdf+ smhhdr
        smvvdr= sumvvdr*pdf+ smvvdr
        smhhdr2= sumhhdr2*pdf+ smhhdr2
        smvvdr2= sumvvdr2*pdf+ smvvdr2
        smvh1 = sumvh1*pdf + smvh1
        smvh3 = sumvh3*pdf + smvh3
        cnti=1.0

    end

    sbhhd = 4.0*pi*smhhd*delph*delth
    sbvhd = 4.0*pi*smvhd*delph*delth
    sbvvd = 4.0*pi*smvvd*delph*delth
    sbhhdr= 4.0*pi*smhhdr*delph*delth
    sbvvdr= 4.0*pi*smvvdr*delph*delth
    sbhhdr2= 4.0*pi*smhhdr2*delph*delth
    sbvvdr2= 4.0*pi*smvvdr2*delph*delth
    sbvh1 = 4.0*pi*smvh1*delph*delth
    sbvh3 = 4.0*pi*smvh3*delph*delth

    return sbhhd, sbvhd, sbvvd, sbhhdr, sbvvdr, sbvh1, sbvh3

end

function scat(nmax, dat_c::data, p1_c::parm1, p2_c::parm2; i=-1)

    k02=dat_c.ak0^2
    cthi=cos(p2_c.thetai)
    cphi=cos(p2_c.phyi)
    cths=cos(p2_c.thetas)
    cphs=cos(p2_c.phys)
    cthi2=cthi^2
    cths2=cths^2

    sthi=sqrt(1.0-cthi2)
    sths=sqrt(1.0-cths2)

    u = dat_c.ak0*p1_c.r0*sqrt(p1_c.epsi-cthi2)
    vi = dat_c.ak0*p1_c.r0*sqrt(1.0-cthi2)
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

        if (p2_c.phyi == p2_c.phys)
            cphn = 1.0
        end 
        sphn=0.0
        π_δ = 0.0001
        if (π - π_δ < p2_c.phys-p2_c.phyi < π + π_δ) 
            cphn=(-1.0)^n
        end

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
    cphi=cos(p2_c.phyi)
    cths=cos(p2_c.thetas)
    cphs=cos(p2_c.phys)

    cthi2=cthi^2
    coef1= sqrt(p1_c.epsi-cthi2)	
    u = dat_c.ak0*p1_c.r0*coef1
    vi = dat_c.ak0*p1_c.r0*sqrt(1.0-cthi2)
    vs = dat_c.ak0*p1_c.r0*sqrt(1.0-(cths^2))
    sthi=sqrt(1.0-cthi2)

    bjnu=besselj(n,u)
    dbjnu=dbessj(n,u)
    bjnp1u=besselj(n+1,u)
    bjnv=besselj(n,vi)
    dbjnv=dbessj(n,vi)
    bjnp1v=besselj(n+1,vi)
    hnv=besselj(n,vi)-p1_c.ej*bessely(n,vi)
    dhnv=dbessj(n,vi)-p1_c.ej*dbessy(n,vi)

    zn,znp1,znm1 = zeta(n, u, vs)

    zn=zn*p1_c.r0^2
    znp1=znp1*p1_c.r0^2
    znm1=znm1*p1_c.r0^2
    an = (znm1-znp1)/(2.0*coef1)
    bn = (znm1+znp1)/(2.0*coef1)

    rn1=0.5*pi*(vi^2)*hnv
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
    bjnp1u=besselj(n+1,u)
    bjnv=besselj(n,vs)
    bjnp1v=besselj(n+1,vs)

    bjnp2v=besselj(n+2,vs)
    bjnp2u=besselj(n+2,u)

    if n == 0
        zn=(u*bjnv*bjnp1u-vs*bjnu*bjnp1v)/(u^2-vs^2)
		znp1=(u*bjnp1v*bjnp2u-vs*bjnp1u*bjnp2v)/(u^2-vs^2)
		znm1= znp1
        return zn, znp1, znm1
    elseif n == 1
        bj = besselj0(u)
        by = bessely0(u)
        bjnm1u = bj

        by = bessely0(vs)
        bj = besselj0(vs)
        bjnm1v = bj
    else 
        bjnm1u=besselj(n-1,u)
		bjnm1v=besselj(n-1,vs)
    end

    zn=(u*bjnv*bjnp1u-vs*bjnu*bjnp1v)/(u^2-vs^2)
    znp1=(u*bjnp1v*bjnp2u-vs*bjnp1u*bjnp2v)/(u^2-vs^2)
    znm1=(u*bjnm1v*bjnu-vs*bjnm1u*bjnv)/(u^2-vs^2)

    return zn, znp1, znm1
end

function dbessj(n, x)

    if n == 0 
        dbessj=-besselj(1,x)
    elseif n == 1
        dbessj = besselj(0,x)-besselj(1,x)/x
    else 
        dbessj=-n/x*besselj(n,x)+besselj(n-1,x)
    end
    return dbessj
end

function dbessy(n, x)

    if n == 0 
        dbessy=-bessely(1,x)
    elseif n == 1
        dbessy = bessely(0,x)-bessely(1,x)/x
    else 
        dbessy=-n/x*bessely(n,x)+bessely(n-1,x)
    end
    return dbessy

end