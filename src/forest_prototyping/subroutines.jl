function afsal(thetir, ai1)
    vol = π * (amaj * bmin) * t / 4.0
    beta = ak0^2 * vol / (4.0 * π)
    xi = epsl-1.0
    xii = xi/(2.0 * (1.0 + xi))
    afhhl = beta * xi * (1.0 - xii * ai1)
    si1 = sin(thetir) ^ 2
    ci2 = cos(thetir) ^ 2
    afvvl = beta * xi * (1.0 - xii * (2.0 * si2 + (ci2 - 2.0 * si2) * ai1))
    return (afhhl, afvvl)
end

function asal(thetir, ntype, parm)
    exi = epsl - 1.0
    ea = ak0^2 * exi / (sqrt(4.0 * π))
    a2 = abs(ea)^2
    eb = exi / (1.0 + exi)
    b2 = abs(eb) ^ 2

    cthi  = cos(thetir)
    sthi  = sin(thetir)
    cthi2 = cthi^2
    sthi2 = sthi^2
    nth1 = nth+1
    nph1 = nph+1
    dth   = pi/nth
    dph   = 2.0*pi/nph
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
    cnst1 = ak0/pi
    cnst2 = 2.0*pi
    cnst3 = cnst2*amaj*bmin/4.0

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
            amaj2=(amaj/2)^2
            bmin2=(bmin/2)^2
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

    t2=t^2
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