"""Independent PROSPECT plate oracle; run with Python + mpmath==1.3.0.

No production Julia formulas are imported. Diffuse attenuation uses angular
quadrature (checked against the integral E3 definition); fractional plate counts
use mpmath's general matrix eigendecomposition, while integers use repeated flux adding.
The CSV is read by Julia tests without a Python runtime dependency.
"""
from pathlib import Path
import csv
import mpmath as mp
mp.mp.dps = 90
D = mp.mpf

def attenuation(k):
    if not k:
        return D(1)
    return 2 * mp.expint(3, k)

def plate(k, talf, t12, nr):
    tau = attenuation(k)
    t21 = t12 / nr**2
    r21 = 1-t21
    denominator = 1-r21**2*tau**2
    t = t12*tau*t21/denominator
    r = 1-t12+r21*tau*t
    ta = talf*tau*t21/denominator
    ra = 1-talf+r21*tau*ta
    return t,r,ta,ra

def budget(k,n,talf,t12,nr):
    t,r,ta,ra = plate(k,talf,t12,nr)
    m = n-1
    if m == int(m) and m <= 8:
        tm,rm = D(1),D(0)
        for _ in range(int(m)):
            tm,rm = t*tm/(1-r*rm),r+t*t*rm/(1-r*rm)
    elif not k:
        tm,rm = t/(t+m*r),m*r/(t+m*r)
    else:
        matrix = mp.matrix([[1,-r],[r,t*t-r*r]])/t
        eigenvalues, vectors = mp.eig(matrix)
        power = vectors * mp.diag([value**m for value in eigenvalues]) * vectors**-1
        tm,rm = 1/power[0,0],power[1,0]/power[0,0]
    trans = ta*tm/(1-rm*r)
    refl = ra+ta*rm*t/(1-rm*r)
    return tuple(mp.re(x) for x in (trans,refl,D(0) if not k else 1-trans-refl))

def fresnel_average(angle,n):
    lower=mp.cos(mp.radians(angle))
    def integrand(mu):
        transmitted=mp.sqrt(1-(1-mu*mu)/(n*n))
        rs=((mu-n*transmitted)/(mu+n*transmitted))**2
        rp=((n*mu-transmitted)/(n*mu+transmitted))**2
        return mu*(2-rs-rp)
    return mp.quad(integrand,[lower,(lower+1)/2,1])/(1-lower*lower)

def main():
    # Validate E3 against angular quadrature without using the recurrence in Julia.
    for k in map(D,['0','1e-20','1e-8','.1','1','10']):
        a=2*mp.quad(lambda mu: mu*(-mp.expm1(-k/mu)),[0,D('.01'),D('.1'),1])
        assert abs((1-attenuation(k))-a) < D('1e-75')
    path=Path(__file__).with_name('prospect_reference.csv')
    with path.open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n')
        w.writerow(['k','N','talf','t12','nr','T','R','A','dTdk','dRdk','dAdk','dTdN','dRdN','dAdN'])
        for ns in ['1','1.2','1.7','2','3','8']:
            n=D(ns)
            for ks in ['0','1e-25','1e-12','1e-8','1e-4','.01','.1','.999999','1','1.000001','5','10','30']:
                k=D(ks); talf=D('.95'); t12=D('.9'); nr=D('1.4')
                values=budget(k,n,talf,t12,nr)
                h=D('1e-30') if k==0 else min(k*D('1e-15'),D('1e-18'))
                dk=tuple((a-b)/h for a,b in zip(budget(k+h,n,talf,t12,nr),values)) if not k else tuple((a-b)/(2*h) for a,b in zip(budget(k+h,n,talf,t12,nr),budget(k-h,n,talf,t12,nr)))
                hn=D('1e-25')
                dn=tuple((a-b)/hn for a,b in zip(budget(k,n+hn,talf,t12,nr),values)) if n==1 else tuple((a-b)/(2*hn) for a,b in zip(budget(k,n+hn,talf,t12,nr),budget(k,n-hn,talf,t12,nr)))
                w.writerow([mp.nstr(x,45) for x in (k,n,talf,t12,nr,*values,*dk,*dn)])
    with Path(__file__).with_name('prospect_fresnel_reference.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['nr','talf','t12'])
        for ns in ['1.25','1.3','1.4','1.5','1.6']:
            n=D(ns);w.writerow([mp.nstr(x,45) for x in (n,fresnel_average(40,n),fresnel_average(90,n))])
    with mp.workdps(750), Path(__file__).with_name('prospect_stack_reference.csv').open('w',newline='') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['t','alpha','m','T','R','A'])
        cases=[('.5','1e-50','1e25'),('.5','1e-30','1e15'),('1e-300','1e-300','.7'),('1e-35','1e-35','.7'),('.5','0','1e30'),('.5','0','1e200'),('.5','1e-4','1000')]
        cases += [(str(D(2)**-149),'.75','.5'),(str(D(2)**-1074),'.75','.5')]
        for ts,aas,ms in cases:
            t,a,m=map(D,(ts,aas,ms));r=1-t-a
            if not a:
                trans=t/(t+m*r);refl=m*r/(t+m*r)
            else:
                matrix=mp.matrix([[1,-r],[r,t*t-r*r]])/t
                eigenvalues,vectors=mp.eig(matrix)
                power=vectors*mp.diag([value**m for value in eigenvalues])*vectors**-1
                trans=1/power[0,0];refl=power[1,0]/power[0,0]
            w.writerow([mp.nstr(mp.re(x),45) for x in (t,a,m,trans,refl,1-trans-refl)])
    print(path)

if __name__=='__main__': main()
