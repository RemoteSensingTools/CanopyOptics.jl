"""90-digit angular attenuation/Fresnel and absorption-source depth quadrature.

Inputs are frozen to Float32/Float64. Unlike the production source, attenuation
and direct loss are independent angular integrals; source strength is a numerical
depth integral. Requires mpmath1.3.0, no production Julia imports.
"""
import csv
import pathlib
import struct
from functools import lru_cache
import mpmath as mp
from generate_prospect_reference import fresnel_average

mp.mp.dps = 90


@lru_cache(None)
def traversal(k):
    k = mp.mpf(k)
    points = [0, mp.mpf('.001'), mp.mpf('.01'), mp.mpf('.1'), 1]
    # Scale out exp(-k), so opaque cases do not rely on an absolute quadrature
    # tolerance to resolve a transmission hundreds of orders below unity.
    tau = 2*mp.exp(-k)*mp.quad(lambda mu: mu*mp.exp(-k*(1/mu-1)), points) if k else mp.mpf(1)
    loss = 2*mp.quad(lambda mu: mu*(-mp.expm1(-k/mu)), points) if k else mp.mpf(0)
    assert abs(tau+loss-1) < mp.mpf('1e-85')
    depth = -mp.log1p(-loss) if loss <= mp.mpf('.5') else -mp.log(tau)
    return tau, loss, depth


@lru_cache(None)
def surfaces(nr):
    nr = mp.mpf(nr)
    return fresnel_average(mp.mpf(90), nr), fresnel_average(mp.mpf(40), nr)


def source_integral(depth, tau, exit_transmission, entry):
    reflection = 1-exit_transmission
    prefactor = depth*entry/(1-(reflection*tau)**2)
    points = [0, mp.mpf('.001'), mp.mpf('.01'), mp.mpf('.1'), mp.mpf('.5'),
              mp.mpf('.9'), mp.mpf('.99'), mp.mpf('.999'), 1]
    return mp.quad(lambda z: prefactor*(mp.exp(-depth*z)+reflection*tau*mp.exp(-depth*(1-z))), points)


def main():
    path = pathlib.Path(__file__).with_name('plate_preparation_reference.csv')
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['precision','k','nr','traversal_transmission','traversal_absorptance',
                         'effective_depth','interface_entry_transmission','interface_exit_transmission',
                         'source_diffuse90','entry40','source_cone40'])
        for precision, fmt in [('Float32','f'),('Float64','d')]:
            freeze = lambda x: struct.unpack(fmt,struct.pack(fmt,x))[0]
            for k in map(freeze, [0,1e-25,1e-12,1e-6,.01,.1,.5,1,5,20,50,90,100,500]):
                tau, loss, depth = traversal(k)
                for nr in map(freeze,[1.25,1.3,1.4,1.5,1.6]):
                    entry, entry40 = surfaces(nr)
                    exit_transmission = entry/mp.mpf(nr)**2
                    source90 = source_integral(depth,tau,exit_transmission,entry)
                    source40 = source_integral(depth,tau,exit_transmission,entry40)
                    writer.writerow([precision,repr(k),repr(nr),*(mp.nstr(x,85) for x in
                        (tau,loss,depth,entry,exit_transmission,source90,entry40,source40))])
    print(f'Wrote 140 optical cases with 280 source-depth integrals to {path}')


if __name__ == '__main__':
    main()
