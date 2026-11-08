"""Independent 90-digit depth quadrature and explicit reflection-path ledger.

Requires mpmath 1.3.0; no Julia, production formulas or optical database. Inputs
are frozen to the requested IEEE precision before high-precision evaluation.
The unprocessed path remainder is bounded, never rescaled into a photon fate.
"""
import csv
import pathlib
import struct
import mpmath as mp

mp.mp.dps = 90
CASES = [
    (0, 0, .6, 1, 0), (0, .01, 1, 1, 0), (0, 2, .6, 1, 0),
    (1, 1, .6, 1, .2), (1, 1, .6, .2, 1), (8, 1e-20, .05, 1, 0),
    (9, 1e-7, .6, 1, 1), (30, .01, .05, 1, .3), (30, 30, .6, 1, 0),
    (40, 2, .6, 1, 0), (2, 40, .6, 1, 0), (1e-7, 20, 1, 1, 0),
    (.1, .1, .6, 1, 1e-20), (.1, .1, .6, 1e-20, 1),
]


def reference(a, b, t, f, g):
    def raw(z):
        return f*mp.exp(-a*z)+g*mp.exp(-a*(1-z))
    points = [0, mp.mpf('.01'), mp.mpf('.1'), mp.mpf('.5'), mp.mpf('.9'), mp.mpf('.99'), 1]
    normalization = mp.quad(raw, points)
    front_arrival = mp.quad(lambda z: raw(z)*mp.exp(-b*z)/2, points)/normalization
    back_arrival = mp.quad(lambda z: raw(z)*mp.exp(-b*(1-z))/2, points)/normalization
    absorbed = mp.quad(lambda z: raw(z)*(-mp.expm1(-b*z)-mp.expm1(-b*(1-z)))/2, points)/normalization
    front, back = mp.mpf(0), mp.mpf(0)
    reflected = 1-t
    traversal_survival = mp.exp(-b)
    traversal_loss = -mp.expm1(-b)
    for paths in range(10000):
        front += t*front_arrival
        back += t*back_arrival
        absorbed += reflected*traversal_loss*(front_arrival+back_arrival)
        front_arrival, back_arrival = (reflected*traversal_survival*back_arrival,
                                       reflected*traversal_survival*front_arrival)
        remainder = front_arrival+back_arrival
        if remainder < mp.mpf('1e-85'):
            assert abs(front+back+absorbed+remainder-1) < mp.mpf('1e-87')
            return front, back, absorbed, remainder, paths+1
    raise RuntimeError('reflection path sum did not converge')


def main():
    path = pathlib.Path(__file__).with_name('plate_emission_reference.csv')
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['precision','excitation_depth','emission_depth','interface_transmission',
                         'forward_amplitude','backward_amplitude','front','back','reabsorbed',
                         'unprocessed_path_bound','path_count'])
        for precision, fmt in [('Float32','f'),('Float64','d')]:
            for values in CASES:
                frozen = [struct.unpack(fmt,struct.pack(fmt,x))[0] for x in values]
                result = reference(*map(mp.mpf,frozen))
                writer.writerow([precision,*map(repr,frozen),
                                 *(mp.nstr(x,85) for x in result[:-1]),result[-1]])
    print(f'Wrote {2*len(CASES)} independent quadrature/path references to {path}')


if __name__ == '__main__':
    main()
