"""Archive completed source-component gates and the exact reviewed implementation."""
import datetime
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(log_directory):
    logs = {}
    for name in ('full-cpu.log', 'batches-gpu.log', 'docs.log'):
        path = log_directory / name
        text = path.read_text()
        if re.search(r'ERROR:|Test Failed|Error During Test|┌ Error:', text):
            raise ValueError('Failed gate: ' + name)
        count = None
        if name == 'docs.log':
            if 'HTMLWriter: rendering HTML pages' not in text:
                raise ValueError('Missing completed docs gate')
        else:
            summaries = re.findall(r'^([^\n|]+)\|\s+(\d+)\s+(\d+)\s+[^\n]+$', text, re.M)
            if not summaries or any(p != n for _, p, n in summaries):
                raise ValueError('Incomplete tests: ' + name)
            count = sum(int(p) for _, p, _ in summaries)
        logs[name] = dict(sha256=digest(path), passed_checks=count, text=text)
    files = ('src/CanopyOptics.jl', 'src/leaf_optics/emission_spectrum.jl',
             'src/leaf_optics/plate_emission.jl', 'test/runtests.jl',
             'test/test_emission_spectrum.jl', 'test/emission_spectrum_batches.jl',
             'test/emission_spectrum_gpu.jl', 'docs/make.jl',
             'docs/src/pages/emission_spectrum.md',
             'test/reference/record_emission_spectrum_validation.py')
    report = dict(generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scope='Host conditional photon-spectrum preparation and shared matrix CPU/CUDA application, Float32 and Float64',
                  source_sha256={f: digest(ROOT / f) for f in files}, logs=logs,
                  references=['BigFloat normalization of retained raw masses',
                              'BigFloat accumulation of prepared operators, independent of CPU/CUDA BLAS'],
                  limitations=['No empirical source shape or calibration is installed.',
                               'No GPU or AD claim for host preparation.',
                               'Existing regular-plate physics is independently validated in its separate certificate.',
                               'Photon and Stokes energy bounds are not complete leaf thermal closure.'])
    output = ROOT / 'docs/src/data/emission_spectrum_validation.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(output)


if __name__ == '__main__':
    main(Path(sys.argv[1]))
