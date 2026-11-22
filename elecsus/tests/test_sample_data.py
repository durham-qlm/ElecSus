"""The bundled scope2elecsus Rb D2 sample spectrum
(elecsus/sample_data/rb_d2_doppler_40C.csv) must be present, well-formed, and
reproducible by ElecSus, so Notebooks/06_fitting fits real published data.
"""
from importlib.resources import files

import numpy as np

from elecsus.elecsus_methods import calculate, fit_data

_P = {'Elem': 'Rb', 'Dline': 'D2', 'lcell': 75e-3, 'Bfield': 0.}


def _load():
    with (files('elecsus') / 'sample_data' / 'rb_d2_doppler_40C.csv').open() as f:
        arr = np.genfromtxt(f, delimiter=',', skip_header=1)
    return arr[:, 0], arr[:, 1]


def test_sample_file_present_and_shaped():
    det, tr = _load()
    assert det.size > 1000 and det.size == tr.size
    assert det.min() < -3000 and det.max() > 3000          # several GHz, in MHz
    assert -0.1 < tr.min() < 0.2 and 0.9 < tr.max() < 1.1   # normalised transmission
    assert det[-1] > det[0]                                 # increasing detuning axis


def test_elecsus_reproduces_sample():
    det, tr = _load()
    [s0] = calculate(det, [1, 0, 0], {**_P, 'T': 40.}, outputs=['S0'])
    rms = float(np.sqrt(np.mean((s0.real - tr) ** 2)))
    assert rms < 0.015, f"forward RMS {rms:.4f} too large for the shipped sample"


def test_fit_recovers_temperature():
    det, tr = _load()
    best, rms, _ = fit_data([det, tr], {**_P, 'T': 60.}, {'T': True},
                            data_type='S0', fit_algorithm='ML')
    assert 37.0 < best['T'] < 43.0, f"recovered T={best['T']:.2f} C (expected ~40)"
    assert rms < 0.015
