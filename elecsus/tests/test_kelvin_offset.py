"""One Kelvin offset: the Boltzmann populations used
T + 273.16 where everything else used T + 273.15. ES3 did the same; the ES3
golden tests set 273.16 back."""
import numpy as np

from elecsus.elecsus_methods import calculate
from elecsus.libs import spectra


def test_the_boltzmann_populations_use_the_same_offset_as_everything_else():
    assert spectra._KELVIN_BOLTZMANN == 273.15


def test_the_es3_offset_is_a_live_switch_that_moves_the_spectrum_by_under_1e_6(monkeypatch):
    X = np.linspace(-4000, 4000, 400)
    p = {'Elem': 'Rb', 'Dline': 'D2', 'T': 80., 'lcell': 2e-3}
    [now] = calculate(X, [1, 0, 0], p, outputs=['S0'])
    monkeypatch.setattr(spectra, '_KELVIN_BOLTZMANN', 273.16)
    [es3] = calculate(X, [1, 0, 0], p, outputs=['S0'])
    diff = np.max(np.abs(np.real(now) - np.real(es3)))
    assert 0 < diff < 1e-6
