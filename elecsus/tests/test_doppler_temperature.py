"""
The Doppler temperature.

- DoppTemp is the temperature of the vapour in the beam. It sets the Doppler
  width and the buffer-gas model (test_buffer_gas_fits.py). With
  Constrain=True, the default, it is T. Constrain=False gives a cell body
  hotter than its stem its own temperature, T then setting the number density.
- A DoppTemp not given, or None (its default), is T. It was 20 C, so
  Constrain=False without DoppTemp put the Doppler width at 20 C whatever T
  was: in calculate(), and in fits, whose model starts DoppTemp at 20 C. That
  setting now warns that it has no effect, and fit results no longer report
  the model's unused DoppTemp.
"""
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import MLFittingRoutine as ML
from elecsus.libs import spectra as SP

X = np.linspace(-4000, 4000, 161)
CELL = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3)


def _S0(p_dict):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.real(calculate(X, [1, 0, 0], p_dict, outputs=['S0'])[0])


class _Stop(Exception):
    """Raised by a stand-in fit routine once it has seen its arguments."""


# --- calculate --------------------------------------------------------------------

def test_dopptemp_sets_the_doppler_width_only_without_constrain():
    assert np.array_equal(_S0(dict(CELL, DoppTemp=150.)), _S0(CELL))
    assert np.array_equal(_S0(dict(CELL, Constrain=False, DoppTemp=80.)), _S0(CELL))
    assert np.max(np.abs(_S0(dict(CELL, Constrain=False, DoppTemp=150.)) - _S0(CELL))) > 1e-3


def test_constrain_false_without_dopptemp_follows_T():
    """The Doppler width was taken at 20 C, whatever T was."""
    unconstrained = _S0(dict(CELL, Constrain=False))
    assert np.array_equal(unconstrained, _S0(CELL))
    assert np.max(np.abs(unconstrained - _S0(dict(CELL, Constrain=False, DoppTemp=20.)))) > 1e-3


def test_constrain_false_without_dopptemp_warns():
    with pytest.warns(UserWarning, match='Constrain=False has no effect without DoppTemp'):
        calculate(X, [1, 0, 0], dict(CELL, Constrain=False), outputs=['S0'])


def test_a_copy_of_the_defaults_follows_T():
    """p_dict_defaults holds DoppTemp=None, which means T: a script that
    starts from a copy of the defaults gets T's Doppler width, not 20 C's."""
    defaults = dict(SP.p_dict_defaults, T=CELL['T'], lcell=CELL['lcell'])
    assert defaults['DoppTemp'] is None
    assert np.array_equal(_S0(defaults), _S0(CELL))
    assert np.array_equal(_S0(dict(defaults, Constrain=False)), _S0(CELL))


# --- fit_data ---------------------------------------------------------------------

def test_a_fit_without_dopptemp_keeps_the_temperatures_together():
    """Constrain=False without DoppTemp: the fit model started DoppTemp at
    20 C and kept the Doppler width there while T varied, so T came out wrong.
    The fit now runs, and reports, the model calculate() would use."""
    data = [X, _S0(CELL)]
    with pytest.warns(UserWarning, match='Constrain=False has no effect without DoppTemp'):
        best, rms, result = fit_data(data, dict(CELL, T=70., Constrain=False), {'T': True})
    assert best['T'] == pytest.approx(80., abs=1e-3)
    assert best['Constrain'] is True
    assert 'DoppTemp' not in best        # not the fit model's unused 20 C
    assert np.allclose(_S0(best), result.best_fit, atol=1e-12)


def test_a_fit_of_dopptemp_without_a_start_value_starts_at_T(monkeypatch):
    """It started at the fit model's default, 20 C."""
    seen = {}

    def stand_in(data, E_in, p_dict, p_dict_bools, **kw):
        seen.update(p_dict)
        raise _Stop

    monkeypatch.setattr(ML, 'ML_fit', stand_in)
    with pytest.raises(_Stop):
        fit_data([X, _S0(CELL)], dict(CELL, Constrain=False), {'DoppTemp': True})
    assert seen.get('DoppTemp') == CELL['T']
    assert seen['Constrain'] is False
