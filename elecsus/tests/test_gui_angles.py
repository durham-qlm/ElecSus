"""
GUI angle-unit convention tests.

The core p_dict convention for Btheta/Bphi is radians (see solve_dielectric.py);
the Dash GUI presents these angles in degrees (matching the ES3 wxPython GUI and
the other angle controls in the panel) and converts at the _p() boundary.
These tests pin that conversion so neither side drifts.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate
from elecsus.libs.spectra import get_spectra
from elecsus import elecsus_dash_gui as gui


def test_p_converts_degrees_to_radians():
    p = gui._p('Rb', 'D2', 20., 0.075, 100., 0., 0., 90., 180., 0.,
               'None', 0., 293.15)
    assert p['Btheta'] == pytest.approx(np.pi / 2)
    assert p['Bphi'] == pytest.approx(np.pi)


def test_gui_voigt_matches_core_radians():
    """A 90-degree GUI Btheta must reproduce the core Voigt (pi/2) spectrum."""
    X = np.linspace(-6000, 6000, 200)
    E_in = np.array([1., 0., 0.])

    p_gui = gui._p('Rb', 'D2', 60., 0.075, 1000., 0., 0., 90., 0., 0.,
                   'None', 0., 293.15)
    p_core = {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 0.075,
              'Bfield': 1000., 'Btheta': np.pi / 2, 'Bphi': 0.}

    S0_gui = calculate(X, E_in, p_gui, outputs=['S0'])[0]
    S0_core = calculate(X, E_in, p_core, outputs=['S0'])[0]
    np.testing.assert_allclose(S0_gui, S0_core, atol=1e-12)

    # and 90 deg is genuinely different physics from 90 rad ≡ its 2π-wrap
    p_bad = dict(p_core, Btheta=90.0)
    S0_bad = calculate(X, E_in, p_bad, outputs=['S0'])[0]
    assert not np.allclose(S0_gui, S0_bad, atol=1e-3)


def test_fit_writeback_returns_degrees():
    """Varying Btheta through the GUI fit callback must write back degrees."""
    np.random.seed(3)
    btheta_true_deg = 45.0
    p_true = {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 0.075,
              'Bfield': 1000., 'Btheta': np.radians(btheta_true_deg), 'Bphi': 0.,
              'GammaBuf': 0., 'shift': 0.}
    X_mhz = np.linspace(-6000, 6000, 150)
    y = get_spectra(X_mhz, np.array([1., 0., 0.]), p_true, outputs=['S0'])[0].real
    y += np.random.randn(len(X_mhz)) * 0.003

    store = {'X': (X_mhz / 1e3).tolist(), 'Y': y.tolist()}
    out = gui.run_fit_single(
        1, store, 'S0', 'ML', ['Btheta'], None, None,
        'Rb', 'D2', 60., 0.075, 1000., 0., 0., 40.0, 0., 0., 'Linear (x)',
        'None', 0., 293.15, 72.17, 0.01, 6.73)
    err, btheta_out = out[2], out[8]

    assert err == ''
    # written-back value is in degrees, close to the 45-degree truth
    assert btheta_out == pytest.approx(btheta_true_deg, abs=1.0)
