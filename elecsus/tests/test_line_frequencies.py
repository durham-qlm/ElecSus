"""
Transition frequencies are no longer truncated to whole MHz.

spectra.FreqStren, as in ES3, truncated every transition frequency (MHz) with
int(). Each line then sat up to 1 MHz too close to the line centre, and the
spectrum was a staircase in Bfield: the line positions moved in 1 MHz jumps.
A fit differentiates by finite differences with a step of ~1e-5 G, so its
Jacobian saw only the change in line strength, and ML fits of Bfield stalled
near their start (Rb D2, true 1000 G: 900 G -> 899.97 G). DE found B, but its
ML round-off reported an error bar ~10x too large.

The ES3-equivalence tests (test_es3_golden.py and the published-figure
anchors) switch the truncation back on with spectra._ES3_WHOLE_MHZ_LINES.
"""
import os
import warnings

import numpy as np
import pytest

from elecsus.elecsus_methods import calculate, fit_data
from elecsus.libs import spectra
from elecsus.libs import AtomConstants as AC
from elecsus.tests.test_es3_golden import ES3_K41_D1_AP

GOLDEN = np.load(os.path.join(os.path.dirname(__file__), 'golden', 'es3_reference.npz'))
CASE_NAMES = [str(n) for n in GOLDEN['case_names']]


def _s0(x, p):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return calculate(x, [1, 0, 0], p, outputs=['S0'])[0].real


def test_es4_differs_from_es3_by_the_truncation_only(monkeypatch):
    """Against ES3's own output, exact frequencies move S0 by at most a few
    1e-3 and the fractional Stokes parameters by a few 1e-2: about what a
    shift of under 1 MHz does to Doppler-broadened lines. (With the
    truncation back on, test_es3_golden.py agrees to ~1e-9.) ES3's K41 P1/2
    constant is put back, as there."""
    monkeypatch.setattr(AC.K41_D1, 'Ap', ES3_K41_D1_AP)
    worst_s0 = worst_stokes = 0.0
    for name in CASE_NAMES:
        X = GOLDEN[f'{name}__X']
        E_in = GOLDEN[f'{name}__E_in']
        p = eval(GOLDEN[f'{name}__p'].tobytes().decode())
        # ES3's geometry is the mirror image of ES4's, and ES3's
        # propagation is no longer the default
        p = dict(p, Btheta=np.pi - float(p.get('Btheta', 0.0)), Bphi=-float(p.get('Bphi', 0.0)),
                 use_non_orthogonal=False)
        g = {k: GOLDEN[f'{name}__{k}'].real for k in ('S0', 'S1', 'S2', 'S3')}
        S0, S1, S2, S3 = (np.asarray(a).squeeze().real
                          for a in calculate(X, E_in, p, outputs=['S0', 'S1', 'S2', 'S3']))
        m = g['S0'] > 1e-2   # where S_i/S0 is not 0/0
        worst_s0 = max(worst_s0, np.max(abs(S0 - g['S0'])))
        for S, key in ((S1, 'S1'), (S2, 'S2'), (S3, 'S3')):
            worst_stokes = max(worst_stokes, np.max(abs(S[m] - g[key][m] / g['S0'][m])))
    assert 1e-4 < worst_s0 < 5e-3
    assert worst_stokes < 3e-2


@pytest.mark.parametrize('elem, dline, bfield', [('Rb', 'D2', 500.), ('K', 'D1', 300.),
                                                 ('Cs', 'D1', 100.)])
def test_spectrum_is_smooth_in_bfield(elem, dline, bfield):
    """dS0/dB by central differences is the same with a step of 1e-4 G (the
    size of a fit's) and of 0.05 G. With whole-MHz lines the small step
    missed ~100 % of it: the lines did not move."""
    x = np.linspace(-8000, 8000, 800)
    p = dict(Elem=elem, Dline=dline, T=60., lcell=10e-3, Bfield=bfield)

    def derivative(h):
        return (_s0(x, dict(p, Bfield=bfield + h)) - _s0(x, dict(p, Bfield=bfield - h))) / (2 * h)

    small, large = derivative(1e-4), derivative(0.05)
    assert np.linalg.norm(small - large) / np.linalg.norm(large) < 1e-2


def _fit(data, p, bools, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return fit_data(data, p, bools, **kwargs)


def test_ml_fit_of_bfield_converges():
    """Rb D2, true 1000 G, ML from 900 G: it stalled at 899.97 G."""
    x = np.linspace(-20000, 20000, 600)
    true = dict(Elem='Rb', Dline='D2', T=50., lcell=75e-3, Bfield=1000.)
    y = _s0(x, true) + np.random.default_rng(0).normal(0, 0.002, len(x))
    best, rms, result = _fit([x, y], dict(true, Bfield=900.), {'Bfield': True})
    stderr = result.params['Bfield'].stderr
    assert abs(best['Bfield'] - 1000.) < 3 * stderr
    assert stderr < 0.5
    assert rms < 0.0022


def test_ml_fit_of_temperature_and_bfield_converges():
    """K D1, true (100 C, 300 G), ML from (90, 250): it stopped at
    (90.4, 210) with an rms 45 times the noise."""
    x = np.linspace(-4000, 4000, 400)
    true = dict(Elem='K', Dline='D1', T=100., lcell=5e-3, Bfield=300.)
    y = _s0(x, true) + np.random.default_rng(0).normal(0, 0.002, len(x))
    best, rms, result = _fit([x, y], dict(true, T=90., Bfield=250.),
                             {'T': True, 'Bfield': True})
    for key, value in (('T', 100.), ('Bfield', 300.)):
        assert abs(best[key] - value) < 3 * result.params[key].stderr
    assert rms < 0.0022


def test_es3_truncation_switch_restores_whole_mhz_lines(monkeypatch):
    """The switch the ES3-equivalence tests use makes the spectrum a
    staircase in B again."""
    x = np.linspace(-8000, 8000, 400)
    p = dict(Elem='Rb', Dline='D2', T=60., lcell=10e-3, Bfield=500.)
    monkeypatch.setattr(spectra, '_ES3_WHOLE_MHZ_LINES', True)
    a, b = _s0(x, p), _s0(x, dict(p, Bfield=500. + 1e-4))
    assert np.max(abs(a - b)) < 1e-6 * np.max(abs(1 - a))
