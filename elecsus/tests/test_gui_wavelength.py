"""The GUI's wavelength axes: the "Plot against" selector under each Detuning
Axis panel, and the conversion line beneath it.

Astronomers give filter profiles in wavelength, as an offset from line centre
in mÅ or as the air wavelength in Å. The detuning boxes stay in GHz; the plots
can be drawn against either wavelength instead, and redrawn without computing
again.
"""
import contextlib
import io
import warnings

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui

# NIST air wavelengths (Å) of the D lines
NIST_AIR = {('Na', 'D1'): 5895.924, ('Na', 'D2'): 5889.951,
            ('K', 'D1'): 7698.964, ('K', 'D2'): 7664.899,
            ('Rb', 'D1'): 7947.603, ('Rb', 'D2'): 7800.268,
            ('Cs', 'D1'): 8943.474, ('Cs', 'D2'): 8521.132}


@pytest.mark.parametrize('line', sorted(NIST_AIR))
def test_line_centre_is_the_nist_air_wavelength(line):
    """Zero detuning in air, from the core's line centre and the IAU air
    index: within 3 mÅ of NIST for every element and line the GUI offers."""
    _, lam0 = gui._line_centre(*line)
    assert lam0 == pytest.approx(NIST_AIR[line], abs=3e-3)


def test_positive_detuning_is_a_shorter_wavelength():
    """K D1: 1 GHz is -19.78 mÅ, so the Oliviero filter's peaks at +-100 mÅ
    are near -+5.06 GHz; the air wavelength is the line centre plus the offset."""
    X = np.array([-5.056, 0., 1., 5.056])
    dlam = gui._x_on(X, 'mA', 'K', 'D1')
    assert dlam[1] == 0.
    assert dlam[2] == pytest.approx(-19.78, abs=0.01)
    np.testing.assert_allclose(dlam[[0, 3]], [100., -100.], atol=0.05)
    _, lam0 = gui._line_centre('K', 'D1')
    np.testing.assert_allclose(gui._x_on(X, 'A', 'K', 'D1'), lam0 + dlam / 1e3, rtol=1e-15)
    np.testing.assert_array_equal(gui._x_on(X, 'GHz', 'K', 'D1'), X)


def test_the_wavelength_inverts_to_the_detuning():
    """lambda = c / (n (nu0 + detuning)), exactly, at the line-centre index."""
    X = np.linspace(-30., 30., 7)
    nu0, lam0 = gui._line_centre('Rb', 'D2')
    n = gui._air_index(gui._C / nu0)
    lam = gui._x_on(X, 'A', 'Rb', 'D2') * 1e-10
    np.testing.assert_allclose((gui._C / (n * lam) - nu0) / 1e9, X, atol=1e-6)


def test_conversion_line_gives_the_factor_and_the_range():
    text = gui._conversion_text('K', 'D1', -20., 20.)
    assert '+1 GHz ≙ -19.78 mÅ' in text
    assert '+100 mÅ ≙ -5.056 GHz' in text and '+395.6 to -395.5 mÅ' in text
    # a blank or unreadable range leaves the factor
    assert 'range' not in gui._conversion_text('K', 'D1', None, 20.)
    assert gui._conversion_text(None, 'D1', -1., 1.) == ''


def test_the_air_wavelength_axis_is_hidden_until_later():
    """The selector offers detuning and the offset in mÅ
    only, and nothing on the page says "air". The code behind the air axis
    stays (the tests above still check it), for later."""
    assert [o['value'] for o in gui.X_UNIT_OPTIONS] == ['GHz', 'mA']
    assert 'air' not in gui.X_UNITS['mA']
    assert 'air' not in gui._conversion_text('K', 'D1', -20., 20.)
    assert 'λ₀' not in gui._conversion_text('Rb', 'D2', -5., 5.)


SINGLE = dict(elem='K', dline='D1', T=120., lcell=0.025, B=1400., gbuf=0., shift=0.,
              btheta=0., bphi=0., theta0=0., pol='Linear (x)', bufgas='None', bufP=0.,
              bufTf=293.15, rb85frac=72.17, K40frac=0.01, K41frac=6.73,
              det_start=-10., det_stop=10., det_N=201, outputs=['S0', 'Iy'])


def _single(xunit):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return gui.compute_single(1, *SINGLE.values(), xunit=xunit)


def test_compute_draws_against_the_chosen_axis():
    fig, store, err, _ = _single('mA')
    assert err == ''
    want = gui._x_on(store['X'], 'mA', 'K', 'D1')
    for trace in fig.data:
        np.testing.assert_allclose(trace.x, want)
    titles = [fig.layout[k].title.text for k in fig.layout if k.startswith('xaxis')]
    assert 'Δλ from line centre (mÅ)' in titles
    # the store keeps GHz, so the export and the fit are unchanged
    assert store['X'][0] == -10. and store['elem'] == 'K' and store['dline'] == 'D1'


def test_the_selector_redraws_without_computing():
    fig, store, _, _ = _single('GHz')
    redrawn = gui._redraw(fig.to_dict(), store, 'mA')
    np.testing.assert_allclose(redrawn.data[0].x, gui._x_on(store['X'], 'mA', 'K', 'D1'))
    assert redrawn.to_dict()['data'][1]['y'] == fig.to_dict()['data'][1]['y']   # same result
    back = gui._redraw(redrawn.to_dict(), store, 'GHz')
    np.testing.assert_allclose(back.data[0].x, store['X'])
    titles = [back.layout[k].title.text for k in back.layout if k.startswith('xaxis')]
    assert 'Detuning (GHz)' in titles
    # nothing computed yet, or the last compute failed: the figure stays
    assert gui._redraw(fig.to_dict(), None, 'mA') is gui.no_update
    assert gui.redraw_single('mA', None, None, None, None) == (gui.no_update,
                                                                gui.no_update)


def test_the_fit_plot_follows_the_selector():
    sample = gui.load_fit_sample(1)[0]
    a = dict(SINGLE, elem='Rb', dline='D2', T=40., lcell=0.075, B=0.)
    for key in ('det_start', 'det_stop', 'det_N', 'outputs'):
        a.pop(key)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        out = gui.run_fit_single(1, sample, 'S0', 'ML', ['T'], None, None,
                                 *a.values(), xunit='mA')
    fig, fit_store = out[0], out[10]
    assert out[2] == ''
    want = gui._x_on(fit_store['X'], 'mA', 'Rb', 'D2')
    for trace in fig.data:
        np.testing.assert_allclose(trace.x, want)
    back = gui._redraw(fig.to_dict(), fit_store, 'GHz')
    np.testing.assert_allclose(back.data[0].x, fit_store['X'])


def test_the_cascade_draws_both_panels_against_the_chosen_axis():
    fig, store, err, _ = gui.compute_cascade(
        1, 'K', 'D1', 120., 0.025, 1400., 400., 0., 0., 0., 0., 'Linear (x)',
        'None', 0., 293.15, 72.17, 0.01, 6.73,
        160., 0.03, 2700., 400., 0., 0., 0., 'None', 0., 293.15,
        ['on'], 90., 'right_cp', -20., 20., 201, ['S0'], xunit='mA')
    assert err == ''
    want = gui._x_on(store['X'], 'mA', 'K', 'D1')
    for trace in fig.data:
        np.testing.assert_allclose(trace.x, want)
    redrawn = gui.redraw_cascade('GHz', fig.to_dict(), store)
    np.testing.assert_allclose(redrawn.data[-1].x, store['X'])
