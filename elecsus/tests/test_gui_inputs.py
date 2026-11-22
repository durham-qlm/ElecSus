"""What the GUI does with what is typed: zero, blank and out-of-range values.

Dash sends None for an empty or non-numeric number field (and, in earlier versions, for
a value outside the field's min/max). The callbacks used to turn None, and also
0, into a default with `float(x or default)`, so a scan from 0 GHz ran from
-5 GHz and Btheta = 200 computed at 0, with nothing on the page to say so
. A required field must now hold a number, and a blank one
is refused with a message naming it.
"""
import warnings

import numpy as np
import pytest

from elecsus import elecsus_dash_gui as gui
from elecsus.elecsus_methods import calculate

# compute_single's arguments, in the order of its State list
ARGS = dict(
    elem='Rb', dline='D2', T=60., lcell=0.075, B=200., gbuf=0., shift=0.,
    btheta=0., bphi=0., theta0=0., pol='Linear (x)',
    bufgas='None', bufP=0., bufTf=293.15,
    rb85frac=72.17, K40frac=0.01, K41frac=6.73,
    det_start=-3., det_stop=3., det_N=40,
)


def _compute(outputs=('S0',), nonortho=gui.NONORTHO_DEFAULT, windows=None, n_window=None, doppT=None,
             **over):
    a = dict(ARGS, **over)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return gui.compute_single(1, *[a[k] for k in ARGS], list(outputs),
                                  nonortho, windows, n_window, doppT)[:3]


def _core(store, E_in=(1., 0., 0.), **p):
    base = dict(Elem='Rb', Dline='D2', T=60., lcell=0.075, Bfield=200.,
                GammaBuf=0., shift=0., Btheta=0., Bphi=0., rb85frac=72.17)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        [S0] = calculate(np.array(store['X']) * 1e3, list(E_in), dict(base, **p),
                         outputs=['S0'])
    return np.real(S0)


# --- zero is a value --------------------------------------------------

def test_the_detuning_axis_starts_and_stops_where_it_is_told():
    np.testing.assert_allclose(gui._X(0., 10., 11) / 1e3, range(11))
    np.testing.assert_allclose(gui._X(-10., 0., 11) / 1e3, range(-10, 1))


def test_a_scan_from_zero_computes_from_zero():
    _, store, err = _compute(det_start=0., det_stop=6.)
    assert err == '' and store['X'][0] == 0. and store['X'][-1] == 6.
    np.testing.assert_allclose(store['results']['S0'], _core(store), rtol=1e-12)


def test_zero_celsius_is_not_twenty():
    assert gui._p('Rb', 'D2', 0., 0.075, 0, 0, 0, 0, 0, 0, 'None', 0, 293.15)['T'] == 0.
    fig, store, err = _compute(T=0.)
    assert err == ''
    np.testing.assert_allclose(store['results']['S0'], _core(store, T=0.), rtol=1e-12)
    assert 'T=0°C' in fig.layout.title.text


def test_a_zero_cell_length_reaches_the_core_which_refuses_it():
    _, _, err = _compute(lcell=0.)
    assert "'lcell'" in err


# --- a blank required field is refused, not defaulted ------------

REQUIRED = [('T', 'T [°C]'), ('lcell', 'l [m]'), ('B', 'B [G]'),
            ('gbuf', 'GammaBuf [MHz]'), ('shift', 'Shift [MHz]'),
            ('btheta', 'Bθ [deg]'), ('bphi', 'Bφ [deg]'), ('theta0', 'θ₀ [deg]'),
            ('rb85frac', '⁸⁵Rb (%)'),
            ('det_start', 'Start [GHz]'), ('det_stop', 'Stop [GHz]'),
            ('det_N', 'N pts')]


@pytest.mark.parametrize('key, label', REQUIRED, ids=[k for k, _ in REQUIRED])
def test_a_blank_required_field_is_refused_by_name(key, label):
    fig, store, err = _compute(**{key: None})
    assert label in err and 'empty' in err
    assert store is None and fig is not None


def test_blank_buffer_gas_fields_are_refused_only_with_a_gas():
    _, _, err = _compute(bufP=None, bufTf=None)
    assert err == ''                                   # no gas: not used
    _, _, err = _compute(bufgas='N2', bufP=None)
    assert 'P [Torr]' in err
    _, _, err = _compute(bufgas='N2', bufP=5., bufTf=None)
    assert 'Fill T [K]' in err


def test_blank_gammabuf_is_fine_while_a_buffer_gas_supplies_it():
    _, _, err = _compute(bufgas='N2', bufP=5., gbuf=None)
    assert err == ''


def test_blank_fields_that_are_not_used_are_ignored():
    _, _, err = _compute(pol='Left CP', theta0=None)               # CP has no angle
    assert err == ''
    _, _, err = _compute(elem='Cs', rb85frac=None, K40frac=None)   # Rb and K only
    assert err == ''
    _, _, err = _compute(elem='K', K40frac=None)
    assert '⁴⁰K' in err


def test_fields_where_blank_means_something_stay_optional():
    _, _, err = _compute(nonortho=['on'], windows=['on'], n_window=None, doppT=None)
    assert err == ''


# --- values outside the old min/max reach the code that can judge them -

@pytest.mark.parametrize('btheta, bphi', [(200., 0.), (-10., 0.), (60., -90.), (60., 400.)])
def test_field_angles_are_used_as_typed(btheta, bphi):
    _, store, err = _compute(B=300., btheta=btheta, bphi=bphi)
    assert err == ''
    np.testing.assert_allclose(
        store['results']['S0'],
        _core(store, Bfield=300., Btheta=np.radians(btheta), Bphi=np.radians(bphi)),
        rtol=1e-12)


@pytest.mark.parametrize('key, value, named', [
    ('rb85frac', 105., 'rb85frac'),
    ('bufP', -1., 'BufferPressure'),
    ('bufTf', 0., 'BufferTemp_fill'),
])
def test_out_of_range_values_are_refused_by_the_core_not_replaced(key, value, named):
    over = {key: value}
    if key != 'rb85frac':
        over['bufgas'] = 'N2'
        over.setdefault('bufP', 5.)
    _, _, err = _compute(**over)
    assert named in err


def test_a_window_index_below_one_is_refused_not_replaced():
    _, _, err = _compute(nonortho=['on'], windows=['on'], n_window=0.5)
    assert 'n_window' in err


@pytest.mark.parametrize('n', [5, 6000, 100.5])
def test_the_point_count_is_checked_and_the_error_says_so(n):
    _, store, err = _compute(det_N=n)
    assert 'N pts' in err and store is None


@pytest.mark.parametrize('n', [10, 5000])
def test_the_point_count_limits_are_inclusive(n):
    assert len(gui._X(-1., 1., n)) == n


# --- the other tabs ----------------------------------------------------------

CASCADE = dict(
    e1='Rb', d1='D2', T1=60., l1=0.075, B1=200., g1=0., sh1=0., bt1=0., bp1=0.,
    th1=0., pol1='Linear (x)', bg1='None', bP1=0., bTf1=293.15,
    rb85=72.17, k40=0.01, k41=6.73,
    T2=80., l2=0.075, B2=50., g2=0., sh2=0., bt2=0., bp2=0.,
    bg2='None', bP2=0., bTf2=293.15,
    use_pol=['on'], pol_angle=90., reset_pol='none',
    det_start=-3., det_stop=3., det_N=40, outputs=['S0'])


def _cascade(**over):
    a = dict(CASCADE, **over)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return gui.compute_cascade(1, *a.values())[:3]


def test_the_cascade_names_the_cell_with_the_blank_field():
    _, _, err = _cascade(T2=None)
    assert err.startswith('Cell 2: T [°C]')
    _, _, err = _cascade(B1=None)
    assert err.startswith('Cell 1: B [G]')


def test_the_cascade_polariser_angle_zero_is_an_angle():
    fig, store, err = _cascade(pol_angle=0.)
    assert err == '' and 'polariser at 0°' in fig.layout.title.text
    _, _, err = _cascade(pol_angle=None)
    assert 'Polariser angle' in err
    _, _, err = _cascade(use_pol=[], pol_angle=None)       # off: the angle is not used
    assert err == ''


def test_a_blank_field_stops_a_fit_with_a_message():
    store = gui.load_fit_sample(1)[0]
    out = gui.run_fit_single(1, store, 'S0', 'ML', ['T'], None, None, 'Rb', 'D2',
                             None, 0.075, 0., 0., 0., 0., 0., 0., 'Linear (x)',
                             'None', 0., 293.15, 72.17, 0.01, 6.73)
    assert 'T [°C]' in out[2]
