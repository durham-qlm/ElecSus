"""
Kramers-Kronig between Re chi and Im chi (a convention-free physics test).

A causal response is analytic in the upper half of the complex frequency
plane. Near a resonance (the rotating-wave form ElecSus computes, chi as a
function of the detuning Delta), and with fields varying as exp(-i omega t),
that gives the Hilbert-transform pair

    Re chi(Delta) = (1/pi) P int Im chi(Delta') / (Delta' - Delta) dDelta',

for every susceptibility component and whatever the lineshape, the field, the
broadening or the shift. The sign is physics, not convention: it makes the
vapour's index exceed 1 just below an absorption line (normal dispersion in
the wings) and is what a causal Lorentz oscillator gives.

The transform is computed here from Im chi alone, by Maclaurin's rule on a
uniform grid (the principal value from the points an odd number of steps
away; Ohta and Ishida, Appl. Spectrosc. 42, 952 (1988)), and compared with
the Re chi the code returns.
"""
import warnings

import numpy as np
import pytest
from scipy.signal import fftconvolve

from elecsus.libs.spectra import calc_chi

STEP = 1.0                                   # MHz
X = np.arange(-100e3, 100e3 + STEP, STEP)    # wide, so the truncated wings
CORE = np.abs(X) < 20e3                      # matter little in the core


def _hilbert_from_absorption(im):
    """(1/pi) P int im(x') / (x' - x) dx' at every grid point, by Maclaurin's
    rule: (2/pi) sum over j - i odd of im_j / (j - i)."""
    n = len(im)
    m = np.arange(-(n - 1), n)
    kernel = np.zeros(len(m))
    odd = m % 2 != 0
    kernel[odd] = 1.0 / m[odd]
    # full[i + n - 1] = sum_j im_j kernel[j - i]
    full = fftconvolve(im, kernel[::-1], mode='full')
    return (2 / np.pi) * full[n - 1:2 * n - 1]


CASES = {
    'Rb D2, zero field': dict(Elem='Rb', Dline='D2', T=50.),
    'Rb D2, 2 kG': dict(Elem='Rb', Dline='D2', T=50., Bfield=2000.),
    'Rb D2, 1 T': dict(Elem='Rb', Dline='D2', T=80., Bfield=10000.),
    'Cs D1, Doppler-free': dict(Elem='Cs', Dline='D1', T=50., Constrain=False,
                                DoppTemp=-272.),
    'K D2, buffer gas (100 Torr Ar)': dict(Elem='K', Dline='D2', T=50., BufferGas='Ar',
                                           BufferPressure=100.),
    'Na D1, GammaBuf and shift': dict(Elem='Na', Dline='D1', T=150., GammaBuf=300.,
                                      shift=-2000.),
}


@pytest.mark.parametrize('case', sorted(CASES))
def test_real_part_is_the_hilbert_transform_of_the_imaginary_part(case):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chis = calc_chi(X, CASES[case])
    for name, chi in zip(('sigma+', 'sigma-', 'pi'), chis):
        scale = np.max(np.abs(chi[CORE]))
        kk = _hilbert_from_absorption(chi.imag)
        error = np.max(np.abs(kk[CORE] - chi.real[CORE])) / scale
        assert error < 1e-4, f'{name}: KK residual {error:.2e} of max |chi|'
        # and the opposite sign (an anti-causal response) is far off
        assert np.max(np.abs(-kk[CORE] - chi.real[CORE])) / scale > 0.5


def test_dispersion_is_normal_in_the_wings():
    """The causal sign, directly: below every line (red detuning) the index
    is above 1, and above every line it is below 1."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chis = calc_chi(np.array([-30e3, 30e3]), dict(Elem='Rb', Dline='D2', T=50.))
    for chi in chis:
        assert chi.real[0] > 0 > chi.real[1]
        assert np.all(chi.imag > 0)                       # and it absorbs
