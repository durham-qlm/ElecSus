"""
The integrated-absorption sum rule (a convention-free physics test).

For a vapour of N atoms per unit volume, the absorption coefficient
alpha = k Im(chi) of an unpolarised ensemble integrates over frequency to

    int alpha dnu = N pi r_e c f,

with r_e the classical electron radius and f the line's absorption oscillator
strength (e.g. Steck, Sec. 2.3; Corney, "Atomic and Laser Spectroscopy",
ch. 4). It holds for any lineshape, so it is independent of the Doppler
temperature and of collisional (buffer-gas) broadening and shifts. A magnetic
field moves strength between the sigma+, sigma- and pi manifolds, and
between lines, but not out of the isotropic part of the susceptibility,
Tr chi / 3 = (chi+ + chi- + chi_pi) / 3: that is the trace of a tensor, which
no rotation and no redistribution among the Zeeman components can change.

The expected value is computed from the literature lifetime and wavelength
(literature.py), never from ElecSus's constants. The number density is the
one input taken from the code, and it is checked against tabulated vapour
pressures separately.

In the optically thin limit the same integral fixes the absorbed fraction,
int (1 - S0) dnu = N L pi r_e c f, for any input polarisation and any field
direction: the integrated susceptibility tensor of the unpolarised ensemble is
isotropic, so every polarisation sees the same integrated absorption.
"""
import warnings

import numpy as np
import pytest
from scipy.constants import c, k as k_B, pi, physical_constants, torr
from scipy.integrate import trapezoid

from elecsus.elecsus_methods import calculate
from elecsus.libs.spectra import calc_chi
from elecsus.libs import numberDensityEqs
from elecsus.tests import literature as lit

R_E = physical_constants['classical electron radius'][0]

# Detuning axis (MHz): 1 MHz steps over +-40 GHz resolve a Doppler-free line
# (Gamma ~ 5 MHz) wherever a field up to 1 T puts it; 20 MHz steps out to
# +-200 GHz carry the Lorentzian wings, of which ~1e-5 is left beyond the axis
# for the far-wing correction in _integral.
X = np.concatenate([np.arange(-200e3, -40e3, 20.), np.arange(-40e3, 40e3, 1.),
                    np.arange(40e3, 200e3 + 1, 20.)])

# FreqStren drops every transition with |<e|q,g>|^2 below 0.0005. The strength
# so lost depends on the field (and, near B = 0, on the arbitrary choice of
# eigenvectors in each degenerate level): up to 3e-4 of the total, measured
# for Rb, Cs, K and Na from 0 to 3 T. With the cut removed, the invariance
# below holds to < 1e-5 up to 1 T. (At 3 T, Na D2 gives 1.5e-4 of its strength
# to D1: fine-structure mixing, of order (muB B / fine-structure splitting)^2.)
CUT = 5e-4

DENSITY = {'Rb': numberDensityEqs.numDenRb, 'Cs': numberDensityEqs.numDenCs,
           'K': numberDensityEqs.numDenK, 'Na': numberDensityEqs.numDenNa}

LINES = sorted(lit.LINES)


def _integral(y, X=X):
    """int y dnu over all frequencies, in Hz times y's unit. Beyond the axis
    the wings are Lorentzian, y ~ C / Delta^2, whose tail from the edge
    Delta_e to infinity is C / |Delta_e| = y(Delta_e) |Delta_e|."""
    return 1e6 * (trapezoid(y, X) + y[0] * abs(X[0]) + y[-1] * abs(X[-1]))


def _chi(p):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return calc_chi(X, p)


def _expected_per_atom(elem, dline):
    """int alpha dnu / N = pi r_e c f, in m^2 Hz."""
    return pi * R_E * c * lit.oscillator_strength(elem, dline)


def _k(elem, dline):
    return 2 * pi / lit.LINES[(elem, dline)]['wavelength']


def _isotropic_integral_per_atom(p):
    """(k / N) int Im(Tr chi / 3) dnu."""
    chi = _chi(p)
    trace = sum(component.imag for component in chi) / 3
    N = DENSITY[p['Elem']](p['T'] + 273.15)
    return _k(p['Elem'], p['Dline']) * _integral(trace) / N


# --- the number density that the sum rule is normalised by --------------------

@pytest.mark.parametrize('elem', sorted(lit.VAPOUR_PRESSURE_25C_TORR))
def test_number_density_matches_the_tabulated_vapour_pressure(elem):
    """N k_B T at 25 C against Steck's tabulated vapour pressure, to its ~5 %."""
    T = 298.15
    pressure = DENSITY[elem](T) * k_B * T / torr
    np.testing.assert_allclose(pressure, lit.VAPOUR_PRESSURE_25C_TORR[elem], rtol=0.05)


# --- the absolute sum rule -------------------------------------------------------

@pytest.mark.parametrize('elem, dline', LINES)
def test_integrated_absorption_matches_the_oscillator_strength(elem, dline):
    """The 1 % tolerance covers the spread of measured lifetimes; an error in
    the dipole moment, the 1/3 of the ground-manifold average or the lineshape
    normalisation would be a factor of order 2 or 3."""
    p = dict(Elem=elem, Dline=dline, T=50.)
    np.testing.assert_allclose(_isotropic_integral_per_atom(p),
                               _expected_per_atom(elem, dline), rtol=0.01)


# --- what the integral must not depend on ----------------------------------------

BASE = dict(Elem='Rb', Dline='D2', T=50.)

VARIATIONS = {
    'B = 100 G': dict(Bfield=100.),
    'B = 1 kG': dict(Bfield=1000.),
    'B = 1 T': dict(Bfield=10000.),
    'B = 1 T, no Boltzmann factor': dict(Bfield=10000., BoltzmannFactor=False),
    'B oblique (angles do not enter chi)': dict(Bfield=1000., Btheta=1.1, Bphi=0.4),
    'GammaBuf = 500 MHz': dict(GammaBuf=500.),
    'shift = 3 GHz': dict(shift=3000.),
    'Doppler-free (DoppTemp = -272 C)': dict(Constrain=False, DoppTemp=-272.),
    'hot Doppler width (DoppTemp = 300 C)': dict(Constrain=False, DoppTemp=300.),
    'buffer gas, 50 Torr N2': dict(BufferGas='N2', BufferPressure=50.),
    'pure Rb87': dict(rb85frac=0.),
    'pure Rb85': dict(rb85frac=100.),
}


@pytest.mark.parametrize('variation', sorted(VARIATIONS))
def test_integrated_absorption_is_invariant(variation):
    """Against the literature value, and against the zero-field,
    natural-abundance integral to the strength the weak-transition cut can
    move (CUT): far tighter than the literature allows, so it catches strength
    lost or created by the change itself."""
    p = dict(BASE, **VARIATIONS[variation])
    value = _isotropic_integral_per_atom(p)
    np.testing.assert_allclose(value, _expected_per_atom('Rb', 'D2'), rtol=0.01)
    np.testing.assert_allclose(value, _isotropic_integral_per_atom(BASE), rtol=CUT)


@pytest.mark.parametrize('elem, dline', [('Rb', 'D2'), ('Cs', 'D1'), ('K', 'D2'),
                                         ('Na', 'D1')])
@pytest.mark.parametrize('B', [0., 30.])
def test_each_manifold_carries_a_third_of_the_strength(elem, dline, B):
    """An unpolarised ensemble (equal ground populations: no Boltzmann
    factor) absorbs sigma+, sigma- and pi light equally in total, whatever
    the field: Tr(|q><q| P_J') is the same for every q because the projector
    on the upper fine-structure level is a scalar. Fine-structure mixing by
    the field breaks this at order muB B / (fine-structure splitting), 1e-4
    for Na at 30 G, so the test stays at low field; the weak-transition cut
    (CUT) sets the tolerance."""
    p = dict(Elem=elem, Dline=dline, T=50., Bfield=B, BoltzmannFactor=False)
    integrals = np.array([_integral(component.imag) for component in _chi(p)])
    np.testing.assert_allclose(integrals / integrals.mean(), 1, rtol=0, atol=CUT)


# --- the optically thin limit ----------------------------------------------------

THIN_CASES = [
    # (element, line, Bfield, Btheta, Bphi, E_in, use_non_orthogonal)
    ('Rb', 'D2', 0., 0., 0., [1, 0, 0], False),
    ('Rb', 'D2', 500., 0., 0., [1, 0, 0], False),
    ('Rb', 'D2', 500., 0., 0., [1, 1j, 0], False),
    ('Rb', 'D2', 500., np.pi / 2, 0.3, [0.6, 0.8, 0], False),
    ('Rb', 'D2', 500., 1.0, -0.7, [0.6, 0.8j, 0], True),
    ('Cs', 'D1', 2000., 2.2, 1.9, [1, np.exp(0.4j), 0], True),
    ('Na', 'D2', 100., 0.4, 0., [0, 1, 0], True),
]


@pytest.mark.parametrize('elem, dline, B, theta, phi, E_in, nonortho', THIN_CASES)
def test_optically_thin_absorption_follows_the_sum_rule(elem, dline, B, theta, phi,
                                                        E_in, nonortho):
    """1 - S0 -> alpha L for a thin sample (peak optical depth ~1e-3 here, so
    the quadratic term is ~1e-3 of the total)."""
    T, L = 20., 1e-4
    p = dict(Elem=elem, Dline=dline, T=T, lcell=L, Bfield=B, Btheta=theta, Bphi=phi,
             use_non_orthogonal=nonortho)
    E_in = np.asarray(E_in, dtype=complex)
    E_in /= np.linalg.norm(E_in)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        [S0] = calculate(X, E_in, p, outputs=['S0'])
    N = DENSITY[elem](T + 273.15)
    absorbed = _integral(1 - np.asarray(S0)) / (N * L)
    np.testing.assert_allclose(absorbed, _expected_per_atom(elem, dline), rtol=0.01)


# Each isotope alone: the oscillator strength is the element's line strength,
# whatever the isotope (the isotope shifts and hyperfine structure only move
# it). The pure-K40 and pure-K41 cells take calc_chi branches that no other
# test reached (found by the coverage report).
ISOTOPE_CELLS = {
    'Rb85': ('Rb', dict(rb85frac=100.)),
    'Rb87': ('Rb', dict(rb85frac=0.)),
    'K39': ('K', dict(K40frac=0., K41frac=0.)),
    'K40': ('K', dict(K40frac=100., K41frac=0.)),
    'K41': ('K', dict(K40frac=0., K41frac=100.)),
}


@pytest.mark.parametrize('dline', ['D1', 'D2'])
@pytest.mark.parametrize('isotope', sorted(ISOTOPE_CELLS))
def test_integrated_absorption_of_an_isotopically_pure_cell(isotope, dline):
    elem, fractions = ISOTOPE_CELLS[isotope]
    p = dict(Elem=elem, Dline=dline, T=50., **fractions)
    np.testing.assert_allclose(_isotropic_integral_per_atom(p),
                               _expected_per_atom(elem, dline), rtol=0.01)
