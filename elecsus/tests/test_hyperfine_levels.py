"""
Atomic energy levels against closed forms and measured constants (a convention-free physics test).

- **Breit-Rabi.** For a J = 1/2 level the hyperfine-Zeeman Hamiltonian
  A I.J + muB (gJ Jz + gI Iz) B splits into 2 x 2 blocks of fixed
  m = mI + mJ, whose eigenvalues are the Breit-Rabi formula (Breit and Rabi,
  Phys. Rev. 38, 2082 (1931); e.g. Steck, Sec. 3.3). The ground levels the
  code diagonalises must follow it at any field, with the measured A and gI.
- **gJ.** The code's ground Hamiltonian uses the free-electron g-factor for
  gJ(S1/2). The measured values differ by up to 1.1e-4 (Cs); a level can then
  move by at most |gJ - g_s| muB B / 2, which is pinned too.
- **Excited levels.** At zero field the hyperfine levels of P1/2 and P3/2 are
  E_F = A K / 2 + B_hf [3K(K+1)/2 - 2I(I+1)J(J+1)] / [2I(2I-1) 2J(2J-1)],
  K = F(F+1) - I(I+1) - J(J+1), with the measured A and B (Steck; Tiecke).
  P1/2 is also a J = 1/2 level, so it follows Breit-Rabi with
  gJ = 4/3 - g_s/3 up to the field's mixing with P3/2.

Levels are compared relative to their centroid, which removes the isotope
shift and the fine-structure offset, neither of which is tested here.
"""
import numpy as np
import pytest

from elecsus.libs import EigenSystem
from elecsus.tests import literature as lit

ISOTOPES = sorted(lit.ISOTOPES)
MU_B = lit.MU_B_MHZ_PER_G

# Fine-structure splittings P3/2 - P1/2 [MHz] (Steck; Tiecke), only to bound
# the field's mixing of P1/2 with P3/2, (muB B)^2 / splitting.
FINE_STRUCTURE = {'Rb': 7.123e6, 'Cs': 16.610e6, 'K': 1.7301e6, 'Na': 0.5156e6}

# P1/2 Breit-Rabi with the free-electron g_s and gL = 1, as in the code.
G_J_P12 = 4 / 3 - lit.G_S / 3


def _element(isotope):
    return isotope.rstrip('0123456789')


def _breit_rabi(I, A, gJ, gI, B):
    """The 2(2I+1) energies of a J = 1/2 level [MHz], sorted, with centroid
    zero. A < 0 (K40) inverts the hyperfine structure; the formula holds."""
    dE = A * (I + 0.5)
    x = (gJ - gI) * MU_B * B / dE
    levels = []
    for m in np.arange(-I - 0.5, I + 1.5):
        if abs(m) == I + 0.5:
            # stretched state: one level, linear in B (the square root's
            # argument is (1 +- x)^2, and its sign matters past x = 1)
            levels.append(dE * I / (2 * I + 1) + np.sign(m) * (gJ / 2 + gI * I) * MU_B * B)
            continue
        root = np.sqrt(1 + 4 * m * x / (2 * I + 1) + x * x)
        for sign in (+1, -1):
            levels.append(-dE / (2 * (2 * I + 1)) + gI * MU_B * m * B + sign * dE / 2 * root)
    return np.sort(levels)


def _zero_field_hyperfine(I, J, A, B_hf):
    """The (2I+1)(2J+1) zero-field energies of a level of angular momentum J
    [MHz], each F repeated 2F+1 times, sorted, centroid zero."""
    levels = []
    for F in np.arange(abs(I - J), I + J + 1):
        K = F * (F + 1) - I * (I + 1) - J * (J + 1)
        E = A * K / 2
        if J > 0.5 and I > 0.5:
            E += B_hf * (1.5 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)) / (
                2 * I * (2 * I - 1) * 2 * J * (2 * J - 1))
        levels += [E] * int(round(2 * F + 1))
    return np.sort(levels)


def _centred(levels):
    levels = np.sort(np.asarray(levels, dtype=float))
    return levels - levels.mean()


def _ground(isotope, B):
    return _centred(EigenSystem.Hamiltonian(isotope, 'D2', 1.0, B).groundEnergies)


def _excited(isotope, J, B):
    """The code's P_J levels: D1 is the lowest (2I+1)(2S+1) of the D1
    Hamiltonian's excited manifold, D2 the rest of the D2 Hamiltonian's (the
    rows FreqStren takes for each line)."""
    H = EigenSystem.Hamiltonian(isotope, 'D1' if J == 0.5 else 'D2', 1.0, B)
    E = np.sort(H.excitedEnergies)
    return _centred(E[:H.ds] if J == 0.5 else E[H.ds:])


# --- the ground state ------------------------------------------------------------

@pytest.mark.parametrize('isotope', ISOTOPES)
@pytest.mark.parametrize('B', [0., 10., 300., 3000., 30000.])
def test_ground_levels_follow_breit_rabi(isotope, B):
    """With the model's g_s, exactly: 1e-6 MHz is round-off on levels of up to
    ~50 GHz. B = 0 is the code's 1e-5 G."""
    c = lit.ISOTOPES[isotope]
    expected = _breit_rabi(c['I'], c['A_S'], lit.G_S, c['gI'], max(B, 1e-5))
    np.testing.assert_allclose(_ground(isotope, B), expected, rtol=0, atol=1e-6)


@pytest.mark.parametrize('isotope', ISOTOPES)
@pytest.mark.parametrize('B', [300., 3000., 30000.])
def test_ground_levels_with_the_measured_gJ(isotope, B):
    """The cost of using g_s for gJ(S1/2): dE/dgJ = muB B <Jz>, |<Jz>| <= 1/2,
    so no level can be off by more than |gJ - g_s| muB B / 2 (1.5 MHz for Cs
    at 1 T). Anything larger is a different error."""
    c = lit.ISOTOPES[isotope]
    expected = _breit_rabi(c['I'], c['A_S'], c['gJ'], c['gI'], B)
    bound = abs(c['gJ'] - lit.G_S) * MU_B * B / 2
    assert np.max(np.abs(_ground(isotope, B) - expected)) <= bound * (1 + 1e-6) + 1e-6


# --- the excited states ----------------------------------------------------------

def _excited_cases():
    # K41 P1/2 failed here and below until it was fixed: its A was the ground state's
    for isotope in ISOTOPES:
        for J in (0.5, 1.5):
            yield pytest.param(isotope, J, id=f'{isotope}-P{int(2 * J)}/2')


@pytest.mark.parametrize('isotope, J', list(_excited_cases()))
def test_zero_field_excited_hyperfine_levels(isotope, J):
    """Within 1 % of the level's hyperfine width, which is wider than the
    spread between published constants (literature.py) and far narrower than
    a wrong constant."""
    c, e = lit.ISOTOPES[isotope], lit.EXCITED[isotope]
    A, B_hf = (e['A_P12'], 0.) if J == 0.5 else (e['A_P32'], e['B_P32'])
    expected = _zero_field_hyperfine(c['I'], J, A, B_hf)
    width = expected.max() - expected.min()
    np.testing.assert_allclose(_excited(isotope, J, 0.), expected, rtol=0, atol=0.01 * width)


@pytest.mark.parametrize('isotope', ISOTOPES)
@pytest.mark.parametrize('B', [10., 100., 1000.])
def test_p12_levels_follow_breit_rabi(isotope, B):
    """With the measured A(P1/2), so to 1 % of the hyperfine width as above,
    plus (muB B)^2 / (fine-structure splitting) for the field's mixing with
    P3/2 (at most 4 MHz, Na at 1 kG). A wrong excited-state g-factor would
    be off by a fraction of muB B, 1.4 GHz at 1 kG."""
    c, e = lit.ISOTOPES[isotope], lit.EXCITED[isotope]
    expected = _breit_rabi(c['I'], e['A_P12'], G_J_P12, c['gI'], B)
    tol = (0.01 * abs(e['A_P12']) * (2 * c['I'] + 1)
           + (MU_B * B) ** 2 / FINE_STRUCTURE[_element(isotope)])
    np.testing.assert_allclose(_excited(isotope, 0.5, B), expected, rtol=0, atol=tol)
