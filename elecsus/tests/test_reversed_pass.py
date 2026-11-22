"""
Reversed pass, wavevector_dirn = -1.

ES3, and ES4 until this was fixed, turned B through pi and also made the
wavenumber negative on a reversed pass, so exp(i n k L) grew wherever the
vapour absorbs: transmission reached ~1e62. The Faraday rotation still came
out right, which is all the old double-pass checks compared.

The reversed pass is now a forward pass in the frame turned by pi about x,
with B given in the same lab frame for both passes. These tests pin the
physics any reversed pass must obey, not that construction:

* it absorbs: transmission stays at or below 1, and matches the forward pass
  where the geometry is reciprocal;
* in Faraday and Voigt geometry, going there and back is one pass through
  twice the length (so the Faraday rotation and the optical depth double);
* Lorentz reciprocity for any field direction, J_back(B) = J_forward(-B)^T
  in lab coordinates;
* the input check accepts only +1 and -1;
* the outputs are in lab components, about +z, on either pass.
"""
import numpy as np
import pytest

from elecsus.elecsus_methods import calculate

X = np.linspace(-6e3, 6e3, 401)
BASE = dict(Elem='Rb', Dline='D2', T=80., lcell=2e-3, Bfield=300.,
            rb85frac=72.17)
PATHS = [False, True]   # use_non_orthogonal


def _field(E_in, p):
    return np.asarray(calculate(X, E_in, p, outputs=['E_out'])[0])


def _jones(p):
    """Transverse Jones matrices, shape (len(X), 2, 2): J[:, out, in]."""
    ex = _field([1, 0, 0], p)
    ey = _field([0, 1, 0], p)
    return np.stack([ex[:2], ey[:2]], axis=1).transpose(2, 0, 1)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta', [0., np.pi / 4, np.pi / 2, 2 * np.pi / 3,
                                   np.pi])
def test_reversed_pass_absorbs(theta, nonortho):
    """75 mm at 80 C is opaque on resonance: the reversed pass must be too.
    Before the fix its maximum S0 here was ~1e62."""
    p = dict(BASE, lcell=75e-3, Btheta=theta, Bphi=0.5,
             use_non_orthogonal=nonortho)
    E_in = [0.6, 0.8, 0]
    outputs = ['S0', 'S0_P'] if nonortho else ['S0']
    back = [np.asarray(s).real for s in
            calculate(X, E_in, dict(p, wavevector_dirn=-1), outputs=outputs)]
    # the physical (flux) transmission: S0_P on the non-orthogonal path,
    # whose |E|^2-based S0 is not bounded by 1
    transmission = back[-1]
    assert np.all(np.isfinite(transmission))
    assert transmission.max() <= 1 + 1e-9
    assert transmission.min() < 1e-6     # it absorbs on resonance
    if theta in (0., np.pi / 2, np.pi):
        # Faraday and Voigt geometry are the same both ways
        fwd = np.asarray(calculate(X, E_in, p, outputs=['S0'])[0]).real
        np.testing.assert_allclose(back[0], fwd, rtol=0, atol=1e-12)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta', [0., np.pi, np.pi / 2, 3 * np.pi / 2])
def test_round_trip_is_one_pass_through_twice_the_length(theta, nonortho):
    """Faraday and Voigt geometry: the Jones matrix is the same both ways, so
    there and back through L equals one pass through 2L. That doubles the
    Faraday rotation (non-reciprocal) and the optical depth."""
    p = dict(BASE, Btheta=theta, Bphi=0.3, use_non_orthogonal=nonortho)
    E0 = np.array([0.6, 0.8j, 0])
    E1 = _field(E0, p)
    E1[2] = 0                                    # outside the cell
    E2 = _field(E1, dict(p, wavevector_dirn=-1))
    E_2L = _field(E0, dict(p, lcell=2 * BASE['lcell']))
    np.testing.assert_allclose(E2[:2], E_2L[:2], rtol=0, atol=1e-12)


@pytest.mark.parametrize('nonortho', PATHS)
@pytest.mark.parametrize('theta, phi', [
    (np.pi / 4, 0.), (np.pi / 3, 0.7), (2 * np.pi / 3, -1.1), (0.2, 2.0),
])
def test_lorentz_reciprocity(theta, phi, nonortho):
    """For a magneto-optical medium, J_back(B) = J_forward(-B)^T when both are
    written in the same lab (x, y) axes. -B has polar angle theta + pi."""
    p = dict(BASE, Btheta=theta, Bphi=phi, use_non_orthogonal=nonortho)
    J_back = _jones(dict(p, wavevector_dirn=-1))
    J_fwd_reversed_B = _jones(dict(p, Btheta=theta + np.pi))
    np.testing.assert_allclose(J_back, J_fwd_reversed_B.transpose(0, 2, 1),
                               rtol=0, atol=1e-9)


@pytest.mark.parametrize('value', [1, -1, -1.0, np.int64(-1)])
def test_wavevector_dirn_accepts_plus_and_minus_one(value):
    calculate(X[:5], [1, 0, 0], dict(BASE, wavevector_dirn=value),
              outputs=['S0'])


@pytest.mark.parametrize('value', [2, 0, 0.5, -2, True, np.True_, '-1',
                                   float('nan')])
def test_wavevector_dirn_rejects_anything_else(value):
    """2 used to double the optical depth in silence."""
    with pytest.raises(ValueError, match='wavevector_dirn'):
        calculate(X[:5], [1, 0, 0], dict(BASE, wavevector_dirn=value),
                  outputs=['S0'])


@pytest.mark.parametrize('dirn', [1, -1])
def test_circular_outputs_are_about_the_lab_z_axis(dirn):
    """[ES4] Il, Ir and S3 are the [1,+i] and [1,-i] components in
    lab x, y, as the documents now say, on a reversed pass too: there,
    [1,-i] (Ir, S3 = +1) is left-circular for the returning beam itself. In
    vacuum (a 1 nm cell far off resonance) the field is unchanged."""
    Xv = np.array([-1e5, 1e5])
    p = dict(Elem='Rb', Dline='D2', T=20., lcell=1e-9, Bfield=0., wavevector_dirn=dirn)
    S3, Il, Ir = calculate(Xv, np.array([1, -1j, 0]) / np.sqrt(2), p,
                           outputs=['S3', 'Il', 'Ir'])
    np.testing.assert_allclose(S3, 1., atol=1e-9)
    np.testing.assert_allclose(Ir, 1., atol=1e-9)
    np.testing.assert_allclose(Il, 0., atol=1e-9)
