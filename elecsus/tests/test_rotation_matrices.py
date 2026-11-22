"""[ES4] Direct tests for libs/RotationMatrices.py.

These functions were converted off the deprecated np.matrix subclass to plain
ndarrays (@ and np.linalg.inv). Two of them, rotate_forward() and
rotate_back(), are called by nothing inside the package -- spectra.py uses only
rotate_around_z() and rotate_around_z2() -- so the ES3 golden regression does
not reach them and they had no coverage at all. They are also where the
conversion did the least mechanical thing: `.I` became np.linalg.inv, and
np.matrix's silent collapse of a trailing length-1 axis (which is how the
stacked (n,3) input path used to work) had to be reproduced explicitly by
_to_2d().

A wrong rotation returns a plausible spectrum rather than an exception, so this
module tests the properties directly -- orthogonality, the round trip, known
axis images -- plus reference values captured from the pre-conversion np.matrix
implementation, which the ndarray version reproduces bit-for-bit (verified over
300 random (phi, theta) with complex input vectors before these were frozen).
"""
import numpy as np
import pytest

from elecsus.libs import RotationMatrices as RM

PHI, THETA = 0.4, 0.9
V = np.array([0.3, -0.7, 0.5])

# Captured from the np.matrix implementation at 6b05dfc^.
REF_FORWARD = np.array([[-0.2051804643325424, -0.8572884174951348,
                         0.23012506661589405]])
REF_BACK = np.array([[0.7607582281295733], [-0.45003352464769425],
                     [0.22071870112504144]])
REF_Z = np.array([[0.6428504158494283], [-0.4083421884182642], [0.5]])


# --- reference values ---------------------------------------------------

def test_rotate_forward_matches_reference():
    np.testing.assert_allclose(RM.rotate_forward(V, PHI, THETA), REF_FORWARD,
                               rtol=0, atol=1e-15)


def test_rotate_back_matches_reference():
    np.testing.assert_allclose(RM.rotate_back(V, PHI, THETA), REF_BACK,
                               rtol=0, atol=1e-15)


def test_rotate_around_z_matches_reference():
    np.testing.assert_allclose(RM.rotate_around_z(V, 0.6), REF_Z,
                               rtol=0, atol=1e-15)


# --- properties ---------------------------------------------------------

def test_forward_then_back_is_the_identity():
    """rotate_back must invert rotate_forward: this is what `.I` was for."""
    out = RM.rotate_back(RM.rotate_forward(V, PHI, THETA)[0], PHI, THETA)
    np.testing.assert_allclose(out.ravel(), V, atol=1e-14)


@pytest.mark.parametrize('phi, theta', [(0., 0.), (0.4, 0.9), (-1.2, 2.7),
                                        (np.pi / 2, np.pi / 3)])
def test_rotations_preserve_length(phi, theta):
    """Rotations are orthogonal; nothing here should change |v|."""
    for out in (RM.rotate_forward(V, phi, theta),
                RM.rotate_back(V, phi, theta),
                RM.rotate_around_z(V, phi)):
        np.testing.assert_allclose(np.linalg.norm(out), np.linalg.norm(V),
                                   rtol=1e-14)


def test_forward_sends_the_field_axis_to_z():
    """The point of rotate_forward: B along (cos(phi)sin(theta), sin(phi),
    cos(phi)cos(theta)) must come out along z' in the rotated frame."""
    phi, theta = 0.4, 0.9
    b = [np.cos(phi) * np.sin(theta), np.sin(phi), np.cos(phi) * np.cos(theta)]
    out = RM.rotate_forward(b, phi, theta).ravel()
    np.testing.assert_allclose(out, [0., 0., 1.], atol=1e-14)


def test_rotate_around_z_is_a_z_rotation():
    x_rotated = RM.rotate_around_z([1., 0., 0.], np.pi / 2).ravel()
    np.testing.assert_allclose(x_rotated, [0., 1., 0.], atol=1e-15)
    z_fixed = RM.rotate_around_z([0., 0., 1.], 1.234).ravel()
    np.testing.assert_allclose(z_fixed, [0., 0., 1.], atol=1e-15)


def test_zero_angle_is_a_no_op():
    np.testing.assert_allclose(RM.rotate_around_z(V, 0.).ravel(), V, atol=1e-15)
    np.testing.assert_allclose(RM.rotate_forward(V, 0., 0.).ravel(), V, atol=1e-15)


# --- shapes, including the stacked path spectra.py uses ------------------

def test_shapes_match_the_np_matrix_behaviour():
    """np.matrix collapsed a trailing length-1 axis; _to_2d reproduces it.
    These shapes are what spectra.get_Efield reshapes, so they are load-bearing."""
    assert RM.rotate_forward(V, PHI, THETA).shape == (1, 3)
    assert RM.rotate_back(V, PHI, THETA).shape == (3, 1)
    assert RM.rotate_around_z(V, PHI).shape == (3, 1)
    stacked = np.ones((7, 3))
    assert RM.rotate_around_z(stacked, PHI).shape == (3, 7)


def test_stacked_input_matches_per_vector_rotation():
    """The (n,3) path must agree with rotating each vector on its own."""
    rng = np.random.default_rng(0)
    vecs = rng.normal(size=(5, 3))
    stacked = RM.rotate_around_z(vecs, 0.7)
    for i, v in enumerate(vecs):
        np.testing.assert_allclose(stacked[:, i],
                                   RM.rotate_around_z(v, 0.7).ravel(), atol=1e-15)


def test_complex_fields_are_supported():
    """The propagated E-field is complex; rotation must be linear over it."""
    v = np.array([0.3 + 0.2j, -0.7 - 0.1j, 0.5 + 0.j])
    out = RM.rotate_around_z(v, 0.6).ravel()
    re = RM.rotate_around_z(v.real, 0.6).ravel()
    im = RM.rotate_around_z(v.imag, 0.6).ravel()
    np.testing.assert_allclose(out, re + 1j * im, atol=1e-15)


def test_rotate_around_z2_agrees_with_rotate_around_z():
    """rotate_around_z2 is the stacked variant used by the fast solver path."""
    rng = np.random.default_rng(1)
    vecs = rng.normal(size=(4, 3))
    a = RM.rotate_around_z2(vecs, 0.55)
    for i, v in enumerate(vecs):
        np.testing.assert_allclose(a[i].ravel(),
                                   RM.rotate_around_z(v, 0.55).ravel(), atol=1e-15)
