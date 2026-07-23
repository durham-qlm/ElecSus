"""
Generate golden reference spectra from genuine ElecSus v3.0.8.

Provenance: run against a git worktree of commit 065fa28 (the final upstream
durham-qlm/ElecSus commit, v3.0.8, 2023-09-15 - the last commit before any
ES4 modification), using an environment with numpy<2 and scipy<1.14 that
ES3 requires (e.g. the `elecsus` conda env: numpy 1.26 / scipy 1.11 / sympy 1.12).

Usage:
    git worktree add /tmp/es3_baseline 065fa28
    conda activate elecsus     # or any numpy<2, scipy<1.14 env
    python elecsus/tests/golden/generate_es3_reference.py /tmp/es3_baseline

Writes es3_reference.npz next to this script. The npz holds, per case,
the detuning axis, E_in, the parameter dict (as a repr string plus
individual keys) and the raw ES3 outputs S0,S1,S2,S3,Ix,Iy.

The consuming test (elecsus/tests/test_es3_golden.py) asserts the
documented ES3->ES4 relationships against these arrays. This script is
committed for provenance/regeneration and is NOT run by pytest.
"""

import sys
import os
import numpy as np

CASES = [
    # name, p_dict, E_in, detuning (MHz)
    # -- API grid: geometries, elements, D-lines --
    ('rb_d2_zerofield',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 20., 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'Bphi': 0.}, [1, 0, 0], (-10e3, 10e3, 400)),
    ('rb_d2_faraday',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 75e-3, 'Bfield': 300.,
      'Btheta': 0., 'Bphi': 0.}, [1, 0, 0], (-10e3, 10e3, 400)),
    ('rb_d2_voigt',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 75e-3, 'Bfield': 300.,
      'Btheta': np.pi / 2, 'Bphi': 0.}, [1, 0, 0], (-10e3, 10e3, 400)),
    ('rb_d2_oblique45',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 75e-3, 'Bfield': 300.,
      'Btheta': np.pi / 4, 'Bphi': 0.}, [1, 0, 0], (-10e3, 10e3, 400)),
    ('rb_d1_faraday',
     {'Elem': 'Rb', 'Dline': 'D1', 'T': 40., 'lcell': 75e-3, 'Bfield': 100.,
      'Btheta': 0., 'Bphi': 0.}, [1, 0, 0], (-6e3, 6e3, 400)),
    ('cs_d1_smallfield',
     {'Elem': 'Cs', 'Dline': 'D1', 'T': 25., 'lcell': 75e-3, 'Bfield': 50.,
      'Btheta': 0., 'Bphi': 0.}, [1, 0, 0], (-8e3, 8e3, 400)),
    ('k_d1_faraday',
     {'Elem': 'K', 'Dline': 'D1', 'T': 80., 'lcell': 75e-3, 'Bfield': 200.,
      'Btheta': 0., 'Bphi': 0.}, [1, 0, 0], (-6e3, 6e3, 400)),
    ('k_d2_voigt',
     {'Elem': 'K', 'Dline': 'D2', 'T': 100., 'lcell': 75e-3, 'Bfield': 500.,
      'Btheta': np.pi / 2, 'Bphi': 0.}, [1, 0, 0], (-8e3, 8e3, 400)),
    ('rb87_d2_hyperfine_pb',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 100., 'lcell': 5e-3, 'Bfield': 6000.,
      'Btheta': 0., 'Bphi': 0., 'rb85frac': 0.}, [1, 0, 0], (-25e3, 25e3, 500)),
    ('rb_d2_gammabuf',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'Bphi': 0., 'GammaBuf': 300., 'shift': -50.},
     [1, 0, 0], (-10e3, 10e3, 400)),
    ('rb_d2_pol45',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 60., 'lcell': 75e-3, 'Bfield': 300.,
      'Btheta': 0., 'Bphi': 0.},
     [1 / np.sqrt(2), 1 / np.sqrt(2), 0], (-10e3, 10e3, 400)),

    # -- published-figure configs (mirror elecsus/tests/test_spectroscopy.py
    #    and test_rotondaro.py so those tests can pin quantitative values) --
    ('siddons_rb_d2_16C',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 16.5, 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'rb85frac': 72.17}, [1, 0, 0], (-4e3, 5e3, 400)),
    ('siddons_rb_d2_25C',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 25., 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'rb85frac': 72.17}, [1, 0, 0], (-4e3, 5e3, 400)),
    ('siddons_rb_d2_36C',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 36.6, 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'rb85frac': 72.17}, [1, 0, 0], (-4e3, 5e3, 400)),
    ('weller_rb_d2_80C',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 80., 'lcell': 75e-3, 'Bfield': 0.,
      'Btheta': 0., 'rb85frac': 72.17}, [1, 0, 0], (-6e3, 6e3, 400)),
    ('zentile_cs_d1_filter',
     {'Elem': 'Cs', 'Dline': 'D1', 'T': 67.8, 'lcell': 75e-3, 'Bfield': 45.7,
      'Btheta': 0.}, [1, 0, 0], (-8e3, 8e3, 400)),
    ('keaveney_rb_d2_filter',
     {'Elem': 'Rb', 'Dline': 'D2', 'T': 90.3, 'lcell': 5e-3, 'Bfield': 248.5,
      'Btheta': 0., 'rb85frac': 72.17}, [1, 0, 0], (-8e3, 8e3, 400)),
    ('rotondaro_cs_d2_faraday',
     {'Elem': 'Cs', 'Dline': 'D2', 'T': 58., 'lcell': 75e-3, 'Bfield': 300.,
      'Btheta': 0.}, [1, 0, 0], (-15e3, 15e3, 500)),
    ('rotondaro_cs_d2_voigt',
     {'Elem': 'Cs', 'Dline': 'D2', 'T': 84., 'lcell': 75e-3, 'Bfield': 700.,
      'Btheta': np.pi / 2}, [1 / np.sqrt(2), 1 / np.sqrt(2), 0], (-15e3, 15e3, 500)),
    ('rotondaro_cs_d2_ob87',
     {'Elem': 'Cs', 'Dline': 'D2', 'T': 100., 'lcell': 75e-3, 'Bfield': 500.,
      'Btheta': 87 * np.pi / 180, 'Bphi': 0.}, [1, 0, 0], (-15e3, 15e3, 300)),
]

OUTPUTS = ['S0', 'S1', 'S2', 'S3', 'Ix', 'Iy']


def main(es3_root):
    es3_root = os.path.abspath(es3_root)
    sys.path.insert(0, es3_root)
    from elecsus.elecsus_methods import calculate  # noqa: E402 (ES3 import)
    import elecsus as es3_pkg
    assert es3_root in os.path.abspath(es3_pkg.__file__), \
        f'imported elecsus from {es3_pkg.__file__}, not the ES3 worktree'

    data = {}
    names = []
    for name, p, E_in, (lo, hi, n) in CASES:
        X = np.linspace(lo, hi, n)
        print(f'{name} ...')
        res = calculate(X, E_in, dict(p), outputs=OUTPUTS)
        data[f'{name}__X'] = X
        data[f'{name}__E_in'] = np.asarray(E_in, dtype=complex)
        data[f'{name}__p'] = np.bytes_(repr(p).encode())
        for key, arr in zip(OUTPUTS, res):
            data[f'{name}__{key}'] = np.asarray(arr).squeeze().astype(complex)
        names.append(name)

    data['case_names'] = np.array(names)
    data['es3_commit'] = np.bytes_(b'065fa28 (durham-qlm/ElecSus v3.0.8, 2023-09-15)')
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'es3_reference.npz')
    np.savez_compressed(out, **data)
    print(f'\nWrote {out} ({os.path.getsize(out)/1024:.0f} kB, {len(names)} cases)')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: python generate_es3_reference.py <path-to-ES3-worktree>')
    main(sys.argv[1])
