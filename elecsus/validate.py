# Copyright 2026 S. A. Wrathmall.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
ElecSus installation validator  [ES4 NEW]
=========================================

One command for any user to confirm their ElecSus installation reproduces
the physics it should:

    python -m elecsus.validate            # writes figures to ./elecsus_validation/
    python -m elecsus.validate -o mydir   # custom output directory

It does two things:

1. Runs the packaged validation test suites (pytest):
     - ES3 golden regression: v4 reproduces genuine ElecSus v3.0.8 output
       (S0/Ix/Iy unchanged; S1/S2/S3 under the documented renormalisation)
       across 20 configurations including oblique fields.
     - Published-figure suites: Siddons 2008, Weller 2011, Zentile 2015,
       Keaveney 2018, Rotondaro 2015 configurations with quantitative anchors.
     - Briscoe 2025 non-orthogonal theory reproduction.

2. Regenerates the classic published-figure spectra with YOUR installation
   and overlays them on the stored ES3 v3.0.8 reference, saving PNGs so you
   can see the agreement, not just trust an exit code.
"""

import argparse
import os
import subprocess
import sys

import numpy as np

TEST_SUITES = [
    ('ES3 v3.0.8 golden regression', 'tests/test_es3_golden.py'),
    ('Published spectra (Siddons/Weller/Zentile/Keaveney)', 'tests/test_spectroscopy.py'),
    ('Rotondaro 2015 filters', 'tests/test_rotondaro.py'),
    ('Briscoe 2025 non-orthogonal theory', 'tests/briscoe_2025_theory_repro.py'),
]

FIGURES = [
    # title, [(case, output, label)], ylabel
    ('Siddons 2008 Fig. 9 — Rb D2 Doppler-broadened transmission',
     [('siddons_rb_d2_16C', 'S0', 'T = 16.5 °C'),
      ('siddons_rb_d2_25C', 'S0', 'T = 25 °C'),
      ('siddons_rb_d2_36C', 'S0', 'T = 36.6 °C')],
     'Transmission $S_0$', 'siddons_2008'),
    ('Zentile 2015 — Cs D1 Faraday filter (45.7 G, 67.8 °C)',
     [('zentile_cs_d1_filter', 'Iy', 'Filter output $I_y$')],
     '$I_y$', 'zentile_2015_filter'),
    ('Keaveney 2018 — Rb D2 5 mm filter (248.5 G, 90.3 °C)',
     [('keaveney_rb_d2_filter', 'Iy', 'Filter output $I_y$')],
     '$I_y$', 'keaveney_2018_filter'),
    ('Rotondaro 2015 Fig. 3 — Cs D2 Faraday filter (300 G, 58 °C)',
     [('rotondaro_cs_d2_faraday', 'Iy', 'Filter output $I_y$')],
     '$I_y$', 'rotondaro_2015_filter'),
]


def run_suites(pkg_dir):
    results = []
    for label, rel_path in TEST_SUITES:
        path = os.path.join(pkg_dir, rel_path)
        proc = subprocess.run(
            [sys.executable, '-m', 'pytest', path, '-q', '--tb=line'],
            capture_output=True, text=True)
        tail = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ''
        results.append((label, proc.returncode == 0, tail))
    return results


def make_figures(pkg_dir, out_dir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from elecsus.elecsus_methods import calculate

    golden = np.load(os.path.join(pkg_dir, 'tests', 'golden', 'es3_reference.npz'))
    written = []
    for title, curves, ylabel, fname in FIGURES:
        fig, ax = plt.subplots(figsize=(7, 4))
        for i, (case, output, label) in enumerate(curves):
            X = golden[f'{case}__X']
            E_in = golden[f'{case}__E_in']
            p = eval(golden[f'{case}__p'].tobytes().decode())
            [arr] = calculate(X, E_in, p, outputs=[output])
            v4 = np.asarray(arr).squeeze().real
            g = golden[f'{case}__{output}'].real
            color = f'C{i}'
            ax.plot(X[::6] / 1e3, g[::6], 'o', ms=3.5, color=color, alpha=0.5,
                    label=f'{label} — ES3 v3.0.8')
            ax.plot(X / 1e3, v4, '-', lw=1.3, color=color,
                    label=f'{label} — this installation')
        ax.set_xlabel('Detuning (GHz)')
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7)
        fig.tight_layout()
        out = os.path.join(out_dir, f'{fname}.png')
        fig.savefig(out, dpi=150)
        plt.close(fig)
        written.append(out)
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog='python -m elecsus.validate',
        description='Validate an ElecSus installation against ES3 reference '
                    'data and published-figure configurations.')
    ap.add_argument('-o', '--output-dir', default='elecsus_validation',
                    help='directory for validation figures (default: ./elecsus_validation)')
    ap.add_argument('--no-figures', action='store_true',
                    help='run the test suites only, skip figure generation')
    args = ap.parse_args(argv)

    pkg_dir = os.path.dirname(os.path.abspath(__file__))

    print('ElecSus v4 installation validation')
    print('=' * 60)

    results = run_suites(pkg_dir)
    all_ok = all(ok for _, ok, _ in results)
    for label, ok, tail in results:
        mark = 'PASS' if ok else 'FAIL'
        print(f'  [{mark}]  {label}')
        if not ok:
            print(f'          {tail}')

    if not args.no_figures:
        os.makedirs(args.output_dir, exist_ok=True)
        print('-' * 60)
        for path in make_figures(pkg_dir, args.output_dir):
            print(f'  figure: {path}')

    print('=' * 60)
    if all_ok:
        print('All validation suites passed: this installation reproduces the')
        print('ES3 v3.0.8 reference physics and the published-figure spectra.')
    else:
        print('VALIDATION FAILED - see suite output above.')
        print('Re-run a single suite for detail, e.g.:')
        print(f'  python -m pytest {pkg_dir}/tests/test_es3_golden.py -v')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
