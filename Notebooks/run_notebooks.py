"""Execute the tutorial notebooks against the installed ElecSus (CI).

    python Notebooks/run_notebooks.py [Notebooks/*.ipynb]

Each notebook runs in a fresh kernel from a scratch directory, so that it
imports the installed package rather than the source tree, and the executed
copies are written there; the originals are not touched. A notebook fails on
the first cell that raises. Needs nbclient and ipykernel.
"""
import glob
import os
import sys
import tempfile
import time

import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))


def main(paths):
    paths = paths or sorted(glob.glob(os.path.join(HERE, '*.ipynb')))
    if not paths:
        sys.exit('no notebooks found')
    failed = []
    with tempfile.TemporaryDirectory() as scratch:
        for path in paths:
            name = os.path.basename(path)
            nb = nbformat.read(path, as_version=4)
            client = NotebookClient(nb, timeout=1200, kernel_name='python3',
                                    resources={'metadata': {'path': scratch}})
            start = time.time()
            try:
                client.execute()
            except Exception as exc:          # CellExecutionError, DeadKernelError, ...
                failed.append(name)
                print(f'FAIL {name} ({time.time() - start:.0f} s)\n{exc}', flush=True)
                continue
            nbformat.write(nb, os.path.join(scratch, name))
            print(f'ok   {name} ({time.time() - start:.0f} s)', flush=True)
    if failed:
        sys.exit(f'{len(failed)} of {len(paths)} notebooks failed: {", ".join(failed)}')
    print(f'all {len(paths)} notebooks ran')


if __name__ == '__main__':
    main(sys.argv[1:])
