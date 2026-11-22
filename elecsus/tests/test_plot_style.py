"""The house plot style draws without warnings where Times New Roman is
not installed, and keeps Times New Roman where it is."""
import logging
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from elecsus.libs import plot_style  # noqa: E402

try:
    # The oldest matplotlib we support (3.7) calls pyparsing by names that
    # pyparsing 3.3 deprecates: matplotlib's warnings, not ours.
    from pyparsing.warnings import PyparsingDeprecationWarning
except ImportError:
    PyparsingDeprecationWarning = None


@pytest.mark.parametrize('size', ['default', 'paper'])
def test_the_style_draws_without_font_warnings(size, caplog):
    with plt.rc_context():
        plot_style.apply_style(size)
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1], label=r'$S_0$')
        ax.set_xlabel('Detuning (GHz)')
        ax.legend()
        with warnings.catch_warnings(), caplog.at_level(logging.WARNING, logger='matplotlib.font_manager'):
            warnings.simplefilter('error')
            if PyparsingDeprecationWarning is not None:
                warnings.simplefilter('ignore', PyparsingDeprecationWarning)
            fig.canvas.draw()
        plt.close(fig)
    assert not [r for r in caplog.records if 'findfont' in r.getMessage()]


def test_times_new_roman_is_used_where_it_is_installed(monkeypatch):
    monkeypatch.setattr(plot_style, '_have_font', lambda name: True)
    with plt.rc_context():
        plot_style.apply_style()
        assert plt.rcParams['font.family'] == ['Times New Roman']
        assert plt.rcParams['mathtext.rm'] == 'Times New Roman'
