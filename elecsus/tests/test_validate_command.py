"""[ES4] `python -m elecsus.validate` says why it failed.

pytest is in the 'test' extra only, and the suites run in a subprocess whose
stderr used to be thrown away: without pytest the command printed
"[FAIL] ES3 golden regression" and an empty line, and the re-run command it
suggested could not work either.
"""
import subprocess
import sys

from elecsus import validate


def test_without_pytest_it_says_so_and_names_the_extra(monkeypatch, capsys):
    monkeypatch.setattr(validate, 'pytest_missing', lambda: True)
    monkeypatch.setattr(sys, 'argv', ['validate', '--no-figures'])
    assert validate.main() == 2
    out = capsys.readouterr().out
    assert 'pytest' in out and 'elecsus4[test]' in out


def test_pytest_is_found_when_it_is_installed():
    assert validate.pytest_missing() is False


def _proc(returncode, stdout='', stderr=''):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


def test_a_suite_that_could_not_start_shows_stderr():
    detail = validate._failure_detail(
        _proc(1, '', '/usr/bin/python: No module named pytest\n'))
    assert 'No module named pytest' in detail


def test_a_failing_suite_shows_pytests_summary_line():
    detail = validate._failure_detail(
        _proc(1, 'F\nFAILED x.py::t\n1 failed, 3 passed in 0.1s\n', 'noise'))
    assert detail == '1 failed, 3 passed in 0.1s'


def test_a_crash_with_output_but_no_summary_shows_stderr():
    detail = validate._failure_detail(
        _proc(2, 'collecting ...\n', 'Traceback\nImportError: boom\n'))
    assert 'ImportError: boom' in detail


def test_run_suites_reports_the_reason(monkeypatch):
    monkeypatch.setattr(
        subprocess, 'run',
        lambda *a, **k: _proc(1, '', 'No module named pytest'))
    results = validate.run_suites('/nowhere')
    assert results and all(not ok and 'No module named pytest' in tail
                           for _, ok, tail in results)
