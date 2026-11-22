"""[ES4] The version has to be the same in all three places that declare it.

`pyproject.toml`, `CITATION.cff` and `.zenodo.json` each carry a version string,
and they have drifted before: at the point this test was written they read
4.0.0, 4.0.0a2 and 4.0.0a1 respectively. Citation tooling reads the latter two
and installers read the first, so a mismatch means the DOI, the citation and the
installed package disagree about which release they are.

Bumping a version is therefore a three-file edit, and this test is what tells
you when you have only done two of them.

Skipped when the repository files are not present (e.g. running against an
installed wheel, which ships none of them).
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PYPROJECT = REPO / 'pyproject.toml'
CITATION = REPO / 'CITATION.cff'
ZENODO = REPO / '.zenodo.json'

pytestmark = pytest.mark.skipif(
    not (PYPROJECT.exists() and CITATION.exists() and ZENODO.exists()),
    reason='release metadata files are not present outside the repository')


def _pyproject_version():
    for line in PYPROJECT.read_text().splitlines():
        match = re.match(r'\s*version\s*=\s*["\']([^"\']+)["\']', line)
        if match:
            return match.group(1)
    raise AssertionError('no version field in pyproject.toml')


def _citation_version():
    for line in CITATION.read_text().splitlines():
        match = re.match(r'\s*version:\s*(\S+)', line)
        if match:
            return match.group(1).strip('"\'')
    raise AssertionError('no version field in CITATION.cff')


def _zenodo_version():
    return json.loads(ZENODO.read_text())['version']


def test_all_three_declared_versions_agree():
    versions = {
        'pyproject.toml': _pyproject_version(),
        'CITATION.cff': _citation_version(),
        '.zenodo.json': _zenodo_version(),
    }
    assert len(set(versions.values())) == 1, (
        'version strings disagree: '
        + ', '.join(f'{name} = {value}' for name, value in sorted(versions.items()))
        + '. Bumping a release means editing all three.')


def test_zenodo_metadata_is_valid_json_with_the_fields_zenodo_needs():
    data = json.loads(ZENODO.read_text())
    for field in ('title', 'version', 'upload_type', 'license', 'creators'):
        assert field in data, f'.zenodo.json is missing {field!r}'
    assert data['creators'], '.zenodo.json lists no creators'
    for creator in data['creators']:
        assert 'name' in creator, f'creator entry without a name: {creator!r}'


def test_citation_cff_declares_the_fields_github_needs():
    text = CITATION.read_text()
    for field in ('cff-version:', 'title:', 'authors:', 'version:'):
        assert field in text, f'CITATION.cff is missing {field}'


def test_the_installed_version_is_exposed_as_dunder_version():
    """`elecsus.__version__` is the installed distribution's version."""
    from importlib.metadata import version
    import elecsus
    assert elecsus.__version__ == version('elecsus4')
    assert re.fullmatch(r'\d+(\.\d+)*(a|b|rc)?\d*(\.dev\d+)?', elecsus.__version__)


def test_the_distribution_is_elecsus4_and_the_package_is_elecsus():
    """PyPI's 'ElecSus' is another project's (two empty 2.0 pre-releases), so
    this one installs as `pip install elecsus4` and is imported as `elecsus`."""
    text = PYPROJECT.read_text()
    assert re.search(r'^name = "elecsus4"', text, re.M)
    assert re.search(r'^include = \["elecsus\*"\]', text, re.M)


def _family_given(name):
    """'S. A. Wrathmall' -> ('Wrathmall', 'S. A.')"""
    given, family = name.rsplit(' ', 1)
    return family, given


def test_the_version_4_authors_agree_everywhere():
    """pyproject.toml's authors are version 4's. CITATION.cff and .zenodo.json
    list them first, in the same order, before the authors of versions 1-3,
    and AUTHORS lists them under "Version 4". The four lists disagreed before
    4.0.0b1."""
    block = re.search(r'^authors = \[(.*?)^\]', PYPROJECT.read_text(), re.M | re.S)
    assert block, 'pyproject.toml has no authors = [...] list'
    v4 = [_family_given(n) for n in re.findall(r'name = "([^"]+)"', block[1])]
    assert v4, 'pyproject.toml lists no authors'
    citation = re.findall(r'^  - family-names: (.+)\n    given-names: (.+)$',
                          CITATION.read_text(), re.M)
    assert citation[:len(v4)] == v4, f'CITATION.cff starts {citation[:len(v4)]}, not {v4}'
    zenodo = [tuple(c['name'].split(', ', 1)) for c in json.loads(ZENODO.read_text())['creators']]
    assert zenodo[:len(v4)] == v4, f'.zenodo.json starts {zenodo[:len(v4)]}, not {v4}'
    authors = REPO / 'AUTHORS'
    if authors.exists():
        listed = re.search(r'^Version 4\b.*\n((?:  \S.*\n)+)', authors.read_text(), re.M)
        assert listed, 'AUTHORS has no "Version 4" list'
        assert [_family_given(n.strip()) for n in listed[1].splitlines()] == v4


def test_the_readme_cites_the_software_by_its_doi():
    """CITATION.cff's DOI (Zenodo's concept DOI, shared by every version) is
    the one the README's Citation section gives."""
    doi = re.search(r'^doi:\s*(\S+)', CITATION.read_text(), re.M)
    readme = REPO / 'README.md'
    if not doi or not readme.exists():
        pytest.skip('no DOI in CITATION.cff, or no README here')
    assert f'https://doi.org/{doi[1]}' in readme.read_text(encoding='utf-8')


def test_the_changelog_opens_with_this_version():
    """The newest entry of CHANGELOG.md is the version in pyproject.toml, or
    'Unreleased' while that is a development version (4.0.0b2.dev0, between
    releases)."""
    changelog = REPO / 'CHANGELOG.md'
    if not changelog.exists():
        pytest.skip('no CHANGELOG.md here')
    newest = re.search(r'^## (\S+)', changelog.read_text(encoding='utf-8'), re.M)
    assert newest, 'CHANGELOG.md has no "## <version>" entry'
    version = _pyproject_version()
    want = 'Unreleased' if '.dev' in version else version
    assert newest[1] == want, (
        f'the newest CHANGELOG entry is {newest[1]!r}; with version {version} it '
        f'should be {want!r}')


def test_no_tracked_file_holds_a_merge_conflict_marker():
    """elecsus/libs/NOTICE shipped '<<<<<<< HEAD' for years."""
    import shutil
    import subprocess
    if not (REPO / '.git').exists() or not shutil.which('git'):
        pytest.skip('needs a git checkout')
    # the pattern is built here so that this file does not match itself
    pattern = '^(' + '<' * 7 + '|' + '>' * 7 + ') '
    done = subprocess.run(['git', 'grep', '-nIE', pattern, '--', '.',
                           ':(exclude)*/reviews/*'],   # review reports quote one
                          cwd=REPO,
                          capture_output=True, text=True)
    assert done.returncode == 1, 'merge-conflict markers:\n' + done.stdout
