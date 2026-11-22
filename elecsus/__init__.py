from importlib.metadata import PackageNotFoundError, version as _version

try:
    # 'elecsus4' on PyPI; 'ElecSus' is the name it had before, which an older
    # editable install may still carry
    __version__ = _version('elecsus4')
except PackageNotFoundError:
    try:
        __version__ = _version('ElecSus')
    except PackageNotFoundError:   # run from a source tree that was never installed
        __version__ = '0+unknown'
