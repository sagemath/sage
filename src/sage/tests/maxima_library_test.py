"""Run special-function doctests with the Pexpect Maxima module unavailable."""

import subprocess
import sys
from pathlib import Path

import pytest


SAGE_SOURCE = Path(__file__).resolve().parents[1]
SOURCES = [
    path.relative_to(SAGE_SOURCE).as_posix()
    for path in sorted((SAGE_SOURCE / "functions").glob("*.py"))
    if path.name not in {"__init__.py", "all.py"}
] + [
    "calculus/calculus.py",
    "calculus/desolvers.py",
    "symbolic/expression_conversions.py",
    "symbolic/integration/external.py",
    "symbolic/maxima_wrapper.py",
]

# Install the guard before importing Sage, including its lazy imports.  Raising
# something other than ImportError prevents optional-import handling from
# silently skipping a dependency on the old interface.  Record attempts as
# well, in case the caller catches a broader exception.
RUN_WITHOUT_PEXPECT_MAXIMA = """
import importlib.abc
import sys

blocked = []

class RejectPexpectMaxima(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sage.interfaces.maxima':
            blocked.append(fullname)
            raise AssertionError('Pexpect Maxima must not be imported')

assert 'sage.interfaces.maxima' not in sys.modules
sys.meta_path.insert(0, RejectPexpectMaxima())

from sage.doctest.control import DocTestController, DocTestDefaults

options = DocTestDefaults(serial=True, long=True, stats_path=sys.argv[2])
result = DocTestController(options, [sys.argv[1]]).run()
assert not blocked, blocked
assert 'sage.interfaces.maxima' not in sys.modules
sys.exit(result)
"""


@pytest.mark.long
@pytest.mark.parametrize("source", SOURCES)
def test_without_pexpect_maxima(source, tmp_path):
    """Each file runs in a fresh interpreter, with all its long doctests."""
    result = subprocess.run(
        [sys.executable, "-c", RUN_WITHOUT_PEXPECT_MAXIMA,
         str(SAGE_SOURCE / source), str(tmp_path / "timings.json")],
        capture_output=True, text=True, timeout=600, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
