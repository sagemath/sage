import os
import subprocess
import sys

import pytest

from sage import config as sage_config


def test_cython_traceback_in_editable_install(tmp_path):
    get_editable_root = getattr(sage_config, "get_editable_root", None)
    if get_editable_root is None or get_editable_root() is None:
        pytest.skip("requires a Meson editable install")

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["DOT_SAGE"] = str(tmp_path / ".sage")
    env["IPYTHONDIR"] = str(tmp_path / "ipython")
    env["MPLCONFIGDIR"] = str(tmp_path / "matplotlib")
    code = """
from sage.repl.interpreter import get_test_shell

shell = get_test_shell()
shell.run_cell("1/0")
shell.quit()
"""

    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    assert 'raise ZeroDivisionError("rational division by zero")' in result.stdout
    assert "Could not get source" not in result.stdout


def test_import_with_config_without_editable_root(tmp_path):
    """The interpreter remains importable with an older generated config.py."""
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    code = """
import sage.config
sage.config.__dict__.pop("get_editable_root", None)
import sage.repl.interpreter
"""
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        check=True,
        timeout=60,
    )
