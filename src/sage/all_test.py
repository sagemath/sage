import os
import subprocess
import sys


def test_import_sage_all_in_fresh_interpreter():
    # Prepare a clean-ish environment, without any sage variables
    env = os.environ.copy()
    env.pop("PYTHONSTARTUP", None)
    for var in list(env):
        if var.startswith("SAGE_"):
            env.pop(var, None)

    proc = subprocess.run(
        [sys.executable, "-c", "import sage.all"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        text=True, check=False,
    )
    assert proc.returncode == 0, (
        "Importing 'sage.all' in a fresh interpreter failed.\n"
        f"Return code: {proc.returncode}\n"
        f"Stdout:\n{proc.stdout}\n"
        f"Stderr:\n{proc.stderr}"
    )
    assert proc.stderr == ""
    assert proc.stdout == ""
