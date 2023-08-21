from pathlib import Path
import subprocess


from sage.env import SAGE_SRC

input_file = Path(SAGE_SRC) / "conftest_inputtest.py"


class TestOldDoctestSageScript:
    """Run `sage --t`."""

    def test_invoke_on_inputtest_file(self):
        result = subprocess.run(
            ["sage", "-t", input_file],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1  # There are failures in the input test
        assert (
            "Failed example:\n"
            "    something()\n"
            "Expected:\n"
            "    44\n"
            "Got:\n"
            "    42\n"
            "" in result.stdout
        )


class TestPytestSageScript:
    """Run `sage --pytest`."""

    def test_invoke_on_inputtest_file(self):
        result = subprocess.run(
            ["sage", "--pytest", input_file],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1  # There are failures in the input test
        assert (
            "004     EXAMPLES::\n"
            "005 \n"
            "006         sage: something()\n"
            "007         42\n"
            "008         sage: something() + 1\n"
            "009         43\n"
            "010 \n"
            "011     TESTS::\n"
            "012 \n"
            "013         sage: something()\n"
            "Expected:\n"
            "    44\n"
            "Got:\n"
            "    42\n"
            "" in result.stdout
        )
