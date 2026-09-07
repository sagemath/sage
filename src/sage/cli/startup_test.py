"""Exercise startup configuration through fresh Python CLI processes."""

import json
import os
import shlex
import stat
import subprocess
import sys

import pytest


@pytest.fixture
def startup_env(tmp_path):
    env = os.environ.copy()
    for name in ("SAGE_RC_FILE", "SAGE_STARTUP_FILE", "SAGE_ENV_SOURCED"):
        env.pop(name, None)
    for name in (
        "DOT_SAGE", "IPYTHONDIR", "MPLCONFIGDIR", "XDG_CACHE_HOME",
        "JUPYTER_CONFIG_DIR",
    ):
        directory = tmp_path / name.lower()
        directory.mkdir()
        env[name] = str(directory)
    return env


def run_cli(env, *args, **kwargs):
    return subprocess.run(
        [sys.executable, "-m", "sage.cli", *args],
        env=env, capture_output=True, text=True, check=False, timeout=60, **kwargs,
    )


@pytest.mark.parametrize("mode", ["command", "file", "interactive"])
def test_default_sagerc(startup_env, tmp_path, mode):
    """The rc affects every entry mode before Sage caches its environment."""
    startup_env["SAGE_STARTUP_TEST_REMOVED"] = "remove me"
    dot_sage = tmp_path / "dot_sage"
    (dot_sage / "sagerc").write_text(
        "export SAGE_STARTUP_TEST_VALUE='value with spaces'\n"
        "export SAGE_DOC_SERVER_URL='https://sagerc.example.invalid/'\n"
        "unset SAGE_STARTUP_TEST_REMOVED\n"
        "printf x >> \"$DOT_SAGE/read-count\"\n"
        "echo sagerc-startup-message\n"
    )
    command = (
        "import json, os, sage.env; "
        "print(json.dumps([os.environ['SAGE_STARTUP_TEST_VALUE'], "
        "os.environ.get('SAGE_STARTUP_TEST_REMOVED'), "
        "sage.env.SAGE_DOC_SERVER_URL]))"
    )
    if mode == "command":
        result = run_cli(startup_env, "-c", command)
    elif mode == "file":
        source = tmp_path / "startup.sage"
        source.write_text(command + "\n")
        result = run_cli(startup_env, str(source))
    else:
        result = run_cli(
            startup_env, "-q", "--simple-prompt", input=command + "\nexit\n",
        )
    assert result.returncode == 0, result.stderr
    assert json.dumps([
        "value with spaces", None, "https://sagerc.example.invalid/",
    ]) in result.stdout
    assert "sagerc-startup-message" not in result.stdout
    assert "sagerc-startup-message" in result.stderr
    assert (dot_sage / "read-count").read_text() == "x"


def test_custom_sagerc_and_arguments(startup_env, tmp_path):
    """Custom paths and shell positional parameters preserve CLI arguments."""
    (tmp_path / "dot_sage" / "sagerc").write_text("exit 99\n")
    rc = tmp_path / "custom 'quoted' rc"
    rc.write_text(
        "set -- 'replacement arguments'\n"
        "export SAGE_STARTUP_TEST_VALUE=custom\n"
    )
    startup_env["SAGE_RC_FILE"] = str(rc)
    source = tmp_path / "script with spaces.py"
    source.write_text(
        "import json, os, sys\n"
        "print(json.dumps([os.environ['SAGE_STARTUP_TEST_VALUE'], sys.argv[1:]]))\n"
    )
    arguments = ["two words", "single ' and double \" quotes", "$(echo wrong)"]
    result = run_cli(startup_env, str(source), *arguments)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == ["custom", arguments]


@pytest.mark.parametrize(("contents", "diagnostic"), [
    pytest.param("return 7\n", True, id="return"),
    pytest.param("if\n", True, id="syntax-error"),
    pytest.param(
        "set -e\nfalse\nexport SAGE_STARTUP_TEST_VALUE=wrong\n",
        False, id="errexit",
    ),
    pytest.param(
        "readonly _SAGE_RC_SOURCED=blocked\n", False, id="restart-marker-failure",
    ),
])
def test_sagerc_failure_stops_startup(startup_env, tmp_path, contents, diagnostic):
    rc = tmp_path / "dot_sage" / "sagerc"
    rc.write_text(contents)
    result = run_cli(startup_env, "-c", "print('payload-ran')")
    assert result.returncode != 0
    assert "payload-ran" not in result.stdout
    if diagnostic:
        assert "Error sourcing" in result.stderr
        assert str(rc) in result.stderr


def test_sagerc_preserves_shell_state(startup_env, tmp_path):
    working_directory = tmp_path / "working directory"
    working_directory.mkdir()
    (tmp_path / "dot_sage" / "sagerc").write_text(
        f"cd {shlex.quote(str(working_directory))}\numask 027\n"
    )
    result = run_cli(
        startup_env, "-c",
        "import os; open('created-file', 'w').close(); print(os.getcwd())",
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(working_directory)
    assert stat.S_IMODE((working_directory / "created-file").stat().st_mode) == 0o640


@pytest.mark.parametrize("script_name", ["console-wrapper.py", "-sage"])
def test_relative_console_script_and_interpreter_options(startup_env, tmp_path, script_name):
    """Restart a console script after cd without losing Python options."""
    working_directory = tmp_path / "working directory"
    working_directory.mkdir()
    (tmp_path / "dot_sage" / "sagerc").write_text(
        f"cd {shlex.quote(str(working_directory))}\n"
    )
    (tmp_path / script_name).write_text(
        "import sys\n"
        "from sage.cli.__main__ import main\n"
        "sys.exit(main())\n"
    )
    command = (
        "import json, os, sys; "
        "print(json.dumps([sys._xoptions.get('utf8'), "
        "'ignore' in sys.warnoptions, os.getcwd()]))"
    )
    result = subprocess.run(
        [sys.executable, "-W", "ignore", "-X", "utf8", "--", script_name,
         "-c", command],
        env=startup_env, cwd=tmp_path, capture_output=True, text=True,
        check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == [True, True, str(working_directory)]


def test_import_does_not_source_sagerc(startup_env, tmp_path):
    marker = tmp_path / "dot_sage" / "imported"
    (tmp_path / "dot_sage" / "sagerc").write_text(
        "touch \"$DOT_SAGE/imported\"\nexit 98\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", "import sage.cli; import sage.cli.__main__"],
        env=startup_env, capture_output=True, text=True, check=False, timeout=60,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


@pytest.mark.parametrize("unreadable", [False, True])
def test_unavailable_custom_sagerc_is_skipped(startup_env, tmp_path, unreadable):
    """An unavailable explicit rc does not fall back to the default file."""
    (tmp_path / "dot_sage" / "sagerc").write_text("exit 99\n")
    rc = tmp_path / "unavailable-rc"
    startup_env["SAGE_RC_FILE"] = str(rc)
    if unreadable:
        rc.write_text("exit 98\n")
        rc.chmod(0)
        if os.access(rc, os.R_OK):
            pytest.skip("the current user can read a file with no read permissions")
    result = run_cli(startup_env, "-c", "print('started')")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "started"


@pytest.mark.parametrize("current", [True, False])
def test_legacy_launcher_marker(startup_env, tmp_path, current):
    """Only the current legacy environment has already loaded the rc."""
    (tmp_path / "dot_sage" / "sagerc").write_text(
        "export SAGE_STARTUP_TEST_VALUE=loaded\n"
    )
    startup_env["SAGE_STARTUP_TEST_VALUE"] = "inherited"
    startup_env["SAGE_ENV_SOURCED"] = (
        "6:" + ":".join(startup_env.get(name, "") for name in (
            "SAGE_LOCAL", "SAGE_VENV", "SAGE_SRC",
        ))
        if current else "6:/old/local:/old/venv:/old/src"
    )
    result = run_cli(
        startup_env, "-c", "import os; print(os.environ['SAGE_STARTUP_TEST_VALUE'])",
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == ("inherited" if current else "loaded")


def test_child_cli_reads_its_own_sagerc(startup_env, tmp_path):
    """A restart guard must not suppress configuration in a child Sage."""
    (tmp_path / "dot_sage" / "sagerc").write_text(
        "export SAGE_STARTUP_TEST_VALUE=parent\n"
    )
    child_directory = tmp_path / "child"
    child_directory.mkdir()
    (child_directory / "sagerc").write_text(
        "export SAGE_STARTUP_TEST_VALUE=child\n"
        "printf x >> \"$DOT_SAGE/read-count\"\n"
    )
    child_command = "import os; print(os.environ['SAGE_STARTUP_TEST_VALUE'])"
    command = (
        "import os, subprocess, sys; "
        "print(os.environ['SAGE_STARTUP_TEST_VALUE'], flush=True); "
        f"env = dict(os.environ, DOT_SAGE={str(child_directory)!r}); "
        f"subprocess.run([sys.executable, '-m', 'sage.cli', '-c', {child_command!r}], "
        "env=env, check=True)"
    )
    result = run_cli(startup_env, "-c", command)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == ["parent", "child"]
    assert (child_directory / "read-count").read_text() == "x"
