"""Tests for initialization of doctest workers."""

import subprocess
import sys
import textwrap

import pytest


@pytest.mark.skipif(sys.platform != "darwin", reason="macOS psutil initialization")
def test_cpu_time_after_fork():
    """The timer can use psutil in a worker forked from a threaded parent."""
    pytest.importorskip("psutil")
    # Start a fresh interpreter: pytest or its plugins may already have
    # imported psutil, which would hide the initialization problem.
    code = textwrap.dedent("""\
        import os
        import resource
        import threading

        from sage.doctest.forker import init_sage
        from sage.doctest.util import Timer

        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        init_sage()

        class RunningInterface:
            def __call__(self):
                return self

            def is_running(self):
                return True

            def pid(self):
                return os.getpid()

        # Make macOS enforce its Objective-C initialization checks after fork.
        threading.Thread(target=threading.Event().wait, daemon=True).start()
        pid = os.fork()
        if pid == 0:
            assert Timer()._quick_cputime([RunningInterface()]) >= 0
            os._exit(0)
        _, status = os.waitpid(pid, 0)
        assert os.waitstatus_to_exitcode(status) == 0, status
        """)
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=60, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
