# sage_setup: distribution = sagemath-repl
"""
Quitting interfaces
"""

################################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of (any version of) the GNU
#  General Public License (GPL). The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
################################################################################

import os
import subprocess
import sys
from typing import TYPE_CHECKING

from sage.env import DOT_SAGE, HOSTNAME
from sage.misc.cachefunc import cached_function

if TYPE_CHECKING:
    from weakref import ReferenceType

    from sage.interfaces.expect import Expect


@cached_function
def sage_spawned_process_file() -> str:
    """
    EXAMPLES::

        sage: from sage.interfaces.quit import sage_spawned_process_file
        sage: len(sage_spawned_process_file()) > 1
        True
    """
    # This is the old value of SAGE_TMP. Until sage-cleaner is
    # completely removed, we need to leave these spawned_processes
    # files where sage-cleaner will look for them.
    d = os.path.join(DOT_SAGE, "temp", HOSTNAME, str(os.getpid()))
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "spawned_processes")


def register_spawned_process(pid: int, cmd: str = "") -> None:
    """
    Write a line to the ``spawned_processes`` file with the given
    ``pid`` and ``cmd``.
    """
    if cmd != "":
        cmd = cmd.strip().split()[0]
    # This is safe, since only this process writes to this file.
    try:
        with open(sage_spawned_process_file(), "a") as file:
            file.write("%s %s\n" % (pid, cmd))
    except OSError:
        pass
    else:
        # If sage is being used as a python library, we need to launch
        # the cleaner ourselves upon being told that there will be
        # something to clean.
        from sage.interfaces.cleaner import start_cleaner

        start_cleaner()


expect_objects: list[ReferenceType[Expect]] = []


def expect_quitall(verbose: bool = False) -> None:
    """
    EXAMPLES::

        sage: sage.interfaces.quit.expect_quitall()
        sage: gp.eval('a=10')                                                           # needs sage.libs.pari
        '10'
        sage: gp('a')                                                                   # needs sage.libs.pari
        10
        sage: sage.interfaces.quit.expect_quitall()
        sage: gp('a')                                                                   # needs sage.libs.pari
        a
        sage: sage.interfaces.quit.expect_quitall(verbose=True)                         # needs sage.libs.pari
        Exiting PARI/GP interpreter with PID ... running .../gp --fast --emacs --quiet --stacksize 10000000
    """
    for reference in expect_objects:
        process = reference()
        if process is not None:
            try:
                process.quit(verbose=verbose)
            except RuntimeError:
                pass
    kill_spawned_jobs()


def kill_spawned_jobs(verbose: bool = False):
    """
    INPUT:

    - ``verbose`` -- boolean (default: ``False``); if ``True``, display a
      message each time a process is sent a kill signal

    EXAMPLES::

        sage: gp.eval('a=10')                                                           # needs sage.libs.pari
        '10'
        sage: sage.interfaces.quit.kill_spawned_jobs(verbose=False)
        sage: sage.interfaces.quit.expect_quitall()
        sage: gp.eval('a=10')                                                           # needs sage.libs.pari
        '10'
        sage: sage.interfaces.quit.kill_spawned_jobs(verbose=True)                      # needs sage.libs.pari
        Killing spawned job ...

    After doing the above, we do the following to avoid confusion in other doctests::

        sage: sage.interfaces.quit.expect_quitall()
    """
    fname = sage_spawned_process_file()
    if not os.path.exists(fname):
        return

    with open(fname) as file:
        for line in file:
            i = line.find(" ")
            pid = line[:i].strip()
            try:
                if verbose:
                    print("Killing spawned job %s" % pid)
                if sys.platform == "win32":
                    # From https://stackoverflow.com/a/47756757/873661
                    subprocess.call(["taskkill", "/F", "/T", "/PID", pid])
                else:
                    os.killpg(int(pid), 9)
            except OSError:
                pass


def is_running(pid: int) -> bool:
    """
    Return ``True`` if and only if there is a process with id pid running.
    """
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def invalidate_all() -> None:
    """
    Invalidate all of the expect interfaces.

    This is used, e.g., by the fork-based ``@parallel`` decorator.

    EXAMPLES::

        sage: # needs sage.libs.pari sage.symbolic
        sage: a = maxima(2); b = gp(3)
        sage: a, b
        (2, 3)
        sage: sage.interfaces.quit.invalidate_all()
        sage: a
        (invalid Maxima object -- The maxima session in which this object was defined is no longer running.)
        sage: b
        (invalid PARI/GP interpreter object -- The pari session in which this object was defined is no longer running.)

    However the maxima and gp sessions should still work out, though with their state reset::

        sage: a = maxima(2); b = gp(3)                                                  # needs sage.libs.pari sage.symbolic
        sage: a, b                                                                      # needs sage.libs.pari sage.symbolic
        (2, 3)
    """
    for reference in expect_objects:
        process = reference()
        if process:
            process.detach()
