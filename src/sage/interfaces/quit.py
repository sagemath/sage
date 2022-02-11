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

expect_objects = []


def expect_quitall(verbose=False):
    """
    EXAMPLES::

        sage: sage.interfaces.quit.expect_quitall()
        sage: gp.eval('a=10')
        '10'
        sage: gp('a')
        10
        sage: sage.interfaces.quit.expect_quitall()
        sage: gp('a')
        a
        sage: sage.interfaces.quit.expect_quitall(verbose=True)
        Exiting PARI/GP interpreter with PID ... running .../gp --fast --emacs --quiet --stacksize 10000000
    """
    for P in expect_objects:
        R = P()
        if R is not None:
            try:
                R.quit(verbose=verbose)
            except RuntimeError:
                pass
    kill_spawned_jobs()


def kill_spawned_jobs(verbose=False):
    """
    INPUT:

        - ``verbose`` -- bool (default: False); if True, display a
          message each time a process is sent a kill signal

    EXAMPLES::

        sage: gp.eval('a=10')
        '10'
        sage: sage.interfaces.quit.kill_spawned_jobs(verbose=False)
        sage: sage.interfaces.quit.expect_quitall()
        sage: gp.eval('a=10')
        '10'
        sage: sage.interfaces.quit.kill_spawned_jobs(verbose=True)
        Killing spawned job ...

    After doing the above, we do the following to avoid confusion in other doctests::

        sage: sage.interfaces.quit.expect_quitall()
    """
    from sage.misc.misc import SAGE_TMP
    file = os.path.join(SAGE_TMP, 'spawned_processes')
    if not os.path.exists(file):
        return
    with open(file) as f:
        for L in f:
            i = L.find(' ')
            pid = L[:i].strip()
            try:
                if verbose:
                    print("Killing spawned job %s" % pid)
                os.killpg(int(pid), 9)
            except OSError:
                pass


def is_running(pid):
    """
    Return True if and only if there is a process with id pid running.
    """
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, ValueError):
        return False


def invalidate_all():
    """
    Invalidate all of the expect interfaces.

    This is used, e.g., by the fork-based @parallel decorator.

    EXAMPLES::

        sage: a = maxima(2); b = gp(3)
        sage: a, b
        (2, 3)
        sage: sage.interfaces.quit.invalidate_all()
        sage: a
        (invalid Maxima object -- The maxima session in which this object was defined is no longer running.)
        sage: b
        (invalid PARI/GP interpreter object -- The pari session in which this object was defined is no longer running.)

    However the maxima and gp sessions should still work out, though with their state reset:

        sage: a = maxima(2); b = gp(3)
        sage: a, b
        (2, 3)
    """
    for I in expect_objects:
        I1 = I()
        if I1:
            I1.detach()
