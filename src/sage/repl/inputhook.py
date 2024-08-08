# sage_setup: distribution = sagemath-repl
"""
The Sage Input Hook

This lets us perform actions while IPython is sitting at the terminal input
prompt. We use it to reload attached files if they have changed.
"""

###########################################################################
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
###########################################################################

import select
import errno
import contextlib
import io

from IPython.core.getipython import get_ipython
from IPython.terminal.pt_inputhooks import register

import sage.repl.attach


TIMEOUT = 0.25   # seconds


def sage_inputhook(context):
    f = context.fileno()
    while True:
        sage.repl.attach.reload_attached_files_if_modified()
        try:
            r, _, _ = select.select([f], [], [], TIMEOUT)
            if f in r:
                return  # IPython signalled us to stop
        except OSError as e:
            if e[0] != errno.EINTR:
                raise


register('sage', sage_inputhook)


def install():
    """
    Install the Sage input hook.

    EXAMPLES:

    Make sure ipython is running so we really test this function::

        sage: from sage.repl.interpreter import get_test_shell
        sage: get_test_shell()
        <sage.repl.interpreter.SageTestShell object at ...>

    Run the function twice, to check it is idempotent (see :issue:`35235`)::

        sage: from sage.repl.inputhook import install
        sage: install()
        sage: install()
    """
    ip = get_ipython()
    if not ip:
        return   # Not running in ipython, e.g. doctests
    if ip._inputhook != sage_inputhook:
        # silence `ip.enable_gui()` useless output
        with contextlib.redirect_stdout(io.StringIO()):
            ip.enable_gui('sage')


def uninstall():
    """
    Uninstall the Sage input hook.

    EXAMPLES::

        sage: from sage.repl.inputhook import uninstall
        sage: uninstall()
    """
    ip = get_ipython()
    if not ip:
        return
    if ip._inputhook == sage_inputhook:
        # silence `ip.enable_gui()` useless output
        with contextlib.redirect_stdout(io.StringIO()):
            ip.enable_gui(None)
