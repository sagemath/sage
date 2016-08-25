"""
The Sage Input Hook

This is a hook into the IPython input prompt and will be called
periodically (every 100ms) while Python is sitting idle. We use it to
reload attached files if they have changed.

IPython has analogous code to set an input hook, but we are not using
their implementation. For once, it unsets signal handlers which will
disable Ctrl-C.
"""

###########################################################################
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include "cysignals/signals.pxi"

cdef extern from 'pythonrun.h':
    int (*PyOS_InputHook)() nogil except -1

cdef extern from 'intrcheck.h':
    int PyOS_InterruptOccurred() nogil

from cpython.exc cimport PyErr_SetInterrupt

import sage.repl.attach


cdef int c_sage_inputhook() nogil except -1:
    """
    This is the C function that is installed as PyOS_InputHook
    """
    if PyOS_InterruptOccurred():   # clears interrupt
        PyErr_SetInterrupt()       # re-set
    else:
        with gil:
            sage_inputhook()
            sig_check()

def install():
    """
    Install the Sage input hook

    EXAMPLES::

        sage: from sage.repl.inputhook import install
        sage: install()
    """
    global PyOS_InputHook
    PyOS_InputHook = c_sage_inputhook

def uninstall():
    """
    Uninstall the Sage input hook

    EXAMPLES::

        sage: from sage.repl.inputhook import uninstall
        sage: uninstall()
    """
    global PyOS_InputHook
    PyOS_InputHook = NULL


def is_installed():
    r"""
    Test whether the Sage input hook is installed

    This is only for doctesting purposes

    EXAMPLES::

        sage: from sage.repl.inputhook import is_installed
        sage: is_installed()
        False

    The Sage input hook is only installed while files are attached::

        sage: tmp = tmp_filename(ext='.py')
        sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
        sage: from sage.repl.attach import attach, detach
        sage: attach(tmp)
        sage: is_installed()
        True
        sage: detach(tmp)
        sage: is_installed()
        False
    """
    global PyOS_InputHook
    return (PyOS_InputHook == c_sage_inputhook)


def sage_inputhook():
    r"""
    The input hook.

    This function will be called every 100ms when IPython is idle at
    the command prompt.

    EXAMPLES::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: tmp = tmp_filename(ext='.py')
        sage: f = open(tmp, 'w'); f.write('a = 2\n'); f.close()
        sage: shell.run_cell('%attach ' + tmp)
        sage: shell.run_cell('a')
        2
        sage: sleep(1)  # filesystem timestamp granularity
        sage: f = open(tmp, 'w'); f.write('a = 3\n'); f.close()

    Note that the doctests are never really at the command prompt, so
    we call the input hook manually::

        sage: shell.run_cell('from sage.repl.inputhook import sage_inputhook')
        sage: shell.run_cell('sage_inputhook()')
        ### reloading attached file tmp_....py modified at ... ###
        0

        sage: shell.run_cell('a')
        3
        sage: shell.run_cell('detach({0})'.format(repr(tmp)))
        sage: shell.run_cell('attached_files()')
        []
        sage: shell.quit()
    """
    sage.repl.attach.reload_attached_files_if_modified()
    return 0





