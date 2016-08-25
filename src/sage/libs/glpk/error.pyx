"""
Error handler for the GLPK library
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
from .env cimport *
from cpython.exc cimport PyErr_SetObject
from sage.numerical.mip import MIPSolverException

class GLPKError(MIPSolverException):
    """
    An error raised by the GLPK library.

    EXAMPLES::

        sage: from sage.libs.glpk.error import GLPKError
        sage: raise GLPKError("trouble!")
        Traceback (most recent call last):
        ...
        GLPKError: trouble!
    """
    pass


# Global error message string
cdef error_message = ""


cdef int sage_glpk_term_hook(void *info, const char *s) with gil:
    """
    A hook to intercept all output written by GLPK.
    """
    global error_message
    if glp_at_error():
        # Save error message and skip normal printing
        error_message += s
        return 1
    else:
        # Normal non-error output: the return value 0 means that GLPK
        # will write the output as usual.
        return 0


cdef void sage_glpk_error_hook(void *info) with gil:
    """
    A hook to intercept GLPK errors.
    """
    global error_message
    PyErr_SetObject(GLPKError, error_message.strip())
    error_message = ""
    sig_error()


def setup_glpk_error_handler():
    r"""
    Setup the GLPK error handler. Called automatically when this module
    is imported at Sage startup.

    TESTS::

        sage: cython('''
        ....: #clib glpk z gmp
        ....: from sage.libs.glpk.env cimport glp_term_out
        ....: include "cysignals/signals.pxi"
        ....:
        ....: sig_on()
        ....: glp_term_out(12345)  # invalid value
        ....: sig_off()
        ....: ''')
        Traceback (most recent call last):
        ...
        GLPKError: glp_term_out: flag = 12345; invalid parameter
        Error detected in file env/stdout.c at line ...

    Check that normal terminal output still works, see :trac:`20832`::

        sage: def verbose_GLPK():
        ....:     from sage.numerical.backends.generic_backend import get_solver
        ....:     s = get_solver(solver = "GLPK")
        ....:     s.set_verbosity(2)
        ....:     return s
        sage: p = MixedIntegerLinearProgram(solver=verbose_GLPK)
        sage: x, y = p['x'], p['y']
        sage: p.add_constraint(2*x + 3*y <= 6)
        sage: p.add_constraint(3*x + 2*y <= 6)
        sage: p.add_constraint(x >= 0)
        sage: p.set_objective(x + y)
        sage: res = p.solve()
              0: obj = ...
        sage: res  # rel tol 1e-15
        2.4
    """
    glp_term_hook(sage_glpk_term_hook, NULL)
    glp_error_hook(sage_glpk_error_hook, NULL)

setup_glpk_error_handler()
