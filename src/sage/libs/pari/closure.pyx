"""
Convert Python functions to PARI closures

AUTHORS:

- Jeroen Demeyer (2015-04-10): initial version, :trac:`18052`.

EXAMPLES::

    sage: def the_answer():
    ....:     return 42
    sage: f = pari(the_answer)
    sage: f()
    42

    sage: cube = pari(lambda i: i^3)
    sage: cube.apply(range(10))
    [0, 1, 8, 27, 64, 125, 216, 343, 512, 729]
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
from __future__ import absolute_import

from cpython.tuple cimport *
from cpython.object cimport PyObject_Call
from cpython.ref cimport Py_INCREF

include "cysignals/signals.pxi"
from .paridecl cimport *

from .pari_instance cimport pari_instance
from .gen cimport objtogen


cdef inline GEN call_python_func_impl "call_python_func"(GEN* args, object py_func) except NULL:
    """
    Call ``py_func(*args)`` where ``py_func`` is a Python function
    and ``args`` is an array of ``GEN``s terminated by ``NULL``.

    The arguments are converted from ``GEN`` to a Sage ``gen`` before
    calling ``py_func``. The result is converted back to a PARI ``GEN``.
    """
    # How many arguments are there?
    cdef Py_ssize_t n = 0
    while args[n] != NULL:
        n += 1

    # Construct a Python tuple for args
    cdef tuple t = PyTuple_New(n)
    cdef Py_ssize_t i
    for i in range(n):
        a = pari_instance.new_gen_noclear(args[i])
        Py_INCREF(a)  # Need to increase refcount because the tuple steals it
        PyTuple_SET_ITEM(t, i, a)

    # Call the Python function
    r = PyObject_Call(py_func, t, <dict>NULL)

    # Convert the result to a GEN and copy it to the PARI stack
    # (with a special case for None)
    if r is None:
        return gnil
    return gcopy(objtogen(r).g)

# We rename this function to be able to call it with a different
# signature. In particular, we want manual exception handling and we
# implicitly convert py_func from a PyObject* to an object.
cdef extern from *:
    GEN call_python_func(GEN* args, PyObject* py_func)


cdef GEN call_python(GEN arg1, GEN arg2, GEN arg3, GEN arg4, GEN arg5, ulong py_func):
    """
    This function, which will be installed in PARI, is a front-end for
    ``call_python_func_impl``.

    It has 5 optional ``GEN``s as argument and one ``ulong``.
    This last argument is actually a Python callable object cast to
    ``ulong``.
    """
    # Convert arguments to a NULL-terminated array. From PARI's point
    # of view, all these arguments are optional: if an argument is not
    # given, PARI will pass NULL as argument and the array will
    # terminate sooner.
    cdef GEN args[6]
    args[0] = arg1
    args[1] = arg2
    args[2] = arg3
    args[3] = arg4
    args[4] = arg5
    args[5] = NULL

    sig_block()
    # Disallow interrupts during the Python code inside
    # call_python_func_impl(). We need to do this because this function
    # is very likely called within sig_on() and interrupting arbitrary
    # Python code is bad.
    cdef GEN r = call_python_func(args, <PyObject*>py_func)
    sig_unblock()
    if not r:  # An exception was raised
        sig_error()
    return r

# Install the function "call_python" for use in the PARI library.
cdef entree* ep_call_python = install(<void*>call_python, "call_python", "DGDGDGDGDGU")


cpdef gen objtoclosure(f):
    """
    Convert a Python function (more generally, any callable) to a PARI
    ``t_CLOSURE``.

    .. NOTE::

        With the current implementation, the function can be called
        with at most 5 arguments.

    .. WARNING::

        The function ``f`` which is called through the closure cannot
        be interrupted. Therefore, it is advised to use this only for
        simple functions which do not take a long time.

    EXAMPLES::

        sage: from sage.libs.pari.closure import objtoclosure
        sage: mul = objtoclosure(lambda i,j: i*j)
        sage: mul
        (v1,v2,v3,v4,v5)->call_python(v1,v2,v3,v4,v5,...)
        sage: mul.type()
        't_CLOSURE'
        sage: mul(6,9)
        54
        sage: def printme(x):
        ....:     print(x)
        sage: objtoclosure(printme)('matid(2)')
        [1, 0; 0, 1]

    Test various kinds of errors::

        sage: mul(4)
        Traceback (most recent call last):
        ...
        TypeError: <lambda>() takes exactly 2 arguments (1 given)
        sage: mul(None, None)
        Traceback (most recent call last):
        ...
        ValueError: Cannot convert None to pari
        sage: mul(*range(100))
        Traceback (most recent call last):
        ...
        PariError: call_python: too many parameters in user-defined function call
        sage: mul([1], [2])
        Traceback (most recent call last):
        ...
        PariError: call_python: forbidden multiplication t_VEC (1 elts) * t_VEC (1 elts)
    """
    sig_on()
    # Convert f to a t_INT containing the address of f
    cdef GEN f_int = utoi(<ulong><PyObject*>f)
    # Create a t_CLOSURE which calls call_python() with py_func equal to f
    cdef gen c = pari_instance.new_gen(snm_closure(ep_call_python, mkvec(f_int)))
    c.refers_to = {0:f}  # c needs to keep a reference to f
    return c
