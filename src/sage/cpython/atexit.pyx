"""
Utilities for interfacing with the standard library's atexit module.
"""

# ****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import atexit


__all__ = ['restore_atexit']


cdef class restore_atexit:
    r"""
    Context manager that restores the state of the atexit module to its
    previous state when exiting the context.

    INPUT:

    - ``run`` -- boolean (default: ``False``); if ``True``, when exiting the
      context (but before restoring the old exit functions), run all
      atexit functions which were added inside the context

    - ``clear`` -- boolean (default: equal to ``run``); if ``True``, clear
      already registered atexit handlers upon entering the context

    .. WARNING::

        The combination ``run=True`` and ``clear=False`` will cause
        already-registered exit functions to be run twice: once when
        exiting the context and again when exiting Python.

    EXAMPLES:

    For this example we will wrap the entire example with
    ``restore_atexit(clear=True)`` so as to start with a fresh atexit
    module state for the sake of the example.

    Note that the function ``atexit._run_exitfuncs()`` runs all registered
    handlers, and then clears the list of handlers, so we can use it to test
    manipulation of the ``atexit`` state::

        sage: import atexit
        sage: from sage.cpython.atexit import restore_atexit
        sage: def handler(*args, **kwargs):
        ....:     import sys
        ....:     # see https://github.com/sagemath/sage/issues/25270#comment:56
        ....:     sys.stdout.write(str((args, kwargs)))
        ....:     sys.stdout.write('\n')
        sage: atexit.register(handler, 1, 2, c=3)
        <function handler at 0x...>
        sage: atexit.register(handler, 4, 5, d=6)
        <function handler at 0x...>
        sage: with restore_atexit(clear=True):
        ....:     atexit._run_exitfuncs()  # Should be none registered
        ....:     atexit.register(handler, 1, 2, c=3)
        ....:     with restore_atexit():
        ....:         atexit._run_exitfuncs()  # Run just registered handler
        ....:     atexit._run_exitfuncs()  # Handler should be run again
        <function handler at 0x...>
        ((1, 2), {'c': 3})
        ((1, 2), {'c': 3})

    We test the ``run`` option::

        sage: with restore_atexit(run=True):
        ....:     # this handler is run when exiting the context
        ....:     _ = atexit.register(handler, 7, 8, e=9)
        ((7, 8), {'e': 9})
        sage: with restore_atexit(clear=False, run=True):
        ....:     # original handlers are run when exiting the context
        ....:     pass
        ((4, 5), {'d': 6})
        ((1, 2), {'c': 3})

    The original handlers are still in place::

        sage: atexit._run_exitfuncs()
        ((4, 5), {'d': 6})
        ((1, 2), {'c': 3})

    TESTS::

        sage: from sage.cpython.atexit import (_get_exithandlers,
        ....:                                  _clear_exithandlers)
        sage: atexit.register(handler, 1, 2, c=3)
        <function handler at 0x...>
        sage: atexit.register(handler, 4, 5, d=6)
        <function handler at 0x...>
        sage: print("Initial exit handlers:\n{}".format(_get_exithandlers()))
        Initial exit handlers:
        [(<function handler at 0x...>, (1, 2), {'c': 3}),
         (<function handler at 0x...>, (4, 5), {'d': 6})]

        sage: with restore_atexit():
        ....:     pass
        sage: print("After restore_atexit:\n{}".format(_get_exithandlers()))
        After restore_atexit:
        [(<function handler at 0x...>, (1, 2), {'c': 3}),
         (<function handler at 0x...>, (4, 5), {'d': 6})]

        sage: with restore_atexit(clear=True):
        ....:     print("Exit handlers in context manager: {}".format(
        ....:           _get_exithandlers()))
        Exit handlers in context manager: []

        sage: print("After restore_atexit with clear=True:\n{}".format(
        ....:       _get_exithandlers()))
        After restore_atexit with clear=True:
        [(<function handler at 0x...>, (1, 2), {'c': 3}),
         (<function handler at 0x...>, (4, 5), {'d': 6})]
        sage: _clear_exithandlers()
        sage: _get_exithandlers()
        []
    """

    cdef list _exithandlers
    cdef bint _clear, _run

    def __init__(self, *, run=False, clear=None):
        self._clear = self._run = run
        if clear is not None:
            self._clear = clear
        self._exithandlers = None

    def __enter__(self):
        self._exithandlers = _get_exithandlers()
        if self._clear:
            _clear_exithandlers()

        return self

    def __exit__(self, *exc):
        if self._run:
            atexit._run_exitfuncs()
        _set_exithandlers(self._exithandlers)

from cpython.ref cimport PyObject
import sys

# Implement a uniform interface for getting atexit callbacks
cdef extern from *:
    """
    #ifndef Py_BUILD_CORE
    #define Py_BUILD_CORE
    #endif
    #undef _PyGC_FINALIZED
    #include "internal/pycore_interp.h"
    #include "internal/pycore_pystate.h"
    
    // Always define this struct for Cython's use
    typedef struct {
        PyObject *func;
        PyObject *args;
        PyObject *kwargs;
    } atexit_callback_struct;
    
    #if PY_VERSION_HEX >= 0x030e0000
    // Python 3.14+: atexit uses a PyList stored in state->callbacks
    // Note: In Python 3.14 the atexit_state struct changed - callbacks is now a PyObject* (PyList)
    
    static PyObject* get_atexit_callbacks_list(PyObject *self) {
        PyInterpreterState *interp = _PyInterpreterState_GET();
        // Access the callbacks list directly from the interpreter state
        // We return a new reference because Cython expects an owned reference
        PyObject *callbacks = interp->atexit.callbacks;
        Py_XINCREF(callbacks);
        return callbacks;
    }
    
    // Dummy function for Python 3.14+ (never called)
    static atexit_callback_struct** get_atexit_callbacks_array(PyObject *self) {
        PyErr_SetString(PyExc_RuntimeError, "Python >= 3.14 has no atexit arrays");
        return NULL;
    }
    #else
    // Python < 3.14: atexit uses C array
    static atexit_callback_struct** get_atexit_callbacks_array(PyObject *self) {
        PyInterpreterState *interp = _PyInterpreterState_GET();
        struct atexit_state state = interp->atexit;
        // Cast from atexit_callback** to our struct type
        return (atexit_callback_struct**)state.callbacks;
    }
    
    // Dummy function for Python < 3.14 (never called)
    static PyObject* get_atexit_callbacks_list(PyObject *self) {
        PyErr_SetString(PyExc_RuntimeError, "Python < 3.14 has no atexit lists");
        return NULL;
    }
    #endif
    """
    # Declare both functions - they exist in all Python versions (one is dummy)
    object get_atexit_callbacks_list(object module)
    
    ctypedef struct atexit_callback_struct:
        PyObject* func
        PyObject* args
        PyObject* kwargs
    atexit_callback_struct** get_atexit_callbacks_array(object module) except NULL


def _get_exithandlers():
    """Return list of exit handlers registered with the atexit module."""
    cdef list exithandlers = []
    cdef atexit_callback_struct ** callbacks
    cdef atexit_callback_struct callback
    cdef int idx
    cdef object kwargs
    
    # Python 3.14+ uses a PyList directly
    if sys.version_info >= (3, 14):
        callbacks_list = get_atexit_callbacks_list(atexit)
        if callbacks_list is None:
            return exithandlers
        # callbacks is a list of tuples: [(func, args, kwargs), ...]
        # Normalize kwargs to ensure it's always a dict (not None)
        # Note: In Python 3.14+, atexit stores callbacks in LIFO order
        # (most recently registered first), but we return them in FIFO
        # order (registration order) for consistency with earlier versions
        for item in reversed(callbacks_list):
            func, args, kwargs = item
            if kwargs is None:
                kwargs = {}
            exithandlers.append((func, args, kwargs))
    else:
        # Python < 3.14 uses C array
        callbacks = get_atexit_callbacks_array(atexit)
        for idx in range(atexit._ncallbacks()):
            callback = callbacks[idx][0]
            if callback.kwargs:
                kwargs = <object>callback.kwargs
            else:
                kwargs = {}
            exithandlers.append((<object>callback.func,
                                 <object>callback.args,
                                 kwargs))
    return exithandlers


def _set_exithandlers(exithandlers):
    """
    Replace the list of exit handlers registered with the atexit module
    with a new list.
    """

    # Clear the existing list
    atexit._clear()

    # We could do this more efficiently by directly rebuilding the array
    # of atexit_callbacks, but this is much simpler
    # Note: exithandlers is in registration order (FIFO).
    # In Python 3.14+, atexit.register prepends to the list (LIFO),
    # so registering in forward order gives us the correct execution order.
    for callback in exithandlers:
        atexit.register(callback[0], *callback[1], **callback[2])


def _clear_exithandlers():
    """Clear the atexit module of all registered exit handlers."""

    atexit._clear()
