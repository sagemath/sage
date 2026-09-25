"""
Standard C helper code for Cython modules
"""
# ***************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from cpython.object cimport Py_TYPE, PyTypeObject, PyObject


cdef inline PY_NEW(type t):
    """
    Return ``t.__new__(t)``.  This works even for types like
    :class:`Integer` where we change ``tp_new`` at runtime (Cython
    optimizations assume that ``tp_new`` doesn't change).

    TESTS:

    ``tp_new`` must be called with a real (empty) argument tuple, not with
    ``NULL``.  Cython's ``__cinit__`` wrapper reads the tuple's size whenever
    ``__cinit__`` takes arguments, and Cython 3.3 also reads it in the
    vectorcall adapter it installs as ``tp_new``; either way a ``NULL``
    argument tuple segfaults.  Check both a base extension type and a
    derived one, which reaches the base ``tp_new`` through the inheritance
    chain::

        sage: # needs sage.misc.cython
        sage: cython(
        ....: '''
        ....: from sage.ext.stdsage cimport PY_NEW
        ....:
        ....: cdef class Base:
        ....:     cdef public tuple stored
        ....:     def __cinit__(self, *args):
        ....:         self.stored = args
        ....:
        ....: cdef class Derived(Base):
        ....:     pass
        ....:
        ....: def new_base():
        ....:     return PY_NEW(Base)
        ....:
        ....: def new_derived():
        ....:     return PY_NEW(Derived)
        ....: ''')
        sage: new_base().stored
        ()
        sage: new_derived().stored
        ()
    """
    # tp_new requires a tuple for positional arguments.  In particular,
    # Cython 3.3's vectorcall wrapper reads its size unconditionally.
    return (<PyTypeObject*>t).tp_new(t, <PyObject*>(), <PyObject*>NULL)


cdef inline void PY_SET_TP_NEW(type dst, type src) noexcept:
    """
    Manually set ``dst.__new__`` to ``src.__new__``.  This is used to
    speed up Cython's boilerplate object construction code by skipping
    irrelevant base class ``tp_new`` methods.
    """
    (<PyTypeObject*>dst).tp_new = (<PyTypeObject*>src).tp_new


cdef inline bint HAS_DICTIONARY(obj) noexcept:
    """
    Test whether the given object has a Python dictionary.
    """
    return Py_TYPE(obj).tp_dictoffset != 0
