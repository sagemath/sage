# sage_setup: distribution = sagemath-objects
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
    """
    return (<PyTypeObject*>t).tp_new(t, <PyObject*>NULL, <PyObject*>NULL)


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
