# sage_setup: distribution = sagemath-objects
"""
Type internals
"""

from cpython.object cimport Py_TYPE, PyTypeObject, Py_TPFLAGS_HEAPTYPE

cdef PyTypeObject* PyInstance_Type = NULL
try:
    from types import InstanceType
    PyInstance_Type = <PyTypeObject*>InstanceType
except ImportError:
    pass


cpdef bint can_assign_class(obj) noexcept:
    """
    Can we assign ``obj.__class__``?

    Note that Python 3.5 has experimented with allowing assigning
    ``__class__`` in more cases but this was mostly reverted. In this
    function, we apply the Python 2.7 semantics.

    EXAMPLES::

        sage: class X: pass
        sage: from sage.cpython.type import can_assign_class
        sage: can_assign_class(X())
        True
        sage: class Y(int): pass
        sage: from sage.cpython.type import can_assign_class
        sage: can_assign_class(Y())
        True
        sage: can_assign_class(1)
        False
    """
    cdef PyTypeObject* tp = Py_TYPE(obj)
    if tp is PyInstance_Type:
        return True
    return (tp.tp_flags & Py_TPFLAGS_HEAPTYPE) != 0
