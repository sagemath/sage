# sage_setup: distribution = sagemath-objects
# Python can be built with extensive debugging support. This file lets
# Sage know about which debugging options are enabled.
#
# To enable debugging, build Sage with the SAGE_DEBUG environment
# variable set to "yes":
#
# [user@localhost:~/sage-x.y.z] export SAGE_DEBUG=yes && make


from cpython.ref cimport PyObject, PyTypeObject

cdef extern from "python_debug.h":

    # This is what is generally meant by "a debug build" of Python.
    # Implies Py_REF_DEBUG, Py_TRACE_REFS, and PYMALLOC_DEBUG (if
    # WITH_PYMALLOC is enabled).  In addition, C assert()s are enabled.
    cdef bint Py_DEBUG "SAGE_Py_DEBUG"

    # Turn on aggregate reference counting.  This arranges that extern
    # _Py_RefTotal hold a count of all references, the sum of
    # ob_refcnt across all objects.  Also checks after every decref to
    # verify that the refcount hasn't gone negative, and causes an
    # immediate fatal error if it has.
    cdef bint Py_REF_DEBUG "SAGE_Py_REF_DEBUG"

    # Heavy reference debugging. Every PyObject grows two more
    # pointers, to maintain a doubly-linked list of all live
    # heap-allocated objects. Implies Py_REF_DEBUG.
    cdef bint Py_TRACE_REFS "SAGE_Py_TRACE_REFS"

    # Conditionally call PyObject_INIT() if Py_TRACE_REFS is
    # enabled. This is necessary to initialize the aforementioned
    # double-linked list.
    void if_Py_TRACE_REFS_then_PyObject_INIT(PyObject *, PyTypeObject *)

    # When this is enabled, calls to the PyObject_ memory routines are
    # handled by Python's own small-object allocator, while calls to
    # the PyMem_ memory routines are directed to the system malloc/
    # realloc/free.
    cdef bint WITH_PYMALLOC "SAGE_WITH_PYMALLOC"

    # If this is also defined in addition to WITH_PYMALLOC, calls to
    # both PyObject_ and PyMem_ memory routines are directed to a
    # special debugging mode of Python's small-object
    # allocator. Requires WITH_PYMALLOC.
    cdef bint PYMALLOC_DEBUG "SAGE_PYMALLOC_DEBUG"
