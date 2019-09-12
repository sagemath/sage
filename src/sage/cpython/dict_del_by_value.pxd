from cpython.object cimport PyObject
cdef extern from "Python.h":
    ctypedef struct PyDictObject

cdef int del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash) except -1

# PyDict_GetItemWithError is underscore-prefixed in Python 2
IF PY_MAJOR_VERSION == 2:
    cdef extern from "Python.h":
        PyObject* PyDict_GetItemWithError "_PyDict_GetItemWithError"(dict op, object key) except? NULL
ELSE:
    cdef extern from "Python.h":
        PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
