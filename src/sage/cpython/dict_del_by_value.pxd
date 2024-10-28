# sage_setup: distribution = sagemath-objects
from cpython.object cimport PyObject
cdef extern from "Python.h":
    ctypedef struct PyDictObject

cdef int del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash) except -1

cdef extern from "Python.h":
    PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
