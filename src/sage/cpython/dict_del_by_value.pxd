from cpython.object cimport PyObject
cdef extern from "Python.h":
    ctypedef struct PyDictObject

cdef del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash)

#This is for compatibility: this routine is available in Py3 and we
#implement it ourselves in Py2.
IF PY_VERSION_HEX<=0x02ffffff:
    cdef PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
ELSE:
    cdef extern from "Python.h":
        PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
