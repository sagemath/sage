# sage.doctest: optional - numpy
r"""
Utility functions for numpy.
"""

cimport numpy as np
import numpy as np
from sage.libs.m4ri cimport *
from libc.stdint cimport uintptr_t


ctypedef fused numpy_integral:
    np.int8_t
    np.int32_t
    np.int64_t


cdef set_mzd_from_numpy_unsafe(mzd_t* entries, np.ndarray[numpy_integral, ndim=1] x):
    """
    Internal function.
    Caller are responsible for checking the two arrays have the same length.
    """
    for i in range(len(x)):
        mzd_write_bit(entries, 0, i, x[i] & 1)


def set_mzd_from_numpy(uintptr_t entries_addr, Py_ssize_t degree, x):
    """
    Set the entries in ``entries`` from numpy array ``x``.

    INPUT:

    - ``entries_addr`` -- must be a ``mzd_t*`` casted to ``uintptr_t``; the casting
      is necessary to pass it through Python boundary because of lazy import

    - ``degree`` -- the length of the array

    - ``x`` -- a numpy array of integers or booleans, or any other object (in which
      case this function will return ``False``)

    OUTPUT: ``True`` if successful, ``False`` otherwise. May throw ``ValueError``.
    """
    cdef Py_ssize_t i
    cdef np.ndarray[np.npy_bool, ndim=1] x_bool
    cdef mzd_t* entries = <mzd_t*>entries_addr
    if isinstance(x, np.ndarray):
        if x.ndim != 1:
            raise ValueError("numpy array must have dimension 1")
        if x.shape[0] != degree:
            raise ValueError("numpy array must have the right length")
        if x.dtype == np.int8:
            set_mzd_from_numpy_unsafe(entries, <np.ndarray[np.int8_t, ndim=1]>x)
            return True
        if x.dtype == np.int32:
            set_mzd_from_numpy_unsafe(entries, <np.ndarray[np.int32_t, ndim=1]>x)
            return True
        if x.dtype == np.int64:
            set_mzd_from_numpy_unsafe(entries, <np.ndarray[np.int64_t, ndim=1]>x)
            return True
        if x.dtype == np.bool_:
            x_bool = x
            for i in range(degree):
                mzd_write_bit(entries, 0, i, x_bool[i])
            return True
    return False
