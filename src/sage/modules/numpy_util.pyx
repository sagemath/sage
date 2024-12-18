# sage.doctest: optional - numpy
# cython: fast_getattr=False
# https://github.com/cython/cython/issues/6442
r"""
Utility functions for numpy.
"""

cimport numpy as np
import numpy as np


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


cpdef int set_mzd_from_numpy(uintptr_t entries_addr, Py_ssize_t degree, x) except -1:
    """
    Set the entries in ``<mzd_t*>entries_addr`` from numpy array ``x``.

    INPUT:

    - ``entries_addr`` -- must be a ``mzd_t*`` casted to ``uintptr_t``; the casting
      is necessary to pass it through Python boundary because of lazy import.
      Do not pass arbitrary integer value here, will crash the program.

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


cpdef int _set_matrix_mod2_from_numpy_helper(Matrix_mod2_dense a, np.ndarray[numpy_integral, ndim=2] b) except -1:
    """
    Internal function, helper for :func:`set_matrix_mod2_from_numpy`.

    TESTS::

        sage: from sage.modules.numpy_util import _set_matrix_mod2_from_numpy_helper
        sage: import numpy as np
        sage: a = matrix(GF(2), 2, 3)
        sage: b = np.array([[1, 0, 1], [0, 1, 0]], dtype=np.int8)
        sage: _set_matrix_mod2_from_numpy_helper(a, b)
        1
        sage: a
        [1 0 1]
        [0 1 0]
        sage: _set_matrix_mod2_from_numpy_helper(a, np.array([[1, 0], [0, 1]], dtype=np.int8))
        Traceback (most recent call last):
        ...
        ValueError: shape mismatch
    """
    if not (a.nrows() == b.shape[0] and a.ncols() == b.shape[1]):
        raise ValueError("shape mismatch")
    for i in range(b.shape[0]):
        for j in range(b.shape[1]):
            a.set_unsafe_int(i, j, b[i, j] & 1)
    return True


cpdef int set_matrix_mod2_from_numpy(Matrix_mod2_dense a, b) except -1:
    """
    Try to set the entries of a matrix from a numpy array.

    INPUT:

    - ``a`` -- the destination matrix
    - ``b`` -- a numpy array, must have dimension 2 and the same shape as ``a``

    OUTPUT: ``True`` (when used as bool) if successful, ``False`` otherwise. May throw ``ValueError``.

    The exact type of the return value is not guaranteed, in the actual current implementation
    it is ``1`` for success and ``0`` for failure.

    TESTS::

        sage: from sage.modules.numpy_util import set_matrix_mod2_from_numpy
        sage: import numpy as np
        sage: a = matrix(GF(2), 2, 3)
        sage: b = np.array([[1, 0, 1], [0, 1, 0]], dtype=np.int8)
        sage: set_matrix_mod2_from_numpy(a, b)
        1
        sage: a
        [1 0 1]
        [0 1 0]
        sage: set_matrix_mod2_from_numpy(a, np.array([[1, 0], [0, 1]], dtype=np.int8))
        Traceback (most recent call last):
        ...
        ValueError: shape mismatch
        sage: # unsupported type (may be supported in the future)
        sage: set_matrix_mod2_from_numpy(a, np.array([[1, 1, 0], [0, 1, 0]], dtype=np.float64))
        0
        sage: set_matrix_mod2_from_numpy(a, np.array([1, 0, 0], dtype=np.int8))  # wrong number of dimensions
        0
    """
    try:
        return (<object>_set_matrix_mod2_from_numpy_helper)(a, b)  # https://github.com/cython/cython/issues/6588
    except TypeError:
        return False


cpdef object mzd_matrix_to_numpy(uintptr_t entries_addr, object dtype):
    """
    Convert ``<mzd_t*>entries_addr`` to a numpy array.

    INPUT:

    - ``entries_addr`` -- must be a ``mzd_t*`` casted to ``uintptr_t``; the casting
      is necessary to pass it through Python boundary because of lazy import.
      Do not pass arbitrary integer value here, will crash the program.

    - ``dtype`` -- numpy dtype. If ``None``, the result will have some convenient dtype.

    OUTPUT: a 2-dimensional array.
    """
    if dtype is not None:
        return mzd_matrix_to_numpy(entries_addr, None).astype(dtype)
    cdef mzd_t* entries = <mzd_t*>entries_addr
    cdef np.ndarray[np.uint8_t, ndim=2] result = np.empty((entries.nrows, entries.ncols), dtype=np.uint8)
    cdef Py_ssize_t i, j
    for i in range(entries.nrows):
        for j in range(entries.ncols):
            result[i, j] = mzd_read_bit(entries, i, j)
    return result
