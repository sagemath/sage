from cpython.longintrepr cimport py_long, digit

cdef extern from "pycore_long.h":
    digit* ob_digit(py_long o)
    bint _PyLong_IsZero(py_long o)
    bint _PyLong_IsNegative(py_long o)
    bint _PyLong_IsPositive(py_long o)
    Py_ssize_t _PyLong_DigitCount(py_long o)
    void _PyLong_SetSignAndDigitCount(py_long o, int sign, Py_ssize_t size)
