/* sage_setup: distribution = sagemath-objects
 */
#include "Python.h"
#include <stdbool.h>

#if PY_VERSION_HEX >= 0x030C00A5
#define ob_digit(o)  (((PyLongObject*)o)->long_value.ob_digit)
#else
#define ob_digit(o)  (((PyLongObject*)o)->ob_digit)
#endif

#if PY_VERSION_HEX >= 0x030C00A7
// taken from cpython:Include/internal/pycore_long.h @ 3.12

/* Long value tag bits:
 * 0-1: Sign bits value = (1-sign), ie. negative=2, positive=0, zero=1.
 * 2: Reserved for immortality bit
 * 3+ Unsigned digit count
 */
#define SIGN_MASK 3
#define SIGN_ZERO 1
#define SIGN_NEGATIVE 2
#define NON_SIZE_BITS 3

static inline bool
_PyLong_IsZero(const PyLongObject *op)
{
    return (op->long_value.lv_tag & SIGN_MASK) == SIGN_ZERO;
}

static inline bool
_PyLong_IsNegative(const PyLongObject *op)
{
    return (op->long_value.lv_tag & SIGN_MASK) == SIGN_NEGATIVE;
}

static inline bool
_PyLong_IsPositive(const PyLongObject *op)
{
    return (op->long_value.lv_tag & SIGN_MASK) == 0;
}

static inline Py_ssize_t
_PyLong_DigitCount(const PyLongObject *op)
{
    assert(PyLong_Check(op));
    return op->long_value.lv_tag >> NON_SIZE_BITS;
}

#define TAG_FROM_SIGN_AND_SIZE(sign, size) ((1 - (sign)) | ((size) << NON_SIZE_BITS))

static inline void
_PyLong_SetSignAndDigitCount(PyLongObject *op, int sign, Py_ssize_t size)
{
    assert(size >= 0);
    assert(-1 <= sign && sign <= 1);
    assert(sign != 0 || size == 0);
    op->long_value.lv_tag = TAG_FROM_SIGN_AND_SIZE(sign, (size_t)size);
}

#else
// fallback for < 3.12

static inline bool
_PyLong_IsZero(const PyLongObject *op)
{
    return Py_SIZE(op) == 0;
}

static inline bool
_PyLong_IsNegative(const PyLongObject *op)
{
    return Py_SIZE(op) < 0;
}

static inline bool
_PyLong_IsPositive(const PyLongObject *op)
{
    return Py_SIZE(op) > 0;
}

static inline Py_ssize_t
_PyLong_DigitCount(const PyLongObject *op)
{
    Py_ssize_t size = Py_SIZE(op);
    return size < 0 ? -size : size;
}

static inline void
_PyLong_SetSignAndDigitCount(PyLongObject *op, int sign, Py_ssize_t size)
{
#if (PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION < 9)
// The function Py_SET_SIZE is defined starting with python 3.9.
    Py_SIZE(op) = size;
#else
    Py_SET_SIZE(op, sign < 0 ? -size : size);
#endif
}

#endif
