from cpython.object cimport PyObject

cdef union double_op_params:
    PyObject* func
    double (*f)(double)
    double (*ff)(double, double)
    double c
    int n

cdef struct fast_double_op:
    char type
    double_op_params params

cdef class FastDoubleFunc:
    cdef readonly int max_height
    cdef readonly int nargs
    cdef readonly int nops
    cdef fast_double_op* ops

    cdef double* argv
    cdef double* stack

    # need to keep this around because structs can't contain (ref-counted) python objects
    cdef py_funcs

    cdef int allocate_stack(FastDoubleFunc self) except -1
    cdef double _call_c(FastDoubleFunc self, double* argv) except? -2
    cpdef bint is_pure_c(self)
    cdef FastDoubleFunc cfunc(FastDoubleFunc self, void* func)
    cdef FastDoubleFunc unop(FastDoubleFunc self, char type)
