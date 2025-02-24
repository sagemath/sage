from sage.misc.binary_tree import BinaryTree
from sage.misc.binary_tree cimport BinaryTree


cdef class generic_pd

cdef class CompiledPolynomialFunction:
    cdef generic_pd _dag
    cdef object _coeffs

    cdef object _parse_structure(CompiledPolynomialFunction)
    cdef generic_pd _get_gap(CompiledPolynomialFunction, BinaryTree, int)
    cdef void _fill_gaps_binary(CompiledPolynomialFunction, BinaryTree) noexcept
    cdef object eval(CompiledPolynomialFunction, object)

cdef class generic_pd:
    cdef object value
    cdef int refs, hits
    cdef int label
    cdef int eval(self, vars, coeffs) except -2
    cdef generic_pd nodummies(generic_pd)
    cdef void reset(self) noexcept

cdef class dummy_pd(generic_pd):
    cdef generic_pd link
    cdef void fill(dummy_pd self, generic_pd link) noexcept

cdef class var_pd(generic_pd):
    cdef int index

cdef class univar_pd(generic_pd):
    pass

cdef class coeff_pd(generic_pd):
    cdef int index


cdef class unary_pd(generic_pd):
    cdef generic_pd operand

cdef class sqr_pd(unary_pd):
    pass

cdef class pow_pd(unary_pd):
    cdef object exponent


cdef class binary_pd(generic_pd):
    cdef generic_pd left, right

cdef class add_pd(binary_pd):
    pass

cdef class mul_pd(binary_pd):
    pass

cdef class abc_pd(binary_pd):
    cdef int index
