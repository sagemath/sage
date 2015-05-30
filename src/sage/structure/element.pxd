from sage_object cimport SageObject
from parent cimport Parent
from cpython.number cimport PyNumber_Check

cdef inline parent_c(x):
    if isinstance(x, Element):
        return (<Element>x)._parent
    # Fast check for "number" types, including int and float
    if PyNumber_Check(x):
        return type(x)
    try:
        p = x.parent
    except AttributeError:
        return type(x)
    else:
        return p()

cdef inline bint have_same_parent_c(left, right):
    """
    Return ``True`` if and only if ``left`` and ``right`` have the
    same parent.
    """
    # We know at least one of the arguments is an Element. So if
    # their types are *equal* (fast to check) then they are both
    # Elements.  Otherwise use the slower test via isinstance().
    if type(left) is type(right):
        return (<Element>left)._parent is (<Element>right)._parent
    if isinstance(right, Element) and isinstance(left, Element):
        return (<Element>left)._parent is (<Element>right)._parent
    return False


cdef str arith_error_message(x, y, op)

cdef class Element(SageObject):
    cdef Parent _parent
    cpdef _richcmp_(left, Element right, int op)
    cpdef int _cmp_(left, Element right) except -2
    cdef _richcmp(self, right, int op)
    cdef int _cmp(self, right) except -2
    cdef _set_parent_c(self, Parent parent)
    cpdef base_extend(self, R)

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef _acted_upon_(self, x, bint self_on_left)

cdef class ElementWithCachedMethod(Element):
    cdef public dict __cached_methods

cdef class ModuleElement(Element)       # forward declaration

cdef class RingElement(ModuleElement)   # forward declaration

cdef class ModuleElement(Element):

    cpdef ModuleElement _add_(self, ModuleElement right)
    cpdef ModuleElement _sub_(self, ModuleElement right)
    cpdef ModuleElement _neg_(self)
    # self._rmul_(x) is x * self
    cpdef ModuleElement _lmul_(self, RingElement right)
    # self._lmul_(x) is self * x, to abide with Python conventions.
    cpdef ModuleElement _rmul_(self, RingElement left)

    cdef ModuleElement _mul_long(self, long n)

    # Inplace operations, override, do *NOT* call directly
    cpdef ModuleElement _iadd_(self, ModuleElement right)
    cpdef ModuleElement _isub_(self, ModuleElement right)
    cpdef ModuleElement _ilmul_(self, RingElement right)

    # Coerce x to the base ring of self and return the result.
    cdef RingElement coerce_to_base_ring(self, x)

cdef class MonoidElement(Element):
    cpdef MonoidElement _mul_(self, MonoidElement right)

cdef class MultiplicativeGroupElement(MonoidElement):
    cpdef MultiplicativeGroupElement _div_(self, MultiplicativeGroupElement right)


cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(ModuleElement):
    cpdef RingElement _mul_(self, RingElement right)
    cpdef RingElement _div_(self, RingElement right)

    # Inplace operations, override, do *NOT* call directly
    cpdef RingElement _imul_(self, RingElement right)
    cpdef RingElement _idiv_(self, RingElement right)

    cdef RingElement _add_long(self, long n)

cdef class CommutativeRingElement(RingElement):
    pass

cdef class IntegralDomainElement(CommutativeRingElement):
    pass

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    pass

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):
    pass

cdef class FieldElement(CommutativeRingElement):
    pass

cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

cdef class CommutativeAlgebra(AlgebraElement):
    pass

cdef class InfinityElement(RingElement):
    pass


cdef class Vector(ModuleElement):
    cdef Py_ssize_t _degree

    # Return the dot product using the simple metric
    # $e_i \cdot e_j = \delta_{ij}$. The first assumes that the parents
    # are equal, both assume that the degrees are equal.
    cpdef Element _dot_product_(Vector left, Vector right)
    cpdef Element _dot_product_coerce_(Vector left, Vector right)

    cpdef Vector _pairwise_product_(Vector left, Vector right) # override, call if parents the same

    cdef bint is_sparse_c(self)
    cdef bint is_dense_c(self)


cdef class Matrix(ModuleElement):
    # All matrix classes must be written in Cython
    cdef Py_ssize_t _nrows
    cdef Py_ssize_t _ncols

    cdef Vector _vector_times_matrix_(matrix_right, Vector vector_left)    # OK to override, AND call directly
    cdef Vector _matrix_times_vector_(matrix_left, Vector vector_right)    # OK to override, AND call directly
    cdef Matrix _matrix_times_matrix_(left, Matrix right)                  # OK to override, AND call directly

    cdef bint is_sparse_c(self)
    cdef bint is_dense_c(self)




cdef class CoercionModel:
    cpdef canonical_coercion(self, x, y)
    cpdef bin_op(self, x, y, op)


cdef generic_power_c(a, nn, one)
