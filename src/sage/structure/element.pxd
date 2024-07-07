# sage_setup: distribution = sagemath-objects
from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.misc.inherit_comparison cimport InheritComparisonMetaclass


cpdef inline parent(x):
    """
    Return the parent of the element ``x``.

    Usually, this means the mathematical object of which ``x`` is an
    element.

    INPUT:

    - ``x`` -- an element

    OUTPUT:

    - If ``x`` is a Sage :class:`Element`, return ``x.parent()``.

    - Otherwise, return ``type(x)``.

    .. SEEALSO::

        `Parents, Conversion and Coercion <http://doc.sagemath.org/html/en/tutorial/tour_coercion.html>`_
        Section in the Sage Tutorial

    EXAMPLES::

        sage: a = 42
        sage: parent(a)
        Integer Ring
        sage: b = 42/1
        sage: parent(b)
        Rational Field
        sage: c = 42.0
        sage: parent(c)                                                                 # needs sage.rings.real_mpfr
        Real Field with 53 bits of precision

    Some more complicated examples::

        sage: x = Partition([3,2,1,1,1])                                                # needs sage.combinat
        sage: parent(x)                                                                 # needs sage.combinat
        Partitions
        sage: v = vector(RDF, [1,2,3])                                                  # needs sage.modules
        sage: parent(v)                                                                 # needs sage.modules
        Vector space of dimension 3 over Real Double Field

    The following are not considered to be elements, so the type is
    returned::

        sage: d = int(42)  # Python int
        sage: parent(d)
        <... 'int'>
        sage: L = list(range(10))
        sage: parent(L)
        <... 'list'>
    """
    if isinstance(x, Element):
        return (<Element>x)._parent
    return type(x)


cdef inline int classify_elements(left, right) noexcept:
    """
    Given two objects, at least one which is an :class:`Element`,
    classify their type and parent. This is a finer version of
    :func:`have_same_parent`.

    OUTPUT: the sum of the following bits:

    - 0o01: left is an Element
    - 0o02: right is an Element
    - 0o04: both are Element
    - 0o10: left and right have the same type
    - 0o20: left and right have the same parent

    These are the possible outcomes:

    - 0o01: left is an Element, right is not
    - 0o02: right is an Element, left is not
    - 0o07: both are Element, different types, different parents
    - 0o17: both are Element, same type, different parents
    - 0o27: both are Element, different types, same parent
    - 0o37: both are Element, same type, same parent
    """
    if type(left) is type(right):
        # We know at least one of the arguments is an Element. So if
        # their types are *equal* (fast to check) then they are both
        # Elements.
        if (<Element>left)._parent is (<Element>right)._parent:
            return 0o37
        else:
            return 0o17
    if not isinstance(right, Element):
        return 0o01
    if not isinstance(left, Element):
        return 0o02
    if (<Element>left)._parent is (<Element>right)._parent:
        return 0o27
    else:
        return 0o07

# Functions to help understand the result of classify_elements()
cdef inline bint BOTH_ARE_ELEMENT(int cl) noexcept:
    return cl & 0o04
cdef inline bint HAVE_SAME_PARENT(int cl) noexcept:
    return cl & 0o20


cpdef inline bint have_same_parent(left, right) noexcept:
    """
    Return ``True`` if and only if ``left`` and ``right`` have the
    same parent.

    .. WARNING::

        This function assumes that at least one of the arguments is a
        Sage :class:`Element`. When in doubt, use the slower
        ``parent(left) is parent(right)`` instead.

    EXAMPLES::

        sage: from sage.structure.element import have_same_parent
        sage: have_same_parent(1, 3)
        True
        sage: have_same_parent(1, 1/2)
        False
        sage: have_same_parent(gap(1), gap(1/2))                                        # needs sage.libs.gap
        True

    These have different types but the same parent::

        sage: a = RLF(2)
        sage: b = exp(a)
        sage: type(a)
        <... 'sage.rings.real_lazy.LazyWrapper'>
        sage: type(b)
        <... 'sage.rings.real_lazy.LazyNamedUnop'>
        sage: have_same_parent(a, b)
        True
    """
    return HAVE_SAME_PARENT(classify_elements(left, right))


cdef unary_op_exception(op, x)
cdef bin_op_exception(op, x, y)


cdef class Element(SageObject):
    cdef Parent _parent
    cpdef _richcmp_(left, right, int op)
    cpdef int _cmp_(left, right) except -2
    cpdef base_extend(self, R)

    cdef getattr_from_category(self, name)

    cpdef _act_on_(self, x, bint self_on_left)
    cpdef _acted_upon_(self, x, bint self_on_left)

    cdef _add_(self, other)
    cdef _sub_(self, other)
    cdef _neg_(self)
    cdef _add_long(self, long n)

    cdef _mul_(self, other)
    cdef _mul_long(self, long n)
    cdef _matmul_(self, other)
    cdef _div_(self, other)
    cdef _floordiv_(self, other)
    cdef _mod_(self, other)

    cdef _pow_(self, other)
    cdef _pow_int(self, n)
    cdef _pow_long(self, long n)


cdef class ElementWithCachedMethod(Element):
    cdef public dict _cached_methods

cdef class ModuleElement(Element)       # forward declaration

cdef class RingElement(ModuleElement)   # forward declaration

cdef class ModuleElement(Element):
    cpdef _add_(self, other)
    cpdef _sub_(self, other)
    cpdef _neg_(self)

    # self._rmul_(x) is x * self
    cpdef _lmul_(self, Element right)
    # self._lmul_(x) is self * x
    cpdef _rmul_(self, Element left)

cdef class ModuleElementWithMutability(ModuleElement):
    cdef bint _is_immutable
    cpdef bint is_immutable(self) noexcept
    cpdef bint is_mutable(self) noexcept

cdef class MonoidElement(Element):
    cpdef _pow_int(self, n)

cdef class MultiplicativeGroupElement(MonoidElement):
    cpdef _div_(self, other)

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(ModuleElement):
    cpdef _mul_(self, other)
    cpdef _div_(self, other)
    cpdef _pow_int(self, n)

cdef class CommutativeRingElement(RingElement):
    pass

cdef class IntegralDomainElement(CommutativeRingElement):
    pass

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    pass

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):
    cpdef _floordiv_(self, other)
    cpdef _mod_(self, other)

cdef class FieldElement(CommutativeRingElement):
    cpdef _floordiv_(self, other)

cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebraElement(CommutativeRingElement):
    pass

cdef class Expression(CommutativeRingElement):
    pass

cdef class InfinityElement(RingElement):
    pass

cdef class Vector(ModuleElementWithMutability):
    cdef Py_ssize_t _degree

    # Return the dot product using the simple metric
    # $e_i \cdot e_j = \delta_{ij}$. The first assumes that the parents
    # are equal, both assume that the degrees are equal.
    cpdef _dot_product_(Vector left, Vector right)
    cpdef _dot_product_coerce_(Vector left, Vector right)

    cpdef _pairwise_product_(Vector left, Vector right) # override, call if parents the same

    cdef bint is_sparse_c(self) noexcept
    cdef bint is_dense_c(self) noexcept


cdef class Matrix(ModuleElement):
    # All matrix classes must be written in Cython
    cdef Py_ssize_t _nrows
    cdef Py_ssize_t _ncols

    cdef _vector_times_matrix_(matrix_right, Vector vector_left)    # OK to override, AND call directly
    cdef _matrix_times_vector_(matrix_left, Vector vector_right)    # OK to override, AND call directly
    cdef _matrix_times_matrix_(left, Matrix right)                  # OK to override, AND call directly

    cdef bint is_sparse_c(self) noexcept
    cdef bint is_dense_c(self) noexcept
