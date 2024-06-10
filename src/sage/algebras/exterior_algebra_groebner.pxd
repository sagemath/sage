"""
Exterior algebras Gröbner bases
"""
from sage.data_structures.bitset cimport FrozenBitset
from sage.rings.integer cimport Integer
from sage.algebras.clifford_algebra_element cimport CliffordAlgebraElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport MonoidElement

cdef long degree(FrozenBitset X) noexcept
cdef CliffordAlgebraElement build_monomial(Parent E, FrozenBitset supp)

cdef class GBElement:
    cdef CliffordAlgebraElement elt
    cdef FrozenBitset ls  # leading support as a bitset
    cdef Integer lsi  # leading support as an Integer

# Grobner basis functions
cdef class GroebnerStrategy:
    cdef Parent E  # the exterior algebra
    cdef int side
    cdef MonoidElement ideal
    cdef bint homogeneous
    cdef Integer rank
    cdef public tuple groebner_basis

    cdef inline GBElement build_elt(self, CliffordAlgebraElement f)
    cdef inline GBElement prod_GB_term(self, GBElement f, FrozenBitset t)
    cdef inline GBElement prod_term_GB(self, FrozenBitset t, GBElement f)
    cdef inline bint build_S_poly(self, GBElement f, GBElement g) noexcept

    cdef inline FrozenBitset leading_support(self, CliffordAlgebraElement f)
    cdef inline partial_S_poly_left(self, GBElement f, GBElement g)
    cdef inline partial_S_poly_right(self, GBElement f, GBElement g)
    cdef set preprocessing(self, list P, list G)
    cdef list reduction(self, list P, list G)

    cpdef CliffordAlgebraElement reduce(self, CliffordAlgebraElement f)
    cdef bint reduce_single(self, CliffordAlgebraElement f, CliffordAlgebraElement g) except -1
    cdef int reduced_gb(self, list G) except -1

    # These are the methods that determine the ordering of the monomials.
    # These must be implemented in subclasses. Declare them as "inline" there.
    cdef Integer bitset_to_int(self, FrozenBitset X)
    cdef FrozenBitset int_to_bitset(self, Integer n)

cdef class GroebnerStrategyNegLex(GroebnerStrategy):
    pass

cdef class GroebnerStrategyDegRevLex(GroebnerStrategy):
    pass

cdef class GroebnerStrategyDegLex(GroebnerStrategy):
    pass
