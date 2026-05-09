from libcpp.memory cimport unique_ptr, shared_ptr, make_shared

from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base, BooleanPolynomialRing_base
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.structure.element cimport MonoidElement

from sage.libs.polybori.decl cimport *


cdef class BooleanPolynomialRing(BooleanPolynomialRing_base):
    cdef PBRing _pbring
    cdef Py_ssize_t* pbind
    cdef public _monom_monoid
    cdef public object __interface
    cdef object _repr

    # it is very important to keep this cached, since otherwise the magma interface will break
    cdef public object __cover_ring

    cdef _convert(self, rhs)

cdef class BooleanPolynomial(MPolynomial):
    cdef PBPoly _pbpoly
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

cdef class BooleSet:
    cdef BooleanPolynomialRing _ring
    cdef PBSet _pbset

cdef class CCuddNavigator:
    cdef PBNavigator _pbnav
    cdef Py_ssize_t* _pbind

cdef class BooleanMonomial(MonoidElement):
    cdef PBMonom _pbmonom
    cdef BooleanPolynomialRing _ring
    cpdef _mul_(self, other)

cdef class BooleanMonomialVariableIterator:
    cdef object parent
    cdef BooleanPolynomialRing _ring
    cdef BooleanMonomial obj
    cdef PBMonomIter _iter
    cdef PBMonomIter _end
    cdef Py_ssize_t* pbind

cdef class BooleanMonomialIterator:
    cdef BooleanMonomial obj
    cdef PBMonomIter _iter
    cdef PBMonomIter _end
    cdef Py_ssize_t* pbind

# Wrap PBPolyIter using pointers because there is no default constructor
cdef class BooleanPolynomialIterator:
    cdef BooleanPolynomial obj
    cdef PBPolyIter* _iter
    cdef PBPolyIter* _end

# Wrap PBSetIter using pointers because there is no default constructor
cdef class BooleSetIterator:
    cdef object _parent
    cdef BooleanPolynomialRing _ring
    cdef PBSetIter* _iter
    cdef PBSetIter* _end
    cdef BooleSet obj

cdef class BooleanPolynomialEntry:
    cdef public BooleanPolynomial p

# Wrap PBRedStrategy using shared_ptr because there is no default
# constructor and because multiple Python objects may point to
# the same PBRedStrategy
cdef class ReductionStrategy:
    cdef shared_ptr[PBRedStrategy] _strat
    cdef BooleanPolynomialRing _parent

# Wrap PBGBStrategy using shared_ptr because there is no default
# constructor and because multiple Python objects may point to
# the same PBGBStrategy
cdef class GroebnerStrategy:
    cdef shared_ptr[PBGBStrategy] _strat
    cdef BooleanPolynomialRing _parent
    cdef public ReductionStrategy reduction_strategy

# Wrap PBFGLMStrategy using unique_ptr for analogy with
# ReductionStrategy and GroebnerStrategy
cdef class FGLMStrategy:
    cdef unique_ptr[PBFGLMStrategy] _strat
    cdef BooleanPolynomialRing _parent

cdef class BooleanPolynomialVector:
    cdef PBPolyVector _vec
    cdef BooleanPolynomialRing _parent

cdef class BooleanPolynomialVectorIterator:
    cdef BooleanPolynomialVector obj
    cdef BooleanPolynomialRing _parent
    cdef PBPolyVectorIter _iter
    cdef PBPolyVectorIter _end

cdef class VariableBlock:
    cdef BooleanPolynomialRing _ring
    cdef PBVarBlock* _block
    cdef public object __name__

cdef class BooleConstant:
    cdef PBConstant _pbconst

cdef class VariableFactory:
    cdef BooleanPolynomialRing _ring
    cdef PBVarFactory _factory

cdef class MonomialFactory:
    cdef BooleanPolynomialRing _ring
    cdef PBMonomFactory _factory

cdef class PolynomialFactory:
    cdef BooleanPolynomialRing _ring
    cdef PBPolyFactory _factory


# Cython doesn't seem to support constructors with additional template
# parameters, so we declare this aliasing constructor as special case
cdef extern from *:
    cdef shared_ptr[T] shared_ptr_alias_PBGBStrategy "std::shared_ptr"[T](shared_ptr[PBGBStrategy]&, T*)
