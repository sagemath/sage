
from sage.structure.parent_base cimport ParentWithBase
from sage.structure.parent_gens cimport ParentWithGens
from sage.rings.polynomial.multi_polynomial_ring_generic cimport \
                                                MPolynomialRing_generic
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.structure.element cimport MonoidElement

from sage.libs.polybori.decl cimport *


cdef class BooleanPolynomialRing(MPolynomialRing_generic):
    cdef PBRing _pbring
    cdef Py_ssize_t* pbind
    cdef public _monom_monoid
    cdef public object __interface
    cdef object _repr

    # it is very important to keep this cached, since otherwise the magma interface will break
    cdef public object __cover_ring

cdef class BooleanPolynomial(MPolynomial):
    cdef PBPoly _pbpoly

cdef class BooleSet:
    cdef BooleanPolynomialRing _ring
    cdef PBSet _pbset

cdef class CCuddNavigator:
    cdef PBNavigator _pbnav
    cdef Py_ssize_t* _pbind

cdef class BooleanMonomial(MonoidElement):
    cdef PBMonom _pbmonom
    cdef BooleanPolynomialRing _ring

cdef class BooleanMonomialVariableIterator:
    cdef object parent
    cdef BooleanPolynomialRing _ring
    cdef BooleanMonomial obj
    cdef PBMonomVarIter _iter
    cdef PBMonomVarIter _end

cdef class BooleanMonomialIterator:
    cdef BooleanMonomial obj
    cdef PBMonomIter _iter
    cdef PBMonomIter _end
    cdef Py_ssize_t* pbind

cdef class BooleanPolynomialIterator:
    cdef BooleanPolynomial obj
    cdef PBPolyIter _iter
    cdef PBPolyIter _end

cdef class BooleSetIterator:
    cdef object _parent
    cdef BooleanPolynomialRing _ring
    cdef PBSetIter _iter
    cdef PBSetIter _end
    cdef BooleSet obj

cdef class BooleanPolynomialEntry:
    cdef public BooleanPolynomial p

cdef class ReductionStrategy:
    cdef PBRedStrategy *_strat
    cdef bint _borrowed
    cdef BooleanPolynomialRing _parent

cdef class GroebnerStrategy:
    cdef PBGBStrategy* _strat
    cdef PBRefCounter _count
    cdef BooleanPolynomialRing _parent
    cdef public ReductionStrategy reduction_strategy

cdef class FGLMStrategy:
    cdef PBFglmStrategy _strat
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
