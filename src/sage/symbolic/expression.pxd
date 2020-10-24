from sage.libs.pynac.pynac cimport GEx
from sage.structure.element cimport CommutativeRingElement

cdef class Expression(CommutativeRingElement):
    cdef GEx _gobj
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _pow_(self, other)
    cdef Expression coerce_in(self, z)
    cpdef object _eval_self(self, R)
    cpdef object _convert(self, R)
    cpdef bint is_polynomial(self, var)
    cpdef bint is_relational(self)
    cpdef bint is_infinity(self)
    cpdef bint is_positive_infinity(self)
    cpdef bint is_negative_infinity(self)
    cpdef object pyobject(self)
    cpdef Expression _subs_expr(self, expr)
    cpdef int _cmp_add(Expression left, Expression right) except -2
    cpdef int _cmp_mul(Expression left, Expression right) except -2
    cdef bint _rel_equal1(Expression self, Expression other) except -1
    cdef bint _rel_equal2(Expression self, Expression other) except -1

cpdef bint is_Expression(x)
cdef Expression new_Expression_from_GEx(parent, GEx juice)
cdef Expression new_Expression_from_pyobject(parent, x)
