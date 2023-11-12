from sage.rings.tate_algebra_element cimport TateAlgebraTerm
from sage.rings.tate_algebra_element cimport TateAlgebraElement

cdef Jpair(p1, p2) noexcept
cdef TateAlgebraElement regular_reduce(sgb, TateAlgebraTerm s, TateAlgebraElement v, stopval) noexcept
cdef TateAlgebraElement reduce(gb, TateAlgebraElement v, stopval) noexcept
