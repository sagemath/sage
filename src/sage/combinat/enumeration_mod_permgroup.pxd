# sage_setup: distribution = sagemath-gap
from sage.structure.list_clone cimport ClonableIntArray

cpdef list all_children(ClonableIntArray v, int max_part)
cpdef int lex_cmp_partial(ClonableIntArray t1, ClonableIntArray t2, int step) noexcept
cpdef int lex_cmp(ClonableIntArray t1, ClonableIntArray t2) noexcept
cpdef bint is_canonical(list sgs, ClonableIntArray v) except -1
cpdef ClonableIntArray canonical_representative_of_orbit_of(list sgs, ClonableIntArray v)
cpdef list canonical_children(list sgs, ClonableIntArray v, int max_part)
cpdef set orbit(list sgs, ClonableIntArray v)
