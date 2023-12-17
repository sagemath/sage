cpdef tuple _sort_uniq(categories) noexcept
cdef class AxiomContainer(dict):
    pass
cpdef tuple canonicalize_axioms(AxiomContainer all_axioms, axioms) noexcept
from sage.misc.classcall_metaclass cimport ClasscallMetaclass
cpdef tuple _flatten_categories(categories, ClasscallMetaclass JoinCategory) noexcept
cpdef tuple join_as_tuple(tuple categories, tuple axioms, tuple ignore_axioms) noexcept
