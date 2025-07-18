# sage_setup: distribution = sagemath-objects
cpdef tuple _sort_uniq(categories)
cdef class AxiomContainer(dict):
    pass
cpdef tuple canonicalize_axioms(AxiomContainer all_axioms, axioms)
from sage.misc.classcall_metaclass cimport ClasscallMetaclass
cpdef tuple _flatten_categories(categories, ClasscallMetaclass JoinCategory)
cpdef tuple join_as_tuple(tuple categories, tuple axioms, tuple ignore_axioms)
