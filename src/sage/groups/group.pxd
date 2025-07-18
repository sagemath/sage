from sage.structure.parent cimport Parent

cdef class Group(Parent):
    pass

cdef class AbelianGroup(Group):
    pass

cdef class FiniteGroup(Group):
    pass

cdef class AlgebraicGroup(Group):
    pass
