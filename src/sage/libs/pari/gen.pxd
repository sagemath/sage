include 'decl.pxi'

cimport sage.structure.element
cimport cython

@cython.final
cdef class gen(sage.structure.element.RingElement):
    cdef GEN g
    cdef object _refers_to
    cdef pari_sp b
    cdef void init(self, GEN g, pari_sp b)

cdef gen objtogen(object s)
