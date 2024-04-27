from sage.rings.ring_extension cimport RingExtension_generic


cpdef backend_parent(R)
cpdef from_backend_parent(R, RingExtension_generic E)

cpdef backend_element(x)
cpdef from_backend_element(x, RingExtension_generic E)

cdef _backend_morphism(f)
cpdef backend_morphism(f, forget=*)
cpdef from_backend_morphism(f, RingExtension_generic E)

cpdef to_backend(arg)
cpdef from_backend(arg, E)


