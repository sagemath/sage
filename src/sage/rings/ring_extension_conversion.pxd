from sage.rings.ring_extension cimport RingExtension_generic


cpdef backend_parent(R) noexcept
cpdef from_backend_parent(R, RingExtension_generic E) noexcept

cpdef backend_element(x) noexcept
cpdef from_backend_element(x, RingExtension_generic E) noexcept

cdef _backend_morphism(f) noexcept
cpdef backend_morphism(f, forget=*) noexcept
cpdef from_backend_morphism(f, RingExtension_generic E) noexcept

cpdef to_backend(arg) noexcept
cpdef from_backend(arg, E) noexcept


