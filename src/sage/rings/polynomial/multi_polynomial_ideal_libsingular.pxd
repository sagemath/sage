from sage.libs.singular.decl cimport ideal, ring

cdef object singular_ideal_to_sage_sequence(ideal *i, ring *r, object parent) noexcept
cdef ideal *sage_ideal_to_singular_ideal(I) except NULL
