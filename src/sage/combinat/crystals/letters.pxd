# sage_setup: distribution = sagemath-combinat
from sage.structure.element cimport Element

cdef class Letter(Element):
    cdef readonly int value

cdef class EmptyLetter(Element):
    cdef readonly str value
    cpdef e(self, int i)
    cpdef f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_A_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_B_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_C_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_D_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_G_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class LetterTuple(Element):
    cdef readonly tuple value
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class Crystal_of_letters_type_E6_element(LetterTuple):
    cpdef LetterTuple e(self, int i)
    cpdef LetterTuple f(self, int i)

cdef class Crystal_of_letters_type_E6_element_dual(LetterTuple):
    cpdef LetterTuple lift(self)
    cpdef LetterTuple retract(self, LetterTuple p)
    cpdef LetterTuple e(self, int i)
    cpdef LetterTuple f(self, int i)

cdef class Crystal_of_letters_type_E7_element(LetterTuple):
    cpdef LetterTuple e(self, int i)
    cpdef LetterTuple f(self, int i)

cdef class BKKLetter(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)

cdef class QueerLetter_element(Letter):
    cpdef Letter e(self, int i)
    cpdef Letter f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept

cdef class LetterWrapped(Element):
    cdef readonly Element value
    cpdef tuple _to_tuple(self)
    cpdef LetterWrapped e(self, int i)
    cpdef LetterWrapped f(self, int i)
    cpdef int epsilon(self, int i) noexcept
    cpdef int phi(self, int i) noexcept
