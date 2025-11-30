"""
Declaration file for degree_sequences.pyx

This file declares the cdef class and its attributes/methods so they can be
accessed from other Cython modules.
"""

cdef class _DegreeSequenceEnumerator:
    cdef int N
    cdef unsigned char * seq

    cdef tuple build_current_seq(self)
