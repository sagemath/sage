"""
Cython data structure for combinatorial faces.
"""
# ****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from memory_allocator                 cimport MemoryAllocator
from sage.data_structures.bitset_base cimport *

ctypedef int simple
ctypedef long standard

cdef struct face_s:
    sparse_bitset_t atoms
    bitset_t coatoms

ctypedef face_s face_t[1]

ctypedef fused algorithm_variant:
    simple
    standard

#############################################################################
# Face Initialization
#############################################################################

cdef inline bint face_init(face_t face, mp_bitcnt_t n_atoms, mp_bitcnt_t n_coatoms) except -1:
    """
    Initialize and clear ``face``.
    """
    if n_coatoms == 0:
        # Special case for trivial polyhedra.
        n_coatoms += 1
    if n_atoms == 0:
        n_atoms += 1
    bitset_init(face.atoms, n_atoms)
    bitset_init(face.coatoms, n_coatoms)

cdef inline void face_free(face_t face) noexcept:
    """
    Free ``face``.
    """
    bitset_free(face.atoms)
    bitset_free(face.coatoms)

cdef inline bint face_check_alignment(face_t face) noexcept:
    """
    Return whether the data is correctly aligned.
    """
    return bitset_check_alignment(face.atoms)

cdef inline void face_clear(face_t face) noexcept:
    """
    Remove all atoms and coatoms from face.
    """
    bitset_clear(face.atoms)
    bitset_clear(face.coatoms)

cdef inline void face_copy(face_t dst, face_t src) noexcept:
    """
    Copy src to dst overwriting dst.

    dst may contain more atoms and coatoms, but not less.
    """
    bitset_copy_flex(dst.atoms, src.atoms)
    bitset_copy_flex(dst.coatoms, src.coatoms)


#############################################################################
# Face Comparison
#############################################################################

cdef inline bint face_isempty(face_t face) noexcept nogil:
    """
    Return whether ``face`` contains no coatoms.
    """
    return bitset_isempty(face.atoms)

cdef inline int face_cmp(face_t a, face_t b) noexcept:
    """
    Return ``0`` if the faces are equal and consistently
    ``-1`` and ``1`` if not.
    """
    return bitset_cmp(a.atoms, b.atoms)

cdef inline bint face_issubset_fused(face_t a, face_t b, algorithm_variant algorithm) noexcept nogil:
    """
    Return whether ``a`` is a subface of ``b``.
    """
    if algorithm_variant is standard:
        return bitset_issubset(a.atoms, b.atoms)
    else:
        return bitset_issuperset(a.coatoms, b.coatoms)

cdef inline bint face_issubset(face_t a, face_t b) noexcept nogil:
    return face_issubset_fused(a, b, <standard> 0)

#############################################################################
# Face Bit Manipulation
#############################################################################

cdef inline bint face_atom_in(face_t face, mp_bitcnt_t n) noexcept:
    """
    Return whether ``n`` is an atom of ``face``.
    """
    return bitset_in(face.atoms, n)

cdef inline void face_add_atom(face_t face, mp_bitcnt_t n) noexcept:
    """
    Add atom `n` to the face.
    """
    bitset_add(face.atoms, n)

cdef inline int face_add_atom_safe(face_t face, mp_bitcnt_t n) except -1:
    """
    Add atom `n` to the face.
    """
    if (n > face.atoms.size):
        raise KeyError(n)
    bitset_add(face.atoms, n)

cdef inline void face_discard_atom(face_t face, mp_bitcnt_t n) noexcept:
    """
    Discard atom `n` of the face.
    """
    bitset_discard(face.atoms, n)

cdef inline void facet_set_coatom(face_t face, mp_bitcnt_t n) noexcept:
    """
    Set the facet to be coatom ``n``.
    """
    bitset_clear(face.coatoms)
    bitset_add(face.coatoms, n)

cdef inline void face_set_first_n_atoms(face_t face, mp_bitcnt_t n) noexcept:
    """
    Set exactly the first ``n`` atoms.
    """
    bitset_set_first_n(face.atoms, n)


#############################################################################
# Face Searching
#############################################################################

cdef inline long face_next_atom(face_t face, mp_bitcnt_t n) noexcept:
    """
    Return the index of the next atom in ``face`` with index >= ``n``.

    In case there are none, return ``-1``.
    """
    return bitset_next(face.atoms, n)

cdef inline long face_first_missing_atom(face_t face) noexcept:
    """
    Return the index of the first atom not in ``face``.

    In case there are none, return ``-1``.
    """
    return bitset_first_in_complement(face.atoms)

cdef inline long face_len_atoms(face_t face) noexcept nogil:
    """
    Calculate the number of atoms in the face.
    """
    return bitset_len(face.atoms)


#############################################################################
# Arithmetic
#############################################################################

cdef inline void face_intersection_fused(face_t dest, face_t A, face_t B, algorithm_variant algorithm) noexcept nogil:
    """
    Set ``dest`` to the intersection of ``A`` and ``B``.
    """
    if algorithm_variant is standard:
        # Also setting the nonzero positions.
        sparse_bitset_intersection(dest.atoms, A.atoms, B.atoms)
    else:
        bitset_intersection(dest.atoms, A.atoms, B.atoms)
        bitset_union(dest.coatoms, A.coatoms, B.coatoms)

cdef inline void face_intersection(face_t dest, face_t A, face_t B) noexcept nogil:
    face_intersection_fused(dest, A, B, <standard> 0)


#############################################################################
# Miscellaneous
#############################################################################

cdef inline void swap_faces(face_t a, face_t b) noexcept nogil:
    cdef face_t tmp
    tmp[0] = a[0]
    a[0] = b[0]
    b[0] = tmp[0]

cdef inline bint faces_are_identical(face_t a, face_t b) noexcept nogil:
    return a.atoms.limbs == b.atoms.limbs
