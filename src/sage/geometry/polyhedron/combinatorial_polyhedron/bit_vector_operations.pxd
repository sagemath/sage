# distutils: language = c++

cimport cython
from libc.stdint                cimport uint64_t

cdef extern from "bit_vector_operations.cc":
    # Any Bit-representation is assumed to be `chunksize`-Bit aligned.
    cdef const size_t chunksize
    cdef void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
                           size_t face_length)
#    Return ``A & ~B == 0``.
#    A is not subset of B, iff there is a vertex in A, which is not in B.
#    ``face_length`` is the length of A and B in terms of uint64_t.

    cdef size_t get_next_level(
        uint64_t **faces, const size_t n_faces, uint64_t **nextfaces,
        uint64_t **nextfaces2, uint64_t **visited_all,
        size_t n_visited_all, size_t face_length)
#        Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
#        that are not contained in a face of ``visited_all``.

#        INPUT:

#        - ``maybe_newfaces`` -- quasi of type ``uint64_t[n_faces -1][face_length]``,
#          needs to be ``chunksize``-Bit aligned
#        - ``newfaces`` -- quasi of type ``*uint64_t[n_faces -1]
#        - ``visited_all`` -- quasi of type ``*uint64_t[n_visited_all]
#        - ``face_length`` -- length of the faces

#        OUTPUT:

#        - return number of ``newfaces``
#        - set ``newfaces`` to point to the new faces

#        ALGORITHM:

#        To get all facets of ``faces[n_faces-1]``, we would have to:
#        - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
#        - Add all the intersection of ``visited_all`` with the last face
#        - Out of both the inclusion-maximal ones are of codimension 1, i.e. facets.

#        As we have visited all faces of ``visited_all``, we alter the algorithm
#        to not revisit:
#        Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
#        Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
#                At least we obtain all of those, that we have not already visited.
#                Maybe, we get some more.
#        Step 3: Only keep those that we have not already visited.
#                We obtain exactly the facets of ``faces[n_faces-1]`` that we have
#                not visited yet.

    cdef size_t count_atoms(uint64_t *A, size_t face_length)
#        Return the number of atoms/vertices in A.
#        This is the number of set bits in A.
#        ``face_length`` is the length of A in terms of uint64_t.

    cdef size_t bit_repr_to_coatom_repr(
            uint64_t *face, uint64_t **coatoms, size_t n_coatoms,
            size_t face_length, size_t *output)
#        Write the coatom-representation of face in output. Return length.
#        ``face_length`` is the length of ``face`` and ``coatoms[i]``
#        in terms of uint64_t.
#        ``n_coatoms`` length of ``coatoms``.
