/*
#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
*/

#include <cstdint>
#include <cstdio>
using namespace std;

/*
The following three functions can be optimized with intrinsics easily.
Enabling the intrinsics is subject of #27103.
The chunksize is the number of vertices/atoms
that can be handled with one operation.
*/

// Any Bit-representation is assumed to be `chunksize`-Bit aligned.
const size_t chunksize = 64;

inline int is_subset(uint64_t *A, uint64_t *B, size_t face_length){
    /*
    Return ``A & ~B == 0``.
    A is not subset of B, iff there is a vertex in A, which is not in B.
    ``face_length`` is the length of A and B in terms of uint64_t.
    */
    for (size_t i = 0; i < face_length; i++){
        if (A[i] & ~B[i]){
            return 0;
        }
    }
    return 1;
}

inline void intersection(uint64_t *dest, uint64_t *A, uint64_t *B, \
                         size_t face_length){
    /*
    Set ``dest = A & B``, i.e. dest is the intersection of A and B.
    ``face_length`` is the length of A, B and dest in terms of uint64_t.
    */
    for (size_t i = 0; i < face_length; i++){
        dest[i] = A[i] & B[i];
    }
}

inline size_t count_atoms(const uint64_t* A, size_t face_length){
    /*
    Return the number of atoms/vertices in A.
    This is the number of set bits in A.
    ``face_length`` is the length of A in terms of uint64_t.
    */
    size_t i;
    unsigned int count = 0;
    for (i=0; i<face_length; i++){
        uint64_t a = A[i];
        a = a - ((a >> 1) & 0x5555555555555555ULL);
        a = (a & 0x3333333333333333ULL) + ((a >> 2) & 0x3333333333333333ULL);
        count += ( ((a + (a >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
    }
    return count;
}

size_t get_next_level(\
        uint64_t **faces, size_t n_faces, uint64_t **maybe_newfaces, \
        uint64_t **newfaces, uint64_t **visited_all, \
        size_t n_visited_all, size_t face_length){
    /*
    Set ``newfaces`` to be the facets of ``faces[n_faces -1]``
    that are not contained in a face of ``visited_all``.

    INPUT:

    - ``maybe_newfaces`` -- quasi of type ``uint64_t[n_faces -1][face_length]``,
      needs to be ``chunksize``-Bit aligned
    - ``newfaces`` -- quasi of type ``*uint64_t[n_faces -1]
    - ``visited_all`` -- quasi of type ``*uint64_t[n_visited_all]
    - ``face_length`` -- length of the faces

    OUTPUT:

    - return number of ``newfaces``
    - set ``newfaces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces[n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    */

    // We keep track, which face in ``maybe_newfaces`` is a new face.
    int is_not_newface[n_faces -1];

    // Step 1:
    for (size_t j = 0; j < n_faces - 1; j++){
        intersection(maybe_newfaces[j], faces[j], faces[n_faces - 1], face_length);
        is_not_newface[j] = 0;
    }


    // For each face we will Step 2 and Step 3.
    for (size_t j = 0; j < n_faces-1; j++){
        // Step 2a:
        for(size_t k = 0; k < j; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j], maybe_newfaces[k],face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_not_newface[j] = 1;
                break;
                }
            }
        if (is_not_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 2b:
        for(size_t k = j+1; k < n_faces-1; k++){
            // Testing if maybe_newfaces[j] is contained in different nextface.
            if(is_subset(maybe_newfaces[j],maybe_newfaces[k], face_length)){
                // If so, it is not inclusion-maximal and hence not of codimension 1.
                is_not_newface[j] = 1;
            break;
            }
        }
        if (is_not_newface[j]) {
            // No further tests needed, if it is not of codimension 1.
            continue;
        }

        // Step 3:
        for (size_t k = 0; k < n_visited_all; k++){
            // Testing if maybe_newfaces[j] is contained in one,
            // we have already completely visited.
            if(is_subset(maybe_newfaces[j], visited_all[k], face_length)){
                // If so, we don't want to revisit.
                is_not_newface[j] = 1;
                break;
            }
        }
    }

    // Set ``newfaces`` to point to the correct ones.
    size_t n_newfaces = 0;  // length of newfaces2
    for (size_t j = 0; j < n_faces -1; j++){
        if (is_not_newface[j]) {
            // Not a new face of codimension 1.
            continue;
        }
        // It is a new face of codimension 1.
        newfaces[n_newfaces] = maybe_newfaces[j];
        n_newfaces++;
    }
    return n_newfaces;
}

size_t bit_rep_to_coatom_rep(uint64_t *face, uint64_t **coatoms, \
                               size_t n_coatoms, size_t face_length, \
                               size_t *output){
    /*
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    */
    size_t count_length = 0;
    for (size_t i = 0; i < n_coatoms; i++){
        if (is_subset(face, coatoms[i], face_length)){
            // ``face`` is contain in ``coatoms[i]``,
            // then ``i`` is an element in the coatom-representation.
            output[count_length] = i;
            count_length++;
        }
    }
    return count_length;
}
