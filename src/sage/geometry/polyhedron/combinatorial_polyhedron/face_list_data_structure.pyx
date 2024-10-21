# sage_setup: distribution = sagemath-polyhedra
"""
Sorting of a list of faces.
"""
# ****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

cdef void sort_faces_list(face_list_t faces) noexcept:
    r"""
    Sorts faces in place.
    """
    cdef MemoryAllocator mem = MemoryAllocator()

    # Merge sort needs a second list of pointers.
    cdef face_t* extra_mem = <face_t*> mem.allocarray(faces.n_faces, sizeof(face_t))

    # Sort the faces using merge sort.
    _sort_faces_loop(faces.faces, faces.faces, extra_mem, faces.n_faces)

cdef void _sort_faces_loop(face_t* inp, face_t* out, face_t* extra_mem, size_t n_faces) noexcept:
    """
    This is merge sort.

    Sorts ``inp`` and returns it in ``out``.

    .. WARNING::

        ``inp`` is the same as ``out`` or ``extra_mem``.

    See :func:`sort_faces`.
    """
    if n_faces == 0:
        # Prevent it from crashing.
        # In this case there is nothing to do anyway.
        return

    if n_faces == 1:
        # The final case, where there is only one element.
        out[0][0] = inp[0][0]
        return

    cdef size_t middle = n_faces//2
    cdef size_t len_upper_half = n_faces - middle

    # Sort the upper and lower half of ``inp`` iteratively into ``output2``.
    _sort_faces_loop(inp, extra_mem, out, middle)
    _sort_faces_loop(inp+middle, extra_mem+middle,
                     out+middle, len_upper_half)

    # Merge lower and upper half into ``output1``.
    cdef size_t i = 0        # index through lower half
    cdef size_t j = middle   # index through upper half
    cdef size_t counter = 0  # counts how many elements have been "merged" already
    cdef int val
    while i < middle and j < n_faces:
        # Compare the lowest elements of lower and upper half.
        val = face_cmp(extra_mem[i], extra_mem[j])
        if val < 0:
            out[counter][0] = extra_mem[i][0]
            i += 1
            counter += 1
        else:
            out[counter][0] = extra_mem[j][0]
            j += 1
            counter += 1
    if i < middle:
        # Add the remaining elements of lower half.
        while i < middle:
            out[counter][0] = extra_mem[i][0]
            i += 1
            counter += 1
    else:
        # Add the remaining elements of upper half.
        while j < n_faces:
            out[counter][0] = extra_mem[j][0]
            j += 1
            counter += 1
