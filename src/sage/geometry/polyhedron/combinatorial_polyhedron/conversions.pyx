r"""
Conversions

This module provides conversions to :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces` from
- an incidence matrix of a polyhedron or
- a tuple of facets (as tuple of vertices each).

Also this module provides a conversion from the data of :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`,
which is a Bit-vector representing incidences of a face,
to a list of entries which are incident.

.. SEEALSO::

    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Obtain the facets of a polyhedron as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_repr_of_facets
    sage: P = polytopes.simplex(4)
    sage: face_list = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
    sage: face_list = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
    sage: face_list.compute_dimension()
    4

Obtain the Vrepresentation of a polyhedron as facet-incidences stored in
:class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_repr_of_Vrepr
    sage: P = polytopes.associahedron(['A',4])
    sage: face_list = incidence_matrix_to_bit_repr_of_Vrepr(P.incidence_matrix())
    sage: face_list.compute_dimension()
    4

Obtain the facets of a polyhedron as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces` from a facet list::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_repr_of_facets
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_repr_of_facets(facets, 4)

Likewise for the Vrepr as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_repr_of_Vrepr
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_repr_of_Vrepr(facets, 4)

AUTHOR:

- Jonathan Kliem (2019-04)
"""

#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__             import absolute_import, division
from sage.structure.element import is_Matrix

from libc.string            cimport memset
from .list_of_faces         cimport ListOfFaces

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef inline uint64_t vertex_to_bit_dictionary(size_t i):
    r"""
    Return an uint64_t with exactly the i-th bit set to 1.

    This "dictionary" helps storing a vector of 64 incidences as ``uint64_t``.
    """
    return (<uint64_t>1) << (64 - i - 1)

cdef int Vrepr_list_to_bit_repr(tuple Vrepr_list, uint64_t *output,
                                 size_t face_length) except -1:
    r"""
    Convert a vertex list into Bit-representation. Store it in ``output``.

    The first bit represents the entry ``0`` and is set to one, iff ``0`` is in
    ``vertex_list``. The second bit represents ``1`` and so on.

    INPUT:

    - ``vertex_list`` -- tuple of pairwise distinct positive integers in
      ``range(face_length*64)``
    - ``output`` -- array of ``uint64_t`` of length ``face_length``
    - ``face_length``

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport Vrepr_list_to_bit_repr
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer cimport smallInteger
        ....:
        ....: def Vrepr_list_to_bit_repr_wrapper(tup):
        ....:     cdef size_t face_length = max(tup)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     Vrepr_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....:
        ....: def Vrepr_list_to_bit_repr_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     Vrepr_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: Vrepr_list_to_bit_repr_wrapper((62, 63))  # long time
        (3,)
        sage: Vrepr_list_to_bit_repr_wrapper((61, 63, 125))  # long time
        (5, 4)
        sage: Vrepr_list_to_bit_repr_wrong_size((62, 70))  # long time
        Traceback (most recent call last):
        ...
        IndexError: output too small to represent 70
        sage: Vrepr_list_to_bit_repr_wrapper((-1, 12))  # long time
        Traceback (most recent call last):
        ...
        OverflowError: can't convert negative value to size_t
        sage: Vrepr_list_to_bit_repr_wrapper((0, 0))  # long time
        Traceback (most recent call last):
        ...
        ValueError: entries of ``tup`` are not distinct
    """
    cdef size_t entry     # will iterate over tup
    cdef size_t position  # determines the position in output of entry
    cdef size_t value     # determines which bit will be set in output[position]

    memset(output, 0, face_length*8)  # initialize output
    if unlikely(len(Vrepr_list) != len(set(Vrepr_list))):
        raise ValueError("entries of ``tup`` are not distinct")
    for entry in Vrepr_list:
        value = entry % 64
        position = entry//64
        if unlikely(position >= face_length):
            raise IndexError("output too small to represent %s"%entry)
        output[position] += vertex_to_bit_dictionary(value)

cdef int incidences_to_bit_repr(tuple incidences, uint64_t *output,
                                size_t face_length) except -1:

    r"""
    Convert a tuple of incidences into Bit-representation.

    Store it in ``output``. Each entry in ``incidences`` represents a bit in
    ``output``. It is set to ``1``, iff the entry in ``incidences`` is non-zero.

    INPUT:

    - ``incidences`` -- tuple of integers representing incidences
      of length at most ``face_length*64``
    - ``output`` -- array of ``uint64_t`` of length ``face_length``
    - ``face_length``

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport incidences_to_bit_repr
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer cimport smallInteger
        ....:
        ....: def incidences_to_bit_reprs_wrapper(tup):
        ....:     cdef size_t face_length = (len(tup)-1)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:          <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....:
        ....: def incidences_to_bit_reprs_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:         <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: incidences_to_bit_reprs_wrapper((0,) * 62 + (1,1))  # long time
        (3,)
        sage: incidences_to_bit_reprs_wrapper((0,) * 61 + (1,0,1) +  # long time
        ....:                                 (0,) * 61 + (1,))
        (5, 4)
        sage: incidences_to_bit_reprs_wrapper((1,) * 64)  # long time
        (-1,)
        sage: incidences_to_bit_reprs_wrong_size((1,) * 70)  # long time
        Traceback (most recent call last):
        ...
        IndexError: output too small to represent all incidences
    """
    cdef size_t entry       # index for the entries in tup
    cdef size_t position    # determines the position in output of entry
    cdef size_t value       # determines which bit will be set in output[position]
    cdef size_t length = len(incidences)

    memset(output, 0, face_length*8)  #initialize
    if unlikely(length > 64*face_length):
        raise IndexError("output too small to represent all incidences")
    for entry in range(length):
        if incidences[entry]:
            # Vrepr ``entry`` is contained in the face, so set the corresponding bit
            value = entry % 64
            position = entry//64
            output[position] += vertex_to_bit_dictionary(value)

def incidence_matrix_to_bit_repr_of_facets(matrix):
    r"""
    Initialize facets in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....: cimport ListOfFaces
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list
        ....: from sage.rings.integer cimport smallInteger
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_Vrepr_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.n_atoms,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_Vrepr_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
        ....: ''') # long time

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import incidence_matrix_to_bit_repr_of_facets
        sage: P = polytopes.permutahedron(4)
        sage: facets = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
        sage: facets.n_faces
        14
        sage: facets.n_atoms
        24
        sage: for i in range(facets.n_faces):  # long time
        ....:     print(bit_repr_to_Vrepr_list_wrapper(facets, i))
        (18, 19, 20, 21, 22, 23)
        (9, 11, 15, 17, 21, 23)
        (16, 17, 22, 23)
        (0, 1, 2, 3, 4, 5)
        (2, 4, 8, 10)
        (0, 1, 6, 7)
        (6, 7, 12, 13, 18, 19)
        (3, 5, 9, 11)
        (1, 3, 7, 9, 13, 15)
        (0, 2, 6, 8, 12, 14)
        (12, 14, 18, 20)
        (4, 5, 10, 11, 16, 17)
        (8, 10, 14, 16, 20, 22)
        (13, 15, 19, 21)
    """

    if unlikely(not is_Matrix(matrix)):
        raise ValueError("input must be matrix")

    # The incidence-matrix of a polyhedron, might contain columns with all 1's.
    # Those correspond to hyperplanes, the polyhedron lies in, not facets.
    # We delete those from the matrix. Also we transpose.
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]

    # Output will be a ``ListOfFaces`` with ``matrix.nrows()`` faces and
    # ``matrix.ncols()`` Vrepr.
    cdef size_t length = matrix.nrows()
    cdef ListOfFaces facets = ListOfFaces(length, matrix.ncols())
    cdef uint64_t **facets_data = facets.data

    cdef size_t i
    for i in range(length):
        # Filling each facet with its Vrepr-incidences, which "is" the
        # "i-th column" of the original matrix (but we have transposed).
        incidences_to_bit_repr(tuple(matrix[i]), facets_data[i], facets.face_length)
    return facets

def incidence_matrix_to_bit_repr_of_Vrepr(matrix):
    r"""
    Initialize Vrepresentatives in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    Each Vrepresenative is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....: cimport ListOfFaces
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list
        ....: from sage.rings.integer cimport smallInteger
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_Vrepr_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.n_atoms,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_Vrepr_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import incidence_matrix_to_bit_repr_of_Vrepr
        sage: P = polytopes.permutahedron(4)
        sage: vertices = incidence_matrix_to_bit_repr_of_Vrepr(P.incidence_matrix())
        sage: vertices.n_faces
        24
        sage: vertices.n_atoms
        14
        sage: for i in range(vertices.n_faces):
        ....:     print(bit_repr_to_Vrepr_list_wrapper(vertices, i))
        (3, 5, 9)
        (3, 5, 8)
        (3, 4, 9)
        (3, 7, 8)
        (3, 4, 11)
        (3, 7, 11)
        (5, 6, 9)
        (5, 6, 8)
        (4, 9, 12)
        (1, 7, 8)
        (4, 11, 12)
        (1, 7, 11)
        (6, 9, 10)
        (6, 8, 13)
        (9, 10, 12)
        (1, 8, 13)
        (2, 11, 12)
        (1, 2, 11)
        (0, 6, 10)
        (0, 6, 13)
        (0, 10, 12)
        (0, 1, 13)
        (0, 2, 12)
        (0, 1, 2)
    """
    if unlikely(not is_Matrix(matrix)):
        raise ValueError("input must be matrix")

    # The incidence-matrix of a polyhedron, might contain columns with all 1's.
    # Those correspond to hyperplanes, the polyhedron lies in, not facets.
    # We delete those from the matrix. Also we transpose.
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]

    # Output will be a ``ListOfFaces`` with ``matrix.ncols()`` faces and
    # ``matrix.nrows()`` Vrepr.
    cdef size_t length = matrix.ncols()
    cdef ListOfFaces Vrepr = ListOfFaces(length, matrix.nrows())
    cdef uint64_t **Vrepr_data = Vrepr.data

    cdef size_t i
    for i in range(length):
        # Filling each facet with its Vrepr-incidences, which "is" the
        # "i-th row" of the original matrix (but we have transposed).
        incidences_to_bit_repr(tuple(matrix.column(i)), Vrepr_data[i], Vrepr.face_length)
    return Vrepr

def facets_tuple_to_bit_repr_of_facets(tuple facets_input, size_t n_Vrepr):
    r"""
    Initializes facets in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of Vrepr,
      Vrepr must be exactly ``range(n_Vrepr)``
    - ``n_Vrepr``

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....: cimport ListOfFaces
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list
        ....: from sage.rings.integer cimport smallInteger
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_Vrepr_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.n_atoms,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_Vrepr_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
        ....: ''') # long time

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_repr_of_facets
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = facets_tuple_to_bit_repr_of_facets(bi_pyr, 6)
        sage: for i in range(8):  # long time
        ....:     print(bit_repr_to_Vrepr_list_wrapper(facets, i))
        (0, 1, 4)
        (1, 2, 4)
        (2, 3, 4)
        (0, 3, 4)
        (0, 1, 5)
        (1, 2, 5)
        (2, 3, 5)
        (0, 3, 5)
    """
    cdef Py_ssize_t i
    cdef ListOfFaces facets = ListOfFaces(len(facets_input), n_Vrepr)
    cdef size_t face_length = facets.face_length
    cdef uint64_t **facets_data = facets.data
    for i in range(len(facets_input)):
        # filling each facet with the the data from the corresponding facet
        Vrepr_list_to_bit_repr(facets_input[i], facets_data[i], face_length)
    return facets

def facets_tuple_to_bit_repr_of_Vrepr(tuple facets_input, size_t n_Vrepr):
    r"""
    Initialize Vrepresentatives in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    Each Vrepresentative is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of Vrepr,
      Vrepr must be exactly ``range(n_Vrepr)``
    - ``n_Vrepr``

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`


    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....: cimport ListOfFaces
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list
        ....: from sage.rings.integer cimport smallInteger
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_Vrepr_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.n_atoms,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_Vrepr_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_repr_of_Vrepr
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: vertices = facets_tuple_to_bit_repr_of_Vrepr(bi_pyr, 6)
        sage: for i in range(6):
        ....:     print(bit_repr_to_Vrepr_list_wrapper(vertices, i))
        (0, 3, 4, 7)
        (0, 1, 4, 5)
        (1, 2, 5, 6)
        (2, 3, 6, 7)
        (0, 1, 2, 3)
        (4, 5, 6, 7)
    """
    cdef size_t n_facets = len(facets_input)

    # Vertices in facet-representation will be a ``ListOfFaces``
    # with number of Vrepr faces and
    # number of facets "Vrepr"/atoms.
    cdef ListOfFaces Vrepr = ListOfFaces(n_Vrepr, n_facets)
    cdef uint64_t **Vrepr_data = Vrepr.data

    # Initializing the data of ListOfFaces.
    cdef size_t face_length = Vrepr.face_length
    cdef size_t i
    for i in range(n_Vrepr):
        memset(Vrepr_data[i], 0, face_length*8)

    cdef size_t input_facet   # will iterate over indices of facets_input
    cdef size_t position      # determines the position in output of entry
    cdef size_t value         # determines which bit will be set in output[position]
    cdef size_t input_Vrepr   # will iterate over vertices in facet ``input_facet``

    for input_facet in range(n_facets):
        value = input_facet % 64
        position = input_facet//64
        for input_Vrepr in facets_input[input_facet]:
            # Iff the input-Vrepr is in the input-facet,
            # then in facet-representation of the Vrepr
            # input-facet is a Vrepr of intput-Vrepr.
            Vrepr_data[input_Vrepr][position] += vertex_to_bit_dictionary(value)
    return Vrepr

cdef inline size_t bit_repr_to_Vrepr_list(uint64_t *face, size_t *output,
                                           size_t face_length) except -1:
    r"""
    Convert a bitrep-representation to a list of Vrepr. Return length of representation.

    Basically this is an inverse to :meth:`Vrepr_list_to_bit_repr`.
    Instead of returning a tuple, it stores the Vrepr in ``output``.

    INPUT:

    - ``face`` -- a Bit-representation of a face
    - ``output`` -- an array of ``size_t`` long enough to contain all Vrepr
      of that face (``face_length*64`` will suffice)
    - ``face_length`` -- the length of ``face``

    OUTPUT:

    - store Vrepr in ``output``
    - return "length" of ``output``

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....: cimport ListOfFaces
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list, Vrepr_list_to_bit_repr
        ....: from sage.rings.integer cimport smallInteger
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_Vrepr_list_wrapper(tup):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef length = len(tup)
        ....:     output = <size_t *> mem.allocarray(length*64,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data
        ....:     data = <uint64_t *> mem.allocarray(length, 8)
        ....:     for i in range(len(tup)):
        ....:         data[i] = tup[i]
        ....:     outputlength = bit_repr_to_Vrepr_list(
        ....:         data, output, length)
        ....:     return tuple(smallInteger(output[i]) for i in range(outputlength))
        ....: ''')  # long time

        sage: bit_repr_to_Vrepr_list_wrapper((17, 31))  # long time
        (59, 63, 123, 124, 125, 126, 127)
        sage: bit_repr_to_Vrepr_list_wrapper((13,))  # long time
        (60, 61, 63)
        sage: bit_repr_to_Vrepr_list_wrapper((0, 61))  # long time
        (122, 123, 124, 125, 127)

    TESTS:

    Testing that :meth`bit_repr_to_Vrepr_list` is the
    inverse to :meth:`Vrepr_list_to_bit_repr`::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....: cimport bit_repr_to_Vrepr_list, Vrepr_list_to_bit_repr
        ....: from sage.misc.prandom import randint
        ....:
        ....: cdef uint64_t[2] face
        ....: cdef size_t length
        ....: cdef size_t[128] output
        ....:
        ....: for _ in range(10):
        ....:     st = set(randint(0,127) for i in range(40))
        ....:     tup = tuple(sorted(tuple(st)))
        ....:     Vrepr_list_to_bit_repr(tup, face, 2)
        ....:     length = bit_repr_to_Vrepr_list(face, output, 2)
        ....:     tup2 = tuple(output[i] for i in range(length))
        ....:     if not tup == tup2:
        ....:         print('``bit_repr_to_Vrepr_list`` does not behave',
        ....:               'as the inverse of ``Vrepr_list_to_bit_repr``')
        ....: ''') # long time
    """
    cdef size_t i, j
    cdef size_t output_length = 0
    cdef uint64_t copy
    for i in range(face_length):
        if face[i]:
            copy = face[i]
            for j in range(64):
                if copy >= vertex_to_bit_dictionary(j):
                    # Then face[i] has the j-th bit set to 1.
                    # This corresponds to face containing Vrepr i*64 + j.
                    output[output_length] = i*64 + j
                    output_length += 1
                    copy -= vertex_to_bit_dictionary(j)
    return output_length
