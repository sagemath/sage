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
    ....:         import incidence_matrix_to_bit_rep_of_facets
    sage: P = polytopes.simplex(4)
    sage: inc = P.incidence_matrix()
    sage: mod_inc = inc.delete_columns([i for i,V in enumerate(P.Hrepresentation()) if V.is_equation()])
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(mod_inc)
    sage: face_list.compute_dimension()
    4

Obtain the Vrepresentation of a polyhedron as facet-incidences stored in
:class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`::

    sage: # needs sage.combinat
    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_rep_of_Vrep
    sage: P = polytopes.associahedron(['A',4])
    sage: face_list = incidence_matrix_to_bit_rep_of_Vrep(P.incidence_matrix())
    sage: face_list.compute_dimension()
    4

Obtain the facets of a polyhedron as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces` from a facet list::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_facets
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)

Likewise for the Vrep as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_Vrep
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_Vrep(facets, 4)

AUTHOR:

- Jonathan Kliem (2019-04)
"""

#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from memory_allocator                 cimport MemoryAllocator

from sage.matrix.matrix_dense cimport Matrix_dense

from .list_of_faces                   cimport ListOfFaces
from .face_data_structure             cimport face_next_atom, face_add_atom_safe, facet_set_coatom, face_clear

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


def _Vrep_list_to_bit_rep_wrapper(tup):
    r"""
    A function to allow doctesting of :func:`Vrep_list_to_bit_rep`.

    TESTS::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _Vrep_list_to_bit_rep_wrapper
        sage: _Vrep_list_to_bit_rep_wrapper((0, 3)).matrix()
        [1 0 0 1]
    """
    cdef ListOfFaces output = ListOfFaces(1, max(tup) + 1, 1)
    Vrep_list_to_bit_rep(tup, output.data.faces[0])
    return output

cdef int Vrep_list_to_bit_rep(tuple Vrep_list, face_t output) except -1:
    r"""
    Convert a vertex list into Bit-representation. Store it in ``output``.

    The first bit represents the entry ``0`` and is set to one, iff ``0`` is in
    ``vertex_list``. The second bit represents ``1`` and so on.

    INPUT:

    - ``vertex_list`` -- tuple of pairwise distinct positive integers that fit into ``output``
    - ``output`` -- an already initialized face

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _Vrep_list_to_bit_rep_wrapper
        sage: _Vrep_list_to_bit_rep_wrapper((0, 1)).matrix()
        [1 1]
        sage: _Vrep_list_to_bit_rep_wrapper((0, 2, 66)).matrix().nonzero_positions_in_row(0)
        [0, 2, 66]
        sage: _Vrep_list_to_bit_rep_wrapper((-1, 12))
        Traceback (most recent call last):
        ...
        OverflowError: can...t convert negative value to size_t
        sage: _Vrep_list_to_bit_rep_wrapper((0, 0))
        Traceback (most recent call last):
        ...
        ValueError: entries of ``tup`` are not distinct
    """
    cdef size_t entry     # will iterate over tup

    face_clear(output)
    if unlikely(len(Vrep_list) != len(set(Vrep_list))):
        raise ValueError("entries of ``tup`` are not distinct")
    for entry in Vrep_list:
        face_add_atom_safe(output, entry)

def _incidences_to_bit_rep_wrapper(tup):
    r"""
    A function to allow doctesting of :func:`incidences_to_bit_rep`.

    TESTS::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _incidences_to_bit_rep_wrapper
        sage: _incidences_to_bit_rep_wrapper((1,0,0,1)).matrix()
        [1 0 0 1]
    """
    cdef ListOfFaces output = ListOfFaces(1, len(tup), 1)
    incidences_to_bit_rep(tup, output.data.faces[0])
    return output

cdef int incidences_to_bit_rep(tuple incidences, face_t output) except -1:
    r"""
    Convert a tuple of incidences into Bit-representation.

    Store it in ``output``. Each entry in ``incidences`` represents a bit in
    ``output``. It is set to ``1``, iff the entry in ``incidences`` is non-zero.

    INPUT:

    - ``incidences`` -- tuple of integers representing incidences that fit into ``output``
    - ``output`` -- an already initialized face

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _incidences_to_bit_rep_wrapper
        sage: _incidences_to_bit_rep_wrapper((1,1)).matrix()
        [1 1]
        sage: _incidences_to_bit_rep_wrapper((1,0,1) + 61*(0,) +
        ....:                                (0,0,1,)).matrix().nonzero_positions_in_row(0)
        [0, 2, 66]
    """
    cdef size_t entry       # index for the entries in tup
    cdef size_t length = len(incidences)

    face_clear(output)
    for entry in range(length):
        if incidences[entry]:
            # Vrep ``entry`` is contained in the face, so set the corresponding bit
            face_add_atom_safe(output, entry)

def incidence_matrix_to_bit_rep_of_facets(Matrix_dense matrix):
    r"""
    Initialize facets in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      with columns corresponding to equations deleted
      of type :class:`sage.matrix.matrix_dense.Matrix_dense`

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import incidence_matrix_to_bit_rep_of_facets, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: P = polytopes.permutahedron(4)
        sage: inc = P.incidence_matrix()
        sage: mod_inc = inc.delete_columns([i for i,V in enumerate(P.Hrepresentation()) if V.is_equation()])
        sage: facets = incidence_matrix_to_bit_rep_of_facets(mod_inc)
        sage: facets.matrix().dimensions()
        (14, 24)
        sage: for row in facets.matrix():
        ....:     row.nonzero_positions()
        [18, 19, 20, 21, 22, 23]
        [3, 5, 8, 10, 12, 17]
        [2, 7, 11, 13, 20, 21]
        [2, 5, 12, 13]
        [4, 6, 14, 15, 19, 23]
        [3, 4, 8, 14]
        [6, 7, 21, 23]
        [2, 3, 4, 5, 6, 7]
        [0, 1, 9, 16, 18, 22]
        [0, 9, 10, 17]
        [1, 11, 20, 22]
        [0, 1, 10, 11, 12, 13]
        [15, 16, 18, 19]
        [8, 9, 14, 15, 16, 17]
    """
    # Output will be a ``ListOfFaces`` with ``matrix.ncols()`` faces and
    # ``matrix.nrows()`` Vrep.
    cdef size_t nrows = matrix._nrows
    cdef size_t ncols = matrix._ncols
    cdef ListOfFaces facets = ListOfFaces(ncols, nrows, ncols)
    cdef face_t output
    cdef size_t entry       # index for the entries in tup

    cdef size_t i
    for i in range(ncols):
        output = facets.data.faces[i]
        facet_set_coatom(output, i)

        # Filling each facet with its Vrep-incidences, which "is" the
        # "i-th column" of the original matrix (but we have transposed).
        for entry in range(nrows):
            if not matrix.get_is_zero_unsafe(entry, i):
                # Vrep ``entry`` is contained in the face, so set the corresponding bit
                face_add_atom_safe(output, entry)
    return facets

def incidence_matrix_to_bit_rep_of_Vrep(Matrix_dense matrix):
    r"""
    Initialize Vrepresentatives in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    Each Vrepresentative is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      with columns corresponding to equations deleted
      of type :class:`sage.matrix.matrix_dense.Matrix_dense`

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import incidence_matrix_to_bit_rep_of_Vrep, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: P = polytopes.permutahedron(4)
        sage: inc = P.incidence_matrix()
        sage: mod_inc = inc.delete_columns([i for i,V in enumerate(P.Hrepresentation()) if V.is_equation()])
        sage: vertices = incidence_matrix_to_bit_rep_of_Vrep(mod_inc)
        sage: vertices.matrix().dimensions()
        (24, 14)
        sage: for row in vertices.matrix():
        ....:     row.nonzero_positions()
        [8, 9, 11]
        [8, 10, 11]
        [2, 3, 7]
        [1, 5, 7]
        [4, 5, 7]
        [1, 3, 7]
        [4, 6, 7]
        [2, 6, 7]
        [1, 5, 13]
        [8, 9, 13]
        [1, 9, 11]
        [2, 10, 11]
        [1, 3, 11]
        [2, 3, 11]
        [4, 5, 13]
        [4, 12, 13]
        [8, 12, 13]
        [1, 9, 13]
        [0, 8, 12]
        [0, 4, 12]
        [0, 2, 10]
        [0, 2, 6]
        [0, 8, 10]
        [0, 4, 6]
    """
    return incidence_matrix_to_bit_rep_of_facets(matrix.transpose())

def facets_tuple_to_bit_rep_of_facets(tuple facets_input, size_t n_Vrep):
    r"""
    Initializes facets in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of Vrep,
      Vrep must be exactly ``range(n_Vrep)``
    - ``n_Vrep``

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_rep_of_facets, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
        sage: for i in range(8):
        ....:     print(_bit_rep_to_Vrep_list_wrapper(facets, i))
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
    cdef ListOfFaces facets = ListOfFaces(len(facets_input), n_Vrep, len(facets_input))
    for i in range(len(facets_input)):
        # filling each facet with the data from the corresponding facet
        Vrep_list_to_bit_rep(facets_input[i], facets.data.faces[i])
        facet_set_coatom(facets.data.faces[i], i)
    return facets

def facets_tuple_to_bit_rep_of_Vrep(tuple facets_input, size_t n_Vrep):
    r"""
    Initialize Vrepresentatives in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    Each Vrepresentative is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of Vrep,
      Vrep must be exactly ``range(n_Vrep)``
    - ``n_Vrep``

    OUTPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`


    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_rep_of_Vrep, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: vertices = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
        sage: for i in range(6):
        ....:     print(_bit_rep_to_Vrep_list_wrapper(vertices, i))
        (0, 3, 4, 7)
        (0, 1, 4, 5)
        (1, 2, 5, 6)
        (2, 3, 6, 7)
        (0, 1, 2, 3)
        (4, 5, 6, 7)
    """
    cdef size_t n_facets = len(facets_input)

    # Vertices in facet-representation will be a ``ListOfFaces``
    # with number of Vrep faces and
    # number of facets "Vrep"/atoms.
    cdef ListOfFaces Vrep = ListOfFaces(n_Vrep, n_facets, n_Vrep)
    cdef face_t* Vrep_data = Vrep.data.faces

    # Initializing the data of ListOfFaces.

    cdef size_t input_facet   # will iterate over indices of facets_input
    cdef size_t input_Vrep    # will iterate over vertices in facet ``input_facet``
    cdef size_t i

    for i in range(n_Vrep):
        facet_set_coatom(Vrep_data[i], i)

    for input_facet in range(n_facets):
        for input_Vrep in facets_input[input_facet]:
            # Iff the input-Vrep is in the input-facet,
            # then in facet-representation of the Vrep
            # input-facet is a Vrep of intput-Vrep.
            face_add_atom_safe(Vrep_data[input_Vrep], input_facet)
    return Vrep

def _bit_rep_to_Vrep_list_wrapper(ListOfFaces faces, index=0):
    r"""
    A function to test :func:`bit_rep_to_Vrep_list`.

    INPUT:

    - ``faces`` -- a :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`
    - ``index`` -- (default: ``0``); the face to obtain

    OUTPUT: The face as tuple of integers.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_rep_of_facets, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: faces = facets_tuple_to_bit_rep_of_facets(((1,5,123,1054),), 1055)
        sage: _bit_rep_to_Vrep_list_wrapper(faces, 0)
        (1, 5, 123, 1054)
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef size_t *output
    output = <size_t *> mem.allocarray(faces.n_atoms(),
                                       sizeof(size_t))

    length = bit_rep_to_Vrep_list(
            faces.data.faces[index], output)
    return tuple(output[i] for i in range(length))

cdef inline size_t bit_rep_to_Vrep_list(face_t face, size_t *output) except -1:
    r"""
    Convert a bitrep-representation to a list of Vrep. Return length of representation.

    Basically this is an inverse to :meth:`Vrep_list_to_bit_rep`.
    Instead of returning a tuple, it stores the Vrep in ``output``.

    INPUT:

    - ``face`` -- a Bit-representation of a face
    - ``output`` -- an array of ``size_t`` long enough to contain all Vrep
      of that face

    OUTPUT:

    - store Vrep in ``output``
    - return "length" of ``output``

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import _bit_rep_to_Vrep_list_wrapper, \
        ....:            _Vrep_list_to_bit_rep_wrapper
        sage: faces = _Vrep_list_to_bit_rep_wrapper((0, 4, 64, 65, 66, 67, 68))
        sage: _bit_rep_to_Vrep_list_wrapper(faces, 0)
        (0, 4, 64, 65, 66, 67, 68)
        sage: faces = _Vrep_list_to_bit_rep_wrapper((0, 2, 3))
        sage: _bit_rep_to_Vrep_list_wrapper(faces, 0)
        (0, 2, 3)
        sage: faces = _Vrep_list_to_bit_rep_wrapper((64, 66, 67, 68, 69))
        sage: _bit_rep_to_Vrep_list_wrapper(faces, 0)
        (64, 66, 67, 68, 69)

    TESTS:

    Testing that :meth`bit_rep_to_Vrep_list` is the
    inverse to :meth:`Vrep_list_to_bit_rep`::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import _bit_rep_to_Vrep_list_wrapper, \
        ....:            _Vrep_list_to_bit_rep_wrapper
        sage: for _ in range(10):
        ....:     st = set(randint(0,127) for i in range(40))
        ....:     tup = tuple(sorted(tuple(st)))
        ....:     faces = _Vrep_list_to_bit_rep_wrapper(tup)
        ....:     output = _bit_rep_to_Vrep_list_wrapper(faces, 0)
        ....:     if not tup == output:
        ....:         print('``bit_rep_to_Vrep_list`` does not behave',
        ....:               'as the inverse of ``Vrep_list_to_bit_rep``')
    """
    cdef size_t output_length = 0
    cdef long j = face_next_atom(face, 0)
    while (j != -1):
        output[output_length] = j
        output_length += 1
        j = face_next_atom(face, j+1)
    return output_length
