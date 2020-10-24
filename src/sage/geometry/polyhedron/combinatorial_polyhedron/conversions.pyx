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
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(mod_inc)
    sage: face_list.compute_dimension()
    4

Obtain the Vrepresentation of a polyhedron as facet-incidences stored in
:class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`::

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
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import is_Matrix

from libc.string                      cimport memset
from .list_of_faces                   cimport ListOfFaces
from sage.misc.superseded              import deprecated_function_alias
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.ext.memory_allocator        cimport MemoryAllocator

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef inline uint64_t vertex_to_bit_dictionary(size_t i):
    r"""
    Return an uint64_t with exactly the i-th bit set to 1.

    This "dictionary" helps storing a vector of 64 incidences as ``uint64_t``.
    """
    return (<uint64_t>1) << (64 - i - 1)

def _Vrep_list_to_bit_rep_wrapper(tup, size_t face_length=-1):
    r"""
    A function to allow doctesting of :func:`Vrep_list_to_bit_rep`.

    TESTS::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _Vrep_list_to_bit_rep_wrapper
        sage: _Vrep_list_to_bit_rep_wrapper((60, 63))
        (9,)
    """
    if face_length == -1:
        face_length = max(tup)//64 + 1
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
    Vrep_list_to_bit_rep(tup, output, face_length)
    return tuple(output[i] for i in range(face_length))

cdef int Vrep_list_to_bit_rep(tuple Vrep_list, uint64_t *output,
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

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _Vrep_list_to_bit_rep_wrapper
        sage: _Vrep_list_to_bit_rep_wrapper((62, 63))
        (3,)
        sage: _Vrep_list_to_bit_rep_wrapper((61, 63, 125))
        (5, 4)
        sage: _Vrep_list_to_bit_rep_wrapper((62, 70), face_length=1)
        Traceback (most recent call last):
        ...
        IndexError: output too small to represent 70
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
    cdef size_t position  # determines the position in output of entry
    cdef size_t value     # determines which bit will be set in output[position]

    memset(output, 0, face_length*8)  # initialize output
    if unlikely(len(Vrep_list) != len(set(Vrep_list))):
        raise ValueError("entries of ``tup`` are not distinct")
    for entry in Vrep_list:
        value = entry % 64
        position = entry//64
        if unlikely(position >= face_length):
            raise IndexError("output too small to represent %s"%entry)
        output[position] += vertex_to_bit_dictionary(value)

def _incidences_to_bit_rep_wrapper(tup, size_t face_length=-1):
    r"""
    A function to allow doctesting of :func:`incidences_to_bit_rep`.

    TESTS::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _incidences_to_bit_rep_wrapper
        sage: _incidences_to_bit_rep_wrapper((0,) * 60 + (1,0,0,1))
        (9,)
    """
    if face_length == -1:
        face_length = (len(tup)-1)//64 + 1
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint64_t *output = \
         <uint64_t *> mem.allocarray(face_length, 8)
    incidences_to_bit_rep(tup, output, face_length)
    return tuple(output[i] for i in range(face_length))

cdef int incidences_to_bit_rep(tuple incidences, uint64_t *output,
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

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:         import _incidences_to_bit_rep_wrapper
        sage: _incidences_to_bit_rep_wrapper((0,) * 62 + (1,1))
        (3,)
        sage: _incidences_to_bit_rep_wrapper((0,) * 61 + (1,0,1) +
        ....:                                 (0,) * 61 + (1,))
        (5, 4)
        sage: _incidences_to_bit_rep_wrapper((1,) * 64) == (2**64-1,)
        True
        sage: _incidences_to_bit_rep_wrapper((1,) * 70, face_length=1)
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
            # Vrep ``entry`` is contained in the face, so set the corresponding bit
            value = entry % 64
            position = entry//64
            output[position] += vertex_to_bit_dictionary(value)

def incidence_matrix_to_bit_rep_of_facets(Matrix_integer_dense matrix):
    r"""
    Initialize facets in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      with columns corresponding to equations deleted
      of type :class:`sage.matrix.matrix_integer_dense.Matrix_integer_dense`

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
        sage: facets.n_faces
        14
        sage: facets.n_atoms
        24
        sage: for i in range(facets.n_faces):
        ....:     print(_bit_rep_to_Vrep_list_wrapper(facets, i))
        (18, 19, 20, 21, 22, 23)
        (3, 5, 8, 10, 12, 17)
        (2, 7, 11, 13, 20, 21)
        (2, 5, 12, 13)
        (4, 6, 14, 15, 19, 23)
        (3, 4, 8, 14)
        (6, 7, 21, 23)
        (2, 3, 4, 5, 6, 7)
        (0, 1, 9, 16, 18, 22)
        (0, 9, 10, 17)
        (1, 11, 20, 22)
        (0, 1, 10, 11, 12, 13)
        (15, 16, 18, 19)
        (8, 9, 14, 15, 16, 17)
    """
    # Output will be a ``ListOfFaces`` with ``matrix.ncols()`` faces and
    # ``matrix.nrows()`` Vrep.
    cdef size_t nrows = matrix._nrows
    cdef size_t ncols = matrix._ncols
    cdef ListOfFaces facets = ListOfFaces(ncols, nrows)
    cdef uint64_t **facets_data = facets.data
    cdef uint64_t *output
    cdef size_t face_length = facets.face_length
    cdef size_t entry       # index for the entries in tup
    cdef size_t position    # determines the position in output of entry
    cdef size_t value       # determines which bit will be set in output[position]

    cdef size_t i
    for i in range(ncols):
        output = facets_data[i]
        memset(output, 0, face_length*8)  #initialize

        # Filling each facet with its Vrep-incidences, which "is" the
        # "i-th column" of the original matrix (but we have transposed).
        for entry in range(nrows):
            if not matrix.get_is_zero_unsafe(entry, i):
                # Vrep ``entry`` is contained in the face, so set the corresponding bit
                value = entry % 64
                position = entry//64
                output[position] += vertex_to_bit_dictionary(value)
    return facets
incidence_matrix_to_bit_repr_of_facets = deprecated_function_alias(28608, incidence_matrix_to_bit_rep_of_facets)

def incidence_matrix_to_bit_rep_of_Vrep(Matrix_integer_dense matrix):
    r"""
    Initialize Vrepresentatives in Bit-representation as :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.

    Each Vrepresentative is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      with columns corresponding to equations deleted
      of type :class:`sage.matrix.matrix_integer_dense.Matrix_integer_dense`

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
        sage: vertices.n_faces
        24
        sage: vertices.n_atoms
        14
        sage: for i in range(vertices.n_faces):
        ....:     print(_bit_rep_to_Vrep_list_wrapper(vertices, i))
        (8, 9, 11)
        (8, 10, 11)
        (2, 3, 7)
        (1, 5, 7)
        (4, 5, 7)
        (1, 3, 7)
        (4, 6, 7)
        (2, 6, 7)
        (1, 5, 13)
        (8, 9, 13)
        (1, 9, 11)
        (2, 10, 11)
        (1, 3, 11)
        (2, 3, 11)
        (4, 5, 13)
        (4, 12, 13)
        (8, 12, 13)
        (1, 9, 13)
        (0, 8, 12)
        (0, 4, 12)
        (0, 2, 10)
        (0, 2, 6)
        (0, 8, 10)
        (0, 4, 6)
    """
    return incidence_matrix_to_bit_rep_of_facets(matrix.transpose())
incidence_matrix_to_bit_repr_of_Vrepr = deprecated_function_alias(28608, incidence_matrix_to_bit_rep_of_Vrep)

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
    cdef ListOfFaces facets = ListOfFaces(len(facets_input), n_Vrep)
    cdef size_t face_length = facets.face_length
    cdef uint64_t **facets_data = facets.data
    for i in range(len(facets_input)):
        # filling each facet with the data from the corresponding facet
        Vrep_list_to_bit_rep(facets_input[i], facets_data[i], face_length)
    return facets
facets_tuple_to_bit_repr_of_facets = deprecated_function_alias(28608, facets_tuple_to_bit_rep_of_facets)

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
    cdef ListOfFaces Vrep = ListOfFaces(n_Vrep, n_facets)
    cdef uint64_t **Vrep_data = Vrep.data

    # Initializing the data of ListOfFaces.
    cdef size_t face_length = Vrep.face_length
    cdef size_t i
    for i in range(n_Vrep):
        memset(Vrep_data[i], 0, face_length*8)

    cdef size_t input_facet   # will iterate over indices of facets_input
    cdef size_t position      # determines the position in output of entry
    cdef size_t value         # determines which bit will be set in output[position]
    cdef size_t input_Vrep    # will iterate over vertices in facet ``input_facet``

    for input_facet in range(n_facets):
        value = input_facet % 64
        position = input_facet//64
        for input_Vrep in facets_input[input_facet]:
            # Iff the input-Vrep is in the input-facet,
            # then in facet-representation of the Vrep
            # input-facet is a Vrep of intput-Vrep.
            Vrep_data[input_Vrep][position] += vertex_to_bit_dictionary(value)
    return Vrep
facets_tuple_to_bit_repr_of_Vrepr = deprecated_function_alias(28608, facets_tuple_to_bit_rep_of_Vrep)

def _bit_rep_to_Vrep_list_wrapper(data, index=0):
    r"""
    A function to test :func:`bit_rep_to_Vrep_list`.

    INPUT:

    - ``data`` -- either a :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`
    or a tuple of integers in ``range(0,2**64)``
    - ``index`` -- ``0`` if ``data`` is a tuple, otherwise the index of the ``face``
      to convert

    OUTPUT: A tuple of integers.

    If the input is a tuple, it will be interpreted as a list of faces over with `64` atoms per element in the
    tuple. Each number in the tuple corresponds to an ``uint64_t``.

    The list of faces is then translated into a tuple of the integers with set bits.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_rep_of_facets, \
        ....:            _bit_rep_to_Vrep_list_wrapper
        sage: _bit_rep_to_Vrep_list_wrapper((1, 1))
        (63, 127)
        sage: faces = facets_tuple_to_bit_rep_of_facets(((1,5,123,1054),), 1055)
        sage: _bit_rep_to_Vrep_list_wrapper(faces, 0)
        (1, 5, 123, 1054)
    """
    cdef ListOfFaces faces
    if isinstance(data, ListOfFaces):
        faces = data
    else:
        assert isinstance(data, tuple)
        faces = ListOfFaces(1, 64*len(data))
        for i in range(len(data)):
            faces.data[0][i] = data[i]

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef size_t *output
    output = <size_t *> mem.allocarray(faces.n_atoms,
                                       sizeof(size_t))
    length = bit_rep_to_Vrep_list(
        faces.data[index], output, faces.face_length)
    return tuple(output[i] for i in range(length))

cdef inline size_t bit_rep_to_Vrep_list(uint64_t *face, size_t *output,
                                           size_t face_length) except -1:
    r"""
    Convert a bitrep-representation to a list of Vrep. Return length of representation.

    Basically this is an inverse to :meth:`Vrep_list_to_bit_rep`.
    Instead of returning a tuple, it stores the Vrep in ``output``.

    INPUT:

    - ``face`` -- a Bit-representation of a face
    - ``output`` -- an array of ``size_t`` long enough to contain all Vrep
      of that face (``face_length*64`` will suffice)
    - ``face_length`` -- the length of ``face``

    OUTPUT:

    - store Vrep in ``output``
    - return "length" of ``output``

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import _bit_rep_to_Vrep_list_wrapper
        sage: _bit_rep_to_Vrep_list_wrapper((17, 31))
        (59, 63, 123, 124, 125, 126, 127)
        sage: _bit_rep_to_Vrep_list_wrapper((13,))
        (60, 61, 63)
        sage: _bit_rep_to_Vrep_list_wrapper((0, 61))
        (122, 123, 124, 125, 127)

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
    cdef size_t i, j
    cdef size_t output_length = 0
    cdef uint64_t copy
    for i in range(face_length):
        if face[i]:
            copy = face[i]
            for j in range(64):
                if copy >= vertex_to_bit_dictionary(j):
                    # Then face[i] has the j-th bit set to 1.
                    # This corresponds to face containing Vrep i*64 + j.
                    output[output_length] = i*64 + j
                    output_length += 1
                    copy -= vertex_to_bit_dictionary(j)
    return output_length
