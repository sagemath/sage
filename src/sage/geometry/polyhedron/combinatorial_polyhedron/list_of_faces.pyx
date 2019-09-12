r"""
List of faces

This module provides a class to store faces of a polyhedron in Bit-representation.

This class allocates memory to store the faces in.
A face will be stored as vertex-incidences, where each Bit represents an incidence.
In :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions` there a methods to actually convert facets of a polyhedron
to bit-representations of vertices stored in :class:`ListOfFaces`.

Moreover, :class:`ListOfFaces` calculates the dimension of a polyhedron, assuming the
faces are the facets of this polyhedron.

Each face is stored over-aligned according to :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.bit_vector_operations.chunktype`.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Provide enough space to store `20` faces as incidences to `60` vertices::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
    ....: import ListOfFaces
    sage: face_list = ListOfFaces(20, 60)

Obtain the facets of a polyhedron::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_repr_of_facets
    sage: P = polytopes.cube()
    sage: face_list = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
    sage: face_list = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the Vrepresentation of a polyhedron as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_repr_of_Vrepr
    sage: P = polytopes.associahedron(['A',3])
    sage: face_list = incidence_matrix_to_bit_repr_of_Vrepr(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the facets of a polyhedron as :class:`ListOfFaces` from a facet list::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_repr_of_facets
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_repr_of_facets(facets, 4)

Likewise for the Vrepresenatives as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_repr_of_Vrepr
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_repr_of_Vrepr(facets, 4)

.. SEEALSO::

    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.base`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_faces_lattice`.

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

from cysignals.signals      cimport sig_on, sig_off
from .bit_vector_operations cimport chunksize, get_next_level, count_atoms

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class ListOfFaces:
    r"""
    A class to store the Bit-representation of faces in.

    This class will allocate the memory for the faces.

    INPUT:

    - ``n_faces`` -- the number of faces to be stored
    - ``n_atoms`` -- the total number of atoms the faces contain

    .. SEEALSO::

        :meth:`incidence_matrix_to_bit_repr_of_facets`,
        :meth:`incidence_matrix_to_bit_repr_of_Vrepr`,
        :meth:`facets_tuple_to_bit_repr_of_facets`,
        :meth:`facets_tuple_to_bit_repr_of_Vrepr`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....:     import ListOfFaces
        sage: facets = ListOfFaces(5, 13)
        sage: facets.face_length in (1, 2, 4)
        True
        sage: facets.n_atoms
        13
        sage: facets.n_faces
        5
    """
    def __init__(self, size_t n_faces, size_t n_atoms):
        r"""
        Initialize :class:`ListOfFaces`.

        See :class:`ListOfFaces`.

        TESTS:

        Checking for correct alignment of the data::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....: cimport ListOfFaces
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef size_t address
            ....: cdef size_t required_alignment
            ....:
            ....: facets = ListOfFaces(10, 13)
            ....: required_alignment = facets.face_length*8
            ....: for i in range(10):
            ....:     address = <size_t> facets.data[i]
            ....:     if not address == address & ~(required_alignment - 1):
            ....:         print('Alignment not correct')
            ....: ''')

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces).run()
        """
        self.n_faces = n_faces
        self.n_atoms = n_atoms
        self._mem = MemoryAllocator()

        # ``data`` will point to the faces as ``*uint64_t``.
        self.data = <uint64_t **> self._mem.allocarray(n_faces, sizeof(uint64_t *))

        # ``face_length`` is the length in terms of ``uint64_t``
        # NOTE: This needs to be divisible by 2, if chunksize is 128
        #       and divisible by 4, if chunksize is 256.
        self.face_length = ((n_atoms - 1)//chunksize + 1)*chunksize//64


        cdef size_t i
        for i in range(n_faces):
            # Allocate the memory for the i-th face.
            # We must allocate the memory for ListOfFaces overaligned:
            # - must be 16-byte aligned if chunksize = 128
            # - must be 32-byte aligned if chunksize = 256
            self.data[i] = <uint64_t *> \
                self._mem.aligned_malloc(chunksize//8, self.face_length*8)

    cpdef int compute_dimension(self) except -2:
        r"""
        Compute the dimension of a polyhedron by its facets.

        This assumes that ``self`` is the list of facets of a polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_repr_of_facets, \
            ....:            facets_tuple_to_bit_repr_of_Vrepr
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_repr_of_facets(bi_pyr, 6)
            sage: Vrepr = facets_tuple_to_bit_repr_of_Vrepr(bi_pyr, 6)
            sage: facets.compute_dimension()
            3
            sage: Vrepr.compute_dimension()
            3

        ALGORITHM:

        This is done by iteration:

        Computes the facets of one of the facets (i.e. the ridges contained in
        one of the facets). Then computes the dimension of the facet, by
        considering its facets.

        Repeats until a face has only one facet. Usually this is a vertex.

        However, in the unbounded case, this might be different. The face with only
        one facet might be a ray or a line. So the correct dimension of a
        polyhedron with one facet is the number of ``[lines, rays, vertices]``
        that the facet contains.

        Hence, we know the dimension of a face, which has only one facet and
        iteratively we know the dimension of entire polyhedron we started from.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import incidence_matrix_to_bit_repr_of_facets, \
            ....:            incidence_matrix_to_bit_repr_of_Vrepr
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: for _ in range(10):
            ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
            ....:                    for _ in range(randint(3,15)))
            ....:     P = Polyhedron(vertices=points)
            ....:     facets = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
            ....:     vertices = incidence_matrix_to_bit_repr_of_Vrepr(P.incidence_matrix())
            ....:     d1 = P.dimension()
            ....:     if d1 == 0:
            ....:         continue
            ....:     d2 = facets.compute_dimension()
            ....:     d3 = vertices.compute_dimension()
            ....:     if not d1 == d2 == d3:
            ....:         print('calculation_dimension() seems to be incorrect')
        """
        if self.n_faces == 0:
            raise TypeError("at least one face needed")
        return self.compute_dimension_loop(self.data, self.n_faces, self.face_length)

    cdef int compute_dimension_loop(self, uint64_t **faces, size_t n_faces,
                                      size_t face_length) except -2:
        r"""
        Compute the dimension of a polyhedron by its facets.

        INPUT:

        - ``faces`` -- facets in Bit-representation
        - ``n_faces`` -- length of facesdata
        - ``face_length`` -- the length of each face in terms of ``uint64_t``

        OUTPUT:

        - dimension of the polyhedron

        .. SEEALSO::

            :meth:`compute_dimension`
        """
        if n_faces == 0:
            raise TypeError("wrong usage of ``compute_dimension_loop``,\n" +
                            "at least one face needed.")

        if n_faces == 1:
            # We expect the face to be the empty polyhedron.
            # Possibly it contains more than one vertex/rays/lines.
            # The dimension of a polyhedron with this face as only facet is
            # the number of atoms it contains.
            return count_atoms(faces[0], face_length)

        # ``maybe_newfaces`` are all intersection of ``faces[n_faces -1]`` with previous faces.
        # It needs to be allcoated to store those faces.
        cdef ListOfFaces maybe_newfaces_mem = ListOfFaces(n_faces, face_length*64)
        cdef uint64_t **maybe_newfaces = maybe_newfaces_mem.data

        # ``newfaces`` point to the actual facets of ``faces[n_faces -1]``.
        cdef MemoryAllocator newfaces_mem = MemoryAllocator()
        cdef uint64_t **newfaces = <uint64_t **> newfaces_mem.allocarray(n_faces, sizeof(uint64_t *))

        # Calculating ``maybe_newfaces`` and ``newfaces``
        # such that ``newfaces`` points to all facets of ``faces[n_faces -1]``.
        cdef size_t new_n_faces
        sig_on()
        new_n_faces = get_next_level(faces, n_faces, maybe_newfaces,
                                      newfaces, NULL, 0, face_length)
        sig_off()

        # compute the dimension of the polyhedron,
        # by calculating dimension of one of its faces.
        return self.compute_dimension_loop(newfaces, new_n_faces, face_length) + 1
