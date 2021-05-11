r"""
List of faces

This module provides a class to store faces of a polyhedron in Bit-representation.

This class allocates memory to store the faces in.
A face will be stored as vertex-incidences, where each Bit represents an incidence.
In :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions` there a methods to actually convert facets of a polyhedron
to bit-representations of vertices stored in :class:`ListOfFaces`.

Moreover, :class:`ListOfFaces` calculates the dimension of a polyhedron, assuming the
faces are the facets of this polyhedron.

Each face is stored over-aligned according to the ``chunktype``.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Provide enough space to store `20` faces as incidences to `60` vertices::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
    ....: import ListOfFaces
    sage: face_list = ListOfFaces(20, 60, 20)
    sage: face_list.matrix().is_zero()
    True

Obtain the facets of a polyhedron::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_rep_of_facets
    sage: P = polytopes.cube()
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(P.incidence_matrix())
    sage: face_list = incidence_matrix_to_bit_rep_of_facets(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the Vrepresentation of a polyhedron as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import incidence_matrix_to_bit_rep_of_Vrep
    sage: P = polytopes.associahedron(['A',3])
    sage: face_list = incidence_matrix_to_bit_rep_of_Vrep(P.incidence_matrix())
    sage: face_list.compute_dimension()
    3

Obtain the facets of a polyhedron as :class:`ListOfFaces` from a facet list::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_facets
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)

Likewise for the Vrepresentatives as facet-incidences::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
    ....:         import facets_tuple_to_bit_rep_of_Vrep
    sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
    sage: face_list = facets_tuple_to_bit_rep_of_Vrep(facets, 4)

Obtain the matrix of a list of faces::

    sage: face_list.matrix()
    [1 1 1 0]
    [1 1 0 1]
    [1 0 1 1]
    [0 1 1 1]

.. SEEALSO::

    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.base`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.conversions`,
    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_faces_lattice`.

AUTHOR:

- Jonathan Kliem (2019-04)
"""

# ****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import is_Matrix

from sage.matrix.matrix_integer_dense  cimport Matrix_integer_dense

from .face_list_data_structure cimport *

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class ListOfFaces:
    r"""
    A class to store the Bit-representation of faces in.

    This class will allocate the memory for the faces.

    INPUT:

    - ``n_faces`` -- the number of faces to be stored
    - ``n_atoms`` -- the total number of atoms the faces contain
    - ``n_coatoms`` -- the total number of coatoms of the polyhedron

    .. SEEALSO::

        :meth:`incidence_matrix_to_bit_rep_of_facets`,
        :meth:`incidence_matrix_to_bit_rep_of_Vrep`,
        :meth:`facets_tuple_to_bit_rep_of_facets`,
        :meth:`facets_tuple_to_bit_rep_of_Vrep`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
        ....:     import ListOfFaces
        sage: facets = ListOfFaces(5, 13, 5)
        sage: facets.matrix().dimensions()
        (5, 13)
    """
    def __init__(self, size_t n_faces, size_t n_atoms, size_t n_coatoms):
        r"""
        Initialize :class:`ListOfFaces`.

        See :class:`ListOfFaces`.

        TESTS::

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces).run()
        """
        self._mem = MemoryAllocator()
        face_list_init(self.data, n_faces, n_atoms, n_coatoms, self._mem)

    def _test_alignment(self):
        r"""
        Check the correct alignment.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:         import ListOfFaces
            sage: facets = ListOfFaces(10, 13, 10)
            sage: facets._test_alignment()
        """
        assert face_list_check_alignment(self.data)

    def __copy__(self):
        r"""
        Return a copy of self.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: facets = ListOfFaces(5, 13, 5)
            sage: copy(facets).matrix().dimensions()
            (5, 13)

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: facets.compute_dimension()
            3
            sage: copy(facets).compute_dimension()
            3
            sage: facets.matrix() == copy(facets).matrix()
            True
            sage: copy(facets) is facets
            False
        """
        cdef ListOfFaces copy = ListOfFaces(self.n_faces(), self.n_atoms(), self.n_coatoms())
        face_list_copy(copy.data, self.data)
        return copy

    cpdef int compute_dimension(self) except -2:
        r"""
        Compute the dimension of a polyhedron by its facets.

        This assumes that ``self`` is the list of facets of a polyhedron.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets, \
            ....:            facets_tuple_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: Vrep = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
            sage: facets.compute_dimension()
            3
            sage: Vrep.compute_dimension()
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
            ....:     import incidence_matrix_to_bit_rep_of_facets, \
            ....:            incidence_matrix_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: for _ in range(10):
            ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
            ....:                    for _ in range(randint(3,15)))
            ....:     P = Polyhedron(vertices=points)
            ....:     inc = P.incidence_matrix()
            ....:     mod_inc = inc.delete_columns([i for i,V in enumerate(P.Hrepresentation()) if V.is_equation()])
            ....:     facets = incidence_matrix_to_bit_rep_of_facets(mod_inc)
            ....:     vertices = incidence_matrix_to_bit_rep_of_Vrep(mod_inc)
            ....:     d1 = P.dimension()
            ....:     if d1 == 0:
            ....:         continue
            ....:     d2 = facets.compute_dimension()
            ....:     d3 = vertices.compute_dimension()
            ....:     if not d1 == d2 == d3:
            ....:         print('calculation_dimension() seems to be incorrect')
        """
        if self.n_faces() == 0:
            raise TypeError("at least one face needed")

        cdef size_t n_faces = self.n_faces()

        cdef face_list_t empty_forbidden
        empty_forbidden.n_faces = 0

        if n_faces == 1:
            # We expect the face to be the empty polyhedron.
            # Possibly it contains more than one vertex/rays/lines.
            # The dimension of a polyhedron with this face as only facet is
            # the number of atoms it contains.
            return face_len_atoms(self.data.faces[0])

        # ``newfaces`` are all intersection of ``faces[n_faces -1]`` with previous faces.
        # It needs to be allocated to store those faces.
        cdef ListOfFaces new_faces = ListOfFaces(self.n_faces(), self.n_atoms(), self.n_coatoms())

        # Calculating ``newfaces``
        # such that ``newfaces`` points to all facets of ``faces[n_faces -1]``.
        cdef size_t new_n_faces = get_next_level(self.data, new_faces.data, empty_forbidden)

        # Undo what ``get_next_level`` does.
        self.data.n_faces += 1

        # compute the dimension of the polyhedron,
        # by calculating dimension of one of its faces.
        return new_faces.compute_dimension() + 1

    cpdef ListOfFaces pyramid(self):
        r"""
        Return the list of faces of the pyramid.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:         import facets_tuple_to_bit_rep_of_facets
            sage: facets = ((0,1,2), (0,1,3), (0,2,3), (1,2,3))
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [1 1 1 0]
            [1 1 0 1]
            [1 0 1 1]
            [0 1 1 1]
            sage: face_list.pyramid().matrix()
            [1 1 1 0 1]
            [1 1 0 1 1]
            [1 0 1 1 1]
            [0 1 1 1 1]
            [1 1 1 1 0]

        Incorrect facets that illustrate how this method works::

            sage: facets = ((0,1,2,3), (0,1,2,3), (0,1,2,3), (0,1,2,3))
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            [1 1 1 1]
            sage: face_list.pyramid().matrix()
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 1]
            [1 1 1 1 0]

        ::

            sage: facets = ((), (), (), ())
            sage: face_list = facets_tuple_to_bit_rep_of_facets(facets, 4)
            sage: face_list.matrix()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: face_list.pyramid().matrix()
            [0 0 0 0 1]
            [0 0 0 0 1]
            [0 0 0 0 1]
            [0 0 0 0 1]
            [1 1 1 1 0]
        """
        cdef ListOfFaces copy
        cdef size_t i
        cdef size_t n_faces = self.n_faces()
        cdef size_t n_atoms = self.n_atoms()

        # ``copy`` has a new atom and a new coatom.
        copy = ListOfFaces(n_faces + 1, n_atoms + 1, n_faces + 1)

        # Note that a pyramid is simple if and only if the base is simple.

        face_list_copy(copy.data, self.data)

        for i in range(n_faces):
            face_add_atom(copy.data.faces[i], n_atoms)
            facet_set_coatom(copy.data.faces[i], i)

        copy.data.n_faces += 1

        # The new coatom contains all atoms, but the new atom.
        face_set_first_n_atoms(copy.data.faces[n_faces], n_atoms)
        facet_set_coatom(copy.data.faces[n_faces], n_faces)

        return copy

    def matrix(self):
        r"""
        Obtain the matrix of self.

        Each row represents a face and each column an atom.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_rep_of_facets, \
            ....:     facets_tuple_to_bit_rep_of_Vrep
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4), (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
            sage: Vrep = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
            sage: facets.matrix()
            [1 1 0 0 1 0]
            [0 1 1 0 1 0]
            [0 0 1 1 1 0]
            [1 0 0 1 1 0]
            [1 1 0 0 0 1]
            [0 1 1 0 0 1]
            [0 0 1 1 0 1]
            [1 0 0 1 0 1]
            sage: facets.matrix().transpose() == Vrep.matrix()
            True
        """
        from sage.rings.all import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_integer_dense M = matrix(
                ZZ, self.n_faces(), self.n_atoms(), 0)

        cdef size_t i
        cdef long j
        for i in range(self.n_faces()):
            j = face_next_atom(self.data.faces[i], 0)
            while j != -1:
                M.set_unsafe_si(i, j, 1)
                j = face_next_atom(self.data.faces[i], j+1)

        M.set_immutable()
        return M
