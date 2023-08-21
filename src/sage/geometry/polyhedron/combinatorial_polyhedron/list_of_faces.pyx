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

    sage: # needs sage.combinat
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
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix_dense  cimport Matrix_dense

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
    def __cinit__(self, size_t n_faces, size_t n_atoms, size_t n_coatoms):
        r"""
        Initialize :class:`ListOfFaces`.

        See :class:`ListOfFaces`.

        TESTS::

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces).run()
        """
        # Note that all values are set to zero at the time ``__cinit__`` is called:
        # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#initialisation-methods
        # In particular, ``__dealloc__`` will not do harm in this case.

        face_list_init(self.data, n_faces, n_atoms, n_coatoms)

    def __dealloc__(self):
        r"""
        TESTS:

        Check that failed ``cinit`` does not give segmentation fault or similar::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces import ListOfFaces
            sage: ListOfFaces(-1, -1, -1)  # indirect doctest
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to size_t

            sage: from memory_allocator.test import TestMemoryAllocator
            sage: t = TestMemoryAllocator()
            sage: m = t.size_t_max()
            sage: ListOfFaces(10, m, 10)
            Traceback (most recent call last):
            ...
            MemoryError: failed to allocate ...
        """
        face_list_free(self.data)

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

    cpdef ListOfFaces __copy__(self):
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

    cdef ListOfFaces delete_atoms_unsafe(self, bint *delete, face_t face):
        r"""
        Return a copy of ``self`` where bits in ``delete`` have been
        removed/contracted.

        The the remaining bits will be shifted to the left.

        If ``delete`` is ``NULL``, keep exactly the bits set in ``face``.

        .. WARNING::

            ``delete`` is assumed to be of length ``self.n_atoms`` or NULL.
            ``face`` is assumed to be of length ``self.face_length`` if ``delete`` is not ``NULL``.
        """

        cdef output_n_atoms
        cdef size_t i, j
        if delete is NULL:
            output_n_atoms = face_len_atoms(face)
        else:
            output_n_atoms = self.n_atoms()
            for i in range(self.n_atoms()):
                if delete[i]:
                    output_n_atoms -= 1

        cdef ListOfFaces output = ListOfFaces(self.n_faces(), output_n_atoms, self.n_coatoms())
        cdef size_t counter = 0
        cdef size_t n_atoms = self.n_atoms()
        cdef size_t n_faces = self.n_faces()
        for i in range(n_atoms):
            if ((delete is NULL and face_atom_in(face, i)) or
                    (delete is not NULL and not delete[i])):
                # The atom will be kept.
                for j in range(n_faces):
                    if face_atom_in(self.data.faces[j], i):
                        face_add_atom(output.data.faces[j], counter)
                counter += 1

        return output

    cdef void delete_faces_unsafe(self, bint *delete, face_t face):
        r"""
        Deletes face ``i`` if and only if ``delete[i]``.

        Alternatively, deletes all faces such that the ``i``-th bit in ``face`` is not set.

        This will modify ``self``.

        .. WARNING::

            ``delete`` is assumed to be of length ``self.n_faces()`` or NULL.
            ``face`` is assumed to contain ``self.n_faces()`` atoms if ``delete`` is not ``NULL``.
        """
        if delete is not NULL:
            face_list_delete_faces_by_array(self.data, delete)
        else:
            face_list_delete_faces_by_face(self.data, face)

    cdef void get_not_inclusion_maximal_unsafe(self, bint *not_inclusion_maximal):
        r"""
        Get all faces that are not inclusion maximal.

        Set ``not_inclusion_maximal[i]`` to one if ``self.data[i]`` is not
        an inclusion-maximal face, otherwise to zero.

        If there are duplicates, all but the last duplicate will be marked as
        not inclusion maximal.

        .. WARNING::

            ``not_inclusion_maximal`` is assumed to be at least of length ``self.n_atoms`` or NULL.
        """
        cdef size_t i
        memset(not_inclusion_maximal, 0, sizeof(bint)*self.n_faces())
        for i in range(self.n_faces()):
            not_inclusion_maximal[i] = is_not_maximal_fused(self.data, i, <standard> 0, not_inclusion_maximal)

    cdef void get_faces_all_set_unsafe(self, bint *all_set):
        r"""
        Get the faces that have all ``bits`` set.

        Set ``all_set[i]`` to one if ``self.data[i]``
        has all bits set, otherwise to zero.

        .. WARNING::

            ``all_set`` is assumed to be at least of length ``self.n_atoms`` or NULL.
        """
        cdef size_t i
        for i in range(self.n_faces()):
            if face_first_missing_atom(self.data.faces[i]) == -1:
                all_set[i] = 1
            else:
                all_set[i] = 0

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
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_dense M = matrix(
                ZZ, self.n_faces(), self.n_atoms(), 0)

        cdef size_t i
        cdef long j
        for i in range(self.n_faces()):
            j = face_next_atom(self.data.faces[i], 0)
            while j != -1:
                M.set_unsafe_int(i, j, 1)
                j = face_next_atom(self.data.faces[i], j+1)

        M.set_immutable()
        return M

cdef tuple face_as_combinatorial_polyhedron(ListOfFaces facets, ListOfFaces Vrep, face_t face, bint dual):
    r"""
    Obtain facets and Vrepresentation of ``face`` as new combinatorial polyhedron.

    INPUT:

    - ``facets`` -- facets of the polyhedron
    - ``Vrep`` -- Vrepresentation of the polyhedron
    - ``face`` -- face in Vrepresentation or ``NULL``
    - ``dual`` -- boolean

    OUTPUT: A tuple of new facets and new Vrepresentation as :class:`ListOfFaces`.
    """
    cdef ListOfFaces new_facets, new_Vrep
    cdef bint* delete
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef size_t i
    cdef face_t null_face

    # Delete all atoms not in the face.
    if not dual:
        new_facets = facets.delete_atoms_unsafe(NULL, face)
        new_Vrep = Vrep.__copy__()
        new_Vrep.delete_faces_unsafe(NULL, face)

        delete = <bint*> mem.allocarray(new_facets.n_faces(), sizeof(bint))
    else:
        delete = <bint*> mem.allocarray(max(facets.n_faces(), facets.n_atoms()), sizeof(bint))

        # Set ``delete[i]`` to one if ``i`` is not an vertex of ``face``.
        for i in range(Vrep.n_faces()):
            if face_issubset(face, Vrep.data.faces[i]):
                delete[i] = 0
            else:
                delete[i] = 1

        new_facets = facets.delete_atoms_unsafe(delete, null_face)
        new_Vrep = Vrep.__copy__()
        new_Vrep.delete_faces_unsafe(delete, null_face)

    # Delete all facets that define the face.
    new_facets.get_faces_all_set_unsafe(delete)
    new_facets.delete_faces_unsafe(delete, null_face)
    new_Vrep = new_Vrep.delete_atoms_unsafe(delete, null_face)

    # Now delete all facets that are not inclusion maximal.
    # the last copy of each duplicate will remain.
    new_facets.get_not_inclusion_maximal_unsafe(delete)
    new_facets.delete_faces_unsafe(delete, null_face)
    new_Vrep = new_Vrep.delete_atoms_unsafe(delete, null_face)

    # Finally set coatoms of the output.
    for i in range(new_facets.n_faces()):
        facet_set_coatom(new_facets.data.faces[i], i)
    for i in range(new_Vrep.n_faces()):
        facet_set_coatom(new_Vrep.data.faces[i], i)

    return (new_facets, new_Vrep)
