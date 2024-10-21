# sage_setup: distribution = sagemath-polyhedra
r"""
PolyhedronFaceLattice

This module provides a class that stores and sorts all faces of the polyhedron.

:class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron` implicitly uses this class to generate
the face lattice of a polyhedron.

Terminology in this module:

- Vrep -- ``[vertices, rays, lines]`` of the polyhedron

- Hrep -- inequalities and equations of the polyhedron

- Facets -- facets of the polyhedron

- Coatoms -- the faces from which all others are constructed in the face
  iterator. This will be facets or Vrep.  In non-dual mode, faces are
  constructed as intersections of the facets. In dual mode, the are constructed
  theoretically as joins of vertices.  The coatoms are repsented as incidences
  with the atoms they contain.

- Atoms -- facets or Vrep depending on application of algorithm.  Atoms are
  repsented as incidences of coatoms they are contained in.

- Vrepresentation -- represents a face by a list of Vrep it contains

- Hrepresentation -- represents a face by a list of Hrep it is contained in

- bit representation -- represents incidences as ``uint64_t``-array, where each
  bit represents one incidence. There might be trailing zeros, to fit alignment
  requirements.  In most instances, faces are represented by the bit
  representation, where each bit corresponds to an atom.

EXAMPLES::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice \
    ....: import PolyhedronFaceLattice
    sage: P = polytopes.octahedron()
    sage: C = CombinatorialPolyhedron(P)
    sage: all_faces = PolyhedronFaceLattice(C)

.. SEEALSO::

    :mod:`~sage.geometry.polyhedron.combinatorial_polyhedron.base`,
    :class:`PolyhedronFaceLattice`.

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

from .conversions import (facets_tuple_to_bit_rep_of_facets,
                          facets_tuple_to_bit_rep_of_Vrep)

from sage.geometry.polyhedron.combinatorial_polyhedron.conversions cimport bit_rep_to_Vrep_list

from sage.geometry.polyhedron.combinatorial_polyhedron.base                     cimport CombinatorialPolyhedron
from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator            cimport FaceIterator
from sage.geometry.polyhedron.combinatorial_polyhedron.face_list_data_structure cimport *


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class PolyhedronFaceLattice:
    r"""
    A class to generate incidences of :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    On initialization all faces of the given :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`
    are added and sorted (except coatoms). The incidences can be used to
    generate the ``face_lattice``.

    Might generate the faces of the dual polyhedron for speed.

    INPUT:

    - :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.baseCombinatorialPolyhedron`

    .. SEEALSO::

        :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron._record_all_faces`,
        :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron._record_all_faces_helper`,
        :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron.face_lattice`,
        :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron._compute_face_lattice_incidences`.

    EXAMPLES::

        sage: P = polytopes.Birkhoff_polytope(3)
        sage: C = CombinatorialPolyhedron(P)
        sage: C._record_all_faces()  # indirect doctests
        sage: C.face_lattice()                                                          # needs sage.combinat
        Finite lattice containing 50 elements

    ALGORITHM:

    The faces are recorded with :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator` in Bit-representation.
    Once created, all level-sets but the coatoms are sorted with merge sort.
    Non-trivial incidences of elements whose rank differs by 1 are determined
    by intersecting with all coatoms. Then each intersection is looked up in
    the sorted level sets.
    """
    def __cinit__(self, CombinatorialPolyhedron C):
        r"""
        Initialize :class:`PolyhedronFaceLattice`.

        See :class:`PolyhedronFaceLattice`.

        TESTS:

        Not initializing the class, does not give segmentation fault::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice import PolyhedronFaceLattice
            sage: P = polytopes.cube()
            sage: F = PolyhedronFaceLattice.__new__(PolyhedronFaceLattice, P.combinatorial_polyhedron())
            sage: F.get_face(2, 3)
            A 2-dimensional face of a 3-dimensional combinatorial polyhedron
        """
        # Note that all values are set to zero at the time ``__cinit__`` is called:
        # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#initialisation-methods
        # In particular, ``__dealloc__`` will not do harm in this case.

        cdef int i
        cdef size_t j

        self.dimension = C.dimension()
        self.dual = False
        if C.bitrep_facets().n_faces() > C.bitrep_Vrep().n_faces():
            self.dual = True
        if not C.is_bounded():
            self.dual = False
        cdef FaceIterator face_iter = C._face_iter(self.dual, -2)
        self._Vrep = C.Vrep()
        self._facet_names = C.facet_names()
        self._equations = C.equations()
        self._bounded = C.is_bounded()

        # copy f_vector for later use
        f_vector = C.f_vector()
        self.f_vector = <size_t *> check_allocarray(self.dimension + 2, sizeof(size_t))
        if self.dual:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[-i-2]
        else:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[i+1]

        # Initialize atoms, coatoms, ``atom_rep`` and ``coatom_rep``.
        if self.dimension == 0:
            # In case of the 0-dimensional polyhedron, we have to fix atoms and coatoms.
            # So far this didn't matter, as we only iterated over proper faces.
            self.atoms = facets_tuple_to_bit_rep_of_Vrep(((),), 1)
            self.coatoms = facets_tuple_to_bit_rep_of_facets(((),), 1)
        else:
            self.atoms = face_iter.atoms
            self.coatoms = face_iter.coatoms

        self.atom_rep = <size_t *> check_allocarray(self.coatoms.n_atoms(), sizeof(size_t))
        self.coatom_rep = <size_t *> check_allocarray(self.coatoms.n_faces(), sizeof(size_t))

        # Setting up a pointer to raw data of ``faces``:
        self.faces = <face_list_t*> check_calloc(self.dimension + 2, sizeof(face_list_t))

        for i in range(self.dimension + 2):
            face_list_init(self.faces[i], self.f_vector[i],
                           self.coatoms.n_atoms(), self.coatoms.n_coatoms())

        # The universe.
        for j in range(self.coatoms.n_atoms()):
            face_add_atom(self.faces[self.dimension+1].faces[0], j)

        # The coatoms.
        if self.dimension > 0:
            # Note that in the other cases, this was fully initialized above.
            face_list_copy(self.faces[self.dimension], self.coatoms.data)

        # Attributes for iterating over the incidences.
        face_init(self.incidence_face, self.coatoms.n_atoms(), self.coatoms.n_coatoms())

        # Adding all faces, using the iterator.
        for i in range(1, self.dimension):
            self.faces[i].n_faces = 0

        cdef int d
        if face_iter.structure.current_dimension != self.dimension:
            # If there are proper faces.
            d = face_iter.next_dimension()
            while (d == self.dimension - 1):
                # We already have the coatoms.
                d = face_iter.next_dimension()
            while (d < self.dimension):
                add_face_deep(self.faces[d+1], face_iter.structure.face)
                d = face_iter.next_dimension()

    def __init__(self, CombinatorialPolyhedron C):
        r"""
        Initialize :class:`PolyhedronFaceLattice`.

        See :class:`PolyhedronFaceLattice`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces() # indirect doctests
            sage: C.face_lattice()                                                      # needs sage.combinat
            Finite lattice containing 28 elements

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice.PolyhedronFaceLattice).run()
        """
        # Sorting the faces, except for coatoms.
        self._sort()

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice import PolyhedronFaceLattice
            sage: PolyhedronFaceLattice()  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: __cinit__() takes exactly 1 positional argument (0 given)
        """
        cdef int i
        sig_free(self.f_vector)
        sig_free(self.atom_rep)
        sig_free(self.coatom_rep)
        if self.faces:
            for i in range(self.dimension + 2):
                face_list_free(self.faces[i])
            sig_free(self.faces)
        face_free(self.incidence_face)

    cdef int _sort(self) except -1:
        r"""
        Sort each list of ``self.faces`` (except for coatoms).

        This method is used on initialization only.
        """
        cdef int dim = self.dimension
        cdef int i
        for i in range(dim + 2):
            if unlikely(self.f_vector[i] != self.faces[i].n_faces):
                raise ValueError("``PolyhedronFaceLattice`` does not contain all faces")

        for i in range(dim - 1):
            # Sort each level set, except for the facets, the full- and empty polyhedron.
            sort_faces_list(self.faces[i+1])

    def _find_face_from_combinatorial_face(self, CombinatorialFace face):
        r"""
        A method to test :meth:`find_face`.

        ``f`` must be a face in dual mode if and only if ``self`` is in dual mode.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice \
            ....:         import PolyhedronFaceLattice
            sage: P = polytopes.hypercube(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: F = PolyhedronFaceLattice(C)
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: F._find_face_from_combinatorial_face(face)
            Traceback (most recent call last):
            ...
            ValueError: cannot find a facet, as those are not sorted

        """
        if not (self.dual == face._dual):
            raise ValueError("iterator and allfaces not in same mode")
        cdef size_t face_index = self.find_face(face.dimension(), face.face)
        if face_index == -1:
            raise ValueError("face is not in the face lattice")
        return face_index

    cdef inline size_t find_face(self, int dimension, face_t face) except -2:
        r"""
        Return the index of ``face``, if it is of dimension ``dimension``.

        Return -1 if the face is not contained.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice \
            ....:         import PolyhedronFaceLattice
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: F = PolyhedronFaceLattice(C)
            sage: it = C.face_generator(dimension=1)
            sage: S = set(F._find_face_from_combinatorial_face(f) for f in it)
            sage: S == set(range(36))
            True
        """
        if unlikely(dimension == self.dimension -1):
            raise ValueError("cannot find a facet, as those are not sorted")
            # of course one can easily add a function to search for a facet as
            # well, but there seems to be no need for that

        if unlikely(dimension < -1 or dimension > self.dimension):
            raise IndexError("dimension out of range")

        return find_face(face, self.faces[dimension+1])

    cpdef CombinatorialFace get_face(self, int dimension, size_t index):
        r"""
        Return the face of dimension ``dimension`` and index ``index``.

        INPUT:

        - ``dimension`` -- dimension of the face
        - ``index`` -- index of the face
        - ``names`` -- if ``True`` returns the names of the ``[vertices, rays, lines]``
          as given on initialization of :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_face_lattice \
            ....:         import PolyhedronFaceLattice
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: F = PolyhedronFaceLattice(C)
            sage: it = C.face_generator(dimension=1)
            sage: face = next(it)
            sage: index = F._find_face_from_combinatorial_face(face)
            sage: F.get_face(face.dimension(), index).ambient_Vrepresentation()
            (A vertex at (2, 1, 4, 3), A vertex at (1, 2, 4, 3))
            sage: face.ambient_Vrepresentation()
            (A vertex at (2, 1, 4, 3), A vertex at (1, 2, 4, 3))
            sage: all(F.get_face(face.dimension(),
            ....:                F._find_face_from_combinatorial_face(face)).ambient_Vrepresentation() ==
            ....:     face.ambient_Vrepresentation() for face in it)
            True

            sage: P = polytopes.twenty_four_cell()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = PolyhedronFaceLattice(C)
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: while (face.dimension() == 3): face = next(it)
            sage: index = F._find_face_from_combinatorial_face(face)
            sage: F.get_face(face.dimension(), index).ambient_Vrepresentation()
            (A vertex at (-1/2, 1/2, -1/2, -1/2),
             A vertex at (-1/2, 1/2, 1/2, -1/2),
             A vertex at (0, 0, 0, -1))
            sage: all(F.get_face(face.dimension(),
            ....:                F._find_face_from_combinatorial_face(face)).ambient_V_indices() ==
            ....:     face.ambient_V_indices() for face in it)
            True
        """
        if self.dual:
            # if dual, the Vrepresentation corresponds to the coatom-representation
            dimension = self.dimension - 1 - dimension  # if dual, the dimensions are reversed
        return CombinatorialFace(self, dimension=dimension, index=index)

    cdef size_t set_coatom_rep(self, int dimension, size_t index) except -1:
        r"""
        Set ``atom_rep`` to be the atom-representation of the face
        of dimension ``dimension`` and index ``index``.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_coatom_rep`
        """
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s" % dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s" % (index, dimension))
        if unlikely(self.coatoms.n_faces() == 0):
            return 0

        cdef face_t face = self.faces[dimension+1].faces[index]
        return bit_rep_to_coatom_rep(face, self.coatoms.data, self.coatom_rep)

    cdef size_t set_atom_rep(self, int dimension, size_t index) except -1:
        r"""
        Set ``atom_rep`` to be the atom-representation of the face
        of dimension ``dimension`` and index ``index``.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_atom_rep`
        """
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s" % dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s" % (index, dimension))

        cdef face_t face = self.faces[dimension+1].faces[index]
        return bit_rep_to_Vrep_list(face, self.atom_rep)

    cdef void incidence_init(self, int dimension_one, int dimension_two) noexcept:
        r"""
        Initialize the :class:`PolyhedronFaceLattice` to give incidences between
        ``dimension_one`` and ``dimension_two``.

        This will enable :meth:`next_incidence` to give all such incidences.

        Currently only ``dimension_one == dimension_two + 1`` and incidences
        with empty and full polyhedron are implemented, which suffices for the
        face-lattice.
        """
        if dimension_one == self.dimension:
            # The full polyhedron is incident to every face.
            if dimension_two < -1:
                raise ValueError("no faces of dimension %s" % dimension_two)
            if dimension_two > self.dimension:
                raise ValueError("no faces of dimension %s" % dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.is_incidence_initialized = 2
            return

        if dimension_two == -1:
            # The empty polyhedron is incident to every face.
            if dimension_one < -1:
                raise ValueError("no faces of dimension %s" % dimension_two)
            if dimension_one > self.dimension:
                raise ValueError("no faces of dimension %s" % dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.is_incidence_initialized = 3
            return

        if dimension_one != dimension_two + 1:
            raise ValueError("``dimension_one = dimension_two`` + 1 must hold")
            # At the moment, this is not implemented.

        if dimension_one > self.dimension:
            raise ValueError("no faces of dimension %s" % dimension_one)
        if dimension_two < -1:
            raise ValueError("no faces of dimension %s" % dimension_two)

        self.incidence_dim_one = dimension_one
        self.incidence_dim_two = dimension_two
        self.incidence_counter_one = 0
        self.incidence_counter_two = 0
        self.is_incidence_initialized = 1

    cdef inline bint next_incidence(self, size_t *one, size_t *two) noexcept:
        r"""
        Set ``one[0]`` and ``two[0]`` to be the next incidence. Return ``True``
        unless there are no more incidences, then return `0`.

        After initialization with :meth:`next_incidence`, this method will give
        all incidences of faces of ``dimension_one`` and ``dimension_two``.
        ``one[0]`` will represent the index of a face in ``dimension_one`` and
        ``two[0]`` will represent the index of a face in ``dimension_two``
        according to their order in :class:`PolyhedronFaceLattice`.

        Use :meth:`Vrep` and :meth:`Hrep` to interpret the output.

        ALGORITHM:

        This is the algorithm for non-trivial cases:
        ``0 < self.dimension_two + 1 == self.dimension_one < self.dimension``

        We intersect each face of dimension ``incidence_dim_one`` with each
        coatom. The result will be looked up in the ``incidence_dim_two`` faces.
        """
        cdef bint result = False
        while ((not result)
                and (self.incidence_counter_one < self.f_vector[self.incidence_dim_one + 1])):
            # Calls next_incidence_loop, until it gives a result or
            # until there are no more incidences.
            result = self.next_incidence_loop(one, two)

        return result

    cdef inline bint next_incidence_loop(self, size_t *one, size_t *two) noexcept:
        r"""
        Set ``one[0]`` and ``two[0]`` to be the next incidence. Return ``True``
        on success and ``False`` otherwise.

        If it returns ``False``, it needs to be called again, unless
        ``self.incidence_counter_one >= self.f_vector[self.incidence_dim_one + 1])``.

        See :meth:`next_incidence`.
        """
        cdef face_list_t coatoms
        coatoms[0] = self.coatoms.data[0]
        cdef face_t dimension_one_face  # depending on the index ``incidence_counter_one``

        cdef size_t location  # the index the intersection has, if of correct dimension
        if self.is_incidence_initialized == 1:
            # The standard case, where
            # ``0 < self.dimension_two + 1 == self.dimension_one < self.dimension``.

            one[0] = self.incidence_counter_one
            dimension_one_face = self.faces[self.incidence_dim_one + 1].faces[self.incidence_counter_one]

            # Get the intersection of ``dimension_one_face`` with the
            # ``self.incidence_counter_two``-th coatom.
            face_intersection(self.incidence_face, dimension_one_face,
                              coatoms.faces[self.incidence_counter_two])

            # Get the location of the intersection and
            # check whether it is correct.
            location = self.find_face(self.incidence_dim_two, self.incidence_face)
            two[0] = location

            # Set counters for next function call.
            self.incidence_counter_two += 1
            if self.incidence_counter_two == self.f_vector[self.dimension]:
                self.incidence_counter_one += 1
                self.incidence_counter_two = 0
            return location != -1

        if self.is_incidence_initialized == 2:
            # the case where ``dimension_one`` is dimension of polyhedron.
            return self.next_trivial_incidence(one, two)

        if self.is_incidence_initialized == 3:
            # the case where ``dimension_two`` is `-1`.
            return self.next_trivial_incidence2(one, two)

        if self.is_incidence_initialized == 0:
            return 0

    cdef inline bint next_trivial_incidence(self, size_t *one, size_t *two) noexcept:
        r"""
        Handling the case where ``dimension_one`` is dimension of polyhedron.

        See :meth:`next_incidence`.
        """
        one[0] = 0
        two[0] = self.incidence_counter_two
        self.incidence_counter_two += 1

        # Once done, raising ``self.incidence_counter_one``, such that
        # :meth:`next_incidence` recognizes that we are done.
        if self.incidence_counter_two >= self.f_vector[self.incidence_dim_two + 1]:
            self.incidence_counter_one += 1

        return (two[0] < self.f_vector[self.incidence_dim_two + 1])

    cdef inline bint next_trivial_incidence2(self, size_t *one, size_t *two) noexcept:
        r"""
        Handling the case where ``dimension_two`` is `-1`.

        See :meth:`next_incidence`.
        """
        two[0] = 0
        one[0] = self.incidence_counter_one
        self.incidence_counter_one += 1
        return (one[0] < self.f_vector[self.incidence_dim_one + 1])
