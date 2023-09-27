# distutils: extra_compile_args = OPENMP_CFLAGS
# distutils: extra_link_args = OPENMP_CFLAGS
r"""
Face iterator for polyhedra

This iterator in principle works on every graded lattice, where
every interval of length two has exactly 4 elements (diamond property).

It also works on unbounded polyhedra, as those satisfy the diamond property,
except for intervals including the empty face.
A (slightly generalized) description of the algorithm can be found in [KS2019]_.

Terminology in this module:

- Coatoms -- the faces from which all others are constructed in the
  face iterator. This will be facets or Vrep.  In non-dual mode, faces
  are constructed as intersections of the facets. In dual mode, they
  are constructed theoretically as joins of vertices.  The coatoms are
  repsented as incidences with the atoms they contain.

- Atoms -- facets or Vrep depending on application of algorithm.  Atoms are
  represented as incidences of coatoms they are contained in.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`.

EXAMPLES:

Construct a face iterator::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator \
    ....:         import FaceIterator
    sage: P = polytopes.octahedron()
    sage: C = CombinatorialPolyhedron(P)

    sage: FaceIterator(C, False)
    Iterator over the proper faces of a 3-dimensional combinatorial polyhedron
    sage: FaceIterator(C, False, output_dimension=2)
    Iterator over the 2-faces of a 3-dimensional combinatorial polyhedron

Iterator in the non-dual mode starts with facets::

    sage: it = FaceIterator(C, False)
    sage: [next(it) for _ in range(9)]
    [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 1-dimensional face of a 3-dimensional combinatorial polyhedron]

Iterator in the dual-mode starts with vertices::

    sage: it = FaceIterator(C, True)
    sage: [next(it) for _ in range(7)]
    [A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
     A 1-dimensional face of a 3-dimensional combinatorial polyhedron]

Obtain the Vrepresentation::

    sage: it = FaceIterator(C, False)
    sage: face = next(it)
    sage: face.ambient_Vrepresentation()
    (A vertex at (0, -1, 0), A vertex at (0, 0, -1), A vertex at (1, 0, 0))
    sage: face.n_ambient_Vrepresentation()
    3

Obtain the facet-representation::

    sage: it = FaceIterator(C, True)
    sage: face = next(it)
    sage: face.ambient_Hrepresentation()
    (An inequality (-1, -1, 1) x + 1 >= 0,
     An inequality (-1, -1, -1) x + 1 >= 0,
      An inequality (-1, 1, -1) x + 1 >= 0,
       An inequality (-1, 1, 1) x + 1 >= 0)
    sage: face.ambient_H_indices()
    (4, 5, 6, 7)
    sage: face.n_ambient_Hrepresentation()
    4

In non-dual mode one can ignore all faces contained in the current face::

    sage: it = FaceIterator(C, False)
    sage: face = next(it)
    sage: face.ambient_H_indices()
    (7,)
    sage: it.ignore_subfaces()
    sage: [face.ambient_H_indices() for face in it]
    [(6,),
    (5,),
    (4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (5, 6),
    (1, 6),
    (0, 1, 5, 6),
    (4, 5),
    (0, 5),
    (0, 3, 4, 5),
    (3, 4),
    (2, 3),
    (0, 3),
    (0, 1, 2, 3),
    (1, 2),
    (0, 1)]

In dual mode one can ignore all faces that contain the current face::

    sage: it = FaceIterator(C, True)
    sage: face = next(it)
    sage: face.ambient_V_indices()
    (5,)
    sage: it.ignore_supfaces()
    sage: [face.ambient_V_indices() for face in it]
    [(4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (3, 4),
    (2, 4),
    (0, 4),
    (0, 3, 4),
    (0, 2, 4),
    (1, 3),
    (0, 3),
    (0, 1, 3),
    (1, 2),
    (0, 2),
    (0, 1, 2),
    (0, 1)]

There is a special face iterator class for geometric polyhedra.
It yields (geometric) polyhedral faces and it also yields trivial faces.
Otherwise, it works exactly the same::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator \
    ....:         import FaceIterator_geom
    sage: P = polytopes.cube()
    sage: it = FaceIterator_geom(P)
    sage: [next(it) for _ in range(5)]
    [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices,
     A -1-dimensional face of a Polyhedron in ZZ^3,
     A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
     A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
     A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices]
    sage: it
    Iterator over the faces of a 3-dimensional polyhedron in ZZ^3

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

from cython.parallel cimport prange, threadid
from cysignals.memory cimport check_allocarray, sig_free
from memory_allocator cimport MemoryAllocator

from cysignals.signals      cimport sig_check
from .conversions           cimport bit_rep_to_Vrep_list
from .base                  cimport CombinatorialPolyhedron

from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face, PolyhedronFace
from .face_list_data_structure cimport *


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


cdef class FaceIterator_base(SageObject):
    r"""
    A base class to iterate over all faces of a polyhedron.

    Construct all proper faces from the facets. In dual mode, construct all proper
    faces from the vertices. Dual will be faster for less vertices than facets.

    See :class:`FaceIterator`.
    """
    def __cinit__(self, P, dual=None, output_dimension=None):
        r"""
        Initialize :class:`FaceIterator_base`.

        See :class:`FaceIterator_base`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator() # indirect doctest

            sage: f_vector = [1, 0, 0, 0, 1]
            sage: for face in it: f_vector[face.dimension()+1] += 1
            sage: print ('f_vector of permutahedron(4): ', f_vector)
            f_vector of permutahedron(4):  [1, 24, 36, 14, 1]

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator).run()
        """
        # Note that all values are set to zero at the time ``__cinit__`` is called:
        # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#initialisation-methods
        # In particular, ``__dealloc__`` will not do harm in this case.

        cdef CombinatorialPolyhedron C

        # Working around that __cinit__ of base and derived class must be the same,
        # as extension classes do not yet have __new__ in Cython 0.29.
        if isinstance(P, CombinatorialPolyhedron):
            C = P
        else:
            C = P.combinatorial_polyhedron()
            if dual is None:
                # Determine the (likely) faster way, to iterate through all faces.
                if not P.is_compact() or P.n_facets() <= P.n_vertices():
                    dual = False
                else:
                    dual = True

            if output_dimension is not None and (output_dimension < 0 or output_dimension >= P.dim()):
                # In those cases the output will be completely handled by :meth:`FaceIterator_geom.__next__`.
                output_dimension = None

        if dual and not C.is_bounded():
            raise ValueError("cannot iterate over dual of unbounded Polyedron")
        cdef int i
        cdef size_t j

        self.dual = dual
        self.structure.dual = dual
        self.structure.dimension = C.dimension()
        self.structure.current_dimension = self.structure.dimension - 1
        self.structure.highest_dimension = self.structure.dimension - 1

        # We will not yield the empty face.
        # If there are `n` lines, than there
        # are no faces below dimension `n`.
        # The dimension of the level-sets in the face lattice jumps from `n` to `-1`.
        self.structure.lowest_dimension = 0

        if output_dimension is not None:
            if output_dimension not in range(self.structure.dimension):
                raise ValueError("``output_dimension`` must be the dimension of proper faces")
            if self.dual:
                # In dual mode, the dimensions are reversed.
                self.structure.output_dimension = self.structure.dimension - 1 - output_dimension
            else:
                self.structure.output_dimension = output_dimension
            self.structure.lowest_dimension = max(0, self.structure.output_dimension)
        else:
            self.structure.output_dimension = -2

        if dual:
            self.atoms = C.bitrep_facets()
            self.coatoms = C.bitrep_Vrep()
        else:
            self.coatoms = C.bitrep_facets()
            self.atoms = C.bitrep_Vrep()
        self._Vrep = C.Vrep()
        self._facet_names = C.facet_names()
        self._n_facets = C.bitrep_facets().n_faces()
        self._equations = C.equations()
        if self._equations:
            self._n_equations = len(self._equations)
        else:
            self._n_equations = 0
        self._bounded = C.is_bounded()
        self._far_face[0] = C._far_face[0]

        self.structure.atom_rep = <size_t *> check_allocarray(self.coatoms.n_atoms(), sizeof(size_t))
        self.structure.coatom_rep = <size_t *> check_allocarray(self.coatoms.n_faces(), sizeof(size_t))

        if self.structure.dimension == 0 or self.coatoms.n_faces() == 0:
            # As we will only yield proper faces,
            # there is nothing to yield in those cases.
            # We have to discontinue initialization,
            # as it assumes ``self.dimension > 0`` and ``self.n_faces > 0``.
            self.structure.current_dimension = self.structure.dimension
            return
        # We may assume ``dimension > 0`` and ``n_faces > 0``.

        # Initialize ``new_faces``.
        self.structure.new_faces = <face_list_t*> check_calloc((self.structure.dimension), sizeof(face_list_t))
        for i in range(self.structure.dimension):
            face_list_init(self.structure.new_faces[i],
                           self.coatoms.n_faces(), self.coatoms.n_atoms(),
                           self.coatoms.n_coatoms())

        face_list_copy(self.structure.new_faces[self.structure.dimension-1], self.coatoms.data)

        # Initialize ``visited_all``.
        self.structure.visited_all = <face_list_t*> check_calloc((self.structure.dimension), sizeof(face_list_t))
        face_list_shallow_init(self.structure.visited_all[self.structure.dimension-1],
                               self.coatoms.n_faces(), self.coatoms.n_atoms(),
                               self.coatoms.n_coatoms())
        self.structure.visited_all[self.structure.dimension-1].n_faces = 0

        if not C.is_bounded():
            # Treating the far face as if we had visited all its elements.
            # Hence we will visit all intersections of facets unless contained in the far face.

            # Regarding the length of ``self.visited_all``:
            # The last facet will not yield any new faces thus the length of ``visited_all``
            # needs to be at most ``n_facets - 1``.
            # Hence it is fine to use the first entry already for the far face,
            # as ``self.visited_all`` holds ``n_facets`` pointers.
            add_face_shallow(self.structure.visited_all[self.structure.dimension-1], self._far_face)

        # Initialize ``first_time``.
        self.structure.first_time = <bint *> check_allocarray(self.structure.dimension, sizeof(bint))
        self.structure.first_time[self.structure.dimension - 1] = True

        self.structure.yet_to_visit = self.coatoms.n_faces()
        self.structure._index = 0

        self.structure.n_coatoms = self.coatoms.n_faces()

        if C.is_bounded() and ((dual and C.is_simplicial()) or (not dual and C.is_simple())):
            # We are in the comfortable situation that for our iterator
            # all intervals not containing the 0 element are boolean.
            # This makes things a lot easier.
            self.structure.new_faces[self.structure.dimension -1].polyhedron_is_simple = True
        else:
            self.structure.new_faces[self.structure.dimension -1].polyhedron_is_simple = False

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator import FaceIterator_base
            sage: FaceIterator_base(2)  # indirect doctest
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.integer.Integer' object has no attribute 'combinatorial_polyhedron'
        """
        cdef int i
        sig_free(self.structure.atom_rep)
        sig_free(self.structure.coatom_rep)
        sig_free(self.structure.first_time)
        if self.structure.visited_all:
            face_list_shallow_free(self.structure.visited_all[self.structure.dimension - 1])
            sig_free(self.structure.visited_all)
        if self.structure.new_faces:
            for i in range(self.structure.dimension):
                face_list_free(self.structure.new_faces[i])
            sig_free(self.structure.new_faces)

    def reset(self):
        r"""
        Reset the iterator.

        The iterator will start with the first face again.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: it = C.face_generator()
            sage: next(it).ambient_V_indices()
            (0, 3, 4, 5)
            sage: it.reset()
            sage: next(it).ambient_V_indices()
            (0, 3, 4, 5)

        TESTS:

        Resetting will fix the order of the coatoms after ``only_subsets``::

            sage: P = polytopes.Birkhoff_polytope(3)
            sage: C = P.combinatorial_polyhedron()
            sage: it = C.face_generator(algorithm='primal')
            sage: face = next(it)
            sage: face.ambient_H_indices(add_equations=False)
            (8,)
            sage: face = next(it)
            sage: face.ambient_H_indices(add_equations=False)
            (7,)
            sage: it.only_subfaces()
            sage: it.reset()
            sage: face = next(it)
            sage: face.ambient_H_indices(add_equations=False)
            (8,)
        """
        if self.structure.dimension == 0 or self.coatoms.n_faces() == 0:
            # As we will only yield proper faces,
            # there is nothing to yield in those cases.
            # We have to discontinue initialization,
            # as it assumes ``self.dimension > 0`` and ``self.n_faces > 0``.
            self.structure.current_dimension = self.structure.dimension
            return
        if self._bounded:
            self.structure.visited_all[self.structure.dimension -1].n_faces = 0
        else:
            self.structure.visited_all[self.structure.dimension -1].n_faces = 1
        self.structure.face_status = FaceStatus.NOT_INITIALIZED
        self.structure.new_faces[self.structure.dimension - 1].n_faces = self.coatoms.n_faces()
        self.structure.current_dimension = self.structure.dimension - 1
        self.structure.highest_dimension = self.structure.dimension - 1
        self.structure.first_time[self.structure.dimension - 1] = True

        self.structure.yet_to_visit = self.coatoms.n_faces()
        self.structure._index = 0

        # ``only_subsets`` might have messed up the coatoms.
        face_list_copy(self.structure.new_faces[self.structure.dimension-1], self.coatoms.data)

    def __next__(self):
        r"""
        Must be implemented by a derived class.

        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator \
            ....:         import FaceIterator_base
            sage: P = polytopes.octahedron()
            sage: C = CombinatorialPolyhedron(P)
            sage: next(FaceIterator_base(C, False))
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        raise NotImplementedError("a derived class must implement this")

    next = __next__

    def current(self):
        r"""
        Retrieve the last value of :meth:`next`.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: it = P.combinatorial_polyhedron().face_generator()
            sage: next(it)
            A 0-dimensional face of a 3-dimensional combinatorial polyhedron
            sage: it.current()
            A 0-dimensional face of a 3-dimensional combinatorial polyhedron
            sage: next(it).ambient_V_indices() == it.current().ambient_V_indices()
            True
        """
        if unlikely(self.structure.face_status == FaceStatus.NOT_INITIALIZED):
            raise ValueError("iterator not set to a face yet")
        return CombinatorialFace(self)

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: [d for d in it]
            [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron]
        """
        return self

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: it1 = loads(it.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def ignore_subfaces(self):
        r"""
        The iterator will not visit any faces of the current face.

        Only possible when not in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(algorithm='primal')
            sage: n_non_simplex_faces = 1
            sage: for face in it:
            ....:     if face.n_ambient_Vrepresentation() > face.dimension() + 1:
            ....:         n_non_simplex_faces += 1
            ....:     else:
            ....:         it.ignore_subfaces()
            ....:
            sage: n_non_simplex_faces
            127

        Face iterator must not be in dual mode::

            sage: it = C.face_generator(algorithm='dual')
            sage: _ = next(it)
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when not in dual mode

        Ignoring the same face as was requested to visit only consumes the iterator::

            sage: it = C.face_generator(algorithm='primal')
            sage: _ = next(it)
            sage: it.only_subfaces()
            sage: it.ignore_subfaces()
            sage: list(it)
            []

        Face iterator must be set to a face first::

            sage: it = C.face_generator(algorithm='primal')
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: iterator not set to a face yet
        """
        if unlikely(self.dual):
            raise ValueError("only possible when not in dual mode")
        self.ignore_subsets()

    def ignore_supfaces(self):
        r"""
        The iterator will not visit any faces containing the current face.

        Only possible when in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(algorithm='dual')
            sage: n_faces_with_non_simplex_quotient = 1
            sage: for face in it:
            ....:     n_facets = face.n_ambient_Hrepresentation(add_equations=False)
            ....:     if n_facets > C.dimension() - face.dimension() + 1:
            ....:         n_faces_with_non_simplex_quotient += 1
            ....:     else:
            ....:         it.ignore_supfaces()
            ....:
            sage: n_faces_with_non_simplex_quotient
            4845

        Face iterator must be in dual mode::

            sage: it = C.face_generator(algorithm='primal')
            sage: _ = next(it)
            sage: it.ignore_supfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when in dual mode
        """
        if unlikely(not self.dual):
            raise ValueError("only possible when in dual mode")
        self.ignore_subsets()

    def meet_of_Hrep(self, *indices):
        r"""
        Construct the meet of the facets indicated by the indices.

        This is the largest face contained in all facets with the given indices.

        The iterator must be reset if not newly initialized.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: it.meet_of_Hrep(1,2)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            sage: it.meet_of_Hrep(1,2).ambient_H_indices()
            (1, 2)
            sage: it.meet_of_Hrep(1,3).ambient_H_indices()
            (1, 3)
            sage: it.meet_of_Hrep(1,5).ambient_H_indices()
            (0, 1, 2, 3, 4, 5)

            sage: P = polytopes.cross_polytope(4)
            sage: it = P.face_generator()
            sage: it.meet_of_Hrep().ambient_H_indices()
            ()
            sage: it.meet_of_Hrep(1,3).ambient_H_indices()
            (1, 2, 3, 4)
            sage: it.meet_of_Hrep(1,2).ambient_H_indices()
            (1, 2)
            sage: it.meet_of_Hrep(1,6).ambient_H_indices()
            (1, 6)
            sage: it.meet_of_Hrep(1,2,6).ambient_H_indices()
            (1, 2, 6, 7)
            sage: it.meet_of_Hrep(1,2,5,6).ambient_H_indices()
            (0, 1, 2, 3, 4, 5, 6, 7)

            sage: s = cones.schur(4)
            sage: C = CombinatorialPolyhedron(s)
            sage: it = C.face_generator()
            sage: it.meet_of_Hrep(1,2).ambient_H_indices()
            (1, 2)
            sage: it.meet_of_Hrep(1,2,3).ambient_H_indices()
            Traceback (most recent call last):
            ...
            IndexError: coatoms out of range

        If the iterator has already been used, it must be reset before::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: it = P.face_generator()
            sage: _ = next(it), next(it)
            sage: next(it).ambient_V_indices()
            (15, 16, 17, 18, 19)
            sage: it.meet_of_Hrep(9,11)
            Traceback (most recent call last):
            ...
            ValueError: please reset the face iterator
            sage: it.reset()
            sage: it.meet_of_Hrep(9,11).ambient_H_indices()
            (9, 11)

        TESTS:

        Check that things work fine, if the face iterator was never properly initialized::

            sage: P = Polyhedron()
            sage: P.meet_of_Hrep()
            A -1-dimensional face of a Polyhedron in ZZ^0
            sage: P = Polyhedron([[0,0]])
            sage: P.meet_of_Hrep()
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: P.meet_of_Hrep(0)
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: P = Polyhedron(lines=[[1]])
            sage: P.meet_of_Hrep()
            A 1-dimensional face of a Polyhedron in ZZ^1 defined as the convex hull of 1 vertex and 1 line
            sage: P = Polyhedron(lines=[[1, 1]])
            sage: P.meet_of_Hrep()
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 line
            sage: P.meet_of_Hrep(0)
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 line
        """
        # Ignore equations.
        indices = [i for i in indices
                   if not (self._n_facets <= i < self._n_facets + self._n_equations)]
        if self.dual:
            return self._join_of_atoms(*indices)
        else:
            return self._meet_of_coatoms(*indices)

    def join_of_Vrep(self, *indices):
        r"""
        Construct the join of the Vrepresentatives indicated by the indices.

        This is the smallest face containing all Vrepresentatives with the given indices.

        The iterator must be reset if not newly initialized.

        .. NOTE::

            In the case of unbounded polyhedra, the smallest face containing given Vrepresentatives
            may not be well defined.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: it.join_of_Vrep(1)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it.join_of_Vrep(1,2).ambient_V_indices()
            (1, 2)
            sage: it.join_of_Vrep(1,3).ambient_V_indices()
            (0, 1, 2, 3)
            sage: it.join_of_Vrep(1,5).ambient_V_indices()
            (0, 1, 5, 6)

            sage: P = polytopes.cross_polytope(4)
            sage: it = P.face_generator()
            sage: it.join_of_Vrep().ambient_V_indices()
            ()
            sage: it.join_of_Vrep(1,3).ambient_V_indices()
            (1, 3)
            sage: it.join_of_Vrep(1,2).ambient_V_indices()
            (1, 2)
            sage: it.join_of_Vrep(1,6).ambient_V_indices()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: it.join_of_Vrep(8)
            Traceback (most recent call last):
            ...
            IndexError: coatoms out of range

        If the iterator has already been used, it must be reset before::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: it = P.face_generator()
            sage: _ = next(it), next(it)
            sage: next(it).ambient_V_indices()
            (15, 16, 17, 18, 19)
            sage: it.join_of_Vrep(1,10)
            Traceback (most recent call last):
            ...
            ValueError: please reset the face iterator
            sage: it.reset()
            sage: it.join_of_Vrep(1,10).ambient_V_indices()
            (1, 10)

        In the case of an unbounded polyhedron, we try to make sense of the input::

            sage: P = polytopes.cube()*Polyhedron(lines=[[1]])
            sage: it = P.face_generator()
            sage: it.join_of_Vrep(1)
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 1 vertex and 1 line
            sage: it.join_of_Vrep(0, 1)
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 1 vertex and 1 line
            sage: it.join_of_Vrep(0)
            Traceback (most recent call last):
            ...
            ValueError: the join is not well-defined

            sage: P = Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]])
            sage: it = P.face_generator()
            sage: it.join_of_Vrep(0)
            A 0-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: it.join_of_Vrep(1)
            A 0-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: it.join_of_Vrep(2)
            Traceback (most recent call last):
            ...
            ValueError: the join is not well-defined
            sage: it.join_of_Vrep(0,2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray

            sage: P = Polyhedron(rays=[[1,0], [0,1]])
            sage: it = P.face_generator()
            sage: it.join_of_Vrep(0)
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: it.join_of_Vrep(1,2)
            A 2-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 rays

        TESTS:

        Check that things work fine, if the face iterator was never properly initialized::

            sage: P = Polyhedron()
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^0
            sage: P = Polyhedron([[0,0]])
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^2
            sage: P.join_of_Vrep(0)
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: P = Polyhedron(lines=[[1]])
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^1
            sage: P.join_of_Vrep(0)
            A 1-dimensional face of a Polyhedron in ZZ^1 defined as the convex hull of 1 vertex and 1 line
            sage: P = Polyhedron(lines=[[1, 1]])
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^2
            sage: P.Vrepresentation()
            (A line in the direction (1, 1), A vertex at (0, 0))
            sage: P.join_of_Vrep(0)
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 line
            sage: P.join_of_Vrep(1)
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 1 line
            sage: P = Polyhedron(lines=[[1, 0], [0, 1]])
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^2
            sage: P.join_of_Vrep(0)
            A 2-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: P.join_of_Vrep(0, 1)
            A 2-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: P.join_of_Vrep(0, 1, 2)
            A 2-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
        """
        if not self.dual:
            return self._join_of_atoms(*indices)
        else:
            return self._meet_of_coatoms(*indices)

    def _meet_of_coatoms(self, *indices):
        r"""
        Construct the meet of the coatoms indicated by the indices.

        The iterator must be reset if not newly initialized.

        .. SEEALSO::

            :meth:`meet_of_Hrep`,
            :meth:`join_of_Vrep`.

        EXAMPLES:

        In non-dual mode we construct the meet of facets::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='primal')
            sage: it._meet_of_coatoms(1,2)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            sage: it._meet_of_coatoms(1,2,3)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it._meet_of_coatoms(1,2,3).ambient_H_indices()
            (1, 2, 3)

        In dual mode we construct the join of vertices/rays::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='dual')
            sage: it._meet_of_coatoms(1,2)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            sage: it._meet_of_coatoms(1,2,3)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: it._meet_of_coatoms(1)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex

        The face iterator must not have the output dimension specified::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: it = P.face_generator(2)
            sage: it._meet_of_coatoms(1,2)
            Traceback (most recent call last):
            ...
            ValueError: face iterator must not have the output dimension specified

        TESTS:

        We prevent a segmentation fault::

            sage: P = polytopes.simplex()
            sage: it = P.face_generator()
            sage: it._meet_of_coatoms(-1)
            Traceback (most recent call last):
            ...
            IndexError: coatoms out of range
            sage: it._meet_of_coatoms(100)
            Traceback (most recent call last):
            ...
            IndexError: coatoms out of range

        The empty face is detected correctly, even with lines or rays::

            sage: P = polytopes.cube()*Polyhedron(lines=[[1]])
            sage: it = P.face_generator()
            sage: it._meet_of_coatoms(1,2,4,5)
            A -1-dimensional face of a Polyhedron in ZZ^4

            sage: P = Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]])
            sage: it = P.face_generator()
            sage: it._meet_of_coatoms(0)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: it._meet_of_coatoms(1)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: it._meet_of_coatoms(2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: it._meet_of_coatoms(1, 2)
            A -1-dimensional face of a Polyhedron in QQ^2
        """
        if unlikely(self.structure.face_status != FaceStatus.NOT_INITIALIZED):
            raise ValueError("please reset the face iterator")
        if unlikely(self.structure.output_dimension != -2):
            raise ValueError("face iterator must not have the output dimension specified")

        cdef size_t n_atoms = self.coatoms.n_atoms()
        cdef size_t n_coatoms = self.coatoms.n_faces()
        cdef ListOfFaces coatoms = self.coatoms

        cdef ListOfFaces face_mem = ListOfFaces(1, n_atoms, n_coatoms)
        cdef face_t face = face_mem.data.faces[0]
        cdef int i
        cdef size_t j

        # Initialize the full polyhedron.
        for j in range(n_atoms):
            face_add_atom(face, j)

        for i in indices:
            if not 0 <= i < n_coatoms:
                raise IndexError("coatoms out of range")
            face_intersection(face, face, coatoms.data.faces[i])

        if not self._bounded and face_issubset(face, self._far_face):
            # The meet is contained in the far face and therefore is the empty face.
            face_clear(face)

        self.find_face(face)
        output = self.current()
        self.reset()
        return output

    def _join_of_atoms(self, *indices):
        r"""
        Construct the join of atoms indicated by the indices.

        The iterator must be reset if not newly initialized.

        .. SEEALSO::

            :meth:`meet_of_Hrep`,
            :meth:`join_of_Vrep`.

        EXAMPLES:

        In dual mode we construct the meet of facets::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='dual')
            sage: it._join_of_atoms(1,2)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            sage: it._join_of_atoms(1,2,3)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it._join_of_atoms(1,2,3).ambient_H_indices()
            (1, 2, 3)

        In non-dual mode we construct the join of vertices/rays::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='primal')
            sage: it._join_of_atoms(1,2)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            sage: it._join_of_atoms(1,2,3)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: it._join_of_atoms(1)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex

        If the iterator has already been used, it must be reset before::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: it = P.face_generator()
            sage: _ = next(it), next(it)
            sage: next(it).ambient_V_indices()
            (15, 16, 17, 18, 19)
            sage: it._join_of_atoms(1,10)
            Traceback (most recent call last):
            ...
            ValueError: please reset the face iterator
            sage: it.reset()
            sage: it._join_of_atoms(1,10).ambient_V_indices()
            (1, 10)

        The face iterator must not have the output dimension specified::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: it = P.face_generator(2)
            sage: it._join_of_atoms(1,2)
            Traceback (most recent call last):
            ...
            ValueError: face iterator must not have the output dimension specified

        TESTS:

        We prevent a segmentation fault::

            sage: P = polytopes.simplex()
            sage: it = P.face_generator()
            sage: it._join_of_atoms(-1)
            Traceback (most recent call last):
            ...
            IndexError: atoms out of range
            sage: it._join_of_atoms(100)
            Traceback (most recent call last):
            ...
            IndexError: atoms out of range
        """
        if unlikely(self.structure.face_status != FaceStatus.NOT_INITIALIZED):
            raise ValueError("please reset the face iterator")
        if unlikely(self.structure.output_dimension != -2):
            raise ValueError("face iterator must not have the output dimension specified")

        cdef size_t n_atoms = self.coatoms.n_atoms()
        cdef size_t n_coatoms = self.coatoms.n_faces()
        cdef ListOfFaces coatoms = self.coatoms

        cdef ListOfFaces face_mem = ListOfFaces(2, n_atoms, n_coatoms)
        cdef face_t face = face_mem.data.faces[0]
        cdef face_t pseudo_face = face_mem.data.faces[1]
        cdef int j
        cdef size_t i

        if not all(0 <= j < n_atoms for j in indices):
            raise IndexError("atoms out of range")

        # Initialize a pseudo_face as indicated by the indices.
        for i in indices:
            face_add_atom(pseudo_face, i)

        # Initialize the full polyhedron.
        for i in range(n_atoms):
            face_add_atom(face, i)

        # Now we intersect all faces that contain our pseudo_face.
        for i in range(n_coatoms):
            if face_issubset(pseudo_face, coatoms.data.faces[i]):
                face_intersection(face, face, coatoms.data.faces[i])

        if not indices:
            # The neutral element of the join.
            face_clear(face)
        elif not self._bounded and face_issubset(face, self._far_face):
            # The join is not well-defined.
            # We allow for unbounded polyhedra to compute the join, even with rays.
            # However, the result is not necesarrily well-defined.
            raise ValueError("the join is not well-defined")

        self.find_face(face)
        output = self.current()
        self.reset()
        return output

    cdef int ignore_subsets(self) except -1:
        r"""
        Ignore sub-/supfaces of the current face.

        In non-dual mode ignores all subfaces of the current face.
        In dual mode ignores all supfaces of the current face.

        See :meth:`FaceIterator_base.ignore_subfaces` and
        :meth:`FaceIterator_base.ignore_supfaces`.
        """
        if unlikely(self.structure.face_status == FaceStatus.NOT_INITIALIZED):
            raise ValueError("iterator not set to a face yet")
        if unlikely(self.structure.face_status == FaceStatus.ONLY_VISIT_SUBSETS):
            # The iterator is consumed, if it was just set to visit only subsets
            # next thing to ignore subsets.
            self.structure.current_dimension = self.structure.dimension
            return 0
        if unlikely(self.structure.face_status == FaceStatus.IGNORE_SUBSETS):
            # Nothing to do.
            return 0
        # The current face is added to ``visited_all``.
        # This will make the iterator skip those faces.
        # Also, this face will not be added a second time to ``visited_all``,
        # as there are no new faces.

        add_face_shallow(self.structure.visited_all[self.structure.current_dimension], self.structure.face)
        self.structure.face_status = FaceStatus.IGNORE_SUBSETS

    def only_subfaces(self):
        r"""
        The iterator will visit all (remaining) subfaces of the current face and then terminate.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: next(it).ambient_H_indices()
            ()
            sage: next(it).ambient_H_indices()
            (0, 1, 2, 3, 4, 5)
            sage: next(it).ambient_H_indices()
            (5,)
            sage: next(it).ambient_H_indices()
            (4,)
            sage: it.only_subfaces()
            sage: list(f.ambient_H_indices() for f in it)
            [(4, 5), (3, 4), (1, 4), (0, 4), (3, 4, 5), (0, 4, 5), (1, 3, 4), (0, 1, 4)]

        ::

            sage: P = polytopes.Birkhoff_polytope(4)
            sage: C = P.combinatorial_polyhedron()
            sage: it = C.face_generator()
            sage: next(it).ambient_H_indices(add_equations=False)
            (15,)
            sage: next(it).ambient_H_indices(add_equations=False)
            (14,)
            sage: it.only_subfaces()
            sage: all(14 in f.ambient_H_indices() for f in it)
            True

        Face iterator needs to be set to a face first::

            sage: it = C.face_generator()
            sage: it.only_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: iterator not set to a face yet

        Face iterator must not be in dual mode::

            sage: it = C.face_generator(algorithm='dual')
            sage: _ = next(it)
            sage: it.only_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when not in dual mode

        Cannot run ``only_subfaces`` after ``ignore_subfaces``::

            sage: it = C.face_generator()
            sage: _ = next(it)
            sage: it.ignore_subfaces()
            sage: it.only_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: cannot only visit subsets after ignoring a face
        """
        if unlikely(self.dual):
            raise ValueError("only possible when not in dual mode")
        self.only_subsets()

    def only_supfaces(self):
        r"""
        The iterator will visit all (remaining) faces
        containing the current face and then terminate.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3)
            sage: it = P.face_generator()
            sage: next(it).ambient_V_indices()
            (0, 1, 2, 3, 4, 5)
            sage: next(it).ambient_V_indices()
            ()
            sage: next(it).ambient_V_indices()
            (5,)
            sage: next(it).ambient_V_indices()
            (4,)
            sage: it.only_supfaces()
            sage: list(f.ambient_V_indices() for f in it)
            [(4, 5), (3, 4), (2, 4), (0, 4), (3, 4, 5), (2, 4, 5), (0, 3, 4), (0, 2, 4)]

        ::

            sage: P = polytopes.Birkhoff_polytope(4)
            sage: C = P.combinatorial_polyhedron()
            sage: it = C.face_generator(algorithm='dual')
            sage: next(it).ambient_V_indices()
            (23,)
            sage: next(it).ambient_V_indices()
            (22,)
            sage: it.only_supfaces()
            sage: all(22 in f.ambient_V_indices() for f in it)
            True
        """
        if unlikely(not self.dual):
            raise ValueError("only possible when in dual mode")
        self.only_subsets()

    cdef int only_subsets(self) except -1:
        r"""
        Only visit sub-/supfaces of the current face and then
        terminate.

        See :meth:`FaceIterator_base.only_subfaces` and
        :meth:`FaceIterator_base.only_supfaces`.
        """
        if unlikely(self.structure.face_status == FaceStatus.NOT_INITIALIZED):
            raise ValueError("iterator not set to a face yet")
        if unlikely(self.structure.face_status == FaceStatus.IGNORE_SUBSETS):
            raise ValueError("cannot only visit subsets after ignoring a face")

        cdef face_list_t* faces = &self.structure.new_faces[self.structure.current_dimension]
        cdef size_t yet_to_visit = self.structure.yet_to_visit

        if unlikely(yet_to_visit >= faces[0].n_faces
                or not faces_are_identical(faces[0].faces[yet_to_visit], self.structure.face)):
            raise ValueError("iterator is not set to the correct face")

        swap_faces(faces[0].faces[yet_to_visit], faces[0].faces[faces[0].n_faces - 1])

        self.structure.face_status = FaceStatus.ONLY_VISIT_SUBSETS
        self.structure.yet_to_visit = 0
        # This will work:
        # ``next_dimension`` will first call ``next_face_loop`` and then check
        # for the dimension. By this time the current dimension has changed.
        self.structure.highest_dimension = self.structure.current_dimension - 1

    cdef inline CombinatorialFace next_face(self):
        r"""
        Set attribute ``face`` to the next face and return it as
        :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace`.
        """
        self.next_dimension()
        if unlikely(self.structure.current_dimension > self.structure.highest_dimension):
            return None
        return CombinatorialFace(self)

    cdef inline int next_dimension(self) except -1:
        r"""
        Set attribute ``face`` to the next face and return the dimension.

        Will return the dimension of the polyhedron on failure.

        The function calls :meth:`FaceIterator_base.next_face_loop` until a new
        face is set or until the iterator is consumed.

        .. NOTE::

            The face_iterator can be prevented from visiting any subfaces
            (or supfaces in dual mode) as in :meth:`FaceIterator_base.ignore_subfaces`
            and :meth`FaceIterator_base.ignore_supfaces`.

            Those methods add the current face to ``visited_all`` before
            visiting sub-/supfaces instead of after. One cannot arbitrarily
            add faces to ``visited_all``, as visited_all has a maximal length.
        """
        return next_dimension(self.structure)

    cdef inline int next_face_loop(self) except -1:
        r"""
        Set attribute ``face`` to the next face. On success return `1`.
        Otherwise `0`. Needs to be recalled then.

        If ``self.current_dimension == self.dimension``, then the iterator is
        consumed.
        """
        return next_face_loop(self.structure)

    cdef inline size_t n_atom_rep(self) except -1:
        r"""
        Compute the number of atoms in the current face by counting the
        number of set bits.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.n_atom_rep`
        """
        return n_atom_rep(self.structure)

    cdef size_t set_coatom_rep(self) except -1:
        r"""
        Set ``coatom_rep`` to be the coatom-representation of the current face.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_coatom_rep`
        """
        return bit_rep_to_coatom_rep(self.structure.face, self.coatoms.data, self.structure.coatom_rep)

    cdef size_t set_atom_rep(self) except -1:
        r"""
        Set ``atom_rep`` to be the atom-representation of the current face.
        Return its length.

        This is a shortcut of :class:`sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace.set_atom_rep`
        """
        return bit_rep_to_Vrep_list(self.structure.face, self.structure.atom_rep)

    cdef int find_face(self, face_t face) except -1:
        """
        Iterate until the current face is ``face``.

        The value can then be obtained with :meth:`current`.

        The iterator is assumed to be newly initialized or reset.
        See :meth:`FaceIterator_base._join_of_atoms` and
        :meth:`FaceIterator_base._meet_of_coatoms`.
        """
        cdef size_t n_atoms = face_len_atoms(face)

        if n_atoms == self.coatoms.n_atoms():
            # The face is the universe.
            self.structure.face[0] = face[0]
            self.structure.face_status = FaceStatus.INITIALIZED
            self.structure.current_dimension = self.structure.dimension
            return 0
        elif n_atoms == 0:
            # The face is the empty face.
            self.structure.face[0] = face[0]
            self.structure.face_status = FaceStatus.INITIALIZED
            self.structure.current_dimension = -1
            return 0

        cdef int d = self.next_dimension()
        while self.structure.current_dimension != self.structure.dimension:
            if face_issubset(face, self.structure.face):
                if face_issubset(self.structure.face, face):
                    # Found our face.
                    return 0
            else:
                # The face is not a subface/supface of the current face.
                self.ignore_subsets()

            d = self.next_dimension()

        raise ValueError("the face appears to be incorrect")


cdef class FaceIterator(FaceIterator_base):
    r"""
    A class to iterate over all combinatorial faces of a polyhedron.

    Construct all proper faces from the facets. In dual mode, construct all proper
    faces from the vertices. Dual will be faster for less vertices than facets.

    INPUT:

    - ``C`` -- a :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`
    - ``dual`` -- if ``True``, then dual polyhedron is used for iteration
      (only possible for bounded Polyhedra)
    - ``output_dimension`` -- if not ``None``, then the face iterator will only yield
      faces of this dimension

    .. SEEALSO::

        :class:`FaceIterator`,
        :class:`FaceIterator_geom`,
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

    EXAMPLES:

    Construct a face iterator::

        sage: P = polytopes.cuboctahedron()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_generator()
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron

    Construct faces by the dual or not::

        sage: it = C.face_generator(algorithm='primal')
        sage: next(it).dimension()
        2

        sage: it = C.face_generator(algorithm='dual')
        sage: next(it).dimension()
        0

    For unbounded polyhedra only non-dual iteration is possible::

        sage: P = Polyhedron(rays=[[0,0,1], [0,1,0], [1,0,0]])
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_generator()
        sage: [face.ambient_Vrepresentation() for face in it]
        [(A vertex at (0, 0, 0),
          A ray in the direction (0, 1, 0),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0),),
         (A vertex at (0, 0, 0), A ray in the direction (0, 0, 1))]
        sage: it = C.face_generator(algorithm='dual')
        Traceback (most recent call last):
        ...
        ValueError: dual algorithm only available for bounded polyhedra

    Construct a face iterator only yielding dimension `2` faces::

        sage: P = polytopes.permutahedron(5)
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_generator(dimension=2)
        sage: counter = 0
        sage: for _ in it: counter += 1
        sage: print ('permutahedron(5) has', counter,
        ....:        'faces of dimension 2')
        permutahedron(5) has 150 faces of dimension 2
        sage: C.f_vector()
        (1, 120, 240, 150, 30, 1)

    In non-dual mode one can ignore all faces contained in the current face::

        sage: P = polytopes.cube()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_generator(algorithm='primal')
        sage: face = next(it)
        sage: face.ambient_H_indices()
        (5,)
        sage: it.ignore_subfaces()
        sage: [face.ambient_H_indices() for face in it]
        [(4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (3, 4),
         (1, 4),
         (0, 4),
         (1, 3, 4),
         (0, 1, 4),
         (2, 3),
         (1, 3),
         (1, 2, 3),
         (1, 2),
         (0, 2),
         (0, 1, 2),
         (0, 1)]

        sage: it = C.face_generator(algorithm='dual')
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: it.ignore_subfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when not in dual mode

    In dual mode one can ignore all faces that contain the current face::

        sage: it = C.face_generator(algorithm='dual')
        sage: next(it)
        A 0-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: face = next(it)
        sage: face.ambient_V_indices()
        (6,)
        sage: [face.ambient_V_indices() for face in it]
        [(5,),
         (4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (6, 7),
         (4, 7),
         (2, 7),
         (4, 5, 6, 7),
         (1, 2, 6, 7),
         (2, 3, 4, 7),
         (5, 6),
         (1, 6),
         (0, 1, 5, 6),
         (4, 5),
         (0, 5),
         (0, 3, 4, 5),
         (3, 4),
         (2, 3),
         (0, 3),
         (0, 1, 2, 3),
         (1, 2),
         (0, 1)]

        sage: it = C.face_generator(algorithm='primal')
        sage: next(it)
        A 2-dimensional face of a 3-dimensional combinatorial polyhedron
        sage: it.ignore_supfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when in dual mode

    ALGORITHM:

    The algorithm to visit all proper faces exactly once is roughly
    equivalent to the following. A (slightly generalized) description of the
    algorithm can be found in [KS2019]_.

    Initialization::

        faces = [set(facet) for facet in P.facets()]
        face_iterator(faces, [])

    The function ``face_iterator`` is defined recursively. It visits all faces of
    the polyhedron `P`, except those contained in any of ``visited_all``.
    It assumes ``faces`` to be exactly those facets of `P`
    that are not contained in any of the ``visited_all``.
    It assumes ``visited_all`` to be some list of faces of
    a polyhedron `P_2`, which contains `P` as one of its faces::

        def face_iterator(faces, visited_all):
            while facets:
                one_face = faces.pop()
                maybe_new_faces = [one_face.intersection(face) for face in faces]
        ...

    At this point we claim that ``maybe_new_faces`` contains all facets of ``one_face``,
    which we have not visited before.

    Proof: Let `F` be a facet of ``one_face``. We have a chain:
    `P \supset{}` ``one_face`` `{}\supset F`.
    By the diamond property, there exists a ``second_face`` with
    `P \supset{}` ``second_face`` `{}\supset F`.

    Now either ``second_face`` is not an element of ``faces``:
    Hence ``second_face`` is contained in one of ``visited_all``.
    In particular, `F` is contained in ``visited_all``.

    Or ``second_face`` is an element of ``faces``:
    Then, intersecting ``one_face`` with ``second_face`` gives ``F``.

    This concludes the proof.

    Moreover, if an element in ``maybe_new_faces`` is inclusion-maximal
    and not contained in any of the ``visited_all``, it is a facet of ``one_face``.
    Any facet in ``maybe_new_faces`` of ``one_face`` is inclusion-maximal.

    Hence, in the following loop, an element ``face1`` in ``maybe_new_faces``
    is a facet of ``one_face`` if and only if it is not contained in another facet::

        ...
                maybe_new_faces2 = []
                for i, face1 in enumerate(maybe_new_faces):
                    if (all(not face1 < face2 for face2 in maybe_new_faces[:i])
                            and all(not face1 <= face2 for face2 in maybe_new_faces[i+1:])):
                        maybe_new_faces2.append(face1)
        ...

    Now ``maybe_new_faces2`` contains only facets of ``one_face``
    and some faces contained in any of ``visited_all``.
    It also contains all the facets not contained in any of ``visited_all``.

    We construct ``new_faces`` as the list of all facets of ``one_face``
    not contained in any of ``visited_all``::

        ...
                new_faces = []
                for face1 in maybe_new_faces2:
                    if all(not face1 < face2 for face2 in visited_all):
                        new_faces.append(face1)
        ...

    By induction we can apply the algorithm, to visit all
    faces of ``one_face`` not contained in ``visited_all``::

        ...
                face_iterator(new_faces, visited_all)
        ...

    Finally we visit ``one_face`` and add it to ``visited_all``::

        ...
                visit(one_face)
                visited_all.append(one_face)

    Note: At this point, we have visited exactly those faces,
    contained in any of the ``visited_all``. The function ends here.


    ALGORITHM for the special case that all intervals of the lattice not
    containing zero are boolean (e.g. when the polyhedron is simple):

    We do not assume any other properties of our lattice in this case.
    Note that intervals of length 2 not containing zero, have exactly 2 elements now.
    But the atom-representation of faces might not be unique.

    We do the following modifications:

    - To check whether an intersection of faces is zero, we check whether the
      atom-representation is zero. Although not unique,
      it works to distinct from zero.

    - The intersection of two (relative) facets has always codimension `1` unless empty.

    - To intersect we now additionally unite the coatom representation.
      This gives the correct representation of the new face
      unless the intersection is zero.

    - To mark a face as visited, we save its coatom representation.

    - To check whether we have seen a face already, we check containment of the coatom representation.
    """
    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.associahedron(['A',3])                                  # needs sage.combinat
            sage: C = CombinatorialPolyhedron(P)                                        # needs sage.combinat
            sage: C.face_generator()                                                    # needs sage.combinat
            Iterator over the proper faces of a 3-dimensional combinatorial polyhedron

            sage: C.face_generator(1)                                                   # needs sage.combinat
            Iterator over the 1-faces of a 3-dimensional combinatorial polyhedron
        """
        if self.structure.output_dimension != -2:
            if self.dual:
                # ouput_dimension is stored with respect to the dual
                intended_dimension = self.structure.dimension - 1 - self.structure.output_dimension
            else:
                intended_dimension = self.structure.output_dimension
            output = "Iterator over the {}-faces".format(intended_dimension)
        else:
            output = "Iterator over the proper faces"
        return output + " of a {}-dimensional combinatorial polyhedron".format(self.structure.dimension)

    def __next__(self):
        r"""
        Return the next face.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: [next(it) for _ in range(7)]
            [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron]
        """
        cdef CombinatorialFace face = self.next_face()
        if unlikely(self.structure.current_dimension > self.structure.highest_dimension):
            raise StopIteration

        return face


cdef class FaceIterator_geom(FaceIterator_base):
    r"""
    A class to iterate over all geometric faces of a polyhedron.

    Construct all faces from the facets. In dual mode, construct all
    faces from the vertices. Dual will be faster for less vertices than facets.

    INPUT:

    - ``P`` -- an instance of :class:`~sage.geometry.polyhedron.base.Polyhedron_base`
    - ``dual`` -- if ``True``, then dual polyhedron is used for iteration
      (only possible for bounded Polyhedra)
    - ``output_dimension`` -- if not ``None``, then the FaceIterator will only yield
      faces of this dimension

    EXAMPLES:

    Construct a geometric face iterator::

        sage: P = polytopes.cuboctahedron()
        sage: it = P.face_generator()
        sage: next(it)
        A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 12 vertices

    Construct faces by the dual or not::

        sage: it = P.face_generator(algorithm='primal')
        sage: _ = next(it), next(it)
        sage: next(it).dim()
        2

        sage: it = P.face_generator(algorithm='dual')
        sage: _ = next(it), next(it)
        sage: next(it).dim()
        0

    For unbounded polyhedra only non-dual iteration is possible::

        sage: P = Polyhedron(rays=[[0,0,1], [0,1,0], [1,0,0]])
        sage: it = P.face_generator()
        sage: [face.ambient_Vrepresentation() for face in it]
        [(A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (0, 1, 0),
          A ray in the direction (1, 0, 0)),
         (),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 1, 0),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0),
          A ray in the direction (0, 0, 1),
          A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (1, 0, 0)),
         (A vertex at (0, 0, 0), A ray in the direction (0, 1, 0)),
         (A vertex at (0, 0, 0),),
         (A vertex at (0, 0, 0), A ray in the direction (0, 0, 1))]
        sage: it = P.face_generator(algorithm='dual')
        Traceback (most recent call last):
        ...
        ValueError: cannot iterate over dual of unbounded Polyedron

    Construct a FaceIterator only yielding dimension `2` faces::

        sage: P = polytopes.permutahedron(5)
        sage: it = P.face_generator(face_dimension=2)
        sage: counter = 0
        sage: for _ in it: counter += 1
        sage: print ('permutahedron(5) has', counter,
        ....:        'faces of dimension 2')
        permutahedron(5) has 150 faces of dimension 2
        sage: P.f_vector()
        (1, 120, 240, 150, 30, 1)

    In non-dual mode one can ignore all faces contained in the current face::

        sage: P = polytopes.cube()
        sage: it = P.face_generator(algorithm='primal')
        sage: _ = next(it), next(it)
        sage: face = next(it)
        sage: face.ambient_H_indices()
        (5,)
        sage: it.ignore_subfaces()
        sage: [face.ambient_H_indices() for face in it]
        [(4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (3, 4),
         (1, 4),
         (0, 4),
         (1, 3, 4),
         (0, 1, 4),
         (2, 3),
         (1, 3),
         (1, 2, 3),
         (1, 2),
         (0, 2),
         (0, 1, 2),
         (0, 1)]

        sage: it = P.face_generator(algorithm='dual')
        sage: _ = next(it), next(it)
        sage: next(it)
        A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
        sage: it.ignore_subfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when not in dual mode

    In dual mode one can ignore all faces that contain the current face::

        sage: P = polytopes.cube()
        sage: it = P.face_generator(algorithm='dual')
        sage: _ = next(it), next(it)
        sage: next(it)
        A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
        sage: face = next(it)
        sage: face.ambient_V_indices()
        (6,)
        sage: [face.ambient_V_indices() for face in it]
        [(5,),
         (4,),
         (3,),
         (2,),
         (1,),
         (0,),
         (6, 7),
         (4, 7),
         (2, 7),
         (4, 5, 6, 7),
         (1, 2, 6, 7),
         (2, 3, 4, 7),
         (5, 6),
         (1, 6),
         (0, 1, 5, 6),
         (4, 5),
         (0, 5),
         (0, 3, 4, 5),
         (3, 4),
         (2, 3),
         (0, 3),
         (0, 1, 2, 3),
         (1, 2),
         (0, 1)]

        sage: it = P.face_generator(algorithm='primal')
        sage: _ = next(it), next(it)
        sage: next(it)
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
        sage: it.ignore_supfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when in dual mode

    .. SEEALSO::

        :class:`FaceIterator_base`.
    """
    def __init__(self, P, dual=None, output_dimension=None):
        r"""
        Initialize :class:`FaceIterator_base`.

        See :class:`FaceIterator_base`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(4)
            sage: it = P.face_generator()  # indirect doctest
            sage: it
            Iterator over the faces of a 3-dimensional polyhedron in ZZ^4
            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator_geom).run()
        """
        self._requested_dim = output_dimension
        self.P = P
        # Base class only has __cinit__ and not __init__
        self.reset()

    def reset(self):
        r"""
        Reset the iterator.

        The iterator will start with the first face again.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: next(it).ambient_V_indices()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: next(it).ambient_V_indices()
            ()
            sage: next(it).ambient_V_indices()
            (0, 3, 4, 5)
            sage: it.reset()
            sage: next(it).ambient_V_indices()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: next(it).ambient_V_indices()
            ()
            sage: next(it).ambient_V_indices()
            (0, 3, 4, 5)
        """
        output_dimension = self._requested_dim
        P = self.P

        if output_dimension is None:
            if P.dim() == -1:
                self._trivial_faces = 1  # the empty polyhedron, only yield the empty face
            else:
                self._trivial_faces = 2  # yield the universe, then the empty face, than all other faces
        elif output_dimension == P.dim():
            self._trivial_faces = 4  # only yield the full-dimensional face and no other faces
        elif output_dimension == -1:
            self._trivial_faces = 3  # only yield the empty face and no other faces
        elif output_dimension < -1 or output_dimension > P.dim():
            self._trivial_faces = -1  # don't yield any faces at all
        else:
            self._trivial_faces = 0  # yield the faces of the requested dimension
        FaceIterator_base.reset(self)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.associahedron(['A',3])                                  # needs sage.combinat
            sage: P.face_generator()                                                    # needs sage.combinat
            Iterator over the faces of a 3-dimensional polyhedron in QQ^3

            sage: P.face_generator(1)                                                   # needs sage.combinat
            Iterator over the 1-faces of a 3-dimensional polyhedron in QQ^3
        """
        if self._requested_dim is not None:
            output = "Iterator over the {}-faces".format(self._requested_dim)
        else:
            output = "Iterator over the faces"
        return output + " of a {}-dimensional polyhedron in {}".format(self.structure.dimension, self.P.parent()._repr_ambient_module())

    def __next__(self):
        r"""
        Return the next face.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: [next(it) for _ in range(7)]
            [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices,
             A -1-dimensional face of a Polyhedron in ZZ^3,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices]
        """
        if unlikely(self._trivial_faces):
            if self._trivial_faces == -1:
                raise StopIteration
            if self._trivial_faces in (2,4):  # Return the polyhedron.
                if self._trivial_faces == 2:
                    self._trivial_faces = 1  # Return the empty face next.
                else:
                    self._trivial_faces = -1  # The iterator is exhausted.
                equations = [eq.index() for eq in self.P.equation_generator()]
                return PolyhedronFace(self.P, range(self.P.n_Vrepresentation()), equations)
            else:  # Return the empty face.
                if self._trivial_faces == 1:
                    self._trivial_faces = 0  # Return the proper faces next.
                else:
                    self._trivial_faces = -1  # The iterator is exhausted.
                return PolyhedronFace(self.P, [], range(self.P.n_Hrepresentation()))

        self.next_dimension()
        if unlikely(self.structure.current_dimension > self.structure.highest_dimension):
            raise StopIteration
        return self.current()

    def current(self):
        r"""
        Retrieve the last value of :meth:`__next__`.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: it = P.face_generator()
            sage: _ = next(it), next(it)
            sage: next(it)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it.current()
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: next(it).ambient_V_indices() == it.current().ambient_V_indices()
            True
        """
        return combinatorial_face_to_polyhedral_face(self.P, FaceIterator_base.current(self))


# Nogil definitions of crucial functions.

cdef inline int next_dimension(iter_t structure, size_t parallelization_depth=0) nogil except -1:
    r"""
    See :meth:`FaceIterator.next_dimension`.

    ``parallelization_depth`` determines when to stop,
    e.g. if it is ``1`` it will stop after having yield all faces of a facet
    """
    cdef int max_dim = structure.highest_dimension - parallelization_depth
    structure.face_status = FaceStatus.NOT_INITIALIZED
    while (not next_face_loop(structure)) and (structure.current_dimension <= max_dim):
        sig_check()
    structure._index += 1
    return structure.current_dimension

cdef inline int next_face_loop(iter_t structure) nogil except -1:
    r"""
    See :meth:`FaceIterator.next_face_loop`.
    """
    if unlikely(structure.current_dimension == structure.dimension):
        # The function is not supposed to be called,
        # just prevent it from crashing.
        # Actually raising an error here results in a bad branch prediction.
        return -1

    # Getting ``[faces, n_faces, n_visited_all]`` according to dimension.
    cdef face_list_t* faces = &structure.new_faces[structure.current_dimension]
    cdef face_list_t* new_faces = &structure.new_faces[structure.current_dimension-1]
    cdef face_list_t* visited_all = &structure.visited_all[structure.current_dimension]
    cdef size_t n_faces = faces[0].n_faces

    if (structure.output_dimension > -2) and (structure.output_dimension != structure.current_dimension):
        # If only a specific dimension was requested (i.e. ``output_dimension > -2``),
        # then we will not yield faces in other dimension.
        structure.yet_to_visit = 0

    if structure.yet_to_visit:
        # Set ``face`` to the next face.
        structure.yet_to_visit -= 1
        structure.face[0] = faces[0].faces[structure.yet_to_visit][0]
        structure.face_status = FaceStatus.INITIALIZED
        return 1

    if structure.current_dimension <= structure.lowest_dimension:
        # We will not yield the empty face.
        # We will not yield below requested dimension.
        structure.current_dimension += 1
        return 0

    if n_faces <= 1:
        # There will be no more faces from intersections.
        structure.current_dimension += 1
        return 0

    if not structure.first_time[structure.current_dimension]:
        # In this case there exists ``faces[0].faces[n_faces]``, of which we
        # have visited all faces, but which was not added to
        # ``visited_all`` yet.

        if not faces[0].polyhedron_is_simple:
            # In case of a simple lattice, this step needs not to be applied:
            # Every element, except the lower bound, has a unique representation of coatoms in this case.
            # Hence, as the face is already removed from faces[0], any subfaces will not be visited.
            # (If we manually ignore subfaces, faces will still be added to visited_all).
            add_face_shallow(visited_all[0], faces[0].faces[n_faces])
    else:
        # Once we have visited all faces of ``faces[n_faces]``, we want
        # to add it to ``visited_all``.
        structure.first_time[structure.current_dimension] = False

    # We will visit the last face now.

    # Get the faces of codimension 1 contained in ``faces[n_faces-1]``,
    # which we have not yet visited.
    cdef size_t new_faces_counter

    new_faces_counter = get_next_level(
        faces[0], new_faces[0], visited_all[0])

    if new_faces_counter:
        # ``faces[n_faces]`` contains new faces.
        # We will visit them on next call, starting with codimension 1.

        # Setting the variables correctly for next call of ``next_face_loop``.
        structure.current_dimension -= 1
        structure.first_time[structure.current_dimension] = True
        structure.visited_all[structure.current_dimension][0] = visited_all[0][0]
        structure.yet_to_visit = new_faces_counter
        return 0
    else:
        # ``faces.faces[n_faces-1]`` contains no new faces.
        # Hence there is no need to add it to ``visited_all``.
        # NOTE:
        #     For the methods ``ignore_subfaces`` and ``ignore_supfaces``
        #     this step needs to be done, as ``faces.faces[n_faces-1]`` might
        #     have been added manually to ``visited_all``.
        #     So this step is required to respect boundaries of ``visited_all``.
        structure.first_time[structure.current_dimension] = True
        return 0

cdef inline size_t n_atom_rep(iter_t structure) nogil except -1:
    r"""
    See :meth:`FaceIterator.n_atom_rep`.
    """
    if structure.face_status:
        return face_len_atoms(structure.face)

    # The face was not initialized properly.
    raise LookupError("``FaceIterator`` does not point to a face")

# Parallel iteration over the faces.
# Currently only the f-vector is implemented, but slight
# modifications would allow collecting other information as well.

cdef struct parallel_f_s:
    # A structure carrying things that each thread should have exclusive access to.
    size_t* f_vector
    size_t* current_job_id

    # Keep track so that we can easily go from one job to the next.
    size_t* original_n_faces
    size_t* original_n_visited_all

ctypedef parallel_f_s parallel_f_t[1]

cdef int parallel_f_vector(iter_t* structures, size_t num_threads, size_t parallelization_depth, size_t *f_vector) except -1:
    """
    Compute the ``f_vector`` in parallel.

    INPUT:

    - ``structures`` -- one structure per thread

    - ``num_threads`` -- the number of threads to use

    - ``parallelization_depth`` -- the codimension at which the threads are released

    - ``f_vector`` -- where the ``f_vector`` is output
    """
    # One job per face of codimension ``parallelization_depth``.
    cdef size_t n_jobs = structures[0].n_coatoms ** parallelization_depth
    cdef size_t i
    cdef int j
    cdef int dim = structures[0].dimension
    f_vector[0] = 1         # Face iterator will only visit proper faces.
    f_vector[dim + 1] = 1   # Face iterator will only visit proper faces.
    if dim <= 0 or structures[0].n_coatoms == 0:
        # Iterator assumes at least one face and at least dimension 1.
        return 0

    if num_threads == 0:
        num_threads = 1

    cdef size_t thread_id
    cdef MemoryAllocator mem = MemoryAllocator()

    # Setting up for each thread some storage space.
    cdef parallel_f_t* parallel_structs = \
            <parallel_f_t*> mem.allocarray(num_threads, sizeof(parallel_f_t))

    for i in range(num_threads):
        # Partial f-vectors.
        parallel_structs[i].f_vector = \
                <size_t*> mem.calloc(dim+2, sizeof(size_t))
        parallel_structs[i].current_job_id = \
                <size_t*> mem.calloc(parallelization_depth+1, sizeof(size_t))

        # Keeping back of the original number of faces allows faster starting the next job.
        parallel_structs[i].original_n_faces = \
                <size_t*> mem.calloc(parallelization_depth+1, sizeof(size_t))
        parallel_structs[i].original_n_faces[0] = \
                structures[0].new_faces[dim - 1].n_faces

        parallel_structs[i].original_n_visited_all = \
                <size_t*> mem.calloc(parallelization_depth+1, sizeof(size_t))
        parallel_structs[i].original_n_visited_all[0] = \
                structures[0].visited_all[dim - 1].n_faces

    for i in prange(n_jobs, schedule='dynamic', chunksize=1,
                    num_threads=num_threads, nogil=True):
        _parallel_f_vector(structures[threadid()],
                           parallelization_depth,
                           parallel_structs[threadid()],
                           i)

    # Gather the results.
    for i in range(num_threads):
        for j in range(structures[0].dimension + 2):
            f_vector[j] += parallel_structs[i].f_vector[j]

cdef int _parallel_f_vector(iter_t structure, size_t parallelization_depth,
                            parallel_f_t parallel_struct, size_t job_id) nogil except -1:
    """
    Set up a job and then visit all faces.
    """
    cdef int max_dimension = structure.dimension - parallelization_depth
    cdef int d
    if prepare_face_iterator_for_partial_job(structure, parallelization_depth,
                                             parallel_struct, job_id):
        d = next_dimension(structure, parallelization_depth)
        while (d < max_dimension):
            parallel_struct.f_vector[d+1] += 1
            d = next_dimension(structure, parallelization_depth)

cdef inline int prepare_face_iterator_for_partial_job(
        iter_t structure, size_t parallelization_depth,
        parallel_f_t parallel_struct, size_t job_id) nogil except -1:
    """
    Set ``structure`` according to ``job_id``.

    ``job_id`` should be thought of as its digits with base ``structure.n_coatoms``
    padded to ``parallelization_depth``.

    The first digit determines which facet to visit.
    The next digit determines which facet of the facet should be visited.

    OUTPUT: ``1`` if the job exists and ``0`` otherwise.

    In addition, the first job treating a face will "visit" this face
    and increase the corresponding entry of the f-vector.
    """
    cdef int d
    cdef size_t current_depth
    if (not structure.first_time[structure.current_dimension]
            and structure.current_dimension == structure.dimension - parallelization_depth):
        # Act as if we had not visited faces in the last depth.
        # Set ``current_job_id[parallelization_depth - 1] = 0``.
        d = structure.dimension - parallelization_depth
        current_depth = parallelization_depth

        # Recover all faces.
        structure.new_faces[d].n_faces = \
                parallel_struct.original_n_faces[current_depth -1]
        structure.visited_all[d].n_faces = \
                parallel_struct.original_n_visited_all[current_depth -1]
        structure.first_time[d] = True
        structure.yet_to_visit = 0  # just to be on the safe side

        parallel_struct.current_job_id[current_depth -1] = 0

        # If the job does not exist, we will set the next value to ``-1``.
        if parallel_struct.original_n_faces[current_depth -1] == 0:
            parallel_struct.current_job_id[current_depth] = -1
        else:
            parallel_struct.current_job_id[current_depth] = 0

    cdef size_t n_coatoms = structure.n_coatoms
    cdef size_t job_id_c
    cdef size_t i
    cdef size_t diff
    cdef size_t new_faces_counter

    for current_depth in range(1, parallelization_depth + 1):
        d = structure.dimension - current_depth

        # Get the corresponding digit of ``job_id``.
        job_id_c = get_digit(job_id, current_depth - 1, parallelization_depth, n_coatoms)

        # Set ``current_job_id[current_depth - 1]`` to ``job_id_c`` if possible.

        if job_id_c != parallel_struct.current_job_id[current_depth - 1]:
            # Set ``current_job_id[current_depth -1] = 0``.
            structure.current_dimension = d
            structure.new_faces[d].n_faces = parallel_struct.original_n_faces[current_depth -1]
            structure.visited_all[d].n_faces = parallel_struct.original_n_visited_all[current_depth -1]
            parallel_struct.current_job_id[current_depth -1] = 0

            # If the job does not exist, we will set the next value to ``-1``.
            if parallel_struct.original_n_faces[current_depth -1] == 0:
                parallel_struct.current_job_id[current_depth] = -1
            else:
                parallel_struct.current_job_id[current_depth] = 0

            structure.first_time[d] = True
            structure.yet_to_visit = 0  # just to be on the safe side

        if parallel_struct.current_job_id[current_depth] == -1:
            # The job does not exist.
            return 0

        if job_id_c > parallel_struct.current_job_id[current_depth -1]:
            if job_id_c >= structure.new_faces[d].n_faces:
                # The job does not exist.
                return 0

            for i in range(job_id_c):
                # Fast forwarding the jobs.

                if not structure.new_faces[d].polyhedron_is_simple:
                    # In case of a simple lattice, this step needs not to be applied:
                    # Every element, except the lower bound, has a unique representation of coatoms in this case.
                    # Hence, as the face is already removed from faces[0], any subfaces will not be visited.
                    # (If we manually ignore subfaces, faces will still be added to visited_all).
                    add_face_shallow(structure.visited_all[d], structure.new_faces[d].faces[structure.new_faces[d].n_faces -1])

                structure.new_faces[d].n_faces -= 1

            parallel_struct.current_job_id[current_depth -1] = job_id_c

        # Apparently the face exists. We add it to the f-vector, if it is the very first job for the face.
        if job_id == 0 or get_digit(job_id -1, current_depth -1, parallelization_depth, n_coatoms) != job_id_c:
            # Visit ``structure.new_faces[d].faces[structure.new_faces[d].n_faces - 1]
            parallel_struct.f_vector[d + 1] += 1

        if structure.current_dimension == d:
            structure.yet_to_visit = 0

            if structure.new_faces[d].n_faces == 0:
                # The job will not exist.
                parallel_struct.current_job_id[current_depth] = -1
                return 0

            new_faces_counter = get_next_level(
                structure.new_faces[d], structure.new_faces[d-1], structure.visited_all[d])

            if new_faces_counter:
                # Setting the variables correctly, for the next call.
                structure.current_dimension -= 1
                structure.first_time[d-1] = True
                structure.visited_all[d-1][0] = structure.visited_all[d][0]
                structure.yet_to_visit = new_faces_counter
                for i in range(current_depth, parallelization_depth + 1):
                    parallel_struct.current_job_id[i] = 0
                parallel_struct.original_n_faces[current_depth] = new_faces_counter
                parallel_struct.original_n_visited_all[current_depth] = structure.visited_all[d-1].n_faces
            else:
                # The job does not exist.
                parallel_struct.current_job_id[current_depth] = -1
                return 0

    if structure.current_dimension != structure.dimension - parallelization_depth - 1:
        return 0

    return 1

cdef inline size_t get_digit(size_t job_id, size_t pos, size_t padto, size_t base) nogil:
    """
    Get the digit ``pos`` of ``job_id`` with base ``base``
    padding the number of digits to ``pad_to``.

    Digits are enumerated started with ``0``.
    """
    # Remove the last ``parallelization_depth - pos - 1`` digits.
    cdef size_t current_output = job_id / base ** (padto - pos - 1)

    # Remove all digits before our current digit, which is digit ``pos``.
    return current_output % base
