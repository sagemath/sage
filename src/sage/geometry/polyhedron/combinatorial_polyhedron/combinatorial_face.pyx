r"""
Combinatorial face of a polyhedron

This module provides the combinatorial type of a polyhedral face.

.. SEEALSO::

    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.base`,
    :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`.

EXAMPLES:

Obtain a face from a face iterator::

    sage: P = polytopes.cube()
    sage: C = CombinatorialPolyhedron(P)
    sage: it = C.face_generator()
    sage: face = next(it); face
    A 2-dimensional face of a 3-dimensional combinatorial polyhedron

Obtain a face from a face lattice index::

    sage: P = polytopes.simplex(2)
    sage: C = CombinatorialPolyhedron(P)
    sage: sorted(C.face_lattice()._elements)                                            # needs sage.combinat
    [0, 1, 2, 3, 4, 5, 6, 7]
    sage: face = C.face_by_face_lattice_index(0); face
    A -1-dimensional face of a 2-dimensional combinatorial polyhedron

Obtain further information regarding a face::

    sage: P = polytopes.octahedron()
    sage: C = CombinatorialPolyhedron(P)
    sage: it = C.face_generator(2)
    sage: face = next(it); face
    A 2-dimensional face of a 3-dimensional combinatorial polyhedron
    sage: face.ambient_Vrepresentation()
    (A vertex at (0, 0, 1), A vertex at (0, 1, 0), A vertex at (1, 0, 0))
    sage: face.n_ambient_Vrepresentation()
    3
    sage: face.ambient_H_indices()
    (5,)
    sage: face.dimension()
    2
    sage: face.ambient_dimension()
    3

.. SEEALSO::

    :class:`sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

AUTHOR:

- Jonathan Kliem (2019-05)
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

from cysignals.memory           cimport check_allocarray, sig_free

import numbers
from sage.rings.integer         cimport smallInteger
from .conversions               cimport bit_rep_to_Vrep_list
from .base                      cimport CombinatorialPolyhedron
from .face_iterator             cimport FaceIterator_base, FaceStatus
from .polyhedron_face_lattice   cimport PolyhedronFaceLattice
from .face_data_structure       cimport face_len_atoms, face_init, face_free, face_copy, face_issubset
from .face_list_data_structure  cimport bit_rep_to_coatom_rep
from .list_of_faces              cimport face_as_combinatorial_polyhedron


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


cdef class CombinatorialFace(SageObject):
    r"""
    A class of the combinatorial type of a polyhedral face.

    EXAMPLES:

    Obtain a combinatorial face from a face iterator::

        sage: P = polytopes.cyclic_polytope(5,8)
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_generator()
        sage: next(it)
        A 0-dimensional face of a 5-dimensional combinatorial polyhedron

    Obtain a combinatorial face from an index of the face lattice::

        sage: F = C.face_lattice()                                                      # needs sage.combinat
        sage: F._elements[3]                                                            # needs sage.combinat
        34
        sage: C.face_by_face_lattice_index(29)
        A 1-dimensional face of a 5-dimensional combinatorial polyhedron

    Obtain the dimension of a combinatorial face::

        sage: face = next(it)
        sage: face.dimension()
        0

    The dimension of the polyhedron::

        sage: face.ambient_dimension()
        5

    The Vrepresentation::

        sage: face.ambient_Vrepresentation()
        (A vertex at (6, 36, 216, 1296, 7776),)
        sage: face.ambient_V_indices()
        (6,)
        sage: face.n_ambient_Vrepresentation()
        1

    The Hrepresentation::

        sage: face.ambient_Hrepresentation()
        (An inequality (60, -112, 65, -14, 1) x + 0 >= 0,
         An inequality (180, -216, 91, -16, 1) x + 0 >= 0,
         An inequality (360, -342, 119, -18, 1) x + 0 >= 0,
         An inequality (840, -638, 179, -22, 1) x + 0 >= 0,
         An inequality (-2754, 1175, -245, 25, -1) x + 2520 >= 0,
         An inequality (504, -450, 145, -20, 1) x + 0 >= 0,
         An inequality (-1692, 853, -203, 23, -1) x + 1260 >= 0,
         An inequality (252, -288, 113, -18, 1) x + 0 >= 0,
         An inequality (-844, 567, -163, 21, -1) x + 420 >= 0,
         An inequality (84, -152, 83, -16, 1) x + 0 >= 0,
         An inequality (-210, 317, -125, 19, -1) x + 0 >= 0)
        sage: face.ambient_H_indices()
        (3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19)
        sage: face.n_ambient_Hrepresentation()
        11
    """
    def __cinit__(self, data, dimension=None, index=None):
        r"""
        Initialize :class:`CombinatorialFace`.

        See :class:`CombinatorialFace`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],
            ....: [0,2,3],[1,2,3]])
            sage: it = C.face_generator()
            sage: next(it)     # indirect doctest
            A 2-dimensional face of a 3-dimensional combinatorial polyhedron

            sage: C.face_by_face_lattice_index(0)
            A -1-dimensional face of a 3-dimensional combinatorial polyhedron

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace).run()
        """
        # Note that all values are set to zero at the time ``__cinit__`` is called:
        # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#initialisation-methods
        # In particular, ``__dealloc__`` will not do harm in this case.

        cdef FaceIterator_base it
        cdef PolyhedronFaceLattice all_faces

        if isinstance(data, FaceIterator_base):
            assert dimension is None and index is None, "dimension and index must be ``None``, when providing a face iterator"

            # Copy data from FaceIterator.
            it = data
            self._dual              = it.dual
            self.atoms              = it.atoms
            self.coatoms            = it.coatoms

            if it.structure.face_status == FaceStatus.NOT_INITIALIZED:
                raise LookupError("face iterator not set to a face")

            face_init(self.face, self.coatoms.n_atoms(), self.coatoms.n_coatoms())
            face_copy(self.face, it.structure.face)

            self._dimension         = it.structure.current_dimension
            self._ambient_dimension = it.structure.dimension
            self._ambient_Vrep      = it._Vrep
            self._ambient_facets    = it._facet_names
            self._n_ambient_facets  = it._n_facets
            self._equations         = it._equations
            self._n_equations       = it._n_equations
            self._hash_index        = it.structure._index
            self._ambient_bounded   = it._bounded

            self._initialized_from_face_lattice = False

        elif isinstance(data, PolyhedronFaceLattice):
            all_faces = data
            assert isinstance(dimension, numbers.Integral), "dimension must be an integer"
            assert isinstance(index, numbers.Integral), "index must be an integer"
            assert -1 <= dimension <= all_faces.dimension, "dimension must be a face dimension of the polyhedron"
            assert 0 <= index < all_faces.f_vector[dimension + 1], "index is out of range"

            # Copy data from PolyhedronFaceLattice.
            self._dual              = all_faces.dual
            self.atoms              = all_faces.atoms
            self.coatoms            = all_faces.coatoms

            face_init(self.face, self.coatoms.n_atoms(), self.coatoms.n_coatoms())
            face_copy(self.face, all_faces.faces[dimension+1].faces[index])

            self._dimension         = dimension
            self._ambient_dimension = all_faces.dimension
            self._ambient_Vrep      = all_faces._Vrep
            self._ambient_facets    = all_faces._facet_names
            self._equations         = all_faces._equations
            self._n_equations       = len(self._equations) if self._equations else 0
            if self._dual:
                self._n_ambient_facets = self.atoms.n_faces()
            else:
                self._n_ambient_facets = self.coatoms.n_faces()
            self._ambient_bounded   = all_faces._bounded

            self._initialized_from_face_lattice = True

            self._hash_index = index
            for i in range(-1,dimension):
                self._hash_index += all_faces.f_vector[i+1]

            # Add the complete ``f-vector`` to the hash index,
            # such that hash values obtained by an iterator or by the face lattice
            # do not collide.
            for i in range(-1,self._ambient_dimension+1):
                self._hash_index += all_faces.f_vector[i+1]
        else:
            raise ValueError("data must be face iterator or a list of all faces")

        if self._dual:
            # Reverse the hash index in dual mode to respect inclusion of faces.
            self._hash_index = -self._hash_index - 1

    def __dealloc__(self):
        r"""
        TESTS::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face import CombinatorialFace
            sage: CombinatorialFace(2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: data must be face iterator or a list of all faces
        """
        face_free(self.face)
        sig_free(self.atom_rep)
        sig_free(self.coatom_rep)

    def _repr_(self):
        r"""
        Return a description of the combinatorial face.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(6)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(dimension=3, algorithm='primal')
            sage: face = next(it)
            sage: face.__repr__()
            'A 3-dimensional face of a 5-dimensional combinatorial polyhedron'
            sage: it = C.face_generator(dimension=3, algorithm='dual')
            sage: face = next(it)
            sage: face.__repr__()
            'A 3-dimensional face of a 5-dimensional combinatorial polyhedron'
        """
        return "A {}-dimensional face of a {}-dimensional combinatorial polyhedron"\
                .format(self.dimension(), self.ambient_dimension())

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: face1 = loads(face.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __hash__(self):
        r"""
        Return an index for the face.

        This is constructed such that for faces `F,G` constructed in the same manner (same face iterator or face lattice)
        it holds that `F` contained in `G` implies ``hash(F) < hash(G)``.

        If the face was constructed from a :class:`sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`,
        then this is the index of the occurrence in the iterator.
        In dual mode this value is then deducted from the maximal value of ``size_t``.

        If the face was constructed from
        :meth:`sage:geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron.face_by_face_lattice_index`,
        then this is the total number of faces plus the index in the level set plus the number of lower dimensional faces
        (or higher dimensional faces in dual mode).
        In dual mode this value is then deducted from the maximal value of ``size_t``.

        EXAMPLES::

            sage: P = polytopes.simplex(2)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: [hash(face) for face in it]
            [1, 2, 3, 4, 5, 6]

        TESTS::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: G = F.relabel(C.face_by_face_lattice_index)

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()                                                  # needs sage.combinat
            sage: G = F.relabel(C.face_by_face_lattice_index)                           # needs sage.combinat
        """
        return self._hash_index

    def __lt__(self, other):
        r"""
        Compare faces of the same polyhedron.

        This is a helper function.
        In order to construct a Hasse diagram (a digraph) with combinatorial faces,
        we must define some order relation that is compatible with the Hasse diagram.

        Any order relation compatible with ordering by dimension is suitable.
        We use :meth:`__hash__` to define the relation.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F1 = C.face_by_face_lattice_index(0)
            sage: F2 = C.face_by_face_lattice_index(1)
            sage: F1 < F2
            True
            sage: for i,j in Combinations(range(28), 2):
            ....:     F1 = C.face_by_face_lattice_index(i)
            ....:     F2 = C.face_by_face_lattice_index(j)
            ....:     if F1.dim() != F2.dim():
            ....:          assert (F1.dim() < F2.dim()) == (F1 < F2)

            sage: P = polytopes.cross_polytope(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: F1 = C.face_by_face_lattice_index(0)
            sage: F2 = C.face_by_face_lattice_index(1)
            sage: F1 < F2
            True
            sage: for i,j in Combinations(range(28), 2):
            ....:     F1 = C.face_by_face_lattice_index(i)
            ....:     F2 = C.face_by_face_lattice_index(j)
            ....:     if F1.dim() != F2.dim():
            ....:          assert (F1.dim() < F2.dim()) == (F1 < F2)
        """
        cdef CombinatorialFace other_face
        if isinstance(other, CombinatorialFace):
            other_face = other
            if (self._initialized_from_face_lattice == other_face._initialized_from_face_lattice and
                    self.atoms is other_face.atoms):
                # They are faces of the same polyhedron obtained in the same way.
                return hash(self) < hash(other)

    def is_subface(self, CombinatorialFace other):
        r"""
        Return whether ``self`` is contained in ``other``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: face.ambient_V_indices()
            (0, 3, 4, 5)
            sage: face2 = next(it)
            sage: face2.ambient_V_indices()
            (0, 1, 5, 6)
            sage: face.is_subface(face2)
            False
            sage: face2.is_subface(face)
            False
            sage: it.only_subfaces()
            sage: face3 = next(it)
            sage: face3.ambient_V_indices()
            (0, 5)
            sage: face3.is_subface(face2)
            True
            sage: face3.is_subface(face)
            True

        Works for faces of the same combinatorial polyhedron;
        also from different iterators::

            sage: it = C.face_generator(algorithm='dual')
            sage: v7 = next(it); v7.ambient_V_indices()
            (7,)
            sage: v6 = next(it); v6.ambient_V_indices()
            (6,)
            sage: v5 = next(it); v5.ambient_V_indices()
            (5,)
            sage: face.ambient_V_indices()
            (0, 3, 4, 5)
            sage: face.is_subface(v7)
            False
            sage: v7.is_subface(face)
            False
            sage: v6.is_subface(face)
            False
            sage: v5.is_subface(face)
            True
            sage: face2.ambient_V_indices()
            (0, 1, 5, 6)
            sage: face2.is_subface(v7)
            False
            sage: v7.is_subface(face2)
            False
            sage: v6.is_subface(face2)
            True
            sage: v5.is_subface(face2)
            True

        Only implemented for faces of the same combinatorial polyhedron::

            sage: P1 = polytopes.cube()
            sage: C1 = P1.combinatorial_polyhedron()
            sage: it = C1.face_generator()
            sage: other_face = next(it)
            sage: other_face.ambient_V_indices()
            (0, 3, 4, 5)
            sage: face.ambient_V_indices()
            (0, 3, 4, 5)
            sage: C is C1
            False
            sage: face.is_subface(other_face)
            Traceback (most recent call last):
            ...
            NotImplementedError: is_subface only implemented for faces of the same polyhedron
        """
        cdef size_t length_self, length_other, counter_self, counter_other
        cdef size_t* self_v_indices
        cdef size_t* other_v_indices

        if self._dual == other._dual:
            if self.atoms is other.atoms:
                if not self._dual:
                    return face_issubset(self.face, other.face)
                else:
                    return face_issubset(other.face, self.face)
            else:
                raise NotImplementedError("is_subface only implemented for faces of the same polyhedron")
        else:
            if self.atoms is other.coatoms:
                if self.dimension() > other.dimension():
                    return False
                if self._dual:
                    length_self = self.set_coatom_rep()
                    self_v_indices = self.coatom_rep
                    length_other = other.set_atom_rep()
                    other_v_indices = other.atom_rep
                else:
                    length_self = self.set_atom_rep()
                    self_v_indices = self.atom_rep
                    length_other = other.set_coatom_rep()
                    other_v_indices = other.coatom_rep
                if length_self > length_other:
                    return False

                # Check if every element in self_v_indices is contained in other_v_indices.
                counter_self = 0
                counter_other = 0
                while counter_self < length_self and counter_other < length_other:
                    if self_v_indices[counter_self] > other_v_indices[counter_other]:
                        counter_other += 1
                    elif self_v_indices[counter_self] == other_v_indices[counter_other]:
                        counter_self += 1
                        counter_other += 1
                    else:
                        return False
                return counter_self == length_self
            else:
                raise NotImplementedError("is_subface only implemented for faces of the same polyhedron")

    cpdef dimension(self):
        r"""
        Return the dimension of the face.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.associahedron(['A', 3])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: face.dimension()
            2

        ``dim`` is an alias::

            sage: face.dim()                                                            # needs sage.combinat
            2
        """
        if self._dual:
            return smallInteger(self._ambient_dimension - self._dimension - 1)
        else:
            return smallInteger(self._dimension)

    dim = dimension

    def ambient_dimension(self):
        r"""
        Return the dimension of the polyhedron.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: face = next(it)
            sage: face.ambient_dimension()
            3
        """
        return smallInteger(self._ambient_dimension)

    def ambient_Vrepresentation(self):
        r"""
        Return the Vrepresentation objects of the ambient polyhedron
        defining the face.

        It consists of the vertices/rays/lines
        that face contains.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(dimension=2)
            sage: face = next(it)
            sage: face.ambient_Vrepresentation()
            (A vertex at (1, 3, 2, 5, 4),
             A vertex at (2, 3, 1, 5, 4),
             A vertex at (3, 1, 2, 5, 4),
             A vertex at (3, 2, 1, 5, 4),
             A vertex at (2, 1, 3, 5, 4),
             A vertex at (1, 2, 3, 5, 4))
            sage: face = next(it)
            sage: face.ambient_Vrepresentation()
            (A vertex at (2, 1, 4, 5, 3),
             A vertex at (3, 2, 4, 5, 1),
             A vertex at (3, 1, 4, 5, 2),
             A vertex at (1, 3, 4, 5, 2),
             A vertex at (1, 2, 4, 5, 3),
             A vertex at (2, 3, 4, 5, 1))

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_generator()
            sage: for face in it: (face.dimension(), face.ambient_Vrepresentation())
            (2, (1, 2, 3))
            (2, (0, 2, 3))
            (2, (0, 1, 3))
            (2, (0, 1, 2))
            (1, (2, 3))
            (1, (1, 3))
            (1, (1, 2))
            (0, (3,))
            (0, (2,))
            (0, (1,))
            (1, (0, 3))
            (1, (0, 2))
            (0, (0,))
            (1, (0, 1))

        .. SEEALSO::

            :meth:`ambient_V_indices`.
        """
        if not self._ambient_Vrep:
            # There are no names, so we return indices instead.
            return self.ambient_V_indices()
        cdef size_t length
        if self._dual:
            # if dual, the Vrepresentation corresponds to the coatom-representation
            length = self.set_coatom_rep()
            return tuple(self._ambient_Vrep[self.coatom_rep[i]]
                         for i in range(length))
        else:
            # if not dual, the Vrepresentation corresponds to the atom-representation
            length = self.set_atom_rep()
            return tuple(self._ambient_Vrep[self.atom_rep[i]]
                         for i in range(length))

    def ambient_V_indices(self):
        r"""
        Return the indices of the Vrepresentation
        objects of the ambient polyhedron defining the face.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(dimension=2)
            sage: face = next(it)
            sage: next(it).ambient_V_indices()
            (32, 91, 92, 93, 94, 95)
            sage: next(it).ambient_V_indices()
            (32, 89, 90, 94)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_generator()
            sage: for face in it: (face.dimension(), face.ambient_V_indices())
            (2, (1, 2, 3))
            (2, (0, 2, 3))
            (2, (0, 1, 3))
            (2, (0, 1, 2))
            (1, (2, 3))
            (1, (1, 3))
            (1, (1, 2))
            (0, (3,))
            (0, (2,))
            (0, (1,))
            (1, (0, 3))
            (1, (0, 2))
            (0, (0,))
            (1, (0, 1))

        .. SEEALSO::

            :meth:`ambient_Vrepresentation`.
        """
        cdef size_t length
        if self._dual:
            # if dual, the Vrepresentation corresponds to the coatom-representation
            length = self.set_coatom_rep()
            return tuple(smallInteger(self.coatom_rep[i])
                         for i in range(length))
        else:
            # if not dual, the Vrepresentation corresponds to the atom-representation
            length = self.set_atom_rep()
            return tuple(smallInteger(self.atom_rep[i])
                         for i in range(length))

    def n_ambient_Vrepresentation(self):
        r"""
        Return the length of the :meth:`CombinatorialFace.ambient_V_indices`.

        Might be faster than using ``len``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: all(face.n_ambient_Vrepresentation() == len(face.ambient_Vrepresentation()) for face in it)
            True

        TESTS::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: face = next(it)
        """
        if self._dual:
            return smallInteger(self.set_coatom_rep())
        else:
            return smallInteger(self.n_atom_rep())

    def ambient_Hrepresentation(self):
        r"""
        Return the Hrepresentation objects of the ambient polyhedron
        defining the face.

        It consists of the facets/inequalities that contain the face
        and the equations defining the ambient polyhedron.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(2)
            sage: next(it).ambient_Hrepresentation()
            (An inequality (1, 1, 1, 0, 0) x - 6 >= 0,
             An inequality (0, 0, 0, -1, 0) x + 5 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: next(it).ambient_Hrepresentation()
            (An inequality (0, 0, -1, -1, 0) x + 9 >= 0,
             An inequality (0, 0, 0, -1, 0) x + 5 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: next(it).ambient_Hrepresentation()
            (An inequality (-20, 29, -10, 1) x + 0 >= 0,
             An inequality (60, -47, 12, -1) x + 0 >= 0,
             An inequality (30, -31, 10, -1) x + 0 >= 0,
             An inequality (10, -17, 8, -1) x + 0 >= 0,
             An inequality (-154, 71, -14, 1) x + 120 >= 0,
             An inequality (-78, 49, -12, 1) x + 40 >= 0)
            sage: next(it).ambient_Hrepresentation()
            (An inequality (-50, 35, -10, 1) x + 24 >= 0,
             An inequality (-12, 19, -8, 1) x + 0 >= 0,
             An inequality (-20, 29, -10, 1) x + 0 >= 0,
             An inequality (60, -47, 12, -1) x + 0 >= 0,
             An inequality (-154, 71, -14, 1) x + 120 >= 0,
             An inequality (-78, 49, -12, 1) x + 40 >= 0)

        .. SEEALSO::

            :meth:`ambient_H_indices`.
        """
        if not self._ambient_facets:
            # There are no names, so we return indices instead.
            return self.ambient_H_indices()
        cdef size_t length
        if not self._dual:
            # if not dual, the facet-representation corresponds to the coatom-representation
            length = self.set_coatom_rep()  # fill self.coatom_repr_face
            return tuple(self._ambient_facets[self.coatom_rep[i]]
                         for i in range(length)) + self._equations
        else:
            # if dual, the facet-representation corresponds to the atom-representation
            length = self.set_atom_rep()  # fill self.atom_repr_face
            return tuple(self._ambient_facets[self.atom_rep[i]]
                         for i in range(length)) + self._equations

    def ambient_H_indices(self, add_equations=True):
        r"""
        Return the indices of the Hrepresentation objects
        of the ambient polyhedron defining the face.

        INPUT:

        - ``add_equations`` -- boolean (default: ``True``); whether or not to include the equations

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(2)
            sage: face = next(it)
            sage: face.ambient_H_indices(add_equations=False)
            (28, 29)
            sage: face2 = next(it)
            sage: face2.ambient_H_indices(add_equations=False)
            (25, 29)

        Add the indices of the equation::

            sage: face.ambient_H_indices(add_equations=True)                            # needs sage.combinat
            (28, 29, 30)
            sage: face2.ambient_H_indices(add_equations=True)                           # needs sage.combinat
            (25, 29, 30)

        Another example::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: _ = next(it); _ = next(it)
            sage: next(it).ambient_H_indices()
            (0, 1, 2, 4, 5, 7)
            sage: next(it).ambient_H_indices()
            (0, 1, 5, 6, 7, 8)
            sage: next(it).ambient_H_indices()
            (0, 1, 2, 3, 6, 8)
            sage: [next(it).dimension() for _ in range(2)]
            [0, 1]
            sage: face = next(it)
            sage: face.ambient_H_indices()
            (4, 5, 7)

        .. SEEALSO::

            :meth:`ambient_Hrepresentation`.
        """
        cdef size_t length, i
        cdef tuple equations

        if add_equations and self._equations:
            equations = tuple(smallInteger(i)
                              for i in range(self._n_ambient_facets,
                                             self._n_ambient_facets + self._n_equations))
        else:
            equations = ()

        if not self._dual:
            # if not dual, the facet-representation corresponds to the coatom-representation
            length = self.set_coatom_rep()  # fill self.coatom_repr_face
            return tuple(smallInteger(self.coatom_rep[i])
                         for i in range(length)) + equations
        else:
            # if dual, the facet-representation corresponds to the atom-representation
            length = self.set_atom_rep()  # fill self.atom_repr_face
            return tuple(smallInteger(self.atom_rep[i])
                         for i in range(length)) + equations

    def n_ambient_Hrepresentation(self, add_equations=True):
        r"""
        Return the length of the :meth:`CombinatorialFace.ambient_H_indices`.

        Might be faster than then using ``len``.

        INPUT:

        - ``add_equations`` -- boolean (default: ``True``); whether or not to count the equations

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: all(face.n_ambient_Hrepresentation() == len(face.ambient_Hrepresentation()) for face in it)
            True

        Specifying whether to count the equations or not::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(2)
            sage: f = next(it)
            sage: f.n_ambient_Hrepresentation(add_equations=True)
            3
            sage: f.n_ambient_Hrepresentation(add_equations=False)
            2

        TESTS::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: face = next(it)
        """
        cdef size_t n_equations = self._n_equations if add_equations else 0
        if not self._dual:
            return smallInteger(self.set_coatom_rep() + n_equations)
        else:
            return smallInteger(self.n_atom_rep() + n_equations)

    def as_combinatorial_polyhedron(self, quotient=False):
        r"""
        Return ``self`` as combinatorial polyhedron.

        If ``quotient`` is ``True``, return the quotient of the
        polyhedron by ``self``.
        Let ``G`` be the face corresponding to ``self`` in the dual/polar polytope.
        The ``quotient`` is the dual/polar of ``G``.

        Let `[\hat{0}, \hat{1}]` be the face lattice of the ambient polyhedron
        and `F` be ``self`` as element of the face lattice.
        The face lattice of ``self`` as polyhedron corresponds to
        `[\hat{0}, F]` and the face lattice of the quotient by ``self``
        corresponds to `[F, \hat{1}]`.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(7,11)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(4)
            sage: f = next(it); f
            A 4-dimensional face of a 7-dimensional combinatorial polyhedron
            sage: F = f.as_combinatorial_polyhedron(); F
            A 4-dimensional combinatorial polyhedron with 5 facets
            sage: F.f_vector()
            (1, 5, 10, 10, 5, 1)
            sage: F_alt = polytopes.cyclic_polytope(4,5).combinatorial_polyhedron()
            sage: F_alt.vertex_facet_graph().is_isomorphic(F.vertex_facet_graph())      # needs sage.graphs
            True

        Obtaining the quotient::

            sage: Q = f.as_combinatorial_polyhedron(quotient=True); Q
            A 2-dimensional combinatorial polyhedron with 6 facets
            sage: Q
            A 2-dimensional combinatorial polyhedron with 6 facets
            sage: Q.f_vector()
            (1, 6, 6, 1)

        The Vrepresentation of the face as polyhedron is given by the
        ambient Vrepresentation of the face in that order::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(2)
            sage: f = next(it)
            sage: F = f.as_combinatorial_polyhedron()
            sage: C.Vrepresentation()
            (A vertex at (1, -1, -1),
            A vertex at (1, 1, -1),
            A vertex at (1, 1, 1),
            A vertex at (1, -1, 1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, -1, -1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1))
            sage: f.ambient_Vrepresentation()
            (A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, -1, -1))
            sage: F.Vrepresentation()
            (0, 1, 2, 3)

        To obtain the facets of the face as polyhedron,
        we compute the meet of each facet with the face.
        The first representative of each element strictly
        contained in the face is kept::

            sage: C.facets(names=False)
            ((0, 1, 2, 3),
             (1, 2, 6, 7),
             (2, 3, 4, 7),
             (4, 5, 6, 7),
             (0, 1, 5, 6),
             (0, 3, 4, 5))
            sage: F.facets(names=False)
            ((0, 1), (1, 2), (2, 3), (0, 3))

        The Hrepresentation of the quotient by the face is given by the
        ambient Hrepresentation of the face in that order::

            sage: it = C.face_generator(1)
            sage: f = next(it)
            sage: Q = f.as_combinatorial_polyhedron(quotient=True)
            sage: C.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0)
            sage: f.ambient_Hrepresentation()
            (An inequality (0, 0, 1) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0)
            sage: Q.Hrepresentation()
            (0, 1)

        To obtain the vertices of the face as polyhedron,
        we compute the join of each vertex with the face.
        The first representative of each element strictly
        containing the face is kept::

            sage: [g.ambient_H_indices() for g in C.face_generator(0)]
            [(3, 4, 5),
            (0, 4, 5),
            (2, 3, 5),
            (0, 2, 5),
            (1, 3, 4),
            (0, 1, 4),
            (1, 2, 3),
            (0, 1, 2)]
            sage: [g.ambient_H_indices() for g in Q.face_generator(0)]
            [(1,), (0,)]

        The method is not implemented for unbounded polyhedra::

            sage: P = Polyhedron(rays=[[0,1]])*polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(2)
            sage: f = next(it)
            sage: f.as_combinatorial_polyhedron()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for bounded polyhedra

        REFERENCES:

            For more information, see Exercise 2.9 of [Zie2007]_.

        .. NOTE::

            This method is tested in
            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base._test_combinatorial_face_as_combinatorial_polyhedron`.
        """
        if not self._ambient_bounded:
            raise NotImplementedError("only implemented for bounded polyhedra")

        cdef ListOfFaces facets = self.atoms if self._dual else self.coatoms
        cdef ListOfFaces Vrep = self.atoms if not self._dual else self.coatoms

        if not quotient:
            return CombinatorialPolyhedron(face_as_combinatorial_polyhedron(facets, Vrep, self.face, self._dual))
        else:
            # We run ``face_as_combinatorial_polyhedron`` for the dual setting.

            # We then interchange the output of it, to obtain the quotient.
            new_Vrep, new_facets = face_as_combinatorial_polyhedron(Vrep, facets, self.face, not self._dual)
            return CombinatorialPolyhedron((new_facets, new_Vrep))

    cdef size_t n_atom_rep(self) except -1:
        r"""
        Compute the number of atoms in the current face by counting the
        number of set bits.
        """
        if self.atom_rep is not NULL:
            return self._n_atom_rep
        return face_len_atoms(self.face)

    cdef size_t set_coatom_rep(self) except -1:
        r"""
        Set ``coatom_rep`` to be the coatom-representation of the current face.
        Return its length.
        """
        if self.coatom_rep is NULL:
            self.coatom_rep = <size_t *> check_allocarray(self.coatoms.n_faces(), sizeof(size_t))
            self._n_coatom_rep = bit_rep_to_coatom_rep(self.face, self.coatoms.data, self.coatom_rep)
        return self._n_coatom_rep

    cdef size_t set_atom_rep(self) except -1:
        r"""
        Set ``atom_rep`` to be the atom-representation of the current face.
        Return its length.
        """
        if self.atom_rep is NULL:
            self.atom_rep = <size_t *> check_allocarray(self.coatoms.n_atoms(), sizeof(size_t))
            self._n_atom_rep = bit_rep_to_Vrep_list(self.face, self.atom_rep)
        return self._n_atom_rep
