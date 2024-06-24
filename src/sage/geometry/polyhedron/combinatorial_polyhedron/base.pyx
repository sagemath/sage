r"""
Combinatorial polyhedron

This module gathers algorithms for polyhedra that only depend on the
vertex-facet incidences and that are called combinatorial polyhedron.
The main class is :class:`CombinatorialPolyhedron`. Most importantly,
this class allows to iterate quickly through the faces (possibly
of given dimension) via the :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator` object. The :class:`CombinatorialPolyhedron`
uses this iterator to quickly generate the f-vector,  the edges,
the ridges and the face lattice.

Terminology used in this module:

- Vrep                  -- ``[vertices, rays, lines]`` of the polyhedron
- Hrep                  -- inequalities and equations of the polyhedron
- Facets                -- facets of the polyhedron
- Vrepresentation       -- represents a face by the list of Vrep it contains
- Hrepresentation       -- represents a face by a list of Hrep it is contained in
- bit representation    -- represents incidences as bitset, where each bit
  represents one incidence. There might be trailing zeros, to fit alignment
  requirements.  In most instances, faces are represented by the bit
  representation, where each bit corresponds to a Vrep or facet. Thus a bit
  representation can either be a Vrep or facet representation depending on
  context.

EXAMPLES:

Construction::

    sage: P = polytopes.hypercube(4)
    sage: C = CombinatorialPolyhedron(P); C
    A 4-dimensional combinatorial polyhedron with 8 facets

Obtaining edges and ridges::

    sage: C.edges()[:2]
    ((A vertex at (1, -1, -1, -1), A vertex at (-1, -1, -1, -1)),
     (A vertex at (-1, -1, -1, 1), A vertex at (-1, -1, -1, -1)))
    sage: C.edges(names=False)[:2]
    ((6, 15), (14, 15))

    sage: C.ridges()[:2]
    ((An inequality (0, 0, 1, 0) x + 1 >= 0,
      An inequality (0, 1, 0, 0) x + 1 >= 0),
     (An inequality (0, 0, 0, 1) x + 1 >= 0,
      An inequality (0, 1, 0, 0) x + 1 >= 0))
    sage: C.ridges(names=False)[:2]
    ((6, 7), (5, 7))

Vertex-graph and facet-graph::

    sage: C.vertex_graph()                                                              # needs sage.graphs
    Graph on 16 vertices
    sage: C.facet_graph()                                                               # needs sage.graphs
    Graph on 8 vertices

Face lattice::

    sage: C.face_lattice()                                                              # needs sage.combinat
    Finite lattice containing 82 elements

Face iterator::

    sage: C.face_generator()
    Iterator over the proper faces of a 4-dimensional combinatorial polyhedron

    sage: C.face_generator(2)
    Iterator over the 2-faces of a 4-dimensional combinatorial polyhedron

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

import numbers
from memory_allocator cimport MemoryAllocator
from cysignals.memory cimport check_calloc, sig_free

from sage.graphs.graph              import Graph
from sage.geometry.polyhedron.base  import Polyhedron_base
from sage.geometry.lattice_polytope import LatticePolytopeClass
from sage.geometry.cone             import ConvexRationalPolyhedralCone
from sage.structure.element         import Matrix
from sage.matrix.matrix_dense      cimport Matrix_dense
from sage.misc.misc                 import is_iterator
from .conversions import (incidence_matrix_to_bit_rep_of_facets,
                          incidence_matrix_to_bit_rep_of_Vrep,
                          facets_tuple_to_bit_rep_of_facets,
                          facets_tuple_to_bit_rep_of_Vrep)
from sage.geometry.polyhedron.combinatorial_polyhedron.conversions cimport Vrep_list_to_bit_rep
from sage.misc.cachefunc            import cached_method

from sage.rings.integer                cimport smallInteger
from cysignals.signals                 cimport sig_check

from sage.geometry.polyhedron.combinatorial_polyhedron.face_data_structure cimport face_len_atoms, face_init, face_free
from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator cimport iter_t, parallel_f_vector


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyhedron, a Polytope.

    INPUT:

    - ``data`` -- an instance of
       * :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
       * or a :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`
       * or a :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone`
       * or an ``incidence_matrix`` as in
         :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
         In this case you should also specify the ``Vrep`` and ``facets`` arguments
       * or list of facets, each facet given as
         a list of ``[vertices, rays, lines]`` if the polyhedron is unbounded,
         then rays and lines and the extra argument ``nr_lines`` are required
         if the polyhedron contains no lines, the rays can be thought of
         as the vertices of the facets deleted from a bounded polyhedron see
         :class:`~sage.geometry.polyhedron.parent.Polyhedron_base` on how to use
         rays and lines
       * or an integer, representing the dimension of a polyhedron equal to its
         affine hull
       * or a tuple consisting of facets and vertices as two
         :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces`.
    - ``Vrep`` -- (optional) when ``data`` is an incidence matrix, it should
      be the list of ``[vertices, rays, lines]``, if the rows in the incidence_matrix
      should correspond to names
    - ``facets`` -- (optional) when ``data`` is an incidence matrix or a list of facets,
      it should be a list of facets that would be used instead of indices (of the columns
      of the incidence matrix).
    - ``unbounded`` -- value will be overwritten if ``data`` is a polyhedron;
      if ``unbounded`` and ``data`` is incidence matrix or a list of facets,
      need to specify ``far_face``
    - ``far_face`` -- (semi-optional); if the polyhedron is unbounded this
      needs to be set to the list of indices of the rays and line unless ``data`` is
      an instance of :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

    EXAMPLES:

    We illustrate all possible input: a polyhedron:

        sage: P = polytopes.cube()
        sage: CombinatorialPolyhedron(P)
        A 3-dimensional combinatorial polyhedron with 6 facets

    a lattice polytope::

        sage: points = [(1,0,0), (0,1,0), (0,0,1),
        ....: (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: L = LatticePolytope(points)
        sage: CombinatorialPolyhedron(L)
        A 3-dimensional combinatorial polyhedron with 8 facets

    a cone::

        sage: M = Cone([(1,0), (0,1)])
        sage: CombinatorialPolyhedron(M)
        A 2-dimensional combinatorial polyhedron with 2 facets

    an incidence matrix::

        sage: P = Polyhedron(rays=[[0,1]])
        sage: data = P.incidence_matrix()
        sage: far_face = [i for i in range(2) if not P.Vrepresentation()[i].is_vertex()]
        sage: CombinatorialPolyhedron(data, unbounded=True, far_face=far_face)
        A 1-dimensional combinatorial polyhedron with 1 facet
        sage: C = CombinatorialPolyhedron(data, Vrep=['myvertex'],
        ....: facets=['myfacet'], unbounded=True, far_face=far_face)
        sage: C.Vrepresentation()
        ('myvertex',)
        sage: C.Hrepresentation()
        ('myfacet',)

    a list of facets::

        sage: CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
        A 3-dimensional combinatorial polyhedron with 4 facets
        sage: facetnames = ['facet0', 'facet1', 'facet2', 'myfacet3']
        sage: facetinc = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        sage: C = CombinatorialPolyhedron(facetinc, facets=facetnames)
        sage: C.Vrepresentation()
        (1, 2, 3, 4)
        sage: C.Hrepresentation()
        ('facet0', 'facet1', 'facet2', 'myfacet3')

    an integer::

        sage: CombinatorialPolyhedron(-1).f_vector()
        (1)
        sage: CombinatorialPolyhedron(0).f_vector()
        (1, 1)
        sage: CombinatorialPolyhedron(5).f_vector()
        (1, 0, 0, 0, 0, 0, 1)

    tuple of ``ListOfFaces``::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
        ....:     import facets_tuple_to_bit_rep_of_facets, \
        ....:            facets_tuple_to_bit_rep_of_Vrep
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = facets_tuple_to_bit_rep_of_facets(bi_pyr, 6)
        sage: Vrep = facets_tuple_to_bit_rep_of_Vrep(bi_pyr, 6)
        sage: C = CombinatorialPolyhedron((facets, Vrep)); C
        A 3-dimensional combinatorial polyhedron with 8 facets
        sage: C.f_vector()
        (1, 6, 12, 8, 1)

    Specifying that a polyhedron is unbounded is important. The following with a
    polyhedron works fine::

        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P)  # this works fine
        sage: C
        A 2-dimensional combinatorial polyhedron with 2 facets

    The following is incorrect, as ``unbounded`` is implicitly set to ``False``::

        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()
        sage: C = CombinatorialPolyhedron(data, Vrep=vert)
        sage: C
        A 2-dimensional combinatorial polyhedron with 2 facets
        sage: C.f_vector()
        Traceback (most recent call last):
        ...
        ValueError: not all vertices are intersections of facets
        sage: C.vertices()
        (A line in the direction (0, 1), A vertex at (1, 0), A vertex at (-1, 0))

    The correct usage is::

        sage: far_face = [i for i in range(3) if not P.Vrepresentation()[i].is_vertex()]
        sage: C = CombinatorialPolyhedron(data, Vrep=vert, unbounded=True, far_face=far_face)
        sage: C
        A 2-dimensional combinatorial polyhedron with 2 facets
        sage: C.f_vector()
        (1, 0, 2, 1)
        sage: C.vertices()
        ()

    TESTS:

    Checking that :issue:`27987` is fixed::

        sage: P1 = Polyhedron(vertices=[[0,1],[1,0]], rays=[[1,1]])
        sage: P2 = Polyhedron(vertices=[[0,1],[1,0],[1,1]])
        sage: P1.incidence_matrix() == P2.incidence_matrix()
        True
        sage: CombinatorialPolyhedron(P1).f_vector()
        (1, 2, 3, 1)
        sage: CombinatorialPolyhedron(P2).f_vector()
        (1, 3, 3, 1)
        sage: P1 = Polyhedron(vertices=[[0,1],[1,0]], rays=[[1,1]])
        sage: P2 = Polyhedron(vertices=[[0,1],[1,0],[1,1]])
        sage: CombinatorialPolyhedron(P1).f_vector()
        (1, 2, 3, 1)
        sage: CombinatorialPolyhedron(P2).f_vector()
        (1, 3, 3, 1)

    Some other tests regarding small polyhedra::

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: C = CombinatorialPolyhedron(P)
        sage: C
        A 2-dimensional combinatorial polyhedron with 2 facets
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0),)
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()
        sage: far_face = [i for i in range(3) if not P.Vrepresentation()[i].is_vertex()]
        sage: C = CombinatorialPolyhedron(data, Vrep=vert, unbounded=True, far_face=far_face)
        sage: C
        A 2-dimensional combinatorial polyhedron with 2 facets
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0),)
        sage: CombinatorialPolyhedron(3r)
        A 3-dimensional combinatorial polyhedron with 0 facets

    Check that on wrong input subsequent calls of ``f_vector`` fail::

        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()
        sage: C = CombinatorialPolyhedron(data, Vrep=vert)
        sage: C.f_vector()
        Traceback (most recent call last):
        ...
        ValueError: not all vertices are intersections of facets
        sage: C.f_vector()
        Traceback (most recent call last):
        ...
        ValueError: not all vertices are intersections of facets

    Check that :issue:`28678` is fixed::

        sage: CombinatorialPolyhedron([])
        A -1-dimensional combinatorial polyhedron with 0 facets
        sage: CombinatorialPolyhedron(LatticePolytope([], lattice=ToricLattice(3)))
        A -1-dimensional combinatorial polyhedron with 0 facets
    """
    def __cinit__(self):
        r"""
        TESTS:

        Not initializing the class, does not give segmentation fault::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base import CombinatorialPolyhedron
            sage: C = CombinatorialPolyhedron.__new__(CombinatorialPolyhedron)
            sage: C.f_vector()
            Traceback (most recent call last):
            ...
            ValueError: the combinatorial polyhedron was not initialized
            sage: C.face_lattice()                                                      # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: the combinatorial polyhedron was not initialized
            sage: C.face_generator()
            Traceback (most recent call last):
            ...
            ValueError: the combinatorial polyhedron was not initialized
        """
        # Note that all values are set to zero at the time ``__cinit__`` is called:
        # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#initialisation-methods
        # In particular, ``__dealloc__`` will not do harm in this case.

        self._dimension = -2  # a "NULL" value
        self._equations = ()
        self._all_faces = None
        self._n_facets = -1

    def __init__(self, data, Vrep=None, facets=None, unbounded=False, far_face=None):
        r"""
        Initialize :class:`CombinatorialPolyhedron`.

        See :class:`CombinatorialPolyhedron`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2], [0,1,3],   # indirect doctest
            ....:                              [0,2,3], [1,2,3]])

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron).run()
        """
        self._equations = ()
        self._far_face_tuple = ()

        if isinstance(data, Polyhedron_base):
            self._init_from_polyhedron(data)
            return
        if isinstance(data, LatticePolytopeClass):
            self._init_from_lattice_polytope(data)
            return
        if isinstance(data, ConvexRationalPolyhedralCone):
            self._init_from_cone(data)
            return

        self._bounded = not unbounded
        if unbounded:
            if not far_face:
                raise ValueError("must specify far face for unbounded polyhedron")
            self._far_face_tuple = tuple(far_face)

        if Vrep:
            self._Vrep = tuple(Vrep)

        self._init_facet_names(facets)

        if data == [] or data == ():
            self._init_as_trivial_polyhedron(-1)
        elif isinstance(data, Matrix):
            self._init_from_incidence_matrix(data)
        elif isinstance(data, numbers.Integral):
            self._init_as_trivial_polyhedron(data)
        elif (isinstance(data, (tuple, list)) and
              len(data) == 2 and
              isinstance(data[0], ListOfFaces) and
              isinstance(data[1], ListOfFaces)):
            self._init_from_ListOfFaces(data[0], data[1])

        else:
            self._init_from_list_of_facets(data)

    cdef _init_from_polyhedron(self, data):
        r'''
        Initialize from :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.
        '''
        self._Vrep = data.Vrepresentation()
        self._facet_names = data.inequalities()
        self._equations = data.equations()
        self._dimension = data.dimension()

        if not data.is_compact():
            self._bounded = False
            self._far_face_tuple = tuple(i for i in range(data.n_Vrepresentation()) if not data.Vrepresentation()[i].is_vertex())
        else:
            self._bounded = True

        return self._init_from_incidence_matrix(data.incidence_matrix())

    cdef _init_from_lattice_polytope(self, data):
        r'''
        Initialize from :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`.
        '''
        self._bounded = True
        self._Vrep = tuple(data.vertices())
        self._facet_names = tuple(data.facet_normals())
        self._dimension = data.dimension()
        return self._init_from_incidence_matrix(data.incidence_matrix())

    cdef _init_from_cone(self, data):
        r'''
        Initialize from :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone`.
        '''
        self._bounded = False
        self._Vrep = tuple(data.rays()) + (data.lattice().zero(),)
        self._facet_names = tuple(data.facet_normals())
        self._far_face_tuple = tuple(i for i in range(len(self._Vrep) - 1))
        self._dimension = data.dim()
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        incidence_matrix = matrix(ZZ, data.incidence_matrix().rows()
                                  + [[ZZ.one() for _ in range(len(data.facet_normals()))]])
        return self._init_from_incidence_matrix(incidence_matrix)

    cdef _init_facet_names(self, facets):
        '''
        Store facet names and compute equations.
        '''
        if facets is not None:
            facets = tuple(facets)

            test = [1] * len(facets)  # 0 if that facet is an equation
            for i in range(len(facets)):
                if hasattr(facets[i], "is_inequality"):
                    # We remove equations.
                    # At the moment only equations with this attribute ``True``
                    # will be detected.
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._facet_names = tuple(facets[i] for i in range(len(facets)) if test[i])

            self._equations = tuple(facets[i] for i in range(len(facets)) if not test[i])
        else:
            self._facet_names = None

    cdef _init_from_incidence_matrix(self, data):
        """
        Initialize from an incidence matrix.
        """
        # Input is incidence-matrix or was converted to it.
        self._n_Hrepresentation = data.ncols()
        self._n_Vrepresentation = data.nrows()

        if not isinstance(data, Matrix_dense):
            from sage.rings.integer_ring import ZZ
            from sage.matrix.constructor import matrix
            data = matrix(ZZ, data, sparse=False)
            assert isinstance(data, Matrix_dense), "conversion to ``Matrix_dense`` didn't work"

        # Store the incidence matrix.
        if not data.is_immutable():
            data = data.__copy__()
            data.set_immutable()
        self.incidence_matrix.set_cache(data)

        # Delete equations.
        data = data.delete_columns(
                [i for i in range(data.ncols())
                 if all(data[j, i] for j in range(data.nrows()))],
                check=False)

        # Initializing the facets in their Bit-representation.
        self._bitrep_facets = incidence_matrix_to_bit_rep_of_facets(data)

        # Initializing the Vrep as their Bit-representation.
        self._bitrep_Vrep = incidence_matrix_to_bit_rep_of_Vrep(data)

        self._n_facets = self.bitrep_facets().n_faces()

        self._initialize_far_face()

    cdef _init_from_list_of_facets(self, data):
        """
        Initialize from a list of facets.

        Tuple and iterator work as well.

        The facets are given by its ``[vertices, rays, lines]``.
        """
        if is_iterator(data):
            data = tuple(data)

        if self._Vrep is None:
            # Get the names of the Vrep.
            Vrep = sorted(set.union(*map(set, data)))
            n_Vrepresentation = len(Vrep)
            if Vrep != range(len(Vrep)):
                self._Vrep = tuple(Vrep)
                Vinv = {v: i for i, v in enumerate(self._Vrep)}
        else:
            # Assuming the user gave as correct names for the vertices
            # and labeled them instead by `0,...,n`.
            n_Vrepresentation = len(self._Vrep)

        self._n_Vrepresentation = n_Vrepresentation

        # Relabel the Vrep to be `0,...,n`.
        if self._Vrep is not None:
            def f(v):
                return Vinv[v]
        else:
            def f(v):
                return int(v)
        facets = tuple(tuple(f(i) for i in j) for j in data)

        self._n_facets = len(facets)
        self._n_Hrepresentation = len(facets)

        # Initializing the facets in their Bit-representation.
        self._bitrep_facets = facets_tuple_to_bit_rep_of_facets(facets, n_Vrepresentation)

        # Initializing the Vrep as their Bit-representation.
        self._bitrep_Vrep = facets_tuple_to_bit_rep_of_Vrep(facets, n_Vrepresentation)

        self._initialize_far_face()

    cdef _init_from_ListOfFaces(self, ListOfFaces facets, ListOfFaces Vrep):
        """
        Initialize ``self`` from two ``ListOfFaces``.
        """
        self._bitrep_facets = facets
        self._bitrep_Vrep = Vrep

        self._n_Hrepresentation = self._bitrep_facets.n_faces()
        self._n_Vrepresentation = self._bitrep_Vrep.n_faces()
        self._n_facets = self._n_Hrepresentation

        self._initialize_far_face()

    cdef _initialize_far_face(self):
        """
        Initialize far_face if unbounded.
        """
        if not self._bounded:
            face_init(self._far_face, self.bitrep_facets().n_atoms(), self._n_facets)
            Vrep_list_to_bit_rep(tuple(self._far_face_tuple), self._far_face)

    cdef _init_as_trivial_polyhedron(self, int dimension):
        """
        Initialize polyhedron equal to its affine hull.
        """
        if dimension < -1:
            raise ValueError("any polyhedron must have dimension at least -1")
        self._dimension = dimension

        if self._dimension == 0:
            self._n_facets = 1
            self._n_Vrepresentation = 1
        else:
            self._n_facets = 0
            self._n_Vrepresentation = 0

        # Initializing the facets in their Bit-representation.
        self._bitrep_facets = facets_tuple_to_bit_rep_of_facets((), 0)

        # Initializing the Vrep as their Bit-representation.
        self._bitrep_Vrep = facets_tuple_to_bit_rep_of_Vrep((), 0)

    def __dealloc__(self):
        """
        TESTS::

            sage: CombinatorialPolyhedron(-2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: any polyhedron must have dimension at least -1
        """
        if not self._bounded:
            face_free(self._far_face)

    def _repr_(self):
        r"""
        Return a description of the combinatorial polyhedron.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'A 3-dimensional combinatorial polyhedron with 4 facets'

            sage: P = Polyhedron(vertices=[])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'A -1-dimensional combinatorial polyhedron with 0 facets'

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'A 0-dimensional combinatorial polyhedron with 0 facets'

            sage: P = Polyhedron(lines=[[0,0,1],[0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'A 2-dimensional combinatorial polyhedron with 0 facets'

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[-1,0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'A 2-dimensional combinatorial polyhedron with 1 facet'
        """
        desc = "A {}-dimensional combinatorial polyhedron with {} facet"\
            .format(self.dimension(), self.n_facets())
        if self.n_facets() != 1:
            desc += "s"
        return desc

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_generator()
            sage: it1 = C1.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(),
            ....:              face.ambient_Hrepresentation()) for face in it)
            sage: tup1 = tuple((face.ambient_Vrepresentation(),
            ....:               face.ambient_Hrepresentation()) for face in it1)
            sage: tup == tup1
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_generator()
            sage: it1 = C1.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it)
            sage: tup1 = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_generator()
            sage: it1 = C1.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it)
            sage: tup1 = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_generator()
            sage: it1 = C1.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it)
            sage: tup1 = tuple((face.ambient_Vrepresentation(), face.ambient_Hrepresentation()) for face in it1)
            sage: tup == tup1
            True
        """
        # Give a constructor by list of facets.
        if not self.is_bounded():
            return (CombinatorialPolyhedron, (self.incidence_matrix(),
                    self.Vrepresentation(), self.Hrepresentation(),
                    True, self.far_face_tuple()))
        else:
            return (CombinatorialPolyhedron, (self.incidence_matrix(),
                    self.Vrepresentation(), self.Hrepresentation()))

    def _test_bitsets(self, tester=None, **options):
        """
        Test if the bitsets are consistent.

        TESTS::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._test_bitsets()
        """
        if tester is None:
            tester = self._tester(**options)

        cdef ListOfFaces facets = self.bitrep_facets()
        cdef ListOfFaces Vrep = self.bitrep_Vrep()

        tester.assertEqual(facets.matrix(), Vrep.matrix().transpose())

    def Vrepresentation(self):
        r"""
        Return a list of names of ``[vertices, rays, lines]``.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0], \
            ....:                      [0,0,1],[0,0,-1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Vrepresentation()
            (A line in the direction (0, 0, 1),
             A ray in the direction (1, 0, 0),
             A vertex at (0, 0, 0),
             A ray in the direction (0, 1, 0))

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....: (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.Vrepresentation()
            (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))

            sage: M = Cone([(1,0), (0,1)])
            sage: CombinatorialPolyhedron(M).Vrepresentation()
            (N(1, 0), N(0, 1), N(0, 0))
        """
        if self.Vrep() is not None:
            return self.Vrep()
        else:
            return tuple(smallInteger(i) for i in range(self.n_Vrepresentation()))

    def Hrepresentation(self):
        r"""
        Return a list of names of facets and possibly some equations.

        EXAMPLES::

            sage: P = polytopes.permutahedron(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Hrepresentation()
            (An inequality (1, 1, 0) x - 3 >= 0,
             An inequality (-1, -1, 0) x + 5 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (-1, 0, 0) x + 3 >= 0,
             An inequality (1, 0, 0) x - 1 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0,
             An equation (1, 1, 1) x - 6 == 0)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....: (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.Hrepresentation()
            (N(1, -1, -1),
             N(1, 1, -1),
             N(1, 1, 1),
             N(1, -1, 1),
             N(-1, -1, 1),
             N(-1, -1, -1),
             N(-1, 1, -1),
             N(-1, 1, 1))

            sage: M = Cone([(1,0), (0,1)])
            sage: CombinatorialPolyhedron(M).Hrepresentation()
            (M(0, 1), M(1, 0))
        """
        if self.facet_names() is not None:
            return self.facet_names() + self.equations()
        else:
            return tuple(smallInteger(i) for i in range(self.n_Hrepresentation()))

    def dimension(self):
        r"""
        Return the dimension of the polyhedron.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron([(1,2,3), (1,2,4),
            ....:                              (1,3,4), (2,3,4)])
            sage: C.dimension()
            3

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: CombinatorialPolyhedron(P).dimension()
            3

        ``dim`` is an alias::

            sage: CombinatorialPolyhedron(P).dim()
            3
        """
        if self._dimension == -2:
            # Dimension not computed yet.
            if self.n_facets() == -1:
                raise ValueError("the combinatorial polyhedron was not initialized")
            elif self.n_facets() == 0:
                # The dimension of a trivial polyhedron is assumed to contain
                # exactly one "vertex" and for each dimension one "line" as in
                # :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
                self._dimension = self.n_Vrepresentation() - 1
            elif not self.is_bounded() or self.n_facets() <= self.n_Vrepresentation():
                self._dimension = self.bitrep_facets().compute_dimension()
            else:
                # If the polyhedron has many facets,
                # calculating the dimension of the dual will be faster.
                # The dual exists, if the polyhedron is bounded.
                self._dimension = self.bitrep_facets().compute_dimension()
        return smallInteger(self._dimension)

    dim = dimension

    @cached_method
    def n_vertices(self):
        r"""
        Return the number of vertices.

        Is equivalent to ``len(self.vertices())``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_vertices()
            8

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_vertices()
            20

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_vertices()
            0

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0]], lines=[[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_vertices()
            0

            sage: C = CombinatorialPolyhedron(4)
            sage: C.f_vector()
            (1, 0, 0, 0, 0, 1)
            sage: C.n_vertices()
            0

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.n_vertices()
            1
        """
        if self.dimension() == 0:
            # This specific trivial polyhedron needs special attention.
            return smallInteger(1)
        if not self.is_bounded():
            # Some elements in the ``Vrep`` might not correspond to actual combinatorial vertices.
            return len(self.vertices())
        else:
            return smallInteger(self.n_Vrepresentation())

    def vertices(self, names=True):
        r"""
        Return the elements in the Vrepresentation that are vertices.

        In case of an unbounded polyhedron, there might be lines and
        rays in the Vrepresentation.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0, 0),)
            sage: C.Vrepresentation()
            (A vertex at (0, 0, 0),
             A ray in the direction (0, 0, 1),
             A ray in the direction (0, 1, 0),
             A ray in the direction (1, 0, 0))
            sage: P = polytopes.cross_polytope(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (-1, 0, 0),
             A vertex at (0, -1, 0),
             A vertex at (0, 0, -1),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 0),
             A vertex at (1, 0, 0))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....:           (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0),)
        """
        if unlikely(self.dimension() == 0):
            # Handling the case of a trivial polyhedron of dimension `0`.
            if names and self.Vrep():
                return (self.Vrep()[0],)
            else:
                return (smallInteger(0),)
        if not self.is_bounded():
            it = self.face_iter(0)
            try:
                # The Polyhedron has at least one vertex.
                # In this case every element in the ``Vrep``
                # that is not contained in the far face
                # is a vertex.
                next(it)
            except StopIteration:
                # The Polyhedron has no vertex.
                return ()
        if names and self.Vrep():
            return tuple(self.Vrep()[i] for i in range(self.n_Vrepresentation()) if i not in self.far_face_tuple())
        else:
            return tuple(smallInteger(i) for i in range(self.n_Vrepresentation()) if i not in self.far_face_tuple())

    def n_facets(self):
        r"""
        Return the number of facets.

        Is equivalent to ``len(self.facets())``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_facets()
            6

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_facets()
            170

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_facets()
            2

            sage: P = Polyhedron(rays=[[1,0], [-1,0], [0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.n_facets()
            1

            sage: C = CombinatorialPolyhedron(-1)
            sage: C.f_vector()
            (1)
            sage: C.n_facets()
            0

        Facets are defined to be the maximal nontrivial faces.
        The ``0``-dimensional polyhedron does not have nontrivial faces::

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.n_facets()
            0
        """
        if unlikely(self._dimension == 0):
            # This trivial polyhedron needs special attention.
            return smallInteger(0)
        return smallInteger(self._n_facets)

    def facets(self, names=True):
        r"""
        Return the facets as lists of ``[vertices, rays, lines]``.

        If ``names`` is ``False``, then the Vrepresentatives in the facets
        are given by their indices in the Vrepresentation.

        The facets are the maximal nontrivial faces.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.facets()
            ((A vertex at (1, -1, -1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1),
              A vertex at (1, -1, 1)),
             (A vertex at (1, 1, -1),
              A vertex at (1, 1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, 1),
              A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (1, -1, -1),
              A vertex at (1, 1, -1),
              A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1)),
             (A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, -1, -1)))
            sage: C.facets(names=False)
            ((0, 1, 2, 3),
             (1, 2, 6, 7),
             (2, 3, 4, 7),
             (4, 5, 6, 7),
             (0, 1, 5, 6),
             (0, 3, 4, 5))

        The empty face is trivial and hence the ``0``-dimensional
        polyhedron does not have facets::

            sage: C = CombinatorialPolyhedron(0)
            sage: C.facets()
            ()
        """
        if unlikely(self.dimension() <= 0):
            # Special attention for this trivial case.
            # Facets are defined to be nontrivial faces of codimension 1.
            # The empty face is trivial.
            return ()

        # It is essential to have the facets in the exact same order as
        # on input, so that pickle/unpickle by :meth:`reduce` works.
        # Every facet knows its index by the facet representation.
        face_iter = self.face_iter(self.dimension() - 1, algorithm='primal')
        facets = [None] * self.n_facets()
        for face in face_iter:
            index = face.ambient_H_indices()[0]
            if names:
                verts = face.ambient_Vrepresentation()
            else:
                verts = face.ambient_V_indices()
            facets[index] = verts

        return tuple(facets)

    @cached_method
    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`Vrepresentation`.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: C.incidence_matrix()
            [1 0 0 0 1 1]
            [1 1 0 0 1 0]
            [1 1 1 0 0 0]
            [1 0 1 0 0 1]
            [0 0 1 1 0 1]
            [0 0 0 1 1 1]
            [0 1 0 1 1 0]
            [0 1 1 1 0 0]

        In this case the incidence matrix is only computed once::

            sage: P.incidence_matrix() is C.incidence_matrix()
            True
            sage: C.incidence_matrix.clear_cache()
            sage: C.incidence_matrix() is P.incidence_matrix()
            False
            sage: C.incidence_matrix() == P.incidence_matrix()
            True

        ::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5, backend='field')
            sage: C = P.combinatorial_polyhedron()
            sage: C.incidence_matrix.clear_cache()
            sage: C.incidence_matrix() == P.incidence_matrix()
            True

        The incidence matrix is consistent with
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`::

            sage: P = Polyhedron([[0,0]])
            sage: P.incidence_matrix()
            [1 1]
            sage: C = P.combinatorial_polyhedron()
            sage: C.incidence_matrix.clear_cache()
            sage: P.combinatorial_polyhedron().incidence_matrix()
            [1 1]

        TESTS:

        Check that :issue:`29455` is fixed::

            sage: C = Polyhedron([[0]]).combinatorial_polyhedron()
            sage: C.incidence_matrix.clear_cache()
            sage: C.incidence_matrix()
            [1]
            sage: C = CombinatorialPolyhedron(-1)
            sage: C.incidence_matrix.clear_cache()
            sage: C.incidence_matrix()
            []

        Check that the base ring is ``ZZ``, see :issue:`29840`::

            sage: C = CombinatorialPolyhedron([[0,1,2], [0,1,3], [0,2,3], [1,2,3]])
            sage: C.incidence_matrix().base_ring()
            Integer Ring
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_dense incidence_matrix = matrix(
                ZZ, self.n_Vrepresentation(), self.n_Hrepresentation(), 0)

        if self.dim() < 1:
            # Small cases.
            if self.dim() == 0:
                # To be consistent with ``Polyhedron_base``,
                for i in range(self.n_Hrepresentation()):
                    incidence_matrix.set_unsafe_int(0, i, 1)
            incidence_matrix.set_immutable()
            return incidence_matrix

        # If equations are present, we add them as last columns.
        n_facets = self.n_facets()
        if self.facet_names() is not None:
            n_equations = len(self.equations())
            for Hindex in range(n_facets, n_facets + n_equations):
                for Vindex in range(self.n_Vrepresentation()):
                    incidence_matrix.set_unsafe_int(Vindex, Hindex, 1)

        facet_iter = self.face_iter(self.dimension() - 1, algorithm='primal')
        for facet in facet_iter:
            Hindex = facet.ambient_H_indices()[0]
            for Vindex in facet.ambient_V_indices():
                incidence_matrix.set_unsafe_int(Vindex, Hindex, 1)

        incidence_matrix.set_immutable()

        return incidence_matrix

    cdef int _algorithm_to_dual(self, algorithm) except -2:
        if algorithm == 'primal':
            return 0
        elif algorithm == 'dual':
            if not self.is_bounded():
                raise ValueError("dual algorithm only available for bounded polyhedra")
            return 1
        elif algorithm is None:
            return -1
        else:
            raise ValueError("algorithm must be 'primal', 'dual' or None")

    def edges(self, names=True, algorithm=None):
        r"""
        Return the edges of the polyhedron, i.e. the rank 1 faces.

        INPUT:

        - ``names`` -- boolean (default: ``True``); if ``False``,
          then the Vrepresentatives in the edges are given by
          their indices in the Vrepresentation

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:
          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        .. NOTE::

            To compute edges and f_vector, first compute the edges.
            This might be faster.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (3, 9, 27), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (4, 16, 64)),
             (A vertex at (1, 1, 1), A vertex at (4, 16, 64)),
             (A vertex at (0, 0, 0), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (3, 9, 27)),
             (A vertex at (0, 0, 0), A vertex at (3, 9, 27)),
             (A vertex at (1, 1, 1), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (1, 1, 1)))

            sage: C.edges(names=False)
            ((3, 4), (2, 4), (1, 4), (0, 4), (2, 3), (0, 3), (1, 2), (0, 2), (0, 1))

            sage: P = Polyhedron(rays=[[-1,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A line in the direction (1, 0), A vertex at (0, 0)),)

            sage: P = Polyhedron(vertices=[[0,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0), A vertex at (1, 0)),)

            sage: from itertools import combinations
            sage: N = combinations(['a','b','c','d','e'], 4)
            sage: C = CombinatorialPolyhedron(N)
            sage: C.edges()
            (('d', 'e'),
             ('c', 'e'),
             ('b', 'e'),
             ('a', 'e'),
             ('c', 'd'),
             ('b', 'd'),
             ('a', 'd'),
             ('b', 'c'),
             ('a', 'c'),
             ('a', 'b'))
        """
        self._compute_edges(self._algorithm_to_dual(algorithm))

        # Mapping the indices of the Vrep to the names, if requested.
        if self.Vrep() is not None and names is True:
            def f(size_t i):
                return self.Vrep()[i]
        else:
            def f(size_t i):
                return smallInteger(i)

        cdef size_t j
        return tuple((f(self._edges.get(j).first),
                      f(self._edges.get(j).second))
                     for j in range(self._edges.length))

    def vertex_graph(self, names=True, algorithm=None):
        r"""
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to bounded rank 1 faces.

        INPUT:

        - ``names`` -- boolean (default: ``True``); if ``False``,
          then the nodes of the graph are labeld by the
          indices of the Vrepresentation

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:
          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: G = C.vertex_graph(); G                                               # needs sage.graphs
            Graph on 5 vertices
            sage: sorted(G.degree())                                                    # needs sage.graphs
            [3, 3, 4, 4, 4]

            sage: P = Polyhedron(rays=[[1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.graph()                                                             # needs sage.graphs
            Graph on 1 vertex
        """
        vertices = self.vertices(names=names)

        # Getting the bounded edges.
        edges = tuple(edge for edge in self.edges(names=names, algorithm=algorithm)
                      if edge[0] in vertices and edge[1] in vertices)

        return Graph([vertices, edges], format='vertices_and_edges')

    graph = vertex_graph

    @cached_method
    def vertex_adjacency_matrix(self, algorithm=None):
        """
        Return the binary matrix of vertex adjacencies.

        INPUT:

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:
          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_adjacency_matrix`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: C.vertex_adjacency_matrix()
            [0 1 0 1 0 1 0 0]
            [1 0 1 0 0 0 1 0]
            [0 1 0 1 0 0 0 1]
            [1 0 1 0 1 0 0 0]
            [0 0 0 1 0 1 0 1]
            [1 0 0 0 1 0 1 0]
            [0 1 0 0 0 1 0 1]
            [0 0 1 0 1 0 1 0]

        TESTS::

            sage: CombinatorialPolyhedron(-1).vertex_adjacency_matrix()
            []
            sage: CombinatorialPolyhedron(0).vertex_adjacency_matrix()
            [0]
            sage: polytopes.cube().vertex_adjacency_matrix().is_immutable()
            True
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_dense adjacency_matrix = matrix(
                ZZ, self.n_Vrepresentation(), self.n_Vrepresentation(), 0)
        cdef size_t i, first, second

        self._compute_edges(self._algorithm_to_dual(algorithm))
        for i in range(self._edges.length):
            first = self._edges.get(i).first
            second = self._edges.get(i).second
            adjacency_matrix.set_unsafe_int(first, second, 1)
            adjacency_matrix.set_unsafe_int(second, first, 1)
        adjacency_matrix.set_immutable()
        return adjacency_matrix

    def ridges(self, add_equations=False, names=True, algorithm=None):
        r"""
        Return the ridges.

        The ridges of a polyhedron are the faces
        contained in exactly two facets.

        To obtain all faces of codimension 1 use
        :meth:`CombinatorialPolyhedron.face_generator` instead.

        The ridges will be given by the facets, they are contained in.

        INPUT:

        - ``add_equations`` -- if ``True``, then equations of the polyhedron
          will be added (only applicable when ``names`` is ``True``)

        - ``names`` -- boolean (default: ``True``);
          if ``False``, then the facets are given by their indices

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:
          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        .. NOTE::

            To compute ridges and f_vector, compute the ridges first.
            This might be faster.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(2)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (1, 0) x - 1 >= 0, An inequality (-1, 0) x + 2 >= 0),)
            sage: C.ridges(add_equations=True)
            (((An inequality (1, 0) x - 1 >= 0, An equation (1, 1) x - 3 == 0),
              (An inequality (-1, 0) x + 2 >= 0, An equation (1, 1) x - 3 == 0)),)

            sage: P = polytopes.cyclic_polytope(4,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (24, -26, 9, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (8, -14, 7, -1) x + 0 >= 0))
            sage: C.ridges(names=False)
            ((3, 4),
             (2, 4),
             (1, 4),
             (0, 4),
             (2, 3),
             (1, 3),
             (0, 3),
             (1, 2),
             (0, 2),
             (0, 1))

            sage: P = Polyhedron(rays=[[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C
            A 1-dimensional combinatorial polyhedron with 1 facet
            sage: C.ridges()
            ()
            sage: it = C.face_generator(0)
            sage: for face in it: face.ambient_Hrepresentation()
            (An inequality (1, 0) x + 0 >= 0, An equation (0, 1) x + 0 == 0)

        TESTS:

        Testing that ``add_equations`` is ignored if ``names`` is ``False``::

            sage: C = CombinatorialPolyhedron(polytopes.simplex())
            sage: C.ridges(names=False, add_equations=True)
            ((2, 3), (1, 3), (0, 3), (1, 2), (0, 2), (0, 1))
        """
        self._compute_ridges(self._algorithm_to_dual(algorithm))
        cdef size_t n_ridges = self._ridges.length

        # Mapping the indices of the Vepr to the names, if requested.
        if self.facet_names() is not None and names is True:
            def f(size_t i):
                return self.facet_names()[i]
        else:
            def f(size_t i):
                return smallInteger(i)

        if add_equations and names:
            return tuple(
                    ((f(self._ridges.get(i).first),) + self.equations(),
                     (f(self._ridges.get(i).second),) + self.equations())
                    for i in range (n_ridges))
        else:
            return tuple(
                    (f(self._ridges.get(i).first),
                     f(self._ridges.get(i).second))
                    for i in range (n_ridges))

    @cached_method
    def facet_adjacency_matrix(self, algorithm=None):
        """
        Return the binary matrix of facet adjacencies.

        INPUT:

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:
          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_adjacency_matrix`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: C.facet_adjacency_matrix()
            [0 1 1 0 1 1]
            [1 0 1 1 1 0]
            [1 1 0 1 0 1]
            [0 1 1 0 1 1]
            [1 1 0 1 0 1]
            [1 0 1 1 1 0]

        TESTS::

            sage: CombinatorialPolyhedron(-1).facet_adjacency_matrix()
            []
            sage: CombinatorialPolyhedron(0).facet_adjacency_matrix()
            []
            sage: polytopes.cube().facet_adjacency_matrix().is_immutable()
            True
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        cdef Matrix_dense adjacency_matrix = matrix(
                ZZ, self.n_facets(), self.n_facets(), 0)
        cdef size_t i

        self._compute_ridges(self._algorithm_to_dual(algorithm))
        for i in range(self._ridges.length):
            first = self._ridges.get(i).first
            second = self._ridges.get(i).second
            adjacency_matrix.set_unsafe_int(first, second, 1)
            adjacency_matrix.set_unsafe_int(second, first, 1)
        adjacency_matrix.set_immutable()
        return adjacency_matrix

    def facet_graph(self, names=True, algorithm=None):
        r"""
        Return the facet graph.

        The facet graph of a polyhedron consists of
        ridges as edges and facets as vertices.

        INPUT:

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        If ``names`` is ``False``, the ``vertices`` of the graph will
        be the indices of the facets in the Hrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.facet_graph()                                                       # needs sage.graphs
            Graph on 9 vertices

        TESTS::

            sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
            sage: CombinatorialPolyhedron(P).facet_graph()                              # needs sage.graphs
            Graph on 2 vertices

        Checking that :issue:`28604` is fixed::

            sage: C = CombinatorialPolyhedron(polytopes.cube()); C
            A 3-dimensional combinatorial polyhedron with 6 facets
            sage: C.facet_graph(names=False)                                            # needs sage.graphs
            Graph on 6 vertices

            sage: C = CombinatorialPolyhedron(polytopes.hypersimplex(5,2)); C
            A 4-dimensional combinatorial polyhedron with 10 facets
            sage: C.facet_graph()                                                       # needs sage.combinat sage.graphs
            Graph on 10 vertices
        """
        face_iter = self.face_iter(self.dimension() - 1, algorithm='primal')
        if names:
            V = list(facet.ambient_Hrepresentation() for facet in face_iter)
        else:
            V = list(facet.ambient_V_indices() for facet in face_iter)
        E = self.ridges(names=names, add_equations=True, algorithm=algorithm)
        if not names:
            # If names is false, the ridges are given as tuple of indices,
            # i.e. (1,2) instead of (('f1',), ('f2',)).
            V = list(v[0] for v in V)
        return Graph([V, E], format='vertices_and_edges')

    @cached_method
    def vertex_facet_graph(self, names=True):
        r"""
        Return the vertex-facet graph.

        This method constructs a directed bipartite graph.
        The nodes of the graph correspond to elements of the Vrepresentation
        and facets. There is a directed edge from Vrepresentation to facets
        for each incidence.

        If ``names`` is set to ``False``, then the vertices (of the graph) are given by
        integers.

        INPUT:

        - ``names`` -- boolean (default: ``True``); if ``True`` label the vertices of the
          graph by the corresponding names of the Vrepresentation resp. Hrepresentation;
          if ``False`` label the vertices of the graph by integers

        EXAMPLES::

            sage: P = polytopes.hypercube(2).pyramid()
            sage: C = CombinatorialPolyhedron(P)
            sage: G = C.vertex_facet_graph(); G                                         # needs sage.graphs
            Digraph on 10 vertices
            sage: C.Vrepresentation()
            (A vertex at (0, -1, -1),
             A vertex at (0, -1, 1),
             A vertex at (0, 1, -1),
             A vertex at (0, 1, 1),
             A vertex at (1, 0, 0))
            sage: sorted(G.neighbors_out(C.Vrepresentation()[4]))                       # needs sage.graphs
            [An inequality (-1, -1, 0) x + 1 >= 0,
             An inequality (-1, 0, -1) x + 1 >= 0,
             An inequality (-1, 0, 1) x + 1 >= 0,
             An inequality (-1, 1, 0) x + 1 >= 0]

        If ``names`` is ``True`` (the default) but the combinatorial polyhedron
        has been initialized without specifying names to
        ``Vrepresentation`` and ``Hrepresentation``,
        then indices of the Vrepresentation and the facets will be used along
        with a string 'H' or 'V'::

            sage: C = CombinatorialPolyhedron(P.incidence_matrix())
            sage: C.vertex_facet_graph().vertices(sort=True)                            # needs sage.graphs
            [('H', 0),
             ('H', 1),
             ('H', 2),
             ('H', 3),
             ('H', 4),
             ('V', 0),
             ('V', 1),
             ('V', 2),
             ('V', 3),
             ('V', 4)]

        If ``names`` is ``False`` then the vertices of the graph are given by integers::

            sage: C.vertex_facet_graph(names=False).vertices(sort=True)                 # needs sage.graphs
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        TESTS:

        Test that :issue:`29898` is fixed::

            sage: Polyhedron().vertex_facet_graph()                                     # needs sage.graphs
            Digraph on 0 vertices
            sage: Polyhedron([[0]]).vertex_facet_graph()                                # needs sage.graphs
            Digraph on 1 vertex
            sage: Polyhedron([[0]]).vertex_facet_graph(False)                           # needs sage.graphs
            Digraph on 1 vertex
        """
        from sage.graphs.digraph import DiGraph
        if self.dimension() == -1:
            return DiGraph()
        if self.dimension() == 0:
            if not names:
                return DiGraph(1)
            else:
                Vrep = self.Vrep()
                if Vrep:
                    v = Vrep[0]
                else:
                    v = ("V", 0)
                return DiGraph([[v], []])

        # The face iterator will iterate through the facets in opposite order.
        facet_iter = self.face_iter(self.dimension() - 1, algorithm='primal')
        n_facets = self.n_facets()
        n_Vrep = self.n_Vrepresentation()

        if not names:
            vertices = [i for i in range(n_facets + n_Vrep)]
            edges = tuple((j, n_Vrep + n_facets - 1 - i) for i, facet in enumerate(facet_iter) for j in facet.ambient_V_indices())
        else:
            facet_names = self.facet_names()
            if facet_names is None:
                # No names were provided at initialisation.
                facet_names = [("H", i) for i in range(n_facets)]

            Vrep = self.Vrep()
            if Vrep is None:
                # No names were provided at initialisation.
                Vrep = [("V", i) for i in range(n_Vrep)]

            vertices = Vrep + facet_names
            edges = tuple((Vrep[j], facet_names[n_facets - 1 - i]) for i, facet in enumerate(facet_iter) for j in facet.ambient_V_indices())
        return DiGraph([vertices, edges], format='vertices_and_edges', immutable=True)

    @cached_method
    def f_vector(self, num_threads=None, parallelization_depth=None, algorithm=None):
        r"""
        Compute the ``f_vector`` of the polyhedron.

        The ``f_vector`` contains the number of faces of dimension `k`
        for each `k` in ``range(-1, self.dimension() + 1)``.

        INPUT:

        - ``num_threads`` -- integer (optional); specify the number of threads

        - ``parallelization_depth`` -- integer (optional); specify
          how deep in the lattice the parallelization is done

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        .. NOTE::

            To obtain edges and/or ridges as well, first do so. This might
            already compute the ``f_vector``.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)

            sage: P = polytopes.cyclic_polytope(6,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 10, 45, 120, 185, 150, 50, 1)

        Using two threads::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector(num_threads=2)
            (1, 120, 240, 150, 30, 1)

        TESTS::

            sage: type(C.f_vector())
            <class 'sage.modules.vector_integer_dense.Vector_integer_dense'>
        """
        if num_threads is None:
            from sage.parallel.ncpus import ncpus
            num_threads = ncpus()

        if parallelization_depth is None:
            # Setting some reasonable defaults.
            if num_threads == 0:
                parallelization_depth = 0
            elif num_threads <= 3:
                parallelization_depth = 1
            elif num_threads <= 8:
                parallelization_depth = 2
            else:
                parallelization_depth = 3

        if not self._f_vector:
            self._compute_f_vector(num_threads, parallelization_depth, self._algorithm_to_dual(algorithm))
        if not self._f_vector:
            raise ValueError("could not determine f_vector")
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ
        f_vector = vector(ZZ, self._f_vector)
        f_vector.set_immutable()
        return f_vector

    def flag_f_vector(self, *args):
        r"""
        Return the flag f-vector.

        For each `-1 < i_0 < \dots < i_n < d` the flag f-vector
        counts the number of flags `F_0 \subset \dots \subset F_n`
        with `F_j` of dimension `i_j` for each `0 \leq j \leq n`,
        where `d` is the dimension of the polyhedron.

        INPUT:

        - ``args`` -- integer (optional); specify an entry of the
          flag-f-vector (must be an increasing sequence of integers)

        OUTPUT:

        - a dictionary, if no arguments were given

        - an integer, if arguments were given

        EXAMPLES:

        Obtain the entire flag-f-vector::

            sage: C = polytopes.hypercube(4).combinatorial_polyhedron()
            sage: C.flag_f_vector()                                                     # needs sage.combinat
                {(-1,): 1,
                 (0,): 16,
                 (0, 1): 64,
                 (0, 1, 2): 192,
                 (0, 1, 2, 3): 384,
                 (0, 1, 3): 192,
                 (0, 2): 96,
                 (0, 2, 3): 192,
                 (0, 3): 64,
                 (1,): 32,
                 (1, 2): 96,
                 (1, 2, 3): 192,
                 (1, 3): 96,
                 (2,): 24,
                 (2, 3): 48,
                 (3,): 8,
                 (4,): 1}

        Specify an entry::

            sage: C.flag_f_vector(0,3)                                                  # needs sage.combinat
            64
            sage: C.flag_f_vector(2)                                                    # needs sage.combinat
            24

        Leading ``-1`` and trailing entry of dimension are allowed::

            sage: C.flag_f_vector(-1,0,3)                                               # needs sage.combinat
            64
            sage: C.flag_f_vector(-1,0,3,4)                                             # needs sage.combinat
            64

        One can get the number of trivial faces::

            sage: C.flag_f_vector(-1)                                                   # needs sage.combinat
            1
            sage: C.flag_f_vector(4)                                                    # needs sage.combinat
            1

        Polyhedra with lines, have ``0`` entries accordingly::

            sage: C = (Polyhedron(lines=[[1]]) * polytopes.hypercube(2)).combinatorial_polyhedron()
            sage: C.flag_f_vector()                                                     # needs sage.combinat
            {(-1,): 1, (0, 1): 0, (0, 2): 0, (0,): 0, (1, 2): 8, (1,): 4, (2,): 4, 3: 1}

        If the arguments are not stricly increasing or out of range, a key error is raised::

            sage: C.flag_f_vector(-1,0,3,5)                                             # needs sage.combinat
            Traceback (most recent call last):
            ...
            KeyError: (0, 3, 5)
            sage: C.flag_f_vector(-1,3,0)                                               # needs sage.combinat
            Traceback (most recent call last):
            ...
            KeyError: (3, 0)
        """
        flag = self._flag_f_vector()
        if len(args) == 0:
            return flag
        elif len(args) == 1:
            return flag[(args[0],)]
        else:
            dim = self.dimension()
            if args[0] == -1:
                args = args[1:]
            if args[-1] == dim:
                args = args[:-1]
            return flag[tuple(args)]

    @cached_method
    def _flag_f_vector(self):
        r"""
        Obtain the flag-f-vector from the flag-f-polynomial from the face lattice.

        See :meth:`flag_f_vector`.

        TESTS::

            sage: C = CombinatorialPolyhedron(3)
            sage: C._flag_f_vector()                                                    # needs sage.combinat
            {(-1,): 1, (0, 1): 0, (0, 2): 0, (0,): 0, (1, 2): 0, (1,): 0, (2,): 0, 3: 1}
        """
        poly = self.face_lattice().flag_f_polynomial()
        variables = poly.variables()
        dim = self.dimension()
        flag = {(smallInteger(-1),): smallInteger(1)}
        for term in poly.monomials():
            index = tuple([variables.index(var) for var in term.variables()[:-1]])
            if index == ():
                flag[(dim,)] = smallInteger(1)
            else:
                flag[index] = poly.monomial_coefficient(term)

        n_lines = sum([1 for x in self.f_vector() if x == 0])
        if n_lines:
            # The polyhedron has lines and we have to account for that.
            # So we basically shift all entries up by the number of lines
            # and add zero entries for the lines.
            from itertools import combinations
            flag_old = flag
            flag = {(smallInteger(-1),): smallInteger(1)}
            ran = [smallInteger(i) for i in range(self.dim())]
            for k in range(1, self.dim()):
                for comb in combinations(ran, self.dim() - k):
                    if comb[0] < n_lines:
                        # There are no faces of dimension 0,...,n_lines.
                        flag[comb] = smallInteger(0)
                    else:
                        # Shift the old entries up by the number of lines.
                        flag[comb] = flag_old[tuple(i - n_lines for i in comb)]

            flag[self.dimension()] = smallInteger(1)

        return flag

    @cached_method
    def neighborliness(self):
        r"""
        Return the largest ``k``, such that the polyhedron is ``k``-neighborly.

        A polyhedron is `k`-neighborly if every set of `n` vertices forms a face
        for `n` up to `k`.

        In case of the `d`-dimensional simplex, it returns `d + 1`.

        .. SEEALSO::

            :meth:`is_neighborly`

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(8,12)
            sage: C = P.combinatorial_polyhedron()
            sage: C.neighborliness()
            4
            sage: P = polytopes.simplex(6)
            sage: C = P.combinatorial_polyhedron()
            sage: C.neighborliness()
            7
            sage: P = polytopes.cyclic_polytope(4,10)
            sage: P = P.join(P)
            sage: C = P.combinatorial_polyhedron()
            sage: C.neighborliness()
            2
        """
        if self.is_simplex():
            return self.dim() + 1
        else:
            from sage.arith.misc import binomial
            k = 1
            while self.f_vector()[k+1] == binomial(self.n_vertices(), k + 1):
                k += 1
            return k

    @cached_method
    def is_neighborly(self, k=None):
        r"""
        Return whether the polyhedron is neighborly.

        If the input `k` is provided, then return whether the polyhedron is `k`-neighborly.

        A polyhedron is neighborly if every set of `n` vertices forms a face
        for `n` up to floor of half the dimension of the polyhedron.
        It is `k`-neighborly if this is true for `n` up to `k`.

        INPUT:

        - ``k`` -- the dimension up to which to check if every set of ``k``
          vertices forms a face. If no ``k`` is provided, check up to floor
          of half the dimension of the polyhedron.

        OUTPUT:

        - ``True`` if the every set of up to ``k`` vertices forms a face,
        - ``False`` otherwise

        .. SEEALSO::

            :meth:`neighborliness`

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(8,12)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_neighborly()
            True
            sage: P = polytopes.simplex(6)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_neighborly()
            True
            sage: P = polytopes.cyclic_polytope(4,10)
            sage: P = P.join(P)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_neighborly()
            False
            sage: C.is_neighborly(k=2)
            True
        """
        from sage.arith.misc import binomial
        if k is None:
            k = self.dim() // 2
        return all(self.f_vector()[i+1] == binomial(self.n_vertices(), i + 1)
                   for i in range(1, k))

    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        A simplex is a bounded polyhedron with `d+1` vertices, where
        `d` is the dimension.

        EXAMPLES::

            sage: CombinatorialPolyhedron(2).is_simplex()
            False
            sage: CombinatorialPolyhedron([[0,1],[0,2],[1,2]]).is_simplex()
            True
        """
        return self.is_bounded() and (self.dim()+1 == self.n_vertices())

    @cached_method
    def is_simplicial(self):
        r"""
        Test whether the polytope is simplicial.

        This method is not implemented for unbounded polyhedra.

        A polytope is simplicial, if each facet contains exactly `d` vertices,
        where `d` is the dimension of the polytope.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_simplicial()
            True
            sage: P = polytopes.hypercube(4)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_simplicial()
            False

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_simplicial()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not self.is_bounded():
            raise NotImplementedError("this function is implemented for polytopes only")

        cdef ListOfFaces facets = self._bitrep_facets
        cdef size_t n_facets = facets.n_faces()
        cdef size_t i
        cdef int dim = self.dimension()

        for i in range(n_facets):
            if face_len_atoms(facets.data.faces[i]) != dim:
                return False
        return True

    @cached_method
    def simpliciality(self):
        r"""
        Return the largest `k` such that the polytope is `k`-simplicial.

        Return the dimension in case of a simplex.

        A polytope is `k`-simplicial, if every `k`-face is a simplex.

        EXAMPLES::

            sage: cyclic = polytopes.cyclic_polytope(10,4)
            sage: CombinatorialPolyhedron(cyclic).simpliciality()
            3

            sage: hypersimplex = polytopes.hypersimplex(5,2)
            sage: CombinatorialPolyhedron(hypersimplex).simpliciality()
            2

            sage: cross = polytopes.cross_polytope(4)
            sage: P = cross.join(cross)
            sage: CombinatorialPolyhedron(P).simpliciality()
            3

            sage: P = polytopes.simplex(3)
            sage: CombinatorialPolyhedron(P).simpliciality()
            3

            sage: P = polytopes.simplex(1)
            sage: CombinatorialPolyhedron(P).simpliciality()
            1

        TESTS::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.simpliciality is C.simpliciality
            True
        """
        if not self.is_bounded():
            raise NotImplementedError("must be bounded")
        cdef FaceIterator face_iter = self._face_iter(False, -2)
        cdef int d
        cdef int dim = self.dimension()

        if self.n_facets() == self.dimension() + 1:
            # A simplex.
            return self.dimension()

        cdef simpliciality = dim - 1

        # For each face in the iterator, check if its a simplex.
        face_iter.structure.lowest_dimension = 2  # every 1-face is a simplex
        d = face_iter.next_dimension()
        while d < dim:
            sig_check()
            if face_iter.n_atom_rep() == d + 1:
                # The current face is a simplex.
                face_iter.ignore_subfaces()
            else:
                # Current face is not a simplex.
                if simpliciality > d - 1:
                    simpliciality = d - 1
            d = face_iter.next_dimension()
            if simpliciality == 1:
                # Every polytope is 1-simplicial.
                d = dim
        return smallInteger(simpliciality)

    @cached_method
    def is_simple(self):
        r"""
        Test whether the polytope is simple.

        If the polyhedron is unbounded, return ``False``.

        A polytope is simple, if each vertex is contained in exactly `d` facets,
        where `d` is the dimension of the polytope.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_simple()
            False
            sage: P = polytopes.hypercube(4)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_simple()
            True

        Return ``False`` for unbounded polyhedra::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_simple()
            False
        """
        if not self.is_bounded():
            return False

        cdef ListOfFaces vertices = self._bitrep_Vrep
        cdef size_t n_vertices = vertices.n_faces()
        cdef size_t i
        cdef int dim = self.dimension()

        for i in range(n_vertices):
            if face_len_atoms(vertices.data.faces[i]) != dim:
                return False
        return True

    @cached_method
    def simplicity(self):
        r"""
        Return the largest `k` such that the polytope is `k`-simple.

        Return the dimension in case of a simplex.

        A polytope `P` is `k`-simple, if every `(d-1-k)`-face
        is contained in exactly `k+1` facets of `P` for `1 \leq k \leq d-1`.

        Equivalently it is `k`-simple if the polar/dual polytope is `k`-simplicial.

        EXAMPLES::

            sage: hyper4 = polytopes.hypersimplex(4,2)
            sage: CombinatorialPolyhedron(hyper4).simplicity()
            1

            sage: hyper5 = polytopes.hypersimplex(5,2)
            sage: CombinatorialPolyhedron(hyper5).simplicity()
            2

            sage: hyper6 = polytopes.hypersimplex(6,2)
            sage: CombinatorialPolyhedron(hyper6).simplicity()
            3

            sage: P = polytopes.simplex(3)
            sage: CombinatorialPolyhedron(P).simplicity()
            3

            sage: P = polytopes.simplex(1)
            sage: CombinatorialPolyhedron(P).simplicity()
            1

        TESTS::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.simplicity is C.simplicity
            True
        """
        if not self.is_bounded():
            raise NotImplementedError("must be bounded")
        cdef FaceIterator coface_iter = self._face_iter(True, -2)
        cdef int d
        cdef int dim = self.dimension()

        if self.n_facets() == self.dimension() + 1:
            # A simplex.
            return self.dimension()

        cdef simplicity = dim - 1

        # For each coface in the iterator, check if its a simplex.
        coface_iter.structure.lowest_dimension = 2  # every coface of dimension 1 is a simplex
        d = coface_iter.next_dimension()
        while d < dim:
            sig_check()
            if coface_iter.n_atom_rep() == d + 1:
                # The current coface is a simplex.
                coface_iter.ignore_supfaces()
            else:
                # Current coface is not a simplex.
                if simplicity > d - 1:
                    simplicity = d - 1
            d = coface_iter.next_dimension()
            if simplicity == 1:
                # Every polytope is 1-simple.
                d = dim
        return smallInteger(simplicity)

    @cached_method
    def is_lawrence_polytope(self):
        r"""
        Return ``True`` if ``self`` is a Lawrence polytope.

        A polytope is called a Lawrence polytope if it has a centrally
        symmetric (normalized) Gale diagram.

        Equivalently, there exists a partition `P_1,\dots,P_k`
        of the vertices `V` such that each part
        `P_i` has size `2` or `1` and for each part there exists
        a facet with vertices exactly `V \setminus P_i`.

        EXAMPLES::

            sage: C = polytopes.simplex(5).combinatorial_polyhedron()
            sage: C.is_lawrence_polytope()
            True
            sage: P = polytopes.hypercube(4).lawrence_polytope()
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_lawrence_polytope()
            True
            sage: P = polytopes.hypercube(4)
            sage: C = P.combinatorial_polyhedron()
            sage: C.is_lawrence_polytope()
            False

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_lawrence_polytope()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only

        AUTHORS:

        - Laith Rastanawi
        - Jonathan Kliem

        REFERENCES:

            For more information, see [BaSt1990]_.
        """
        if not self.is_compact():
            raise NotImplementedError("this function is implemented for polytopes only")
        if self.n_Vrepresentation() <= 2:
            return True

        cdef FaceIterator facet_iterator = self._face_iter(False, self.dimension()-1)
        cdef CombinatorialFace facet
        cdef size_t n_vertices = self.n_Vrepresentation()
        cdef size_t one, two, length, counter
        cdef list vertices = [1 for _ in range(n_vertices)]

        for facet in facet_iterator:
            length = facet.n_atom_rep()
            if length >= n_vertices - 2:
                # The facet has at most two non-vertices and corresponds to
                # two symmetric vertices or a vertex at the origin
                # in the Gale transform.
                facet.set_atom_rep()
                counter = 0
                while counter < length:
                    if facet.atom_rep[counter] != counter:
                        # We have found our first non-vertex.
                        one = counter
                        break
                    counter += 1
                else:
                    # The facet contains the first ``length`` vertices.
                    one = length

                if length == n_vertices - 1:
                    # The facet corresponds to a vertex at the origin
                    # of the Gale transform.
                    vertices[one] = 0
                else:
                    # The facet corresponds to two symmetric vertices
                    # of the Gale transform.
                    while counter < length:
                        if facet.atom_rep[counter] != counter + 1:
                            # We have found our second non-vertex.
                            two = counter + 1
                            break
                        counter += 1
                    else:
                        # The second non-vertex is the very last vertex.
                        two = length + 1

                    if vertices[one] == vertices[two]:
                        # Possibly the Gale transform contains duplicates,
                        # we must make sure that the mulitplicites are symmetric as well.
                        # (And not two vertices are symmetric to just one).
                        vertices[one] = 0
                        vertices[two] = 0

        return not any(vertices)

    @cached_method
    def is_pyramid(self, certificate=False):
        r"""
        Test whether the polytope is a pyramid over one of its facets.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return a vertex of the polytope which is the apex of a pyramid,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. The apex of the pyramid or ``None``.

        If ``certificate`` is ``False`` returns a boolean.

        AUTHORS:

        - Laith Rastanawi
        - Jonathan Kliem

        EXAMPLES::

            sage: C = polytopes.cross_polytope(4).combinatorial_polyhedron()
            sage: C.is_pyramid()
            False
            sage: C.is_pyramid(certificate=True)
            (False, None)
            sage: C = polytopes.cross_polytope(4).pyramid().combinatorial_polyhedron()
            sage: C.is_pyramid()
            True
            sage: C.is_pyramid(certificate=True)
            (True, A vertex at (1, 0, 0, 0, 0))
            sage: C = polytopes.simplex(5).combinatorial_polyhedron()
            sage: C.is_pyramid(certificate=True)
            (True, A vertex at (1, 0, 0, 0, 0, 0))

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_pyramid()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact

        TESTS::

            sage: CombinatorialPolyhedron(-1).is_pyramid()
            False
            sage: CombinatorialPolyhedron(-1).is_pyramid(True)
            (False, None)
            sage: CombinatorialPolyhedron(0).is_pyramid()
            True
            sage: CombinatorialPolyhedron(0).is_pyramid(True)
            (True, 0)

        Check that :issue:`30292` is fixed::

            sage: Polyhedron([[0, -1, -1], [0, -1, 1], [0, 1, -1], [0, 1, 1], [1, 0, 0]]).is_pyramid(certificate=True)
            (True, A vertex at (1, 0, 0))
        """
        if not self.is_bounded():
            raise ValueError("polyhedron has to be compact")

        if self.dim() == -1:
            if certificate:
                return (False, None)
            return False

        if self.dim() == 0:
            if certificate:
                return (True, self.Vrepresentation()[0])
            return True

        # Find a vertex that is incident to all elements in Hrepresentation but one.
        vertex_iter = self._face_iter(True, 0)
        n_facets = self.n_facets()
        for vertex in vertex_iter:
            if vertex.n_ambient_Hrepresentation(add_equations=False) == n_facets - 1:
                if certificate:
                    return (True, vertex.ambient_Vrepresentation()[0])
                return True

        if certificate:
            return (False, None)
        return False

    @cached_method
    def is_bipyramid(self, certificate=False):
        r"""
        Test whether the polytope is a bipyramid over some other polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return a vertex of the polytope which is the apex of a pyramid,
          if found

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two vertices of the polytope which are the apices of a
          bipyramid, if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. The first apex.
            b. The second apex.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: C = polytopes.hypercube(4).combinatorial_polyhedron()
            sage: C.is_bipyramid()
            False
            sage: C.is_bipyramid(certificate=True)
            (False, None)
            sage: C = polytopes.cross_polytope(4).combinatorial_polyhedron()
            sage: C.is_bipyramid()
            True
            sage: C.is_bipyramid(certificate=True)
            (True, [A vertex at (1, 0, 0, 0), A vertex at (-1, 0, 0, 0)])

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_pyramid()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact

        TESTS::

            sage: CombinatorialPolyhedron(-1).is_bipyramid()
            False
            sage: CombinatorialPolyhedron(-1).is_bipyramid(True)
            (False, None)
            sage: C = polytopes.cross_polytope(1)
            sage: C.is_bipyramid()
            True
            sage: C.is_bipyramid(True)
            (True, [A vertex at (1), A vertex at (-1)])

        Check that bug analog to :issue:`30292` is avoided::

            sage: Polyhedron([[0, 1, 0], [0, 0, 1], [0, -1, -1], [1, 0, 0], [-1, 0, 0]]).is_bipyramid(certificate=True)
            (True, [A vertex at (1, 0, 0), A vertex at (-1, 0, 0)])

        ALGORITHM:

        Assume all faces of a polyhedron to be given as lists of vertices.

        A polytope is a bipyramid with apexes `v`, `w` if and only if for each
        proper face `v \in F` there exists a face `G` with
        `G \setminus \{w\} = F \setminus \{v\}`
        and vice versa (for each proper face
        `w \in F` there exists ...).

        To check this property it suffices to check for all facets of the polyhedron.
        """
        if not self.is_compact():

            raise ValueError("polyhedron has to be compact")

        n_facets = self.n_facets()
        if n_facets % 2 or self.dim() < 1:
            if certificate:
                return (False, None)
            return False

        facets_incidences = [set(f) for f in self.facets(names=False)]
        verts_incidences = dict()
        for v in self.face_iter(0):
            verts_incidences[v.ambient_V_indices()[0]] = set(v.ambient_H_indices(add_equations=False))

        # Find two vertices ``vert1`` and ``vert2`` such that one of them
        # lies on exactly half of the facets, and the other one lies on
        # exactly the other half.
        from itertools import combinations
        for index1, index2 in combinations(verts_incidences, 2):
            vert1_incidences = verts_incidences[index1]
            vert2_incidences = verts_incidences[index2]
            vert1and2 = vert1_incidences.union(vert2_incidences)
            if len(vert1and2) == n_facets:
                # We have found two candidates for apexes.
                # Remove from each facet ``index1`` resp. ``index2``.
                test_facets = set(frozenset(facet_inc.difference({index1, index2}))
                                  for facet_inc in facets_incidences)
                if len(test_facets) == n_facets/2:
                    # For each `F` containing `index1` there is
                    # `G` containing `index2` such that
                    # `F \setminus \{index1\} =  G \setminus \{index2\}
                    # and vice versa.
                    if certificate:
                        V = self.vertices()
                        return (True, [V[index1], V[index2]])
                    return True

        if certificate:
            return (False, None)
        return False

    @cached_method
    def is_prism(self, certificate=False):
        r"""
        Test whether the polytope is a prism of some polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two facets of the polytope which are the bases of a prism,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. List of the vertices of the first base facet.
            b. List of the vertices of the second base facet.

        If ``certificate`` is ``False`` returns a boolean.

        TESTS::

            sage: CombinatorialPolyhedron(-1).is_prism()
            False
            sage: CombinatorialPolyhedron(1).is_prism()
            False
            sage: C = polytopes.cross_polytope(3).prism().combinatorial_polyhedron()
            sage: C.is_prism(certificate=True)
            (True,
            [(A vertex at (0, 0, 1, 0),
            A vertex at (0, 1, 0, 0),
            A vertex at (0, 0, 0, -1),
            A vertex at (0, 0, -1, 0),
            A vertex at (0, -1, 0, 0),
            A vertex at (0, 0, 0, 1)),
            (A vertex at (1, 1, 0, 0),
            A vertex at (1, 0, 0, -1),
            A vertex at (1, 0, -1, 0),
            A vertex at (1, -1, 0, 0),
            A vertex at (1, 0, 0, 1),
            A vertex at (1, 0, 1, 0))])
            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_prism()
            Traceback (most recent call last):
            ...
            ValueError: self must be bounded
        """
        if not certificate:
            return self.dual().is_bipyramid()

        val, cert = self.dual().is_bipyramid(True)
        if val:
            facets = self.facets()
            return (True, [facets[cert[0]], facets[cert[1]]])

        return (False, None)

    def join_of_Vrep(self, *indices):
        r"""
        Return the smallest face containing all Vrepresentatives indicated by the indices.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator_base.join_of_Vrep`.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.join_of_Vrep(0,1)
            A 1-dimensional face of a 3-dimensional combinatorial polyhedron
            sage: C.join_of_Vrep(0,11).ambient_V_indices()
            (0, 1, 10, 11, 12, 13)
            sage: C.join_of_Vrep(8).ambient_V_indices()
            (8,)
            sage: C.join_of_Vrep().ambient_V_indices()
            ()
        """
        return self.face_generator().join_of_Vrep(*indices)

    def meet_of_Hrep(self, *indices):
        r"""
        Return the largest face contained in all facets indicated by the indices.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator_base.meet_of_Hrep`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.meet_of_Hrep(0)
            A 2-dimensional face of a 3-dimensional combinatorial polyhedron
            sage: C.meet_of_Hrep(0).ambient_H_indices()
            (0,)
            sage: C.meet_of_Hrep(0,1).ambient_H_indices()
            (0, 1)
            sage: C.meet_of_Hrep(0,2).ambient_H_indices()
            (0, 2)
            sage: C.meet_of_Hrep(0,2,3).ambient_H_indices()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
            sage: C.meet_of_Hrep().ambient_H_indices()
            ()
        """
        return self.face_generator().meet_of_Hrep(*indices)

    def face_generator(self, dimension=None, algorithm=None):
        r"""
        Iterator over all proper faces of specified dimension.

        INPUT:

        - ``dimension`` -- if specified, then iterate over only this dimension

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        OUTPUT: :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`

        .. NOTE::

            :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`
            can ignore subfaces or supfaces of the current face.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(dimension=2)
            sage: face = next(it); face
            A 2-dimensional face of a 4-dimensional combinatorial polyhedron
            sage: face.ambient_Vrepresentation()
            (A vertex at (1, 3, 2, 5, 4),
             A vertex at (2, 3, 1, 5, 4),
             A vertex at (3, 1, 2, 5, 4),
             A vertex at (3, 2, 1, 5, 4),
             A vertex at (2, 1, 3, 5, 4),
             A vertex at (1, 2, 3, 5, 4))
            sage: face = next(it); face
            A 2-dimensional face of a 4-dimensional combinatorial polyhedron
            sage: face.ambient_Vrepresentation()
            (A vertex at (2, 1, 4, 5, 3),
             A vertex at (3, 2, 4, 5, 1),
             A vertex at (3, 1, 4, 5, 2),
             A vertex at (1, 3, 4, 5, 2),
             A vertex at (1, 2, 4, 5, 3),
             A vertex at (2, 3, 4, 5, 1))
            sage: face.ambient_Hrepresentation()
            (An inequality (0, 0, -1, -1, 0) x + 9 >= 0,
             An inequality (0, 0, 0, -1, 0) x + 5 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: face.ambient_H_indices()
            (25, 29, 30)
            sage: face = next(it); face
            A 2-dimensional face of a 4-dimensional combinatorial polyhedron
            sage: face.ambient_H_indices()
            (24, 29, 30)
            sage: face.ambient_V_indices()
            (32, 89, 90, 94)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_generator()
            sage: for face in it: face.ambient_Vrepresentation()
            (1, 2, 3)
            (0, 2, 3)
            (0, 1, 3)
            (0, 1, 2)
            (2, 3)
            (1, 3)
            (1, 2)
            (3,)
            (2,)
            (1,)
            (0, 3)
            (0, 2)
            (0,)
            (0, 1)

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator(1)
            sage: for face in it: face.ambient_Vrepresentation()
            (A vertex at (0, 1), A vertex at (1, 0))
            (A ray in the direction (1, 0), A vertex at (1, 0))
            (A ray in the direction (0, 1), A vertex at (0, 1))

        .. SEEALSO::

            :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`,
            :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace`.
        """
        cdef int dual

        dual = self._algorithm_to_dual(algorithm)

        if dual == -1:
            # Determine the faster way, to iterate through all faces.
            if not self.is_bounded() or self.n_facets() <= self.n_Vrepresentation():
                dual = 0
            else:
                dual = 1

        return FaceIterator(self, dual, output_dimension=dimension)

    face_iter = face_generator

    cdef FaceIterator _face_iter(self, bint dual, int dimension):
        r"""
        A method to obtain the FaceIterator as Cython object.

        ``dimension`` is the ``output_dimension`` of
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`.
        If ``dimension == -2`` this will indicate no ``output_dimension``.

        See :meth:`CombinatorialPolyhedron.face_iter`
        """
        if dual and not self.is_bounded():
            raise ValueError("cannot iterate over dual of unbounded polyhedron")
        if dimension == -2:
            return FaceIterator(self, dual)
        else:
            return FaceIterator(self, dual, output_dimension=dimension)

    def face_lattice(self):
        r"""
        Generate the face-lattice.

        OUTPUT: :class:`~sage.combinat.posets.lattices.FiniteLatticePoset`

        .. NOTE::

            Use :meth:`CombinatorialPolyhedron.face_by_face_lattice_index` to get
            the face for each index.

        .. WARNING::

            The labeling of the face lattice might depend on architecture
            and implementation. Relabeling the face lattice with
            :meth:`CombinatorialPolyhedron.face_by_face_lattice_index` or
            the properties obtained from this face will be platform independent.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()                                                      # needs sage.combinat
            Finite lattice containing 5 elements

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: P1 = Polyhedron(rays=[[1,0], [-1,0]])
            sage: C1 = CombinatorialPolyhedron(P1)
            sage: C.face_lattice().is_isomorphic(C1.face_lattice())                     # needs sage.combinat
            True

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()                                                      # needs sage.combinat
            Finite lattice containing 542 elements

        TESTS::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())                      # needs sage.combinat
            True

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())                      # needs sage.combinat
            True
        """
        from sage.combinat.posets.lattices import FiniteLatticePoset
        return FiniteLatticePoset(self.hasse_diagram())

    @cached_method
    def hasse_diagram(self):
        r"""
        Return the Hasse diagram of ``self``.

        This is the Hasse diagram of the poset of the faces of ``self``:
        A directed graph consisting of a vertex for each face
        and an edge for each minimal inclusion of faces.

        .. NOTE::

            The vertices of the Hasse diagram are given by indices.
            Use :meth:`CombinatorialPolyhedron.face_by_face_lattice_index`
            to relabel.

        .. WARNING::

            The indices of the Hasse diagram might depend on architecture
            and implementation. Relabeling the face lattice with
            :meth:`CombinatorialPolyhedron.face_by_face_lattice_index` or
            the properties obtained from this face will be platform independent

        EXAMPLES::

            sage: # needs sage.graphs sage.rings.number_field
            sage: P = polytopes.regular_polygon(4).pyramid()
            sage: C = CombinatorialPolyhedron(P)
            sage: D = C.hasse_diagram(); D
            Digraph on 20 vertices
            sage: D.average_degree()
            21/5
            sage: D.relabel(C.face_by_face_lattice_index)
            sage: dim_0_vert = D.vertices(sort=True)[1:6]; dim_0_vert
            [A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron]
            sage: sorted(D.out_degree(vertices=dim_0_vert))
            [3, 3, 3, 3, 4]
        """
        if not self._face_lattice_incidences:
            # compute all incidences.
            self._compute_face_lattice_incidences()
        if self._face_lattice_incidences is None:
            raise TypeError("could not determine face lattice")

        # Edges of the face-lattice/Hasse diagram.
        cdef size_t j
        cdef size_t n_incidences = self._face_lattice_incidences.length
        edges = tuple(self._face_lattice_incidences[j] for j in range(n_incidences))

        V = tuple(smallInteger(i) for i in range(sum(self._f_vector)))

        from sage.graphs.digraph import DiGraph
        D = DiGraph([V, edges], format='vertices_and_edges', vertex_labels=False)
        return D

    def _face_lattice_dimension(self, index):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its dimension.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()                                                  # needs sage.combinat
            sage: def f(i):
            ....:     return (i, C._face_lattice_dimension(i))
            ....:
            sage: G = F.relabel(f)                                                      # needs sage.combinat
            sage: set(G._elements)                                                      # needs sage.combinat
            {(0, -1),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 1),
             (10, 1),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 1),
             (17, 1),
             (18, 1),
             (19, 1),
             (20, 1),
             (21, 2),
             (22, 2),
             (23, 2),
             (24, 2),
             (25, 2),
             (26, 2),
             (27, 3)}
        """
        f_vector = self.f_vector()
        dim = self.dimension()

        # Getting the dimension, by considering the following:
        # The level-set of dimension `d` will have indices `k, k+1, ..., k+n-1`,
        # where `n` is the number of faces of dimension `d` ( ``n = f_vector[d + 1]``)
        # and `k` is the number of face of dimension up to `d`, i.e.
        # ``k = sum(f_vector[:d])``.
        return max(d for d in range(dim+2) if sum(f_vector[:d]) <= index) - 1

    def face_by_face_lattice_index(self, index):
        r"""
        Return the element of :meth:`CombinatorialPolyhedron.face_lattice` with corresponding index.

        The element will be returned as
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace`.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: F
            Finite lattice containing 28 elements
            sage: G = F.relabel(C.face_by_face_lattice_index)
            sage: G.level_sets()[0]
            [A -1-dimensional face of a 3-dimensional combinatorial polyhedron]
            sage: G.level_sets()[3]
            [A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron]

            sage: P = Polyhedron(rays=[[0,1], [1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()                                                  # needs sage.combinat
            sage: G = F.relabel(C.face_by_face_lattice_index)                           # needs sage.combinat
            sage: G._elements                                                           # needs sage.combinat
            (A -1-dimensional face of a 2-dimensional combinatorial polyhedron,
              A 0-dimensional face of a 2-dimensional combinatorial polyhedron,
              A 1-dimensional face of a 2-dimensional combinatorial polyhedron,
              A 1-dimensional face of a 2-dimensional combinatorial polyhedron,
              A 2-dimensional face of a 2-dimensional combinatorial polyhedron)

            sage: def f(i): return C.face_by_face_lattice_index(i).ambient_V_indices()
            sage: G = F.relabel(f)                                                      # needs sage.combinat
            sage: G._elements                                                           # needs sage.combinat
            ((), (0,), (0, 1), (0, 2), (0, 1, 2))
        """
        self._record_all_faces()                            # Initialize ``_all_faces``, if not done yet.
        dim = self._face_lattice_dimension(index)           # Determine dimension to that index.
        newindex = index - sum(self._f_vector[:dim + 1])    # Index in that level-set.

        # Let ``_all_faces`` determine Vrepresentation.
        return self._all_faces.get_face(dim, newindex)

    def a_maximal_chain(self, Vindex=None, Hindex=None):
        r"""
        Return a maximal chain of the face lattice in increasing order
        without empty face and whole polyhedron/maximal face.

        INPUT:

        - ``Vindex`` -- integer (default: ``None``); prescribe the index of the vertex in the chain
        - ``Hindex`` -- integer (default: ``None``); prescribe the index of the facet in the chain

        Each face is given as
        :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.combinatorial_face.CombinatorialFace`.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(4)
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain(); chain
            [A 0-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 3-dimensional face of a 4-dimensional combinatorial polyhedron]
            sage: [face.ambient_V_indices() for face in chain]
            [(7,), (6, 7), (5, 6, 7), (4, 5, 6, 7)]

            sage: P = polytopes.hypercube(4)
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain(); chain
            [A 0-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 4-dimensional combinatorial polyhedron,
             A 3-dimensional face of a 4-dimensional combinatorial polyhedron]
            sage: [face.ambient_V_indices() for face in chain]
            [(15,), (6, 15), (5, 6, 14, 15), (0, 5, 6, 7, 8, 9, 14, 15)]

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(4)
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain(); chain
            [A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron]
            sage: [face.ambient_V_indices() for face in chain]
            [(16,), (15, 16), (8, 9, 14, 15, 16, 17)]

            sage: P = Polyhedron(rays=[[1,0]], lines=[[0,1]])
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain()
            sage: [face.ambient_V_indices() for face in chain]
            [(0, 1)]

            sage: P = Polyhedron(rays=[[1,0,0],[0,0,1]], lines=[[0,1,0]])
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain()
            sage: [face.ambient_V_indices() for face in chain]
            [(0, 1), (0, 1, 3)]

            sage: P = Polyhedron(rays=[[1,0,0]], lines=[[0,1,0],[0,0,1]])
            sage: C = P.combinatorial_polyhedron()
            sage: chain = C.a_maximal_chain()
            sage: [face.ambient_V_indices() for face in chain]
            [(0, 1, 2)]

        Specify an index for the vertex of the chain::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: [face.ambient_V_indices() for face in C.a_maximal_chain()]
            [(5,), (0, 5), (0, 3, 4, 5)]
            sage: [face.ambient_V_indices() for face in C.a_maximal_chain(Vindex=2)]
            [(2,), (2, 7), (2, 3, 4, 7)]

        Specify an index for the facet of the chain::

            sage: [face.ambient_H_indices() for face in C.a_maximal_chain()]
            [(3, 4, 5), (4, 5), (5,)]
            sage: [face.ambient_H_indices() for face in C.a_maximal_chain(Hindex=3)]
            [(3, 4, 5), (3, 4), (3,)]
            sage: [face.ambient_H_indices() for face in C.a_maximal_chain(Hindex=2)]
            [(2, 3, 5), (2, 3), (2,)]

        If the specified vertex is not contained in the specified facet an error is raised::

            sage: C.a_maximal_chain(Vindex=0, Hindex=3)
            Traceback (most recent call last):
            ...
            ValueError: the given Vindex is not compatible with the given Hindex

        An error is raised, if the specified index does not correspond to a facet::

            sage: C.a_maximal_chain(Hindex=40)
            Traceback (most recent call last):
            ...
            ValueError: the given Hindex does not correspond to a facet

        An error is raised, if the specified index does not correspond to a vertex::

            sage: C.a_maximal_chain(Vindex=40)
            Traceback (most recent call last):
            ...
            ValueError: the given Vindex does not correspond to a vertex

        ::

            sage: P = Polyhedron(rays=[[1,0,0],[0,0,1]], lines=[[0,1,0]])
            sage: C = P.combinatorial_polyhedron()
            sage: C.a_maximal_chain(Vindex=0)
            Traceback (most recent call last):
            ...
            ValueError: the given Vindex does not correspond to a vertex

        ::

            sage: P = Polyhedron(rays=[[1,0,0],[0,0,1]])
            sage: C = P.combinatorial_polyhedron()
            sage: C.a_maximal_chain(Vindex=0)
            [A 0-dimensional face of a 2-dimensional combinatorial polyhedron,
            A 1-dimensional face of a 2-dimensional combinatorial polyhedron]
            sage: C.a_maximal_chain(Vindex=1)
            Traceback (most recent call last):
            ...
            ValueError: the given Vindex does not correspond to a vertex
        """
        if self.n_facets() == 0 or self.dimension() == 0:
            return []

        # We take a face iterator and do one depth-search.
        # Depending on whether it is dual or not,
        # the search will be from the top or bottom.
        cdef FaceIterator it = self.face_generator()
        chain = [None]*(self.dimension())
        dual = it.dual
        final_dim = 0 if not dual else self.dimension()-1

        cdef bint found_Vindex = Vindex is None
        cdef bint found_Hindex = Hindex is None

        # For each dimension we save the first face we see.
        # This is the face whose sub-/supfaces we visit in the next step.
        current_dim = self.dimension()
        for face in it:
            if not found_Hindex:
                if Hindex not in face.ambient_H_indices():
                    continue
                if face.dimension() == self.dimension() - 1:
                    found_Hindex = True
                    if not found_Vindex and Vindex not in face.ambient_V_indices():
                        raise ValueError("the given Vindex is not compatible with the given Hindex")
            if not found_Vindex:
                if Vindex not in face.ambient_V_indices():
                    continue
                if face.dimension() == 0:
                    found_Vindex = True
                    if not found_Hindex and Hindex not in face.ambient_H_indices():
                        raise ValueError("the given Vindex is not compatible with the given Hindex")

            it.only_subsets()
            current_dim = face.dimension()
            chain[current_dim] = face

        if found_Vindex is False:
            raise ValueError("the given Vindex does not correspond to a vertex")
        if found_Hindex is False:
            raise ValueError("the given Hindex does not correspond to a facet")

        if current_dim != final_dim:
            # The polyhedron contains lines.
            # Note that the iterator was always not dual
            # in this case.
            return chain[current_dim:]
        return chain

    def _test_a_maximal_chain(self, tester=None, **options):
        """
        Run tests on the method :meth:`.a_maximal_chain`.

        TESTS::

            sage: polytopes.cross_polytope(3).combinatorial_polyhedron()._test_a_maximal_chain()
        """
        if tester is None:
            tester = self._tester(**options)

        def test_a_chain(b):
            for i in range(len(b) - 1):
                tester.assertTrue(b[i].is_subface(b[i+1]))

        if self.is_bounded():
            b = self.a_maximal_chain()
            test_a_chain(b)
            if not self.n_vertices():
                return

            from sage.misc.prandom import randrange

            if self.n_vertices():
                # We obtain a chain containing a random vertex.
                i = randrange(self.n_vertices())
                b = self.a_maximal_chain(Vindex=i)
                test_a_chain(b)
                tester.assertTrue(all(i in f.ambient_V_indices() for f in b))

            if self.n_facets():
                # We obtain a chain containing a random facet.
                i = randrange(self.n_facets())
                b = self.a_maximal_chain(Hindex=i)
                test_a_chain(b)
                tester.assertTrue(all(i in f.ambient_H_indices() for f in b))

                # We obtain a chain containing that facet
                # and a random vertex contained in it.
                facet = self.facets(names=False)[i]
                j = facet[randrange(len(facet))]
                b = self.a_maximal_chain(Vindex=j, Hindex=i)
                test_a_chain(b)
                tester.assertTrue(all(j in f.ambient_V_indices() for f in b))
                tester.assertTrue(all(i in f.ambient_H_indices() for f in b))

    cdef tuple Vrep(self):
        r"""
        Return the names of the Vrepresentation, if they exist. Else return ``None``.
        """
        return self._Vrep

    cdef tuple facet_names(self):
        r"""
        Return the names Hrepresentatives, which are facets.

        If not given, return ``None``.
        """
        return self._facet_names

    cdef tuple equations(self):
        r"""
        Return the names of the equations.

        If not equations are given, return ``None``.
        """
        return self._equations

    cdef unsigned int n_Vrepresentation(self) noexcept:
        r"""
        Return the number of elements in the Vrepresentation.
        """
        return self._n_Vrepresentation

    cdef unsigned int n_Hrepresentation(self) noexcept:
        r"""
        Return the number of elements in the Hrepresentation.
        """
        return self._n_Hrepresentation

    def is_compact(self):
        r"""
        Return whether the polyhedron is compact.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.is_compact()
            False
            sage: C = CombinatorialPolyhedron([[0,1], [0,2], [1,2]])
            sage: C.is_compact()
            True
            sage: P = polytopes.simplex()
            sage: P.combinatorial_polyhedron().is_compact()
            True
            sage: P = Polyhedron(rays=P.vertices())
            sage: P.combinatorial_polyhedron().is_compact()
            False
        """
        return self.is_bounded()

    cdef bint is_bounded(self) noexcept:
        r"""
        Return whether the polyhedron is bounded.
        """
        return self._bounded

    cdef ListOfFaces bitrep_facets(self):
        r"""
        Return the facets in bit representation.
        """
        return self._bitrep_facets

    cdef ListOfFaces bitrep_Vrep(self):
        r"""
        Return the Vrepresentations in bit representation.
        """
        return self._bitrep_Vrep

    cdef tuple far_face_tuple(self):
        r"""
        Return the far face as it was given on initialization.
        """
        return self._far_face_tuple

    def __eq__(self, other):
        r"""
        Return whether ``self`` and ``other`` are equal.
        """
        if not isinstance(other, CombinatorialPolyhedron):
            return False
        cdef CombinatorialPolyhedron other_C = other
        return (self.n_facets() == other.n_facets()
                and self.Vrepresentation() == other.Vrepresentation()
                and self.facet_names() == other_C.facet_names()
                and self.equations() == other_C.equations()
                and self.dimension() == other.dimension()
                and self.far_face_tuple() == other_C.far_face_tuple()
                and self.incidence_matrix() == other.incidence_matrix())

    # Methods to obtain a different combinatorial polyhedron.

    cpdef CombinatorialPolyhedron dual(self):
        r"""
        Return the dual/polar of ``self``.

        Only defined for bounded polyhedra.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.polar`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = P.combinatorial_polyhedron()
            sage: D = C.dual()
            sage: D.f_vector()
            (1, 6, 12, 8, 1)
            sage: D1 = P.polar().combinatorial_polyhedron()
            sage: D1.face_lattice().is_isomorphic(D.face_lattice())                     # needs sage.combinat
            True

        Polar is an alias to be consistent with :class:`~sage.geometry.polyhedron.base.Polyhedron_base`::

            sage: C.polar().f_vector()
            (1, 6, 12, 8, 1)

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.dual()
            Traceback (most recent call last):
            ...
            ValueError: self must be bounded
        """
        if not self.is_bounded():
            raise ValueError("self must be bounded")
        cdef ListOfFaces new_facets = self.bitrep_Vrep().__copy__()
        cdef ListOfFaces new_Vrep = self.bitrep_facets().__copy__()

        return CombinatorialPolyhedron((new_facets, new_Vrep))

    polar = dual

    cpdef CombinatorialPolyhedron pyramid(self, new_vertex=None, new_facet=None):
        r"""
        Return the pyramid of ``self``.

        INPUT:

        - ``new_vertex`` -- (optional); specify a new vertex name to set up
          the pyramid with vertex names
        - ``new_facet`` -- (optional); specify a new facet name to set up
          the pyramid with facet names

        EXAMPLES::

            sage: C = CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
            sage: C1 = C.pyramid()
            sage: C1.facets()
            ((0, 1, 2, 4), (0, 1, 3, 4), (0, 2, 3, 4), (1, 2, 3, 4), (0, 1, 2, 3))

        ::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = C.pyramid()
            sage: P1 = P.pyramid()
            sage: C2 = P1.combinatorial_polyhedron()
            sage: C2.vertex_facet_graph().is_isomorphic(C1.vertex_facet_graph())        # needs sage.combinat
            True

        One can specify a name for the new vertex::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = P.combinatorial_polyhedron()
            sage: C1 = C.pyramid(new_vertex='apex')
            sage: C1.is_pyramid(certificate=True)
            (True, 'apex')
            sage: C1.facets()[0]
            (A vertex at (0, 0, 0, 0),
             A vertex at (1, 1, 1, 1),
             A vertex at (2, 4, 8, 16),
             A vertex at (3, 9, 27, 81),
             'apex')

        One can specify a name for the new facets::

            sage: # needs sage.rings.number_field
            sage: P = polytopes.regular_polygon(4)
            sage: C = P.combinatorial_polyhedron()
            sage: C1 = C.pyramid(new_facet='base')
            sage: C1.Hrepresentation()
            (An inequality (-1/2, 1/2) x + 1/2 >= 0,
             An inequality (-1/2, -1/2) x + 1/2 >= 0,
             An inequality (1/2, 0.50000000000000000?) x + 1/2 >= 0,
             An inequality (1/2, -1/2) x + 1/2 >= 0,
             'base')

        For unbounded polyhedra, an error is raised::

            sage: C = CombinatorialPolyhedron([[0,1], [0,2]], far_face=[1,2], unbounded=True)
            sage: C.pyramid()
            Traceback (most recent call last):
            ...
            ValueError: self must be bounded
        """
        if not self.is_bounded():
            raise ValueError("self must be bounded")
        cdef ListOfFaces new_facets = self.bitrep_facets().pyramid()
        cdef ListOfFaces new_Vrep = self.bitrep_Vrep().pyramid()

        if new_vertex is not None:
            new_Vrep_names = self.Vrepresentation() + (new_vertex,)
        else:
            new_Vrep_names = None

        if new_facet is not None:
            if self.facet_names() is not None:
                new_facet_names = self.facet_names() + (new_facet,)
            else:
                # Closures inside cpdef functions not yet supported
                new_facet_names = self.Hrepresentation()[:self.n_facets()] + (new_facet,)
        else:
            new_facet_names = None

        return CombinatorialPolyhedron((new_facets, new_Vrep), Vrep=new_Vrep_names, facets=new_facet_names)

    # Internal methods.

    cdef int _compute_f_vector(self, size_t num_threads, size_t parallelization_depth, int dual) except -1:
        r"""
        Compute the ``f_vector`` of the polyhedron.

        See :meth:`f_vector`.
        """
        if self._f_vector:
            return 0  # There is no need to recompute the f_vector.

        cdef int dim = self.dimension()
        cdef MemoryAllocator mem = MemoryAllocator()

        if num_threads == 0:
            # No need to complain.
            num_threads = 1

        if parallelization_depth > dim - 1:
            # Is a very bad choice anyway, but prevent segmenation faults.
            parallelization_depth = dim - 1

        if dual == -1:
            if not self.is_bounded() or self.n_facets() <= self.n_Vrepresentation():
                # In this case the non-dual approach is faster.
                dual = 0
            else:
                # In this case the dual approach is faster.
                dual = 1

        cdef FaceIterator face_iter
        cdef iter_t* structs = <iter_t*> mem.allocarray(num_threads, sizeof(iter_t))
        cdef size_t i

        # For each thread an independent structure.
        face_iters = [self._face_iter(dual, -2) for _ in range(num_threads)]
        for i in range(num_threads):
            face_iter = face_iters[i]
            structs[i][0] = face_iter.structure[0]

        # Initialize ``f_vector``.
        cdef size_t *f_vector = <size_t *> mem.calloc((dim + 2), sizeof(size_t))

        parallel_f_vector(structs, num_threads, parallelization_depth, f_vector)

        self._persist_f_vector(f_vector, dual)

    cdef int _persist_f_vector(self, size_t* input_f_vector, bint input_is_reversed) except -1:
        cdef int dim = self.dimension()

        if input_is_reversed:
            f_vector = \
                tuple(smallInteger(input_f_vector[dim + 1 - i]) for i in range(dim + 2))
        else:
            f_vector = \
                tuple(smallInteger(input_f_vector[i]) for i in range(dim + 2))

        # Sanity checks.
        if dim > 1:
            if f_vector[-2] < self.n_facets():
                raise ValueError("not all facets are joins of vertices")
            if self.is_bounded() and f_vector[1] < self.n_Vrepresentation():
                raise ValueError("not all vertices are intersections of facets")

        self._f_vector = f_vector

    cdef int _compute_edges_or_ridges(self, int dual, bint do_edges) except -1:
        r"""
        Compute the edges of the polyhedron if ``edges`` else the ridges.

        If ``dual``, use the face iterator in dual mode, else in non-dual.
        If ``dual`` is ``-1`` determine this automatically.

        If the ``f_vector`` is unkown computes it as well if computing the edges
        in non-dual mode or the ridges in dual-mode.

        See :meth:`CombinatorialPolyhedron.edges` and :meth:`CombinatorialPolyhedron.ridges`.
        """
        if (self._edges is not None and do_edges) or (self._ridges is not None and not do_edges):
            return 0  # There is no need to recompute.

        if dual == -1:
            # Determine whether to use dual mode or not.
            if not self.is_bounded():
                dual = 0
            else:
                algorithm = self.choose_algorithm_to_compute_edges_or_ridges("edges" if do_edges else "ridges")
                dual = self._algorithm_to_dual(algorithm)

        cdef FaceIterator face_iter
        cdef int dim = self.dimension()

        cdef ListOfPairs edges = ListOfPairs()
        cdef int output_dim_init = 1 if do_edges else dim - 2

        cdef size_t* f_vector = NULL

        try:
            if dim == 1 and (do_edges or self.n_facets() > 1):
                # In this case there is an edge/ridge, but its not a proper face.
                edges.add(0, 1)

            elif dim <= 1 or self.n_facets() == 0:
                # There is no edge/ridge.
                # Prevent an error when calling the face iterator.
                pass

            else:
                if not self._f_vector and ((dual ^ do_edges)):
                    # While doing edges in non-dual mode or ridges in dual-mode
                    # one might as well do the f-vector.
                    f_vector = <size_t *> check_calloc((dim + 2), sizeof(size_t))
                    f_vector[0] = 1
                    f_vector[dim + 1] = 1
                    face_iter = self._face_iter(dual, -2)
                else:
                    face_iter = self._face_iter(dual, output_dim_init)
                self._compute_edges_or_ridges_with_iterator(face_iter, (dual ^ do_edges),
                                                            edges, f_vector)

            # Success, persist the data.
            if f_vector is not NULL:
                self._persist_f_vector(f_vector, dual)

            if do_edges:
                self._edges = edges
            else:
                self._ridges = edges
        finally:
            sig_free(f_vector)

        if do_edges and self._edges is None:
            raise ValueError('could not determine edges')
        elif not do_edges and self._ridges is None:
            raise ValueError('could not determine ridges')

    def choose_algorithm_to_compute_edges_or_ridges(self, edges_or_ridges):
        """
        Use some heuristics to pick primal or dual algorithm for
        computation of edges resp. ridges.

        We estimate how long it takes to compute a face using the primal
        and the dual algorithm. This may differ significantly, so that e.g.
        visiting all faces with the primal algorithm is faster than using
        the dual algorithm to just visit vertices and edges.

        We guess the number of edges and ridges and do a wild estimate on
        the total number of faces.

        INPUT:

        - ``edges_or_ridges`` -- string; one of:
          * ``'edges'``
          * ``'ridges'``

        OUTPUT: either ``'primal'`` or ``'dual'``

        EXAMPLES::

            sage: C = polytopes.permutahedron(5).combinatorial_polyhedron()
            sage: C.choose_algorithm_to_compute_edges_or_ridges("edges")
            'primal'
            sage: C.choose_algorithm_to_compute_edges_or_ridges("ridges")
            'primal'

        ::

            sage: C = polytopes.cross_polytope(5).combinatorial_polyhedron()
            sage: C.choose_algorithm_to_compute_edges_or_ridges("edges")
            'dual'
            sage: C.choose_algorithm_to_compute_edges_or_ridges("ridges")
            'dual'


        ::

            sage: C = polytopes.Birkhoff_polytope(5).combinatorial_polyhedron()
            sage: C.choose_algorithm_to_compute_edges_or_ridges("edges")
            'dual'
            sage: C.choose_algorithm_to_compute_edges_or_ridges("ridges")
            'primal'
            sage: C.choose_algorithm_to_compute_edges_or_ridges("something_else")
            Traceback (most recent call last):
            ...
            ValueError: unknown computation goal something_else
        """
        if self.is_simple():
            per_face_primal = self.n_Vrepresentation() * self.n_facets()
        else:
            per_face_primal = self.n_Vrepresentation() * self.n_facets() ** 2

        if self.is_simplicial():
            per_face_dual = self.n_Vrepresentation() * self.n_facets()
        else:
            per_face_dual = self.n_Vrepresentation() ** 2 * self.n_facets()

        from sage.arith.misc import binomial
        estimate_n_faces = self.dimension() * binomial(min(self.n_facets(), self.n_Vrepresentation()),
                                                       self.dimension() // 2)

        # Note that the runtime per face already computes the coatoms of the next level, i.e.
        # the runtime for each facet suffices to compute all ridges in primal,
        # the runtime for each vertex suffices to compute all edges in dual.
        if edges_or_ridges == "edges":
            estimate_primal = estimate_n_faces * per_face_primal
            estimate_dual = self.n_Vrepresentation() * per_face_dual
        elif edges_or_ridges == "ridges":
            estimate_primal = self.n_facets() * per_face_primal
            estimate_dual = estimate_n_faces * per_face_dual
        else:
            raise ValueError(f"unknown computation goal {edges_or_ridges}")

        return 'dual' if (estimate_dual < estimate_primal) else 'primal'

    cdef size_t _compute_edges_or_ridges_with_iterator(
            self, FaceIterator face_iter, const bint do_atom_rep,
            ListOfPairs edges, size_t* f_vector) except -1:
        r"""
        See :meth:`CombinatorialPolyhedron._compute_edges`.
        """
        cdef size_t a, b                # facets of an edge
        cdef int dim = self.dimension()
        cdef bint do_f_vector = f_vector is not NULL

        # The dimension in which to record the edges or ridges.
        cdef output_dimension = 1 if do_atom_rep else dim - 2

        cdef int d = face_iter.next_dimension()
        while d < dim:
            sig_check()
            if do_f_vector:
                f_vector[d + 1] += 1

            # If ``not do_f_vector`` the iterator is set up
            # for ``output_dimension`` and
            # ``d < dim`` implies
            # ``d == ouput_dimension``.
            if not do_f_vector or d == output_dimension:
                if do_atom_rep:
                    # Set up face_iter.atom_rep
                    face_iter.set_atom_rep()

                    # Copy the information.
                    a = face_iter.structure.atom_rep[0]
                    b = face_iter.structure.atom_rep[1]
                else:
                    # Set up face_iter.coatom_rep
                    face_iter.set_coatom_rep()

                    # Copy the information.
                    a = face_iter.structure.coatom_rep[0]
                    b = face_iter.structure.coatom_rep[1]
                edges.add(a, b)
            d = face_iter.next_dimension()

    cdef int _compute_face_lattice_incidences(self) except -1:
        r"""
        Compute all incidences for the face lattice.

        See :meth:`face_lattice`.
        """
        if self._face_lattice_incidences:
            return 1  # There is no need to recompute the incidences.

        cdef int dim = self.dimension()
        f_vector = self.f_vector()
        self._record_all_faces()  # set up ``self._all_faces``
        cdef PolyhedronFaceLattice all_faces = self._all_faces

        # ``all_faces`` will store its incidences in ``first`` and ``second``.
        cdef size_t first = 0, second = 0

        # ``dimension_one`` and ``dimension_two`` will be the dimensions of the
        # incidences, we currently obtain from ``all_faces``.
        # Almost always ``dimension_two = dimension_one - 1``.
        cdef int dimension_one, dimension_two
        cdef int j  # an index for ``range(dimension_two + 1)``

        # The indices of the incidences in ``all_faces`` are levelwise.
        # Hence, we have to add to each index dependent on dimension:

        # For ``dimension_two`` we add:
        cdef size_t already_seen       # = sum(f_vector[j] for j in range(dimension_two + 1))

        # For ``dimension_one`` we add:
        cdef size_t already_seen_next  # = sum(f_vector[j] for j in range(dimension_two + 2))

        cdef ListOfPairs incidences = ListOfPairs()

        if all_faces is None:
            raise ValueError("could not determine a list of all faces")

        dimension_one = 0
        if dim > -1:
            while f_vector[dimension_one + 1] == 0:
                # Taking care of cases, where there might be no faces
                # of dimension 0, 1, etc (``n_lines > 0``).
                dimension_one += 1
            dimension_two = -1

        while dimension_one < dim + 1:
            already_seen = sum(f_vector[j] for j in range(dimension_two + 1))
            already_seen_next = already_seen + f_vector[dimension_two + 1]

            if all_faces.dual:
                # If ``dual``, then ``all_faces`` has the dimensions reversed.
                all_faces.incidence_init(dim - 1 - dimension_two, dim - 1 - dimension_one)
            else:
                all_faces.incidence_init(dimension_one, dimension_two)

            # Get all incidences for fixed ``[dimension_one, dimension_two]``.
            while all_faces.next_incidence(&second, &first):
                if all_faces.dual:
                    # If ``dual``, then ``second`` and ``first are flipped.
                    second += already_seen
                    first += already_seen_next
                    incidences.add(second, first)
                else:
                    second += already_seen_next
                    first += already_seen
                    incidences.add(first, second)

                sig_check()

            # Increase dimensions.
            dimension_one += 1
            dimension_two = dimension_one - 1

        # Success, persist the data.
        self._face_lattice_incidences = incidences

    def _record_all_faces(self):
        r"""
        Initialize :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.polyhedron_faces_lattice.PolyhedronFaceLattice` for the polyhedron.

        Record and sort all faces of the polyhedron in that class.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces()

        TESTS::

            sage: # needs sage.combinat
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(),
            ....:              face.ambient_Hrepresentation()) for face in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple(
            ....:     (C.face_by_face_lattice_index(i).ambient_Vrepresentation(),
            ....:      C.face_by_face_lattice_index(i).ambient_Hrepresentation())
            ....:     for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(),face.ambient_Hrepresentation()) for face in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_by_face_lattice_index(i).ambient_Vrepresentation(),
            ....:               C.face_by_face_lattice_index(i).ambient_Hrepresentation()) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(),face.ambient_Hrepresentation()) for face in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_by_face_lattice_index(i).ambient_Vrepresentation(),
            ....:               C.face_by_face_lattice_index(i).ambient_Hrepresentation()) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_generator()
            sage: tup = tuple((face.ambient_Vrepresentation(),face.ambient_Hrepresentation()) for face in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_by_face_lattice_index(i).ambient_Vrepresentation(),
            ....:               C.face_by_face_lattice_index(i).ambient_Hrepresentation()) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True
        """
        if self._all_faces:
            return  # Have recorded all faces already.

        self._all_faces = PolyhedronFaceLattice(self)
        if self._all_faces is None:
            raise RuntimeError("could not determine a list of all faces")
