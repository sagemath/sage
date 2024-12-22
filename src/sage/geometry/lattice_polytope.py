r"""
Lattice and reflexive polytopes

This module provides tools for work with lattice and reflexive
polytopes. A *convex polytope* is the convex hull of finitely many
points in `\RR^n`. The dimension `n` of a
polytope is the smallest `n` such that the polytope can be
embedded in `\RR^n`.

A *lattice polytope* is a polytope whose vertices all have integer
coordinates.

If `L` is a lattice polytope, the dual polytope of
`L` is

.. MATH::

       \{y \in \ZZ^n :   x\cdot y \geq -1 \text{ all } x \in L\}


A *reflexive polytope* is a lattice polytope, such that its polar
is also a lattice polytope, i.e. it is bounded and has vertices with
integer coordinates.

This Sage module uses Package for Analyzing Lattice Polytopes
(PALP), which is a program written in C by Maximilian Kreuzer and
Harald Skarke, which is freely available under the GNU license
terms at http://hep.itp.tuwien.ac.at/~kreuzer/CY/. Moreover, PALP is
included standard with Sage.

PALP is described in the paper :arxiv:`math.SC/0204356`. Its distribution
also contains the application ``nef.x``, which was created by Erwin
Riegler and computes nef-partitions and Hodge data for toric
complete intersections.

ACKNOWLEDGMENT: polytope.py module written by William Stein was
used as an example of organizing an interface between an external
program and Sage. William Stein also helped Andrey Novoseltsev with
debugging and tuning of this module.

Robert Bradshaw helped Andrey Novoseltsev to realize plot3d
function.

.. NOTE::

   IMPORTANT: PALP requires some parameters to be determined during
   compilation time, i.e., the maximum dimension of polytopes, the
   maximum number of points, etc. These limitations may lead to errors
   during calls to different functions of these module.  Currently, a
   :exc:`ValueError` exception will be raised if the output of ``poly.x``
   or ``nef.x`` is empty or contains the exclamation mark. The error
   message will contain the exact command that caused an error, the
   description and vertices of the polytope, and the obtained output.

Data obtained from PALP and some other data is cached and most
returned values are immutable. In particular, you cannot change the
vertices of the polytope or their order after creation of the
polytope.

If you are going to work with large sets of data, take a look at
``all_*`` functions in this module. They precompute different data
for sequences of polynomials with a few runs of external programs.
This can significantly affect the time of future computations. You
can also use dump/load, but not all data will be stored (currently
only faces and the number of their internal and boundary points are
stored, in addition to polytope vertices and its polar).

AUTHORS:

- Andrey Novoseltsev (2007-01-11): initial version

- Andrey Novoseltsev (2007-01-15): ``all_*`` functions

- Andrey Novoseltsev (2008-04-01): second version, including:

    - dual nef-partitions and necessary convex_hull and minkowski_sum

    - built-in sequences of 2- and 3-dimensional reflexive polytopes

    - plot3d, skeleton_show

- Andrey Novoseltsev (2009-08-26): dropped maximal dimension requirement

- Andrey Novoseltsev (2010-12-15): new version of nef-partitions

- Andrey Novoseltsev (2013-09-30): switch to PointCollection.

- Maximilian Kreuzer and Harald Skarke: authors of PALP (which was
  also used to obtain the list of 3-dimensional reflexive polytopes)

- Erwin Riegler: the author of nef.x
"""

# ****************************************************************************
#       Copyright (C) 2007-2017 Andrey Novoseltsev <novoselt@gmail.com>
#                     2007-2013 William Stein <wstein@gmail.com>
#                     2009      Mike Hansen
#                     2009-2020 John H. Palmieri
#                     2010-2014 Volker Braun
#                     2012      Samuel Gonshaw
#                     2013      Jan Keitel
#                     2014      André Apitzsch
#                     2014      Wilfried Luebbe
#                     2015-2022 Frédéric Chapoton
#                     2015      Ursula Whitcher
#                     2016      Jori Mäntysalo
#                     2017      Travis Scrimshaw
#                     2018      Christian Stump
#                     2018      Vincent Klein
#                     2019      Vincent Delecroix
#                     2019      Jonathan Kliem
#                     2020      Samuel Lelièvre
#                     2021-2023 Matthias Koeppe
#                     2022      David Coudert
#                     2023      Luze Xu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.posets.posets', 'FinitePoset')
from sage.arith.misc import GCD as gcd
from sage.features.databases import DatabaseReflexivePolytopes
from sage.geometry.cone import _ambient_space_point, integral_length
lazy_import('sage.geometry.hasse_diagram', 'lattice_from_incidences')
from sage.geometry.point_collection import (PointCollection,
                                            read_palp_point_collection)
from sage.geometry.toric_lattice import ToricLattice, ToricLattice_generic
lazy_import('sage.groups.perm_gps.permgroup_named', 'SymmetricGroup')

from sage.features import PythonModule
from sage.features.palp import PalpExecutable
lazy_import('ppl', ['C_Polyhedron', 'Generator_System', 'Linear_Expression'],
                    feature=PythonModule("ppl", spkg='pplpy', type='standard'))
lazy_import('ppl', 'point', as_='PPL_point',
                    feature=PythonModule("ppl", spkg='pplpy', type='standard'))

from sage.matrix.constructor import matrix
from sage.structure.element import Matrix
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.temporary_file import tmp_filename
from sage.modules.free_module_element import vector
lazy_import('sage.numerical.mip', 'MixedIntegerLinearProgram')
lazy_import("sage.plot.plot3d.index_face_set", "IndexFaceSet")
lazy_import("sage.plot.plot3d.all", ["line3d", "point3d"])
lazy_import("sage.plot.plot3d.shapes2", "text3d")
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.set import Set_generic
from sage.structure.all import Sequence
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp
from sage.geometry.convex_set import ConvexSet_compact
import sage.geometry.abc

from copy import copy
from collections.abc import Hashable
from copyreg import constructor as copyreg_constructor
import os
import shlex
from subprocess import Popen, PIPE
from warnings import warn
from functools import reduce
from io import IOBase, StringIO


class SetOfAllLatticePolytopesClass(Set_generic):
    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: lattice_polytope.SetOfAllLatticePolytopesClass()._repr_()
            'Set of all Lattice Polytopes'
        """
        return "Set of all Lattice Polytopes"

    def __call__(self, x):
        r"""
        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: lattice_polytope.SetOfAllLatticePolytopesClass().__call__(o)
            3-d reflexive polytope in 3-d lattice M
        """
        if isinstance(x, LatticePolytopeClass):
            return x
        raise TypeError


SetOfAllLatticePolytopes = SetOfAllLatticePolytopesClass()


def LatticePolytope(data, compute_vertices=True, n=0, lattice=None):
    r"""
    Construct a lattice polytope.

    INPUT:

    - ``data`` -- points spanning the lattice polytope, specified as one of:

        * a :class:`point collection
          <sage.geometry.point_collection.PointCollection>` (this is the
          preferred input and it is the quickest and the most memory efficient
          one);

        * an iterable of iterables (for example, a list of vectors)
          defining the point coordinates;

        * a file with matrix data, opened for reading, or

        * a filename of such a file, see
          :func:`~sage.geometry.point_collection.read_palp_point_collection`
          for the file format;

    - ``compute_vertices`` -- boolean (default: ``True``); if ``True``, the
      convex hull of the given points will be computed for determining
      vertices. Otherwise, the given points must be vertices.

    - ``n`` -- integer (default: 0); if ``data`` is a name of a file,
      that contains data blocks for several polytopes, the ``n``-th block
      will be used

    - ``lattice`` -- the ambient lattice of the polytope. If not given, a
      suitable lattice will be determined automatically, most likely the
      :class:`toric lattice <sage.geometry.toric_lattice.ToricLatticeFactory>`
      `M` of the appropriate dimension.

    OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

    EXAMPLES::

        sage: points = [(1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: p = LatticePolytope(points)
        sage: p
        3-d reflexive polytope in 3-d lattice M
        sage: p.vertices()
        M( 1,  0,  0),
        M( 0,  1,  0),
        M( 0,  0,  1),
        M(-1,  0,  0),
        M( 0, -1,  0),
        M( 0,  0, -1)
        in 3-d lattice M

    We draw a pretty picture of the polytope in 3-dimensional space::

        sage: p.plot3d().show()                                                         # needs palp sage.plot

    Now we add an extra point, which is in the interior of the
    polytope...

    ::

        sage: points.append((0,0,0))
        sage: p = LatticePolytope(points)
        sage: p.nvertices()
        6

    You can suppress vertex computation for speed but this can lead to
    mistakes::

        sage: p = LatticePolytope(points, compute_vertices=False)
        ...
        sage: p.nvertices()
        7

    Given points must be in the lattice::

        sage: LatticePolytope([[1/2], [3/2]])
        Traceback (most recent call last):
        ...
        ValueError: points
        [[1/2], [3/2]]
        are not in 1-d lattice M!

    But it is OK to create polytopes of non-maximal dimension::


        sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,0),
        ....:       (-1,0,0), (0,-1,0), (0,0,0), (0,0,0)])
        sage: p
        2-d lattice polytope in 3-d lattice M
        sage: p.vertices()
        M(-1,  0, 0),
        M( 0, -1, 0),
        M( 1,  0, 0),
        M( 0,  1, 0)
        in 3-d lattice M

    An empty lattice polytope can be considered as well::

        sage: p = LatticePolytope([], lattice=ToricLattice(3).dual()); p
        -1-d lattice polytope in 3-d lattice M
        sage: p.lattice_dim()
        3
        sage: p.npoints()
        0
        sage: p.nfacets()
        0
        sage: p.points()
        Empty collection
        in 3-d lattice M
        sage: p.faces()                                                                 # needs sage.graphs
        ((-1-d lattice polytope in 3-d lattice M,),)
    """
    if isinstance(data, LatticePolytopeClass):
        data = data._vertices
        compute_vertices = False
    if (isinstance(data, PointCollection) and
        (lattice is None or lattice is data.module())):
        return LatticePolytopeClass(data, compute_vertices)
    if isinstance(data, str):
        f = open(data)
        skip_palp_matrix(f, n)
        data = read_palp_point_collection(data)
        f.close()
    if isinstance(data, (IOBase, StringIO)):
        data = read_palp_point_collection(data)
    if not isinstance(data, PointCollection) and not isinstance(data, (list, tuple)):
        try:
            data = list(data)
        except TypeError:
            raise TypeError("cannot construct a polytope from\n%s" % data)
    if lattice is None:
        if not data:
            raise ValueError("lattice must be given explicitly for "
                             "empty polytopes!")
        try:
            if isinstance(data[0].parent(), ToricLattice_generic):
                lattice = data[0].parent()
        except AttributeError:
            pass
    if lattice is None:
        try:
            lattice = ToricLattice(len(data[0])).dual()
        except TypeError:
            raise TypeError("cannot construct a polytope from\n%s" % data)
    try:
        data = tuple(map(lattice, data))
    except TypeError:
        raise ValueError("points\n%s\nare not in %s!" % (data, lattice))
    for p in data:
        p.set_immutable()
    data = PointCollection(data, lattice)
    return LatticePolytopeClass(data, compute_vertices)


copyreg_constructor(LatticePolytope)   # "safe for unpickling"


def ReflexivePolytope(dim, n):
    r"""
    Return the `n`-th 2- or 3-dimensional reflexive polytope.

    .. NOTE::

       #. Numeration starts with zero: `0 \leq n \leq 15` for `{\rm dim} = 2`
          and `0 \leq n \leq 4318` for `{\rm dim} = 3`.

       #. During the first call, all reflexive polytopes of requested
          dimension are loaded and cached for future use, so the first
          call for 3-dimensional polytopes can take several seconds,
          but all consecutive calls are fast.

       #. Equivalent to ``ReflexivePolytopes(dim)[n]`` but checks bounds
          first.

    EXAMPLES:

    The 3rd 2-dimensional polytope is "the diamond"::

        sage: ReflexivePolytope(2, 3)
        2-d reflexive polytope #3 in 2-d lattice M
        sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
        M( 1,  0),
        M( 0,  1),
        M( 0, -1),
        M(-1,  0)
        in 2-d lattice M

    There are 16 reflexive polygons and numeration starts with 0::

        sage: ReflexivePolytope(2,16)
        Traceback (most recent call last):
        ...
        ValueError: there are only 16 reflexive polygons!

    It is not possible to load a 4-dimensional polytope in this way::

        sage: ReflexivePolytope(4,16)
        Traceback (most recent call last):
        ...
        NotImplementedError: only 2- and 3-dimensional reflexive polytopes are available!
    """
    if dim == 2:
        if n > 15:
            raise ValueError("there are only 16 reflexive polygons!")
        return ReflexivePolytopes(2)[n]
    elif dim == 3:
        if n > 4318:
            raise ValueError("there are only 4319 reflexive 3-polytopes!")
        return ReflexivePolytopes(3)[n]
    else:
        raise NotImplementedError("only 2- and 3-dimensional reflexive polytopes are available!")


# Sequences of reflexive polytopes
_rp = [None] * 4


def ReflexivePolytopes(dim):
    r"""
    Return the sequence of all 2- or 3-dimensional reflexive polytopes.

    .. NOTE::

       During the first call the database is loaded and cached for
       future use, so repetitive calls will return the same object in
       memory.

    INPUT:

    - ``dim`` -- integer (2 or 3); dimension of required reflexive polytopes

    OUTPUT: list of lattice polytopes

    EXAMPLES:

    There are 16 reflexive polygons::

        sage: len(ReflexivePolytopes(2))
        16

    It is not possible to load 4-dimensional polytopes in this way::

        sage: ReflexivePolytopes(4)
        Traceback (most recent call last):
        ...
        NotImplementedError: only 2- and 3-dimensional reflexive polytopes are available!
    """
    global _rp
    if dim not in [2, 3]:
        raise NotImplementedError("only 2- and 3-dimensional reflexive polytopes are available!")
    if _rp[dim] is None:
        db = DatabaseReflexivePolytopes()
        rp = read_all_polytopes(
                os.path.join(os.path.dirname(db.absolute_filename()),
                             f'reflexive_polytopes_{dim}d'))
        for n, p in enumerate(rp):
            # Data files have normal form of reflexive polytopes
            p.normal_form.set_cache(p._vertices)
            p.index.set_cache(n)
            # Prevents dimension computation later
            p.dim.set_cache(dim)
        _rp[dim] = rp
    return _rp[dim]


def is_LatticePolytope(x):
    r"""
    Check if ``x`` is a lattice polytope.

    INPUT:

    - ``x`` -- anything

    OUTPUT:

    - ``True`` if ``x`` is a :class:`lattice polytope <LatticePolytopeClass>`,
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.lattice_polytope import is_LatticePolytope
        sage: is_LatticePolytope(1)
        doctest:warning...
        DeprecationWarning: is_LatticePolytope is deprecated, use isinstance instead
        See https://github.com/sagemath/sage/issues/34307 for details.
        False
        sage: p = LatticePolytope([(1,0), (0,1), (-1,-1)])
        sage: p                                                                         # needs palp
        2-d reflexive polytope #0 in 2-d lattice M
        sage: is_LatticePolytope(p)
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(34307, "is_LatticePolytope is deprecated, use isinstance instead")
    return isinstance(x, LatticePolytopeClass)


@richcmp_method
class LatticePolytopeClass(ConvexSet_compact, Hashable, sage.geometry.abc.LatticePolytope):
    r"""
    Create a lattice polytope.

    .. WARNING::

        This class does not perform any checks of correctness of input nor
        does it convert input into the standard representation. Use
        :func:`LatticePolytope` to construct lattice polytopes.

    Lattice polytopes are immutable, but they cache most of the returned values.

    INPUT:

    The input can be either:

    - ``points`` -- :class:`~sage.geometry.point_collection.PointCollection`

    - ``compute_vertices`` -- boolean

    or (these parameters must be given as keywords):

    - ``ambient`` -- ambient structure, this polytope *must be a face of*
      ``ambient``

    - ``ambient_vertex_indices`` -- increasing list or tuple of integers,
      indices of vertices of ``ambient`` generating this polytope

    - ``ambient_facet_indices`` -- increasing list or tuple of integers,
      indices of facets of ``ambient`` generating this polytope

    OUTPUT: lattice polytope

    .. NOTE::

        Every polytope has an ambient structure. If it was not specified, it is
        this polytope itself.
    """

    def __init__(self, points=None, compute_vertices=None,
                 ambient=None, ambient_vertex_indices=None,
                 ambient_facet_indices=None):
        r"""
        Construct a lattice polytope.

        See :func:`LatticePolytope` for documentation.

        TESTS::

            sage: LatticePolytope([(1,2,3), (4,5,6)]) # indirect test
            1-d lattice polytope in 3-d lattice M
            sage: TestSuite(_).run()
        """
        if ambient is None:
            self._ambient = self
            if compute_vertices:
                P = C_Polyhedron(Generator_System(
                    [PPL_point(Linear_Expression(p, 0)) for p in points]))
                self._PPL.set_cache(P)
                vertices = P.minimized_generators()
                if len(vertices) != len(points):
                    M = points.module()
                    points = tuple(M(v.coefficients()) for v in vertices)
                    for p in points:
                        p.set_immutable()
                    points = PointCollection(points, M)
            self._vertices = points
            self._ambient_vertex_indices = tuple(range(self.nvertices()))
            self._ambient_facet_indices = ()
        else:
            self._ambient = ambient
            self._ambient_vertex_indices = tuple(ambient_vertex_indices)
            self._ambient_facet_indices = tuple(ambient_facet_indices)
            self._vertices = ambient.vertices(self._ambient_vertex_indices)

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: sage_input(p, verify=True)
            # Verified
            LatticePolytope(sage.geometry.point_collection.PointCollection((vector(ZZ, [1, 0]),
                                                                            vector(ZZ, [0, 1]),
                                                                            vector(ZZ, [-1, 0]),
                                                                            vector(ZZ, [0, -1]))),
                            compute_vertices=False)
        """
        if self._ambient is not self:
            raise NotImplementedError
        return sib.name('LatticePolytope')(sib(self._vertices), compute_vertices=False)

    def __contains__(self, point):
        r"""
        Check if ``point`` is contained in ``self``.

        See :meth:`_contains` (which is called by this function) for
        documentation.

        TESTS::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: (1,0) in p
            True
            sage: [1,0] in p
            True
            sage: (-2,0) in p
            False
        """
        return self._contains(point)

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything

        .. NOTE::

            Two lattice polytopes are equal if they have the same vertices
            listed in the same order.

        TESTS::

            sage: p1 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p2 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p3 = LatticePolytope([(0,1), (1,0), (-1,-1)])
            sage: p1 == p1
            True
            sage: p1 == p2
            True
            sage: p1 is p2
            False
            sage: p1 == p3
            False
            sage: p1 == 0
            False
            sage: p1 < p2
            False
            sage: p2 < p1
            False
            sage: p1 < p3
            False
            sage: p3 < p1
            True
            sage: p1 <= p2
            True
            sage: p2 <= p1
            True
            sage: p1 <= p3
            False
            sage: p3 <= p1
            True
            sage: p1 > p2
            False
            sage: p2 > p1
            False
            sage: p1 > p3
            True
            sage: p3 > p1
            False
            sage: p1 >= p2
            True
            sage: p2 >= p1
            True
            sage: p1 >= p3
            True
            sage: p3 >= p1
            False
        """
        if not isinstance(other, LatticePolytopeClass):
            return NotImplemented
        return richcmp(self._vertices, other._vertices, op)

    @cached_method
    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT: integer

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: hash(o) == hash(o)
            True
        """
        # FIXME: take into account other things that may be preset?..
        return hash(self._vertices)

    def __reduce__(self):
        r"""
        Reduction function. Does not store data that can be relatively fast
        recomputed.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices() == loads(o.dumps()).vertices()
            True
        """
        state = self.__dict__.copy()
        state.pop('_vertices')
        state.pop('_distances', None)
        state.pop('_skeleton', None)
        try:
            state['_npoints'] = len(state['_points'])
            state.pop('_points')
        except KeyError:
            pass
        return (LatticePolytope, (self._vertices, None, False), state)

    def __setstate__(self, state):
        r"""
        Restore the state of pickled polytope.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices() == loads(o.dumps()).vertices()
            True
        """
        self.__dict__.update(state)

    def _compute_embedding(self):
        r"""
        Compute embedding data for this polytope.

        Useful only if the dimension of this polytope is not equal to its
        ambient dimension.

        TESTS::

            sage: p = LatticePolytope(([1], [2], [3]))
            sage: hasattr(p, "_sublattice_polytope")
            False
            sage: p._compute_embedding()
            sage: p._sublattice_polytope
            1-d lattice polytope in 1-d lattice M
        """
        if hasattr(self, "_sublattice_polytope"):
            return
        points = self._vertices
        if not points:  # the empty lattice polytope
            return
        p0 = self._shift_vector = points[0]
        points = [point - p0 for point in points]
        H = self._sublattice = self.lattice().submodule(points).saturation()
        self._sublattice_polytope = LatticePolytope([H.coordinates(point)
                                                     for point in points])
        M = self._embedding_matrix = H.basis_matrix().transpose()
        # In order to use facet normals obtained from subpolytopes, we
        # need the following (see Issue #9188).
        # Basis for the ambient space with spanned subspace in front
        basis = M.columns() + M.integer_kernel().basis()
        # Let's represent it as columns of a matrix
        basis = matrix(basis).transpose()
        # Absolute value helps to keep normals "inner"
        self._dual_embedding_scale = abs(basis.det())
        dualbasis = matrix(ZZ, self._dual_embedding_scale * basis.inverse())
        self._dual_embedding_matrix = dualbasis.submatrix(0,0,M.ncols())

    def _compute_facets(self):
        r"""
        Compute and cache equations of facets of ``self``.

        TESTS::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p._compute_facets()
            sage: p._facet_normals
            N( 1,  1, 0),
            N( 1, -1, 0),
            N(-1, -1, 0),
            N(-1,  1, 0)
            in 3-d lattice N

        Check that :issue:`28741` is fixed::

            sage: # needs sage.graphs
            sage: p = LatticePolytope([], lattice=ToricLattice(3).dual()); p
            -1-d lattice polytope in 3-d lattice M
            sage: a = p.faces()[0][0]
            sage: p = LatticePolytope([], lattice=ToricLattice(3).dual()); p
            -1-d lattice polytope in 3-d lattice M
            sage: a = p.faces()[0][0]; a
            -1-d lattice polytope in 3-d lattice M
            sage: a.facet_normals()
            Empty collection
            in 3-d lattice N
            sage: a
            -1-d lattice polytope in 3-d lattice M
        """
        assert not hasattr(self, "_facet_normals")
        N = self.dual_lattice()
        normals = []
        constants = []
        for c in self._PPL().minimized_constraints():
            if c.is_inequality():
                n = N.element_class(N, [Integer(mpz) for mpz in c.coefficients()])
                n.set_immutable()
                normals.append(n)
                constants.append(Integer(c.inhomogeneous_term()))
        # Sort normals if facets are vertices
        if (self.dim() == 1
            and normals[0] * self.vertex(0) + constants[0] != 0):
            normals = (normals[1], normals[0])
            constants = (constants[1], constants[0])
        self._facet_normals = PointCollection(normals, N)
        # vector(ZZ, constants) is slow
        self._facet_constants = (ZZ**len(constants))(constants)
        self._facet_constants.set_immutable()
        self.is_reflexive.set_cache(self.dim() == self.lattice_dim() and
                                    all(c == 1 for c in constants))
        if self.is_reflexive():
            polar = LatticePolytope(
                self._facet_normals, compute_vertices=False)
            polar.dim.set_cache(self.dim())
            polar.is_reflexive.set_cache(True)
            polar._polar = self
            self._polar = polar
            polar._facet_normals = self._vertices
            ones = [1] * self.nvertices()
            ones = (ZZ**len(ones))(ones)
            ones.set_immutable()
            polar._facet_constants = ones

    def _compute_hodge_numbers(self):
        r"""
        Compute Hodge numbers for the current nef_partitions.

        This function (currently) always raises an exception directing to
        use another way for computing Hodge numbers.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._compute_hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
        """
        raise NotImplementedError("use nef_partitions(hodge_numbers=True)!")

    def _contains(self, point, region='whole polytope'):
        r"""
        Check if ``point`` is contained in ``self``.

        This function is called by :meth:`__contains__` and :meth:`contains`
        to ensure the same call depth for warning messages.

        INPUT:

        - ``point`` -- an attempt will be made to convert it into a
          single element of the ambient space of ``self``; if it fails,
          ``False`` is returned

        - ``region`` -- string; can be either ``'whole polytope'`` (default),
          ``'interior'``, or ``'relative interior'``

        OUTPUT:

        - ``True`` if ``point`` is contained in the specified ``region`` of
          ``self``, ``False`` otherwise

        TESTS::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: p._contains((1,0))
            True
        """
        try:
            point = _ambient_space_point(self, point)
        except TypeError as ex:
            if str(ex).endswith("have incompatible lattices"):
                warn("you have checked if a cone contains a point "
                     "from an incompatible lattice, this is False",
                     stacklevel=3)
            return False

        if region not in ("whole polytope", "relative interior", "interior"):
            raise ValueError("%s is an unknown region" % region)
        if region == "interior" and self.dim() < self.lattice_dim():
            return False
        need_strict = region.endswith("interior")
        N = self.dual_lattice()
        for c in self._PPL().minimized_constraints():
            pr = N([Integer(mpz) for mpz in c.coefficients()]) * point + Integer(c.inhomogeneous_term())
            if c.is_equality():
                if pr != 0:
                    return False
            elif pr < 0 or need_strict and pr == 0:
                return False
        return True

    def _embed(self, data):
        r"""
        Embed given point(s) into the ambient space of this polytope.

        INPUT:

        - ``data`` -- point or matrix of points (as columns) in the affine
          subspace spanned by this polytope

        OUTPUT:

        The same point(s) in the coordinates of the ambient space of
        this polytope.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._embed(o.vertices()) == o.vertices()
            True
            sage: m = matrix(ZZ, 3)
            sage: m[0, 0] = 1
            sage: m[1, 1] = 1
            sage: p = o.affine_transform(m)
            sage: p._embed((0,0))
            M(-1, 0, 0)
        """
        if self.lattice_dim() == self.dim():
            return data
        self._compute_embedding()
        M = self.lattice()
        if isinstance(data, PointCollection):
            r = [M(self._embedding_matrix * point + self._shift_vector)
                 for point in data]
            for point in r:
                point.set_immutable()
            return PointCollection(r, M)
        elif isinstance(data, Matrix):
            r = self._embedding_matrix * data
            for i, col in enumerate(r.columns(copy=False)):
                r.set_column(i, col + self._shift_vector)
            return r
        else:
            return M(self._embedding_matrix * vector(QQ, data) +
                     self._shift_vector)

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        OUTPUT: string

        EXAMPLES:

        Arbitrary lattice polytopes are printed as `\Delta^d`, where `d` is
        the (actual) dimension of the polytope::

            sage: LatticePolytope([(1,1), (0,0)])._latex_()
            '\\Delta^{1}'

        For 2- and 3-d reflexive polytopes the index in the internal database
        appears as a subscript::

            sage: print(ReflexivePolytope(2, 3)._latex_())
            \Delta^{2}_{3}
        """
        result = r"\Delta^{%d}" % self.dim()
        if self.dim() in (2, 3) and self.is_reflexive():
            result += "_{%d}" % self.index()
        return result

    def _palp(self, command, reduce_dimension=False):
        r"""
        Run ``command`` on vertices of this polytope.

        Returns the output of ``command`` as a string.

        .. NOTE::

          PALP cannot be called for polytopes that do not span the ambient space.
          If you specify ``reduce_dimension=True`` argument, PALP will be
          called for vertices of this polytope in some basis of the affine space
          it spans.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._palp("poly.x -f")                                                  # needs palp
            'M:7 6 N:27 8 Pic:17 Cor:0\n'
            sage: print(o._palp("nef.x -f -N -p"))  # random time information           # needs palp
            M:27 8 N:7 6  codim=2 #part=5
            H:[0] P:0 V:2 4 5       0sec  0cpu
            H:[0] P:2 V:3 4 5       0sec  0cpu
            H:[0] P:3 V:4 5       0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu

            sage: p = LatticePolytope([[1]])
            sage: p._palp("poly.x -f")                                                  # needs palp
            Traceback (most recent call last):
            ...
            ValueError: Cannot run "poly.x -f" for the zero-dimensional polytope!
            Polytope: 0-d lattice polytope in 1-d lattice M

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p._palp("poly.x -f")                                                  # needs palp
            Traceback (most recent call last):
            ...
            ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!
            sage: p._palp("poly.x -f", reduce_dimension=True)                           # needs palp
            'M:5 4 F:4\n'
        """
        if self.dim() <= 0:
            raise ValueError(("Cannot run \"%s\" for the zero-dimensional "
                + "polytope!\nPolytope: %s") % (command, self))
        if self.dim() < self.lattice_dim() and not reduce_dimension:
            raise ValueError(("Cannot run PALP for a %d-dimensional polytope " +
            "in a %d-dimensional space!") % (self.dim(), self.lattice_dim()))
        fn = _palp(command, [self], reduce_dimension)
        with open(fn) as f:
            result = f.read()
        os.remove(fn)
        if (not result or
            "!" in result or
            "failed." in result or
            "increase" in result or
            "Unable" in result):
            lines = ["Error executing '%s' for the given polytope!" % command,
                     "Output:", result]
            raise ValueError("\n".join(lines))
        return result

    @cached_method
    def _PPL(self):
        r"""
        Return the Parma Polyhedra Library (PPL) representation of ``self``.

        OUTPUT: :class:`~ppl.polyhedron.C_Polyhedron`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._PPL()
            A 3-dimensional polyhedron in QQ^3
            defined as the convex hull of 6 points
            sage: o._PPL().minimized_generators()
            Generator_System {point(-1/1, 0/1, 0/1),
                              point(0/1, -1/1, 0/1),
                              point(0/1, 0/1, -1/1),
                              point(0/1, 0/1, 1/1),
                              point(0/1, 1/1, 0/1),
                              point(1/1, 0/1, 0/1)}
            sage: o._PPL().minimized_constraints()
            Constraint_System {x0-x1-x2+1>=0,
                               x0+x1-x2+1>=0,
                               x0+x1+x2+1>=0,
                               x0-x1+x2+1>=0,
                               -x0-x1+x2+1>=0,
                               -x0-x1-x2+1>=0,
                               -x0+x1-x2+1>=0,
                               -x0+x1+x2+1>=0}
        """
        P = C_Polyhedron(Generator_System(
            [PPL_point(Linear_Expression(v, 0)) for v in self.vertices()]))
        return P

    def _pullback(self, data):
        r"""
        Pull back given point(s) to the affine subspace spanned by this polytope.

        INPUT:

        - ``data`` -- rational point or matrix of points (as columns) in the
          ambient space

        OUTPUT:

        The same point(s) in the coordinates of the affine subspace
        space spanned by this polytope.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._pullback(o.vertices().column_matrix()) == o.vertices().column_matrix()
            True
            sage: m = matrix(ZZ, 3)
            sage: m[0, 0] = 1
            sage: m[1, 1] = 1
            sage: p = o.affine_transform(m)
            sage: p._pullback((0, 0, 0))
            [1, 0]
        """
        if self.lattice_dim() == self.dim():
            return data
        self._compute_embedding()
        if data is self._vertices:
            return self._sublattice_polytope._vertices
        if isinstance(data, PointCollection):
            r = [self._pullback(point) for point in data]
            for point in r:
                point.set_immutable()
            return PointCollection(r, self._sublattice)
        if isinstance(data, Matrix):
            r = matrix([self._pullback(col)
                    for col in data.columns(copy=False)]).transpose()
            return r
        data = vector(QQ, data)
        return self._sublattice.coordinates(data - self._shift_vector)

    def _read_equations(self, data):
        r"""
        Read equations of facets/vertices of polar polytope from string or
        file.

        TESTS:

        For a reflexive polytope construct the polar polytope::

            sage: # needs palp
            sage: p = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1, -1)
            in 2-d lattice M
            sage: s = p.poly_x("e")
            sage: print(s)
            3 2  Vertices of P-dual <-> Equations of P
               2  -1
              -1   2
              -1  -1
            sage: "_polar" in p.__dict__
            False
            sage: p._read_equations(s)
            sage: p._polar._vertices
            N( 2, -1),
            N(-1,  2),
            N(-1, -1)
            in 2-d lattice N

        For a non-reflexive polytope cache facet equations::

            sage: # needs palp
            sage: p = LatticePolytope([(1,0), (0,2), (-1,-3 )])
            sage: p.vertices()
            M( 1,  0),
            M( 0,  2),
            M(-1, -3)
            in 2-d lattice M
            sage: "_facet_normals" in p.__dict__
            False
            sage: "_facet_constants" in p.__dict__
            False
            sage: s = p.poly_x("e")
            sage: print(s)
            3 2  Equations of P
               5  -1     2
              -2  -1     2
              -3   2     3
            sage: p._read_equations(s)
            sage: p._facet_normals
            N( 5, -1),
            N(-2, -1),
            N(-3,  2)
            in 2-d lattice N
            sage: p._facet_constants
            (2, 2, 3)
        """
        if isinstance(data, str):
            f = StringIO(data)
            self._read_equations(f)
            f.close()
            return
        if self.is_reflexive.cache is not None:
            # If it is already known that this polytope is reflexive, its
            # polar (whose vertices are equations of facets of this one)
            # is already computed and there is no need to read equations
            # of facets of this polytope. Moreover, doing so can corrupt
            # data if this polytope was constructed as polar. Skip input.
            skip_palp_matrix(data)
            return
        pos = data.tell()
        line = data.readline()
        self.is_reflexive.set_cache(line.find("Vertices of P-dual") != -1)
        N = self.dual_lattice()
        if self.is_reflexive():
            data.seek(pos)
            polar = LatticePolytope(
                read_palp_point_collection(data, N), compute_vertices=False)
            polar.dim.set_cache(self.dim())
            polar.is_reflexive.set_cache(True)
            polar._constructed_as_polar = True
            polar._polar = self
            self._polar = polar
            self._facet_normals = polar._vertices
            polar._facet_normals = self._vertices
            ones = [1] * polar.nvertices()
            ones = (ZZ**len(ones))(ones)
            ones.set_immutable()
            self._facet_constants = ones
            ones = [1] * self.nvertices()
            ones = (ZZ**len(ones))(ones)
            ones.set_immutable()
            polar._facet_constants = ones
        else:
            normals = []
            constants = []
            for i in range(int(line.split()[0])):
                line = data.readline()
                numbers = [int(number) for number in line.split()]
                constants.append(numbers.pop())
                normals.append(N(numbers))
                normals[-1].set_immutable()
            self._facet_normals = PointCollection(normals, N)
            self._facet_constants = vector(ZZ, constants)
            self._facet_constants.set_immutable()

    def _read_nef_partitions(self, data):
        r"""
        Read nef-partitions of ``self`` from ``data``.

        INPUT:

        - ``data`` -- string or file

        OUTPUT: none

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: s = o.nef_x("-p -N -Lv")                                              # needs palp
            sage: print(s)  # random time values                                        # needs palp
            M:27 8 N:7 6  codim=2 #part=5
            3 6 Vertices in N-lattice:
                1    0    0   -1    0    0
                0    1    0    0   -1    0
                0    0    1    0    0   -1
            ------------------------------
                1    0    0    1    0    0  d=2  codim=2
                0    1    0    0    1    0  d=2  codim=2
                0    0    1    0    0    1  d=2  codim=2
             P:0 V:2 4 5   (0 2) (1 1) (2 0)     0sec  0cpu
             P:2 V:3 4 5   (1 1) (1 1) (1 1)     0sec  0cpu
             P:3 V:4 5   (0 2) (1 1) (1 1)     0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu

        We make a copy of the octahedron since direct use of this function may
        destroy cache integrity and lead so strange effects in other doctests::

            sage: o_copy = LatticePolytope(o.vertices())
            sage: "_nef_partitions" in o_copy.__dict__
            False
            sage: o_copy._read_nef_partitions(s)                                        # needs palp
            sage: o_copy._nef_partitions                                                # needs palp
            [Nef-partition {0, 1, 3} ⊔ {2, 4, 5},
             Nef-partition {0, 1, 2} ⊔ {3, 4, 5},
             Nef-partition {0, 1, 2, 3} ⊔ {4, 5}]
        """
        if isinstance(data, str):
            f = StringIO(data)
            self._read_nef_partitions(f)
            f.close()
            return
        nvertices = self.nvertices()
        data.readline() # Skip M/N information
        nef_vertices = read_palp_point_collection(data, self.lattice())
        if self.vertices() != nef_vertices:
            raise RuntimeError("nef.x changed the order of vertices!")
        line = data.readline()
        if line == "":
            raise ValueError("more data expected!")
        partitions = Sequence([], cr=True)
        while line and line.find("np=") == -1:
            if line.find("V:") == -1:
                line = data.readline()
                continue
            start = line.find("V:") + 2
            end = line.find("  ", start)  # Find DOUBLE space
            pvertices = Sequence(line[start:end].split(),int)
            partition = [0] * nvertices
            for v in pvertices:
                partition[v] = 1
            partition = NefPartition(partition, self)
            partition._is_product = line.find(" D ") != -1
            partition._is_projection = line.find(" DP ") != -1
            # Set the stuff
            start = line.find("H:")
            if start != -1:
                start += 2
                end = line.find("[", start)
                partition._hodge_numbers = tuple(int(h)
                                            for h in line[start:end].split())
            partitions.append(partition)
            line = data.readline()
        start = line.find("np=")
        if start == -1:
            raise ValueError("""Wrong data format, cannot find "np="!""")
        partitions.set_immutable()
        self._nef_partitions = partitions

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o._repr_()
            '3-d reflexive polytope in 3-d lattice M'
        """
        parts = ["%d-d" % self.dim()]
        if self.ambient() is self:
            parts.extend(["lattice", "polytope", "in"])
            try:
                if self.is_reflexive():
                    parts[1] = "reflexive"
                    if self.dim() == 2 or self.index.is_in_cache():
                        parts.insert(-1, "#%d" % self.index())
            except ValueError:
                pass
            if isinstance(self.lattice(), ToricLattice_generic):
                parts.append(str(self.lattice()))
            else:
                parts.append("%d-d lattice" % self.lattice_dim())
        else:
            parts.extend(["face of", str(self.ambient())])
        return " ".join(parts)

    def _sort_faces(self, faces):
        r"""
        Return sorted (if necessary) ``faces`` as a tuple.

        This function ensures that zero-dimensional faces are listed in
        agreement with the order of corresponding vertices and facets with
        facet normals.

        INPUT:

        - ``faces`` -- iterable of :class:`lattice polytopes
          <LatticePolytopeClass>`

        OUTPUT: :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: # indirect doctest
            sage: for i, face in enumerate(o.faces(0)):                                 # needs sage.graphs
            ....:     if face.vertex(0) != o.vertex(i):
            ....:         print("Wrong order!")
        """
        faces = tuple(faces)
        if len(faces) > 1: # Otherwise there is nothing to sort
            if faces[0].nvertices() == 1:
                faces = tuple(sorted(faces,
                                     key=lambda f: f._ambient_vertex_indices))
            elif faces[0].dim() == self.dim() - 1 and \
                    hasattr(self, "_facet_normals"):
                # If we already have facet normals, sort according to them
                faces = set(faces)
                sorted_faces = [None] * len(faces)
                for i, n in enumerate(self.facet_normals()):
                    for f in faces:
                        if set(n * f.vertices()) == set([- self.facet_constant(i)]):
                            sorted_faces[i] = f
                            faces.remove(f)
                            break
                faces = tuple(sorted_faces)
        return faces

    @cached_method
    def adjacent(self):
        r"""
        Return faces adjacent to ``self`` in the ambient face lattice.

        Two *distinct* faces `F_1` and `F_2` of the same face lattice are
        **adjacent** if all of the following conditions hold:

        * `F_1` and `F_2` have the same dimension `d`;

        * `F_1` and `F_2` share a facet of dimension `d-1`;

        * `F_1` and `F_2` are facets of some face of dimension `d+1`, unless
          `d` is the dimension of the ambient structure.

        OUTPUT: :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.adjacent()                                                          # needs sage.graphs
            ()
            sage: face = o.faces(1)[0]                                                  # needs sage.graphs
            sage: face.adjacent()                                                       # needs sage.graphs
            (1-d face of 3-d reflexive polytope in 3-d lattice M,
             1-d face of 3-d reflexive polytope in 3-d lattice M,
             1-d face of 3-d reflexive polytope in 3-d lattice M,
             1-d face of 3-d reflexive polytope in 3-d lattice M)
        """
        L = self._ambient.face_lattice()
        adjacent = set()
        for superface in self.facet_of():
            for facet in self.facets():
                adjacent.update(L.open_interval(facet, superface))
        adjacent.discard(self)
        return self._sort_faces(adjacent)

    def affine_transform(self, a=1, b=0):
        r"""
        Return a*P+b, where P is this lattice polytope.

        .. NOTE::

          #. While ``a`` and ``b`` may be rational, the final result must be a
             lattice polytope, i.e. all vertices must be integral.

          #. If the transform (restricted to this polytope) is bijective, facial
             structure will be preserved, e.g. the first facet of the image will
             be spanned by the images of vertices which span the first facet of
             the original polytope.

        INPUT:

        - ``a`` -- (default: 1) rational scalar or matrix

        - ``b`` -- (default: 0) rational scalar or vector, scalars are
          interpreted as vectors with the same components

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1,  0),
            M( 0, -1)
            in 2-d lattice M
            sage: o.affine_transform(2).vertices()
            M( 2,  0),
            M( 0,  2),
            M(-2,  0),
            M( 0, -2)
            in 2-d lattice M
            sage: o.affine_transform(1,1).vertices()
            M(2, 1),
            M(1, 2),
            M(0, 1),
            M(1, 0)
            in 2-d lattice M
            sage: o.affine_transform(b=1).vertices()
            M(2, 1),
            M(1, 2),
            M(0, 1),
            M(1, 0)
            in 2-d lattice M
            sage: o.affine_transform(b=(1, 0)).vertices()
            M(2,  0),
            M(1,  1),
            M(0,  0),
            M(1, -1)
            in 2-d lattice M
            sage: a = matrix(QQ, 2, [1/2, 0, 0, 3/2])
            sage: o.polar().vertices()
            N( 1,  1),
            N( 1, -1),
            N(-1, -1),
            N(-1,  1)
            in 2-d lattice N
            sage: o.polar().affine_transform(a, (1/2, -1/2)).vertices()
            M(1,  1),
            M(1, -2),
            M(0, -2),
            M(0,  1)
            in 2-d lattice M

        While you can use rational transformation, the result must be integer::

            sage: o.affine_transform(a)
            Traceback (most recent call last):
            ...
            ValueError: points
            [(1/2, 0), (0, 3/2), (-1/2, 0), (0, -3/2)]
            are not in 2-d lattice M!
        """
        new_vertices = self.vertices() * a
        if b in QQ:
            b = vector(QQ, [b]*new_vertices.ncols())
        else:
            b = vector(QQ, b)
        new_vertices = [c + b for c in new_vertices]
        r = LatticePolytope(new_vertices)
        if (a in QQ and a != 0) or r.dim() == self.dim():
            r._constructed_as_affine_transform = True
            if hasattr(self, "_constructed_as_affine_transform"):
                # Prevent long chains of "original-transform"
                r._original = self._original
            else:
                r._original = self
        return r

    def ambient(self):
        r"""
        Return the ambient structure of ``self``.

        OUTPUT: lattice polytope containing ``self`` as a face

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.ambient()
            3-d reflexive polytope in 3-d lattice M
            sage: o.ambient() is o
            True

            sage: # needs sage.graphs
            sage: face = o.faces(1)[0]
            sage: face
            1-d face of 3-d reflexive polytope in 3-d lattice M
            sage: face.ambient()
            3-d reflexive polytope in 3-d lattice M
            sage: face.ambient() is o
            True
        """
        return self._ambient

    def ambient_facet_indices(self):
        r"""
        Return indices of facets of the ambient polytope containing ``self``.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES:

        The polytope itself is not contained in any of its facets::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.ambient_facet_indices()
            ()

        But each of its other faces is contained in one or more facets::

            sage: # needs sage.graphs
            sage: face = o.faces(1)[0]
            sage: face.ambient_facet_indices()
            (4, 5)
            sage: face.vertices()
            M(1, 0, 0),
            M(0, 1, 0)
            in 3-d lattice M
            sage: o.facets()[face.ambient_facet_indices()[0]].vertices()
            M(1, 0,  0),
            M(0, 1,  0),
            M(0, 0, -1)
            in 3-d lattice M
        """
        return self._ambient_facet_indices

    @cached_method
    def ambient_point_indices(self):
        r"""
        Return indices of points of the ambient polytope contained in this one.

        OUTPUT:

        - :class:`tuple` of integers, the order corresponds to the order of
          points of this polytope.

        EXAMPLES::

            sage: cube = lattice_polytope.cross_polytope(3).polar()
            sage: face = cube.facets()[0]                                               # needs sage.graphs
            sage: face.ambient_point_indices()                                          # needs palp sage.graphs
            (4, 5, 6, 7, 8, 9, 10, 11, 12)
            sage: cube.points(face.ambient_point_indices()) == face.points()            # needs palp sage.graphs
            True
        """
        if self._ambient is self:
            return tuple(range(self.npoints()))
        points = self._ambient.points()
        point_to_index = {p: i for i, p in enumerate(points)}
        return tuple(point_to_index[p] for p in self.points())

    @cached_method
    def ambient_ordered_point_indices(self):
        r"""
        Return indices of points of the ambient polytope contained in this one.

        OUTPUT:

        - :class:`tuple` of integers such that ambient points in this order are
          geometrically ordered, e.g. for an edge points will appear from one
          end point to the other.

        EXAMPLES::

            sage: cube = lattice_polytope.cross_polytope(3).polar()
            sage: face = cube.facets()[0]                                               # needs sage.graphs
            sage: face.ambient_ordered_point_indices()                                  # needs palp sage.graphs
            (5, 8, 4, 9, 10, 11, 6, 12, 7)
            sage: cube.points(face.ambient_ordered_point_indices())                     # needs palp sage.graphs
            N(-1, -1, -1),
            N(-1, -1,  0),
            N(-1, -1,  1),
            N(-1,  0, -1),
            N(-1,  0,  0),
            N(-1,  0,  1),
            N(-1,  1, -1),
            N(-1,  1,  0),
            N(-1,  1,  1)
            in 3-d lattice N
        """
        if self._ambient is self:
            return tuple(range(self.npoints()))
        points = self._ambient.points()
        point_to_index = {p: i for i, p in enumerate(points)}
        return tuple(point_to_index[p] for p in sorted(self.points()))

    def ambient_vertex_indices(self):
        r"""
        Return indices of vertices of the ambient structure generating ``self``.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.ambient_vertex_indices()
            (0, 1, 2, 3, 4, 5)
            sage: face = o.faces(1)[0]                                                  # needs sage.graphs
            sage: face.ambient_vertex_indices()                                         # needs sage.graphs
            (0, 1)
        """
        return self._ambient_vertex_indices

    @cached_method
    def boundary_point_indices(self):
        r"""
        Return indices of (relative) boundary lattice points of this polytope.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES:

        All points but the origin are on the boundary of this square::

            sage: square = lattice_polytope.cross_polytope(2).polar()
            sage: square.points()                                                       # needs palp
            N( 1,  1),
            N( 1, -1),
            N(-1, -1),
            N(-1,  1),
            N(-1,  0),
            N( 0, -1),
            N( 0,  0),
            N( 0,  1),
            N( 1,  0)
            in 2-d lattice N
            sage: square.boundary_point_indices()                                       # needs palp
            (0, 1, 2, 3, 4, 5, 7, 8)

        For an edge the boundary is formed by the end points::

            sage: face = square.edges()[0]                                              # needs sage.graphs
            sage: face.points()                                                         # needs sage.graphs
            N(-1, -1),
            N(-1,  1),
            N(-1,  0)
            in 2-d lattice N
            sage: face.boundary_point_indices()                                         # needs sage.graphs
            (0, 1)
        """
        return tuple(i
                     for i, c in enumerate(self.distances().columns(copy=False))
                     if len(c.nonzero_positions()) < self.nfacets())

    def boundary_points(self):
        r"""
        Return (relative) boundary lattice points of this polytope.

        OUTPUT: a :class:`point collection <PointCollection>`

        EXAMPLES:

        All points but the origin are on the boundary of this square::

            sage: square = lattice_polytope.cross_polytope(2).polar()
            sage: square.boundary_points()                                              # needs palp
            N( 1,  1),
            N( 1, -1),
            N(-1, -1),
            N(-1,  1),
            N(-1,  0),
            N( 0, -1),
            N( 0,  1),
            N( 1,  0)
            in 2-d lattice N

        For an edge the boundary is formed by the end points::

            sage: face = square.edges()[0]                                              # needs sage.graphs
            sage: face.boundary_points()                                                # needs sage.graphs
            N(-1, -1),
            N(-1,  1)
            in 2-d lattice N
        """
        return self.points(self.boundary_point_indices())

    def contains(self, *args):
        r"""
        Check if a given point is contained in ``self``.

        INPUT:

        - an attempt will be made to convert all arguments into a
          single element of the ambient space of ``self``; if it fails,
          ``False`` will be returned

        OUTPUT:

        - ``True`` if the given point is contained in ``self``, ``False``
          otherwise

        EXAMPLES::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: p.contains(p.lattice()(1,0))
            True
            sage: p.contains((1,0))
            True
            sage: p.contains(1,0)
            True
            sage: p.contains((2,0))
            False
        """
        point = flatten(args)
        if len(point) == 1:
            point = point[0]
        return self._contains(point)

    @cached_method
    def dim(self):
        r"""
        Return the dimension of this polytope.

        EXAMPLES:

        We create a 3-dimensional octahedron and check its dimension::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.dim()
            3

        Now we create a 2-dimensional diamond in a 3-dimensional space::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.dim()
            2
            sage: p.lattice_dim()
            3
        """
        nv = self.nvertices()
        return self._PPL().affine_dimension() if nv > 3 else nv - 1

    def distances(self, point=None):
        r"""
        Return the matrix of distances for this polytope or distances for
        the given point.

        The matrix of distances m gives distances m[i,j] between the `i`-th
        facet (which is also the `i`-th vertex of the polar polytope in the
        reflexive case) and `j`-th point of this polytope.

        If point is specified, integral distances from the point to all
        facets of this polytope will be computed.

        EXAMPLES: The matrix of distances for a 3-dimensional octahedron::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.distances()                                                         # needs palp
            [2 0 0 0 2 2 1]
            [2 2 0 0 0 2 1]
            [2 2 2 0 0 0 1]
            [2 0 2 0 2 0 1]
            [0 0 2 2 2 0 1]
            [0 0 0 2 2 2 1]
            [0 2 0 2 0 2 1]
            [0 2 2 2 0 0 1]

        Distances from facets to the point (1,2,3)::

            sage: o.distances([1,2,3])
            (-3, 1, 7, 3, 1, -5, -1, 5)

        It is OK to use RATIONAL coordinates::

            sage: o.distances([1,2,3/2])
            (-3/2, 5/2, 11/2, 3/2, -1/2, -7/2, 1/2, 7/2)
            sage: o.distances([1,2,sqrt(2)])                                            # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(2) to an element of Rational Field

        Now we create a non-spanning polytope::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.distances()                                                         # needs palp
            [2 2 0 0 1]
            [2 0 0 2 1]
            [0 0 2 2 1]
            [0 2 2 0 1]
            sage: p.distances((1/2, 3, 0))                                              # needs palp
            (9/2, -3/2, -5/2, 7/2)

        This point is not even in the affine subspace of the polytope::

            sage: p.distances((1, 1, 1))                                                # needs palp
            (3, 1, -1, 1)
        """
        if point is not None:
            return (vector(QQ, point) * self.facet_normals() +
                    self.facet_constants())
        try:
            return self._distances
        except AttributeError:
            P = self.points()
            n = self.npoints()
            self._distances = matrix(ZZ, [F * P + vector(ZZ, [c]*n)
                for F, c in zip(self.facet_normals(), self.facet_constants())])
            self._distances.set_immutable()
            return self._distances

    @cached_method
    def dual(self):
        r"""
        Return the dual face under face duality of polar reflexive polytopes.

        This duality extends the correspondence between vertices and facets.

        OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

        EXAMPLES::

            sage: # needs sage.graphs
            sage: o = lattice_polytope.cross_polytope(4)
            sage: e = o.edges()[0]; e
            1-d face of 4-d reflexive polytope in 4-d lattice M
            sage: ed = e.dual(); ed
            2-d face of 4-d reflexive polytope in 4-d lattice N
            sage: ed.ambient() is e.ambient().polar()
            True
            sage: e.ambient_vertex_indices() == ed.ambient_facet_indices()
            True
            sage: e.ambient_facet_indices() == ed.ambient_vertex_indices()
            True
        """
        for f in self._ambient.polar().faces(codim=self.dim() + 1):
            if f._ambient_vertex_indices == self._ambient_facet_indices:
                f.dual.set_cache(self)
                return f

    @cached_method
    def dual_lattice(self):
        r"""
        Return the dual of the ambient lattice of ``self``.

        OUTPUT:

        - a lattice. If possible (that is, if :meth:`lattice` has a
          ``dual()`` method), the dual lattice is returned. Otherwise,
          `\ZZ^n` is returned, where `n` is the dimension of ``self``.

        EXAMPLES::

            sage: LatticePolytope([(1,0)]).dual_lattice()
            2-d lattice N
            sage: LatticePolytope([], lattice=ZZ^3).dual_lattice()
            Ambient free module of rank 3
            over the principal ideal domain Integer Ring
        """
        try:
            return self.lattice().dual()
        except AttributeError:
            return ZZ**self.lattice_dim()

    def edges(self):
        r"""
        Return edges (faces of dimension 1) of ``self``.

        OUTPUT: :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.edges()                                                             # needs sage.graphs
            (1-d face of 3-d reflexive polytope in 3-d lattice M,
            ...
             1-d face of 3-d reflexive polytope in 3-d lattice M)
            sage: len(o.edges())                                                        # needs sage.graphs
            12
        """
        return self.faces(dim=1)

    @cached_method
    def face_lattice(self):
        r"""
        Return the face lattice of ``self``.

        This lattice will have the empty polytope as the bottom and this
        polytope itself as the top.

        OUTPUT:

        - :class:`finite poset <sage.combinat.posets.posets.FinitePoset>` of
          :class:`lattice polytopes <LatticePolytopeClass>`.

        EXAMPLES:

        Let's take a look at the face lattice of a square::

            sage: square = LatticePolytope([(0,0), (1,0), (1,1), (0,1)])
            sage: L = square.face_lattice(); L                                          # needs sage.graphs
            Finite lattice containing 10 elements with distinguished linear extension

        To see all faces arranged by dimension, you can do this::

            sage: for level in L.level_sets(): print(level)                             # needs sage.graphs
            [-1-d face of 2-d lattice polytope in 2-d lattice M]
            [0-d face of 2-d lattice polytope in 2-d lattice M,
             0-d face of 2-d lattice polytope in 2-d lattice M,
             0-d face of 2-d lattice polytope in 2-d lattice M,
             0-d face of 2-d lattice polytope in 2-d lattice M]
            [1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M]
            [2-d lattice polytope in 2-d lattice M]

        For a particular face you can look at its actual vertices... ::

            sage: face = L.level_sets()[1][0]                                           # needs sage.graphs
            sage: face.vertices()                                                       # needs sage.graphs
            M(0, 0)
            in 2-d lattice M

        ... or you can see the index of the vertex of the original polytope that
        corresponds to the above one::

            sage: face.ambient_vertex_indices()                                         # needs sage.graphs
            (0,)
            sage: square.vertex(0)
            M(0, 0)

        An alternative to extracting faces from the face lattice is to use
        :meth:`faces` method::

            sage: face is square.faces(dim=0)[0]                                        # needs sage.graphs
            True

        The advantage of working with the face lattice directly is that you
        can (relatively easily) get faces that are related to the given one::

            sage: face = L.level_sets()[1][0]                                           # needs sage.graphs
            sage: D = L.hasse_diagram()                                                 # needs sage.graphs
            sage: sorted(D.neighbors(face))                                             # needs sage.graphs
            [-1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M]

        However, you can achieve some of this functionality using
        :meth:`facets`, :meth:`facet_of`, and :meth:`adjacent` methods::

            sage: # needs sage.graphs
            sage: face = square.faces(0)[0]
            sage: face
            0-d face of 2-d lattice polytope in 2-d lattice M
            sage: face.vertices()
            M(0, 0)
            in 2-d lattice M
            sage: face.facets()
            (-1-d face of 2-d lattice polytope in 2-d lattice M,)
            sage: face.facet_of()
            (1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M)
            sage: face.adjacent()
            (0-d face of 2-d lattice polytope in 2-d lattice M,
             0-d face of 2-d lattice polytope in 2-d lattice M)
            sage: face.adjacent()[0].vertices()
            M(1, 0)
            in 2-d lattice M

        Note that if ``p`` is a face of ``superp``, then the face
        lattice of ``p`` consists of (appropriate) faces of ``superp``::

            sage: # needs sage.graphs
            sage: superp = LatticePolytope([(1,2,3,4), (5,6,7,8),
            ....:                           (1,2,4,8), (1,3,9,7)])
            sage: superp.face_lattice()
            Finite lattice containing 16 elements with distinguished linear extension
            sage: superp.face_lattice().top()
            3-d lattice polytope in 4-d lattice M
            sage: p = superp.facets()[0]
            sage: p
            2-d face of 3-d lattice polytope in 4-d lattice M
            sage: p.face_lattice()
            Finite poset containing 8 elements with distinguished linear extension
            sage: p.face_lattice().bottom()
            -1-d face of 3-d lattice polytope in 4-d lattice M
            sage: p.face_lattice().top()
            2-d face of 3-d lattice polytope in 4-d lattice M
            sage: p.face_lattice().top() is p
            True
        """
        if self._ambient is self:
            # We need to compute face lattice on our own.
            vertex_to_facets = [row.nonzero_positions()
                                for row in self.incidence_matrix().rows()]
            facet_to_vertices = [column.nonzero_positions()
                                 for column in self.incidence_matrix().columns()]

            def LPFace(vertices, facets):
                if not facets:
                    return self
                return LatticePolytopeClass(ambient=self,
                                            ambient_vertex_indices=vertices,
                                            ambient_facet_indices=facets)

            return lattice_from_incidences(
                vertex_to_facets, facet_to_vertices, LPFace, key=id(self))
        else:
            # Get face lattice as a sublattice of the ambient one
            allowed_indices = frozenset(self._ambient_vertex_indices)
            from sage.graphs.digraph import DiGraph
            L = DiGraph()
            empty = self._ambient.face_lattice().bottom()
            L.add_vertex(0) # In case it is the only one
            dfaces = [empty]
            faces = [empty]
            face_to_index = {empty:0}
            next_index = 1
            next_d = 0 # Dimension of faces to be considered next.
            while next_d < self.dim():
                ndfaces = []
                for face in dfaces:
                    face_index = face_to_index[face]
                    for new_face in face.facet_of():
                        if not allowed_indices.issuperset(
                                        new_face._ambient_vertex_indices):
                            continue
                        if new_face in ndfaces:
                            new_face_index = face_to_index[new_face]
                        else:
                            ndfaces.append(new_face)
                            face_to_index[new_face] = next_index
                            new_face_index = next_index
                            next_index += 1
                        L.add_edge(face_index, new_face_index)
                faces.extend(ndfaces)
                dfaces = ndfaces
                next_d += 1
            if self.dim() > 0:
                # Last level is very easy to build, so we do it separately
                # even though the above cycle could do it too.
                faces.append(self)
                for face in dfaces:
                    L.add_edge(face_to_index[face], next_index)
            D = dict(enumerate(faces))
            L.relabel(D)
            return FinitePoset(L, faces, key=id(self))

    def faces(self, dim=None, codim=None):
        r"""
        Return faces of ``self`` of specified (co)dimension.

        INPUT:

        - ``dim`` -- integer; dimension of the requested faces

        - ``codim`` -- integer; codimension of the requested faces

        .. NOTE::

            You can specify at most one parameter. If you don't give any, then
            all faces will be returned.

        OUTPUT:

        - if either ``dim`` or ``codim`` is given, the output will be a
          :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`;

        - if neither ``dim`` nor ``codim`` is given, the output will be the
          :class:`tuple` of tuples as above, giving faces of all existing
          dimensions. If you care about inclusion relations between faces,
          consider using :meth:`face_lattice` or :meth:`adjacent`,
          :meth:`facet_of`, and :meth:`facets`.

        EXAMPLES:

        Let's take a look at the faces of a square::

            sage: square = LatticePolytope([(0,0), (1,0), (1,1), (0,1)])
            sage: square.faces()                                                        # needs sage.graphs
            ((-1-d face of 2-d lattice polytope in 2-d lattice M,),
             (0-d face of 2-d lattice polytope in 2-d lattice M,
              0-d face of 2-d lattice polytope in 2-d lattice M,
              0-d face of 2-d lattice polytope in 2-d lattice M,
              0-d face of 2-d lattice polytope in 2-d lattice M),
             (1-d face of 2-d lattice polytope in 2-d lattice M,
              1-d face of 2-d lattice polytope in 2-d lattice M,
              1-d face of 2-d lattice polytope in 2-d lattice M,
              1-d face of 2-d lattice polytope in 2-d lattice M),
             (2-d lattice polytope in 2-d lattice M,))

        Its faces of dimension one (i.e., edges)::

            sage: square.faces(dim=1)                                                   # needs sage.graphs
            (1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M,
             1-d face of 2-d lattice polytope in 2-d lattice M)

        Its faces of codimension one are the same (also edges)::

            sage: square.faces(codim=1) is square.faces(dim=1)                          # needs sage.graphs
            True

        Let's pick a particular face::

            sage: face = square.faces(dim=1)[0]                                         # needs sage.graphs

        Now you can look at the actual vertices of this face... ::

            sage: face.vertices()                                                       # needs sage.graphs
            M(0, 0),
            M(0, 1)
            in 2-d lattice M

        ... or you can see indices of the vertices of the original polytope that
        correspond to the above ones::

            sage: face.ambient_vertex_indices()                                         # needs sage.graphs
            (0, 3)
            sage: square.vertices(face.ambient_vertex_indices())                        # needs sage.graphs
            M(0, 0),
            M(0, 1)
            in 2-d lattice M
        """
        if dim is not None and codim is not None:
            raise ValueError(
                    "dimension and codimension cannot be specified together!")
        dim = self.dim() - codim if codim is not None else dim
        if "_faces" not in self.__dict__:
            self._faces = tuple(map(self._sort_faces,
                                    self.face_lattice().level_sets()))
        if dim is None:
            return self._faces
        else:
            return self._faces[dim + 1] if -1 <= dim <= self.dim() else ()

    def facet_constant(self, i):
        r"""
        Return the constant in the `i`-th facet inequality of this polytope.

        This is equivalent to ``facet_constants()[i]``.

        INPUT:

        - ``i`` -- integer; the index of the facet

        OUTPUT: integer; the constant in the `i`-th facet inequality

        .. SEEALSO::

            :meth:`facet_constants`,
            :meth:`facet_normal`,
            :meth:`facet_normals`,
            :meth:`facets`.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.facet_constant(0)
            1
            sage: o.facet_constant(0) == o.facet_constants()[0]
            True
        """
        return self.facet_constants()[i]

    def facet_constants(self):
        r"""
        Return facet constants of ``self``.

        Facet inequalities have form `n \cdot x + c \geq 0` where `n` is the
        inner normal and `c` is a constant.

        OUTPUT: integer vector

        .. SEEALSO::

            :meth:`facet_constant`,
            :meth:`facet_normal`,
            :meth:`facet_normals`,
            :meth:`facets`.

        EXAMPLES:

        For reflexive polytopes all constants are 1::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: o.facet_constants()
            (1, 1, 1, 1, 1, 1, 1, 1)

        Here is an example of a 3-dimensional polytope in a 4-dimensional
        space with 3 facets containing the origin::

            sage: p = LatticePolytope([(0,0,0,0), (1,1,1,3),
            ....:                      (1,-1,1,3), (-1,-1,1,3)])
            sage: p.vertices()
            M( 0,  0, 0, 0),
            M( 1,  1, 1, 3),
            M( 1, -1, 1, 3),
            M(-1, -1, 1, 3)
            in 4-d lattice M
            sage: p.facet_constants()
            (0, 0, 3, 0)
        """
        try:
            return self._facet_constants
        except AttributeError:
            self._compute_facets()
            return self._facet_constants

    def facet_normal(self, i):
        r"""
        Return the inner normal to the ``i``-th facet of this polytope.

        This is equivalent to ``facet_normals()[i]``.

        INPUT:

        - ``i`` -- integer; the index of the facet

        OUTPUT: a vector

        .. SEEALSO::

            :meth:`facet_constant`,
            :meth:`facet_constants`,
            :meth:`facet_normals`,
            :meth:`facets`.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.facet_normal(0)
            N(1, -1, -1)
            sage: o.facet_normal(0) is o.facet_normals()[0]
            True
        """
        return self.facet_normals()[i]

    def facet_normals(self):
        r"""
        Return inner normals to the facets of ``self``.

        If this polytope is not full-dimensional, facet normals will define
        this polytope in the affine subspace spanned by it.

        OUTPUT:

        - a :class:`point collection <PointCollection>` in the
          :meth:`dual_lattice` of ``self``.

        .. SEEALSO::

            :meth:`facet_constant`,
            :meth:`facet_constants`,
            :meth:`facet_normal`,
            :meth:`facets`.

        EXAMPLES:

        Normals to facets of an octahedron are vertices of a cube::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: o.facet_normals()
            N( 1, -1, -1),
            N( 1,  1, -1),
            N( 1,  1,  1),
            N( 1, -1,  1),
            N(-1, -1,  1),
            N(-1, -1, -1),
            N(-1,  1, -1),
            N(-1,  1,  1)
            in 3-d lattice N

        Here is an example of a 3-dimensional polytope in a 4-dimensional
        space::

            sage: p = LatticePolytope([(0,0,0,0), (1,1,1,3),
            ....:                      (1,-1,1,3), (-1,-1,1,3)])
            sage: p.vertices()
            M( 0,  0, 0, 0),
            M( 1,  1, 1, 3),
            M( 1, -1, 1, 3),
            M(-1, -1, 1, 3)
            in 4-d lattice M
            sage: p.facet_normals()
            N( 0,  3, 0,  1),
            N( 1, -1, 0,  0),
            N( 0,  0, 0, -1),
            N(-3,  0, 0,  1)
            in 4-d lattice N
            sage: p.facet_constants()
            (0, 0, 3, 0)

        Now we manually compute the distance matrix of this polytope. Since it
        is a simplex, each line (corresponding to a facet) should consist of
        zeros (indicating generating vertices of the corresponding facet) and
        a single positive number (since our normals are inner)::

            sage: matrix([[n * v + c for v in p.vertices()]
            ....:     for n, c in zip(p.facet_normals(), p.facet_constants())])
            [0 6 0 0]
            [0 0 2 0]
            [3 0 0 0]
            [0 0 0 6]
        """
        try:
            return self._facet_normals
        except AttributeError:
            self._compute_facets()
            return self._facet_normals

    @cached_method
    def facet_of(self):
        r"""
        Return elements of the ambient face lattice having ``self`` as a facet.

        OUTPUT: :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`

        EXAMPLES::

            sage: # needs sage.graphs
            sage: square = LatticePolytope([(0,0), (1,0), (1,1), (0,1)])
            sage: square.facet_of()
            ()
            sage: face = square.faces(0)[0]
            sage: len(face.facet_of())
            2
            sage: face.facet_of()[1]
            1-d face of 2-d lattice polytope in 2-d lattice M
        """
        L = self._ambient.face_lattice()
        H = L.hasse_diagram()
        return self._sort_faces(f for f in H.neighbors_out(L(self)))

    def facets(self):
        r"""
        Return facets (faces of codimension 1) of ``self``.

        OUTPUT: :class:`tuple` of :class:`lattice polytopes <LatticePolytopeClass>`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.facets()                                                            # needs sage.graphs
            (2-d face of 3-d reflexive polytope in 3-d lattice M,
            ...
             2-d face of 3-d reflexive polytope in 3-d lattice M)
            sage: len(o.facets())                                                       # needs sage.graphs
            8
        """
        return self.faces(codim=1)

    # Dictionaries of normal forms
    _rp_dict = [None] * 4

    @cached_method
    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

           The columns correspond to facets/facet normals
           in the order of :meth:`facet_normals`, the rows
           correspond to the vertices in the order of
           :meth:`vertices`.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.incidence_matrix()
            [0 0 1 1]
            [0 1 1 0]
            [1 1 0 0]
            [1 0 0 1]
            sage: o.faces(1)[0].incidence_matrix()                                      # needs sage.graphs
            [1 0]
            [0 1]

            sage: o = lattice_polytope.cross_polytope(4)
            sage: o.incidence_matrix().column(3).nonzero_positions()
            [3, 4, 5, 6]
            sage: o.facets()[3].ambient_vertex_indices()                                # needs sage.graphs
            (3, 4, 5, 6)

        TESTS::

            sage: o.incidence_matrix().is_immutable()
            True

        Check that the base ring is ``ZZ``, see :issue:`29840`::

            sage: o.incidence_matrix().base_ring()
            Integer Ring
        """
        incidence_matrix = matrix(ZZ, self.nvertices(),
                                  self.nfacets(), 0)

        for Hindex, normal in enumerate(self.facet_normals()):
            facet_constant = self.facet_constant(Hindex)
            for Vindex, vertex in enumerate(self.vertices()):
                if normal*vertex + facet_constant == 0:
                    incidence_matrix[Vindex, Hindex] = 1

        incidence_matrix.set_immutable()
        return incidence_matrix

    @cached_method
    def index(self):
        r"""
        Return the index of this polytope in the internal database of 2- or
        3-dimensional reflexive polytopes. Databases are stored in the
        directory of the package.

        .. NOTE::

           The first call to this function for each dimension can take
           a few seconds while the dictionary of all polytopes is
           constructed, but after that it is cached and fast.

        :rtype: integer

        EXAMPLES: We check what is the index of the "diamond" in the
        database::

            sage: d = lattice_polytope.cross_polytope(2)
            sage: d.index()                                                             # needs palp
            3

        Note that polytopes with the same index are not necessarily the
        same::

            sage: d.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1,  0),
            M( 0, -1)
            in 2-d lattice M
            sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M

        But they are in the same `GL(\ZZ^n)` orbit and have the same
        normal form::

            sage: d.normal_form()                                                       # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M
            sage: lattice_polytope.ReflexivePolytope(2,3).normal_form()                 # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M
        """
        if not self.is_reflexive():
            raise NotImplementedError("only reflexive polytopes can be indexed!")
        dim = self.dim()
        if dim not in [2, 3]:
            raise NotImplementedError("only 2- and 3-dimensional polytopes can be indexed!")
        if LatticePolytopeClass._rp_dict[dim] is None:
            rp_dict = dict()
            for n, p in enumerate(ReflexivePolytopes(dim)):
                rp_dict[p.normal_form().matrix()] = n
            LatticePolytopeClass._rp_dict[dim] = rp_dict
        return LatticePolytopeClass._rp_dict[dim][self.normal_form().matrix()]

    @cached_method
    def interior_point_indices(self):
        r"""
        Return indices of (relative) interior lattice points of this polytope.

        OUTPUT: increasing :class:`tuple` of integers

        EXAMPLES:

        The origin is the only interior point of this square::

            sage: square = lattice_polytope.cross_polytope(2).polar()
            sage: square.points()                                                       # needs palp
            N( 1,  1),
            N( 1, -1),
            N(-1, -1),
            N(-1,  1),
            N(-1,  0),
            N( 0, -1),
            N( 0,  0),
            N( 0,  1),
            N( 1,  0)
            in 2-d lattice N
            sage: square.interior_point_indices()                                       # needs palp
            (6,)

        Its edges also have a single interior point each::

            sage: face = square.edges()[0]                                              # needs sage.graphs
            sage: face.points()                                                         # needs sage.graphs
            N(-1, -1),
            N(-1,  1),
            N(-1,  0)
            in 2-d lattice N
            sage: face.interior_point_indices()                                         # needs sage.graphs
            (2,)
        """
        return tuple(i
                     for i, c in enumerate(self.distances().columns(copy=False))
                     if len(c.nonzero_positions()) == self.nfacets())

    def interior_points(self):
        r"""
        Return (relative) boundary lattice points of this polytope.

        OUTPUT: a :class:`point collection <PointCollection>`

        EXAMPLES:

        The origin is the only interior point of this square::

            sage: square = lattice_polytope.cross_polytope(2).polar()
            sage: square.interior_points()                                              # needs palp
            N(0, 0)
            in 2-d lattice N

        Its edges also have a single interior point each::

            sage: face = square.edges()[0]                                              # needs sage.graphs
            sage: face.interior_points()                                                # needs sage.graphs
            N(-1, 0)
            in 2-d lattice N
        """
        return self.points(self.interior_point_indices())

    @cached_method
    def is_reflexive(self):
        r"""
        Return ``True`` if this polytope is reflexive.

        EXAMPLES: The 3-dimensional octahedron is reflexive (and 4319 other
        3-polytopes)::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.is_reflexive()
            True

        But not all polytopes are reflexive::

            sage: p = LatticePolytope([(1,0,0), (0,1,17), (-1,0,0), (0,-1,0)])
            sage: p.is_reflexive()
            False

        Only full-dimensional polytopes can be reflexive (otherwise the polar
        set is not a polytope at all, since it is unbounded)::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.is_reflexive()
            False
        """
        return self.dim() == self.lattice_dim() and \
                all(c == 1 for c in self.facet_constants())

    def lattice(self):
        r"""
        Return the ambient lattice of ``self``.

        OUTPUT: a lattice

        EXAMPLES::

            sage: lattice_polytope.cross_polytope(3).lattice()
            3-d lattice M
        """
        return self._vertices.module()

    def lattice_dim(self):
        r"""
        Return the dimension of the ambient lattice of ``self``.

        An alias is :meth:`ambient_dim`.

        OUTPUT: integer

        EXAMPLES::

            sage: p = LatticePolytope([(1,0)])
            sage: p.lattice_dim()
            2
            sage: p.dim()
            0
        """
        return self.lattice().dimension()

    ambient_dim = lattice_dim

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        It is the ambient lattice (:meth:`lattice`) tensored with a field.

        INPUT:

        - ``base_field`` -- (default: the rationals) a field

        EXAMPLES::

            sage: p = LatticePolytope([(1,0)])
            sage: p.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
            sage: p.ambient_vector_space(AA)                                            # needs sage.rings.number_field
            Vector space of dimension 2 over Algebraic Real Field
        """
        return self.lattice().vector_space(base_field=base_field)

    def linearly_independent_vertices(self):
        r"""
        Return a maximal set of linearly independent vertices.

        OUTPUT: a tuple of vertex indices

        EXAMPLES::

            sage: L = LatticePolytope([[0, 0], [-1, 1], [-1, -1]])
            sage: L.linearly_independent_vertices()
            (1, 2)
            sage: L = LatticePolytope([[0, 0, 0]])
            sage: L.linearly_independent_vertices()
            ()
            sage: L = LatticePolytope([[0, 1, 0]])
            sage: L.linearly_independent_vertices()
            (0,)
        """
        return self.vertices().matrix().pivot_rows()

    def nef_partitions(self, keep_symmetric=False, keep_products=True,
        keep_projections=True, hodge_numbers=False):
        r"""
        Return 2-part nef-partitions of ``self``.

        INPUT:

        - ``keep_symmetric`` -- boolean (default: ``False``); if ``True``, "-s" option
          will be passed to ``nef.x`` in order to keep symmetric partitions,
          i.e. partitions related by lattice automorphisms preserving ``self``

        - ``keep_products`` -- boolean (default: ``True``); if ``True``, "-D" option
          will be passed to ``nef.x`` in order to keep product partitions,
          with corresponding complete intersections being direct products

        - ``keep_projections`` -- boolean (default: ``True``); if ``True``, "-P" option
          will be passed to ``nef.x`` in order to keep projection partitions,
          i.e. partitions with one of the parts consisting of a single vertex

        - ``hodge_numbers`` -- boolean (default: ``False``); if ``False``, "-p" option
          will be passed to ``nef.x`` in order to skip Hodge numbers
          computation, which takes a lot of time

        OUTPUT: a sequence of :class:`nef-partitions <NefPartition>`

        Type ``NefPartition?`` for definitions and notation.

        EXAMPLES:

        Nef-partitions of the 4-dimensional cross-polytope::

            sage: p = lattice_polytope.cross_polytope(4)
            sage: p.nef_partitions()                                                    # needs palp
            [Nef-partition {0, 1, 4, 5} ⊔ {2, 3, 6, 7} (direct product),
             Nef-partition {0, 1, 2, 4} ⊔ {3, 5, 6, 7},
             Nef-partition {0, 1, 2, 4, 5} ⊔ {3, 6, 7},
             Nef-partition {0, 1, 2, 4, 5, 6} ⊔ {3, 7} (direct product),
             Nef-partition {0, 1, 2, 3} ⊔ {4, 5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4} ⊔ {5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4, 5} ⊔ {6, 7},
             Nef-partition {0, 1, 2, 3, 4, 5, 6} ⊔ {7} (projection)]

        Now we omit projections::

            sage: p.nef_partitions(keep_projections=False)                              # needs palp
            [Nef-partition {0, 1, 4, 5} ⊔ {2, 3, 6, 7} (direct product),
             Nef-partition {0, 1, 2, 4} ⊔ {3, 5, 6, 7},
             Nef-partition {0, 1, 2, 4, 5} ⊔ {3, 6, 7},
             Nef-partition {0, 1, 2, 4, 5, 6} ⊔ {3, 7} (direct product),
             Nef-partition {0, 1, 2, 3} ⊔ {4, 5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4} ⊔ {5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4, 5} ⊔ {6, 7}]

        Currently Hodge numbers cannot be computed for a given nef-partition::

            sage: p.nef_partitions()[1].hodge_numbers()                                 # needs palp
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!

        But they can be obtained from ``nef.x`` for all nef-partitions at once.
        Partitions will be exactly the same::

            sage: p.nef_partitions(hodge_numbers=True)  # long time (2s on sage.math, 2011), needs palp
            [Nef-partition {0, 1, 4, 5} ⊔ {2, 3, 6, 7} (direct product),
             Nef-partition {0, 1, 2, 4} ⊔ {3, 5, 6, 7},
             Nef-partition {0, 1, 2, 4, 5} ⊔ {3, 6, 7},
             Nef-partition {0, 1, 2, 4, 5, 6} ⊔ {3, 7} (direct product),
             Nef-partition {0, 1, 2, 3} ⊔ {4, 5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4} ⊔ {5, 6, 7},
             Nef-partition {0, 1, 2, 3, 4, 5} ⊔ {6, 7},
             Nef-partition {0, 1, 2, 3, 4, 5, 6} ⊔ {7} (projection)]

        Now it is possible to get Hodge numbers::

            sage: p.nef_partitions(hodge_numbers=True)[1].hodge_numbers()               # needs palp
            (20,)

        Since nef-partitions are cached, their Hodge numbers are accessible
        after the first request, even if you do not specify
        ``hodge_numbers=True`` anymore::

            sage: p.nef_partitions()[1].hodge_numbers()                                 # needs palp
            (20,)

        We illustrate removal of symmetric partitions on a diamond::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: p.nef_partitions()                                                    # needs palp
            [Nef-partition {0, 2} ⊔ {1, 3} (direct product),
             Nef-partition {0, 1} ⊔ {2, 3},
             Nef-partition {0, 1, 2} ⊔ {3} (projection)]
            sage: p.nef_partitions(keep_symmetric=True)                                 # needs palp
            [Nef-partition {0, 1, 3} ⊔ {2} (projection),
             Nef-partition {0, 2, 3} ⊔ {1} (projection),
             Nef-partition {0, 3} ⊔ {1, 2},
             Nef-partition {1, 2, 3} ⊔ {0} (projection),
             Nef-partition {1, 3} ⊔ {0, 2} (direct product),
             Nef-partition {2, 3} ⊔ {0, 1},
             Nef-partition {0, 1, 2} ⊔ {3} (projection)]

        Nef-partitions can be computed only for reflexive polytopes::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,2),
            ....:                      (-1,0,0), (0,-1,0), (0,0,-1)])
            sage: p.nef_partitions()                                                    # needs palp
            Traceback (most recent call last):
            ...
            ValueError: The given polytope is not reflexive!
            Polytope: 3-d lattice polytope in 3-d lattice M
        """
        if not self.is_reflexive():
            raise ValueError(("The given polytope is not reflexive!\n"
                                + "Polytope: %s") % self)
        keys = "-N -V"
        if keep_symmetric:
            keys += " -s"
        if keep_products:
            keys += " -D"
        if keep_projections:
            keys += " -P"
        if not hodge_numbers:
            keys += " -p"
        if hasattr(self, "_npkeys"):
            oldkeys = self._npkeys
            if oldkeys == keys:
                return self._nef_partitions
            if not (hodge_numbers and oldkeys.find("-p") != -1
                or keep_symmetric and oldkeys.find("-s") == -1
                or not keep_symmetric and oldkeys.find("-s") != -1
                or keep_projections and oldkeys.find("-P") == -1
                or keep_products and oldkeys.find("-D") == -1):
                # Select only necessary partitions
                return Sequence([p for p in self._nef_partitions
                                 if (keep_projections or not p._is_projection)
                                     and (keep_products or not p._is_product)],
                                cr=True, check=False)
        self._read_nef_partitions(self.nef_x(keys))
        self._npkeys = keys
        return self._nef_partitions

    def nef_x(self, keys):
        r"""
        Run ``nef.x`` with given ``keys`` on vertices of this
        polytope.

        INPUT:

        - ``keys`` -- string of options passed to ``nef.x``; the
          key "-f" is added automatically

        OUTPUT: the output of ``nef.x`` as a string

        EXAMPLES: This call is used internally for computing
        nef-partitions::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: s = o.nef_x("-N -V -p")                                               # needs palp
            sage: s                      # output contains random time                  # needs palp
            M:27 8 N:7 6  codim=2 #part=5
            3 6  Vertices of P:
                1    0    0   -1    0    0
                0    1    0    0   -1    0
                0    0    1    0    0   -1
             P:0 V:2 4 5       0sec  0cpu
             P:2 V:3 4 5       0sec  0cpu
             P:3 V:4 5       0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu
        """
        return self._palp("nef.x -f " + keys)

    def nfacets(self):
        r"""
        Return the number of facets of this polytope.

        EXAMPLES: The number of facets of the 3-dimensional octahedron::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.nfacets()
            8

        The number of facets of an interval is 2::

            sage: LatticePolytope(([1],[2])).nfacets()
            2

        Now consider a 2-dimensional diamond in a 3-dimensional space::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.nfacets()
            4
        """
        return len(self.facet_normals()) if self.dim() > 0 else 0

    @cached_method
    def normal_form(self, algorithm='palp_native', permutation=False):
        r"""
        Return the normal form of vertices of ``self``.

        Two full-dimensional lattice polytopes are in the same
        `GL(\ZZ^n)`-orbit if and only if their normal forms are the
        same. Normal form is not defined and thus cannot be used for polytopes
        whose dimension is smaller than the dimension of the ambient space.

        The original algorithm was presented in [KS1998]_ and implemented
        in PALP. A modified version of the PALP algorithm is discussed in
        [GK2013]_ and available here as ``'palp_modified'``.

        INPUT:

        - ``algorithm`` -- (default: ``'palp_native'``) the algorithm which is used
          to compute the normal form. Options are:

          * ``'palp'`` -- run external PALP code, usually the fastest option
            when it works; but reproducible crashes have been observed in dimension
            5 and higher.

          * ``'palp_native'`` -- the original PALP algorithm implemented
            in sage. Currently competitive with PALP in many cases.

          * ``'palp_modified'`` -- a modified version of the PALP
            algorithm which determines the maximal vertex-facet
            pairing matrix first and then computes its
            automorphisms, while the PALP algorithm does both things
            concurrently.

        - ``permutation`` -- boolean (default: ``False``); if ``True``, the permutation
          applied to vertices to obtain the normal form is returned as well.
          Note that the different algorithms may return different results
          that nevertheless lead to the same normal form.

        OUTPUT:

        - a :class:`point collection <PointCollection>` in the :meth:`lattice`
          of ``self`` or a tuple of it and a permutation.

        EXAMPLES:

        We compute the normal form of the "diamond"::

            sage: d = LatticePolytope([(1,0), (0,1), (-1,0), (0,-1)])
            sage: d.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1,  0),
            M( 0, -1)
            in 2-d lattice M
            sage: d.normal_form()                                                       # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M

        The diamond is the 3rd polytope in the internal database::

            sage: d.index()                                                             # needs palp
            3
            sage: d                                                                     # needs palp
            2-d reflexive polytope #3 in 2-d lattice M

        You can get it in its normal form (in the default lattice) as ::

            sage: lattice_polytope.ReflexivePolytope(2, 3).vertices()
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M

        It is not possible to compute normal forms for polytopes which do not
        span the space::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.normal_form()
            Traceback (most recent call last):
            ...
            ValueError: normal form is not defined for
            2-d lattice polytope in 3-d lattice M

        We can perform the same examples using other algorithms::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.normal_form(algorithm='palp_native')                                # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.normal_form(algorithm='palp_modified')                              # needs sage.groups
            Traceback (most recent call last):
            ...
            in 2-d lattice M

        The following examples demonstrate the speed of the available algorithms.
        In low dimensions, the default algorithm, ``'palp_native'``, is the fastest.
        As the dimension increases, ``'palp'`` is relatively faster than ``'palp_native'``.
        ``'palp_native'`` is usually much faster than ``'palp_modified'``.
        In some cases when the polytope has high symmetry, however, ``'palp_native'`` is slower::

            sage: # not tested
            sage: o = lattice_polytope.cross_polytope(2)
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp")
            625 loops, best of 3: 3.07 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_native")
            625 loops, best of 3: 0.445 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_modified")
            625 loops, best of 3: 5.01 ms per loop
            sage: o = lattice_polytope.cross_polytope(3)
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp")
            625 loops, best of 3: 3.22 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_native")
            625 loops, best of 3: 2.73 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_modified")
            625 loops, best of 3: 20.7 ms per loop
            sage: o = lattice_polytope.cross_polytope(4)
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp")
            625 loops, best of 3: 4.84 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_native")
            625 loops, best of 3: 55.6 ms per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_modified")
            625 loops, best of 3: 129 ms per loop
            sage: o = lattice_polytope.cross_polytope(5)
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp")
            10 loops, best of 3: 0.0364 s per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_native")
            10 loops, best of 3: 1.68 s per loop
            sage: %timeit o.normal_form.clear_cache(); o.normal_form("palp_modified")
            10 loops, best of 3: 0.858 s per loop

        Note that the algorithm ``'palp'`` may crash for higher dimensions because of
        the overflow errors as mentioned in :issue:`13525#comment:9`.
        Then use ``'palp_native'`` instead, which is usually faster than ``'palp_modified'``.
        Below is an example where ``'palp'`` fails and
        ``'palp_native'`` is much faster than ``'palp_modified'``::

            sage: P = LatticePolytope([[-3, -3, -6, -6, -1], [3, 3, 6, 6, 1], [-3, -3, -6, -6, 1],
            ....:                      [-3, -3, -3, -6, 0], [-3, -3, -3, 0, 0], [-3, -3, 0, 0, 0],
            ....:                      [-3, 0, -6, -6, 0], [-3, 0, -3, -6, 0], [-3, 0, -3, 0, 0],
            ....:                      [-3, 0, 0, 0, -1], [3, 3, 6, 6, -1], [-3, 0, 0, 0, 1],
            ....:                      [0, -3, -6, -6, 0], [0, -3, -3, -6, 0], [0, -3, -3, 0, 0],
            ....:                      [0, -3, 0, 0, -1], [3, 3, 3, 6, 0], [0, -3, 0, 0, 1],
            ....:                      [0, 0, -6, -6, 0], [0, 0, -3, -6, -1], [3, 3, 3, 0, 0],
            ....:                      [0, 0, -3, -6, 1], [0, 0, -3, 0, -1], [3, 3, 0, 0, 0],
            ....:                      [0, 0, -3, 0, 1], [0, 0, 3, 0, -1], [3, 0, 6, 6, 0],
            ....:                      [0, 0, 3, 0, 1], [0, 0, 3, 6, -1], [3, 0, 3, 6, 0],
            ....:                      [0, 0, 3, 6, 1], [0, 0, 6, 6, 0], [0, 3, 0, 0, -1],
            ....:                      [3, 0, 3, 0, 0], [0, 3, 0, 0, 1], [0, 3, 3, 0, 0],
            ....:                      [0, 3, 3, 6, 0], [0, 3, 6, 6, 0], [3, 0,0, 0, -1], [3, 0, 0, 0, 1]])
            sage: P.normal_form(algorithm='palp')  # not tested
            Traceback (most recent call last):
            ...
            RuntimeError: Error executing ... for a polytope sequence!
            Output:
            b'*** stack smashing detected ***: terminated\nAborted\n'
            sage: P.normal_form(algorithm='palp_native')                                # needs sage.groups
            M(  6,  0,  0,  0,  0),
            M( -6,  0,  0,  0,  0),
            M(  0,  1,  0,  0,  0),
            M(  0,  0,  3,  0,  0),
            M(  0,  1,  0,  3,  0),
            M(  0,  0,  0,  0,  3),
            M( -6,  1,  6,  3, -6),
            M( -6,  0,  6,  0, -3),
            M(-12,  1,  6,  3, -3),
            M( -6,  1,  0,  3,  0),
            M( -6,  0,  3,  3,  0),
            M(  6,  0, -6, -3,  6),
            M(-12,  1,  6,  3, -6),
            M(-12,  0,  9,  3, -6),
            M(  0,  0,  0, -3,  0),
            M(-12,  1,  6,  6, -6),
            M(-12,  0,  6,  3, -3),
            M(  0,  1, -3,  0,  0),
            M(  0,  0, -3, -3,  3),
            M(  0,  1,  0,  3, -3),
            M(  0, -1,  0, -3,  3),
            M(  0,  0,  3,  3, -3),
            M(  0, -1,  3,  0,  0),
            M( 12,  0, -6, -3,  3),
            M( 12, -1, -6, -6,  6),
            M(  0,  0,  0,  3,  0),
            M( 12,  0, -9, -3,  6),
            M( 12, -1, -6, -3,  6),
            M( -6,  0,  6,  3, -6),
            M(  6,  0, -3, -3,  0),
            M(  6, -1,  0, -3,  0),
            M(-12,  1,  9,  6, -6),
            M(  6,  0, -6,  0,  3),
            M(  6, -1, -6, -3,  6),
            M(  0,  0,  0,  0, -3),
            M(  0, -1,  0, -3,  0),
            M(  0,  0, -3,  0,  0),
            M(  0, -1,  0,  0,  0),
            M( 12, -1, -9, -6,  6),
            M( 12, -1, -6, -3,  3)
            in 5-d lattice M
            sage: P.normal_form(algorithm='palp_modified')      # not tested (22s; MemoryError on 32 bit), needs sage.groups
            M(  6,  0,  0,  0,  0),
            M( -6,  0,  0,  0,  0),
            M(  0,  1,  0,  0,  0),
            M(  0,  0,  3,  0,  0),
            M(  0,  1,  0,  3,  0),
            M(  0,  0,  0,  0,  3),
            M( -6,  1,  6,  3, -6),
            M( -6,  0,  6,  0, -3),
            M(-12,  1,  6,  3, -3),
            M( -6,  1,  0,  3,  0),
            M( -6,  0,  3,  3,  0),
            M(  6,  0, -6, -3,  6),
            M(-12,  1,  6,  3, -6),
            M(-12,  0,  9,  3, -6),
            M(  0,  0,  0, -3,  0),
            M(-12,  1,  6,  6, -6),
            M(-12,  0,  6,  3, -3),
            M(  0,  1, -3,  0,  0),
            M(  0,  0, -3, -3,  3),
            M(  0,  1,  0,  3, -3),
            M(  0, -1,  0, -3,  3),
            M(  0,  0,  3,  3, -3),
            M(  0, -1,  3,  0,  0),
            M( 12,  0, -6, -3,  3),
            M( 12, -1, -6, -6,  6),
            M(  0,  0,  0,  3,  0),
            M( 12,  0, -9, -3,  6),
            M( 12, -1, -6, -3,  6),
            M( -6,  0,  6,  3, -6),
            M(  6,  0, -3, -3,  0),
            M(  6, -1,  0, -3,  0),
            M(-12,  1,  9,  6, -6),
            M(  6,  0, -6,  0,  3),
            M(  6, -1, -6, -3,  6),
            M(  0,  0,  0,  0, -3),
            M(  0, -1,  0, -3,  0),
            M(  0,  0, -3,  0,  0),
            M(  0, -1,  0,  0,  0),
            M( 12, -1, -9, -6,  6),
            M( 12, -1, -6, -3,  3)
            in 5-d lattice M
            sage: %timeit P.normal_form.clear_cache(); P.normal_form("palp_native")    # not tested
            10 loops, best of 3: 0.137 s per loop
            sage: %timeit P.normal_form.clear_cache(); P.normal_form("palp_modified")  # not tested
            10 loops, best of 3:  22.2 s per loop

        TESTS::

            sage: d.normal_form("palp_fiction")
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'palp', 'palp_native', or 'palp_modified'
        """
        if self.dim() < self.lattice_dim():
            raise ValueError("normal form is not defined for %s" % self)
        M = self.lattice()
        if algorithm == "palp":
            result = read_palp_point_collection(
                StringIO(self.poly_x("N")), M, permutation=permutation)
        elif algorithm == "palp_native":
            result = self._palp_native_normal_form(permutation=permutation)
        elif algorithm == "palp_modified":
            result = self._palp_modified_normal_form(permutation=permutation)
        else:
            raise ValueError("algorithm must be 'palp', 'palp_native', or 'palp_modified'")
        if permutation:
            vertices, perm = result
        else:
            vertices = result
        if algorithm != "palp":
            vertices = [M(_) for _ in vertices]
            for v in vertices:
                v.set_immutable()
            vertices = PointCollection(vertices, M)
        return (vertices, perm) if permutation else vertices

    def _palp_modified_normal_form(self, permutation=False):
        r"""
        Return the normal form of ``self`` using the modified PALP algorithm.

        This is a helper function for :meth:`normal_form` and should not
        be called directly. The modified PALP algorithm can be faster than the
        native algorithm in case the automorphism group of the
        vertex-facet pairing matrix is large.

        INPUT:

        - ``permutation`` -- boolean (default: ``False``); whether to return
          the permutation of the order of the vertices that was applied to
          obtain this matrix

        OUTPUT: a matrix or a tuple of a matrix and a permutation

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1,  0),
            M( 0, -1)
            in 2-d lattice M
            sage: o._palp_modified_normal_form()                                        # needs sage.graphs sage.groups
            Traceback (most recent call last):
            ...
            in 2-d lattice M
            sage: o._palp_modified_normal_form(permutation=True)                        # needs sage.graphs sage.groups
            Traceback (most recent call last):
            ...
            (3,4))
        """
        PM = self.vertex_facet_pairing_matrix()
        PM_max = PM.permutation_normal_form()
        perm = PM.is_permutation_of(PM_max, check=True)[1]
        permutations = PM.automorphisms_of_rows_and_columns()
        permutations = {k:[(perm[0])*p[0], (perm[1])*p[1]]
                        for k, p in enumerate(permutations)}
        out = _palp_canonical_order(self.vertices(), PM_max, permutations)
        if permutation:
            return out
        else:
            return out[0]

    def _palp_native_normal_form(self, permutation=False):
        r"""
        Return the normal form of ``self`` using the native PALP algorithm
        implemented in Sage.

        This is a helper function for :meth:`normal_form` and should not
        be called directly.

        INPUT:

        -   ``permutation`` -- a Boolean, whether to return the permutation
            of the order of the vertices that was applied to obtain this
            matrix.

        OUTPUT: a matrix or a tuple of a matrix and a permutation

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.vertices()
            M( 1,  0),
            M( 0,  1),
            M(-1,  0),
            M( 0, -1)
            in 2-d lattice M
            sage: o._palp_native_normal_form()                                          # needs sage.groups
            M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M
            sage: o._palp_native_normal_form(permutation=True)                          # needs sage.groups
            (M( 1,  0),
            M( 0,  1),
            M( 0, -1),
            M(-1,  0)
            in 2-d lattice M, (1,3,2,4))
        """
        PM_max, permutations = self._palp_PM_max(check=True)
        out = _palp_canonical_order(self.vertices(), PM_max, permutations)
        if permutation:
            return out
        else:
            return out[0]

    def _palp_PM_max(self, check=False):
        r"""
        Compute the permutation normal form of the vertex facet pairing
        matrix .

        The permutation normal form of a matrix is defined as the lexicographic
        maximum under all permutations of its rows and columns. For more
        more detail, see also
        :meth:`~sage.matrix.matrix2.Matrix.permutation_normal_form`.

        Instead of using the generic method for computing the permutation
        normal form, this method uses the PALP algorithm to compute
        the permutation normal form and its automorphisms concurrently.

        INPUT:

        - ``check`` -- boolean (default: ``False``); whether to return
          the permutations leaving the maximal vertex-facet pairing
          matrix invariant

        OUTPUT:

        A matrix or a tuple of a matrix and a dict whose values are the
        permutation group elements corresponding to the permutations
        that permute :meth:`vertices` such that the vertex-facet pairing
        matrix is maximal.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: PM = o.vertex_facet_pairing_matrix()
            sage: PM_max = PM.permutation_normal_form()                                 # needs sage.graphs
            Traceback (most recent call last):
            ...
            <BLANKLINE>
            sage: PM_max == o._palp_PM_max()                                            # needs sage.graphs sage.groups
            Traceback (most recent call last):
            ...
            True
            sage: P2 = ReflexivePolytope(2, 0)
            sage: PM_max, permutations = P2._palp_PM_max(check=True)                    # needs sage.groups
            sage: PM_max                                                                # needs sage.graphs sage.groups
            [3 0 0]
            [0 3 0]
            [0 0 3]
            sage: list(permutations.values())                                           # needs sage.groups
            [[(1,2,3), (1,2,3)],
             [(1,3,2), (1,3,2)],
             [(1,3), (1,3)],
             [(1,2), (1,2)],
             [(), ()],
             [(2,3), (2,3)]]
            sage: PM_max.automorphisms_of_rows_and_columns()                            # needs sage.graphs sage.groups
            [((), ()),
             ((1,2,3), (1,2,3)),
             ((1,3,2), (1,3,2)),
             ((2,3), (2,3)),
             ((1,2), (1,2)),
             ((1,3), (1,3))]
            sage: PMs = ( i._palp_PM_max(check=True)
            ....:         for i in ReflexivePolytopes(2) )
            sage: results = ( len(i) == len(j.automorphisms_of_rows_and_columns())
            ....:             for j, i in PMs )
            sage: all(results)  # long time
            True

        TESTS:

        Check that a bug introduced in :issue:`35997` is fixed::

            sage: P = LatticePolytope([(-4,-6),(-4,-5),(0,0),(1,0),(5,6)])
            sage: P._palp_PM_max()
            [9 5 4 0 0]
            [6 0 6 5 0]
            [1 5 0 0 4]
            [0 6 0 1 6]
            [0 0 3 5 3]
        """
        from .palp_normal_form import _palp_PM_max
        return _palp_PM_max(self.vertex_facet_pairing_matrix(), check)

    def npoints(self):
        r"""
        Return the number of lattice points of this polytope.

        EXAMPLES: The number of lattice points of the 3-dimensional
        octahedron and its polar cube::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.npoints()                                                           # needs palp
            7
            sage: cube = o.polar()
            sage: cube.npoints()                                                        # needs palp
            27
        """
        try:
            return self._npoints
        except AttributeError:
            return len(self.points())

    def nvertices(self):
        r"""
        Return the number of vertices of this polytope.

        EXAMPLES: The number of vertices of the 3-dimensional octahedron
        and its polar cube::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.nvertices()
            6
            sage: cube = o.polar()
            sage: cube.nvertices()
            8
        """
        return len(self._vertices)

    @cached_method
    def origin(self):
        r"""
        Return the index of the origin in the list of points of ``self``.

        OUTPUT: integer if the origin belongs to this polytope, ``None`` otherwise

        EXAMPLES::

            sage: p = lattice_polytope.cross_polytope(2)
            sage: p.origin()                                                            # needs palp
            4
            sage: p.point(p.origin())                                                   # needs palp
            M(0, 0)

            sage: p = LatticePolytope(([1],[2]))
            sage: p.points()
            M(1),
            M(2)
            in 1-d lattice M
            sage: print(p.origin())
            None

        Now we make sure that the origin of non-full-dimensional polytopes can
        be identified correctly (:issue:`10661`)::

            sage: LatticePolytope([(1,0,0), (-1,0,0)]).origin()
            2
        """
        origin = self.lattice().zero()
        try:
            return self.points().index(origin)
        except ValueError:
            pass

    def parent(self):
        """
        Return the set of all lattice polytopes.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.parent()
            Set of all Lattice Polytopes
        """
        return SetOfAllLatticePolytopes

    def plot3d(self,
            show_facets=True, facet_opacity=0.5, facet_color=(0,1,0),
            facet_colors=None,
            show_edges=True, edge_thickness=3, edge_color=(0.5,0.5,0.5),
            show_vertices=True, vertex_size=10, vertex_color=(1,0,0),
            show_points=True, point_size=10, point_color=(0,0,1),
            show_vindices=None, vindex_color=(0,0,0),
            vlabels=None,
            show_pindices=None, pindex_color=(0,0,0),
            index_shift=1.1):
        r"""
        Return a 3d-plot of this polytope.

        Polytopes with ambient dimension 1 and 2 will be plotted along x-axis
        or in xy-plane respectively. Polytopes of dimension 3 and less with
        ambient dimension 4 and greater will be plotted in some basis of the
        spanned space.

        By default, everything is shown with more or less pretty
        combination of size and color parameters.

        INPUT:

        Most of the parameters are self-explanatory:

        - ``show_facets`` -- (default: ``True``)

        - ``facet_opacity`` -- (default:0.5)

        - ``facet_color`` -- (default:(0,1,0))

        - ``facet_colors`` -- (default:None) if specified, must be a list of
          colors for each facet separately, used instead of ``facet_color``

        - ``show_edges`` -- boolean (default: ``True``); whether to draw
          edges as lines

        - ``edge_thickness`` -- (default:3)

        - ``edge_color`` -- (default:(0.5,0.5,0.5))

        - ``show_vertices`` -- boolean (default: ``True``); whether to draw
          vertices as balls

        - ``vertex_size`` -- (default:10)

        - ``vertex_color`` -- (default:(1,0,0))

        - ``show_points`` -- boolean (default: ``True``); whether to draw
          other points as balls

        - ``point_size`` -- (default:10)

        - ``point_color`` -- (default:(0,0,1))

        - ``show_vindices`` -- (default: same as
          ``show_vertices``) whether to show indices of vertices

        - ``vindex_color`` -- (default:(0,0,0)) color for
          vertex labels

        - ``vlabels`` -- (default:None) if specified, must be a list of labels
          for each vertex, default labels are vertex indices

        - ``show_pindices`` -- (default: same as ``show_points``)
          whether to show indices of other points

        - ``pindex_color`` -- (default:(0,0,0)) color for
          point labels

        - ``index_shift`` -- (default:1.1)) if 1, labels are
          placed exactly at the corresponding points. Otherwise the label
          position is computed as a multiple of the point position vector.

        EXAMPLES: The default plot of a cube::

            sage: c = lattice_polytope.cross_polytope(3).polar()
            sage: c.plot3d()                                                            # needs palp sage.plot
            Graphics3d Object

        Plot without facets and points, shown without the frame::

            sage: c.plot3d(show_facets=false,                                           # needs palp sage.plot
            ....:          show_points=false).show(frame=False)

        Plot with facets of different colors::

            sage: c.plot3d(facet_colors=rainbow(c.nfacets(), 'rgbtuple'))               # needs palp sage.plot
            Graphics3d Object

        It is also possible to plot lower dimensional polytops in 3D (let's
        also change labels of vertices)::

            sage: c2 = lattice_polytope.cross_polytope(2)
            sage: c2.plot3d(vlabels=["A", "B", "C", "D"])                               # needs palp sage.plot
            Graphics3d Object

        TESTS::

            sage: p = LatticePolytope([[0,0,0],[0,1,1],[1,0,1],[1,1,0]])
            sage: p.plot3d()                                                            # needs palp sage.plot
            Graphics3d Object
        """
        dim = self.dim()
        amb_dim = self.lattice_dim()
        if dim > 3:
            raise ValueError("%d-dimensional polytopes cannot be plotted in 3D!" % self.dim())
        elif amb_dim > 3:
            return self._sublattice_polytope.plot3d(
                show_facets, facet_opacity, facet_color,
                facet_colors,
                show_edges, edge_thickness, edge_color,
                show_vertices, vertex_size, vertex_color,
                show_points, point_size, point_color,
                show_vindices, vindex_color,
                vlabels,
                show_pindices, pindex_color,
                index_shift)
        elif dim == 3:
            vertices = self.vertices()
            if show_points or show_pindices:
                points = self.points()[self.nvertices():]
        else:
            vertices = [vector(ZZ, list(self.vertex(i))+[0]*(3-amb_dim))
                        for i in range(self.nvertices())]
            if show_points or show_pindices:
                points = [vector(ZZ, list(self.point(i))+[0]*(3-amb_dim))
                        for i in range(self.nvertices(), self.npoints())]
        pplot = 0
        if show_facets:
            if dim == 2:
                pplot += IndexFaceSet([self.traverse_boundary()],
                        vertices, opacity=facet_opacity, rgbcolor=facet_color)
            elif dim == 3:
                if facet_colors is None:
                    facet_colors = [facet_color] * self.nfacets()
                vertex_to_index = {v: i for i, v in enumerate(self.vertices())}
                for f, c in zip(self.facets(), facet_colors):
                    pplot += IndexFaceSet([[vertex_to_index[v] for v in f.vertices(f.traverse_boundary())]],
                        vertices, opacity=facet_opacity, rgbcolor=c)
        if show_edges:
            if dim == 1:
                pplot += line3d(vertices, thickness=edge_thickness, rgbcolor=edge_color)
            else:
                for e in self.edges():
                    start, end = e.ambient_vertex_indices()
                    pplot += line3d([vertices[start], vertices[end]],
                            thickness=edge_thickness, rgbcolor=edge_color)
        if show_vertices:
            pplot += point3d(vertices, size=vertex_size, rgbcolor=vertex_color)
        if show_vindices is None:
            show_vindices = show_vertices
        if show_pindices is None:
            show_pindices = show_points
        if show_vindices or show_pindices:
            # Compute the barycenter and shift text of labels away from it
            bc = 1/Integer(len(vertices)) * vector(QQ, sum(vertices))
        if show_vindices:
            if vlabels is None:
                vlabels = list(range(len(vertices)))
            for i, v in enumerate(vertices):
                pplot += text3d(vlabels[i], bc+index_shift*(v-bc), rgbcolor=vindex_color)
        if show_points and len(points):
            pplot += point3d(points, size=point_size, rgbcolor=point_color)
        if show_pindices:
            for i, p in enumerate(points):
                pplot += text3d(i+self.nvertices(), bc+index_shift*(p-bc), rgbcolor=pindex_color)
        return pplot

    def polyhedron(self, **kwds):
        r"""
        Return the Polyhedron object determined by this polytope's vertices.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(2)
            sage: o.polyhedron()
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(vertices=[list(v) for v in self._vertices], **kwds)

    def show3d(self):
        """
        Show a 3d picture of the polytope with default settings and without axes or frame.

        See self.plot3d? for more details.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.show3d()                                                            # needs palp sage.plot
        """
        self.plot3d().show(axis=False, frame=False)

    def point(self, i):
        r"""
        Return the `i`-th point of this polytope, i.e. the `i`-th column of the
        matrix returned by ``points()``.

        EXAMPLES: First few points are actually vertices::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: o.point(1)                                                            # needs palp
            M(0, 1, 0)

        The only other point in the octahedron is the origin::

            sage: o.point(6)                                                            # needs palp
            M(0, 0, 0)
            sage: o.points()                                                            # needs palp
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1),
            M( 0,  0,  0)
            in 3-d lattice M
        """
        return self.points()[i]

    def points(self, *args, **kwds):
        r"""
        Return all lattice points of ``self``.

        INPUT:

        - any arguments given will be passed on to the returned object.

        OUTPUT: a :class:`point collection <PointCollection>`

        EXAMPLES:

        Lattice points of the octahedron and its polar cube::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.points()                                                            # needs palp
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1),
            M( 0,  0,  0)
            in 3-d lattice M
            sage: cube = o.polar()
            sage: cube.points()                                                         # needs palp
            N( 1, -1, -1),
            N( 1,  1, -1),
            N( 1,  1,  1),
            N( 1, -1,  1),
            N(-1, -1,  1),
            N(-1, -1, -1),
            N(-1,  1, -1),
            N(-1,  1,  1),
            N(-1, -1,  0),
            N(-1,  0, -1),
            N(-1,  0,  0),
            N(-1,  0,  1),
            N(-1,  1,  0),
            N( 0, -1, -1),
            N( 0, -1,  0),
            N( 0, -1,  1),
            N( 0,  0, -1),
            N( 0,  0,  0),
            N( 0,  0,  1),
            N( 0,  1, -1),
            N( 0,  1,  0),
            N( 0,  1,  1),
            N( 1, -1,  0),
            N( 1,  0, -1),
            N( 1,  0,  0),
            N( 1,  0,  1),
            N( 1,  1,  0)
            in 3-d lattice N

        Lattice points of a 2-dimensional diamond in a 3-dimensional space::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.points()                                                            # needs palp
            M( 1,  0, 0),
            M( 0,  1, 0),
            M(-1,  0, 0),
            M( 0, -1, 0),
            M( 0,  0, 0)
            in 3-d lattice M

        Only two of the above points::

            sage: p.points(1, 3)                                                        # needs palp
            M(0,  1, 0),
            M(0, -1, 0)
            in 3-d lattice M

        We check that points of a zero-dimensional polytope can be computed::

            sage: p = LatticePolytope([[1]])
            sage: p.points()
            M(1)
            in 1-d lattice M
        """
        if not hasattr(self, "_points"):
            M = self.lattice()
            nv = self.nvertices()
            self._points = points = self._vertices
            if self.dim() == 1:
                v = points[1] - points[0]
                l = gcd(v)
                if l > 1:
                    v = M(v.base_extend(QQ) / l)
                    points = list(points)
                    current = points[0]
                    for i in range(l - 1):
                        current += v
                        current.set_immutable()
                        points.append(current)
            if self.dim() > 1:
                result = self.poly_x("p", reduce_dimension=True)
                if self.dim() == self.lattice_dim():
                    points = read_palp_point_collection(StringIO(result), M)
                else:
                    m = self._embed(read_palp_matrix(result))
                    if m.ncols() > nv:
                        points = list(points)
                        for j in range(nv, m.ncols()):
                            current = M.element_class(
                                M, [m[i, j] for i in range(M.rank())])
                            current.set_immutable()
                            points.append(current)
            if len(points) > nv:
                self._points = PointCollection(points, M)
        if args or kwds:
            return self._points(*args, **kwds)
        else:
            return self._points

    def _some_elements_(self):
        r"""
        Generate some points of ``self`` as a convex polytope.

        In contrast to :meth:`points`, these are not necessarily lattice points.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.some_elements()  # indirect doctest
            [(1, 0, 0),
            (1/2, 1/2, 0),
            (1/4, 1/4, 1/2),
            (-3/8, 1/8, 1/4),
            (-3/16, -7/16, 1/8),
            (-3/32, -7/32, -7/16)]
        """
        if not self._vertices:
            return
        V = self.ambient_vector_space()
        v_iter = iter(self._vertices)
        p = V(next(v_iter))
        yield p
        for i in range(5):
            try:
                p = (p + next(v_iter)) / 2
            except StopIteration:
                return
            yield p

    def polar(self):
        r"""
        Return the polar polytope, if this polytope is reflexive.

        EXAMPLES: The polar polytope to the 3-dimensional octahedron::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: cube = o.polar()
            sage: cube
            3-d reflexive polytope in 3-d lattice N

        The polar polytope "remembers" the original one::

            sage: cube.polar()
            3-d reflexive polytope in 3-d lattice M
            sage: cube.polar().polar() is cube
            True

        Only reflexive polytopes have polars::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,2),
            ....:                      (-1,0,0), (0,-1,0), (0,0,-1)])
            sage: p.polar()
            Traceback (most recent call last):
            ...
            ValueError: The given polytope is not reflexive!
            Polytope: 3-d lattice polytope in 3-d lattice M
        """
        if self.is_reflexive():
            return self._polar
        else:
            raise ValueError(("The given polytope is not reflexive!\n"
                                + "Polytope: %s") % self)

    def poly_x(self, keys, reduce_dimension=False):
        r"""
        Run ``poly.x`` with given ``keys`` on vertices of this
        polytope.

        INPUT:

        - ``keys`` -- string of options passed to ``poly.x``. The
          key "f" is added automatically

        - ``reduce_dimension`` -- boolean (default: ``False``); if ``True`` and this
          polytope is not full-dimensional, ``poly.x`` will be called for the
          vertices of this polytope in some basis of the spanned affine space

        OUTPUT: the output of ``poly.x`` as a string

        EXAMPLES: This call is used for determining if a polytope is
        reflexive or not::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: print(o.poly_x("e"))                                                  # needs palp
            8 3  Vertices of P-dual <-> Equations of P
              -1  -1   1
               1  -1   1
              -1   1   1
               1   1   1
              -1  -1  -1
               1  -1  -1
              -1   1  -1
               1   1  -1

        Since PALP has limits on different parameters determined during
        compilation, the following code is likely to fail, unless you
        change default settings of PALP::

            sage: BIG = lattice_polytope.cross_polytope(7)
            sage: BIG
            7-d reflexive polytope in 7-d lattice M
            sage: BIG.poly_x("e")                                                       # needs palp
            Traceback (most recent call last):
            ...
            ValueError: Error executing 'poly.x -fe' for the given polytope!
            Output:
            Please increase POLY_Dmax to at least 7

        You cannot call ``poly.x`` for polytopes that don't span the space (if you
        could, it would crush anyway)::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
            sage: p.poly_x("e")                                                         # needs palp
            Traceback (most recent call last):
            ...
            ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!

        But if you know what you are doing, you can call it for the polytope in
        some basis of the spanned space::

            sage: print(p.poly_x("e", reduce_dimension=True))                           # needs palp
            4 2  Equations of P
              -1   1     0
               1   1     2
              -1  -1     0
               1  -1     2
        """
        return self._palp("poly.x -f" + keys, reduce_dimension)

    @cached_method
    def skeleton(self):
        r"""
        Return the graph of the one-skeleton of this polytope.

        EXAMPLES::

            sage: d = lattice_polytope.cross_polytope(2)
            sage: g = d.skeleton(); g                                                   # needs palp sage.graphs
            Graph on 4 vertices
            sage: g.edges(sort=True)                                                    # needs palp sage.graphs
            [(0, 1, None), (0, 3, None), (1, 2, None), (2, 3, None)]
        """
        from sage.graphs.graph import Graph
        skeleton = Graph()
        skeleton.add_vertices(self.skeleton_points(1))
        for edge in self.edges():
            points = edge.ambient_ordered_point_indices()
            for i in range(len(points) - 1):
                skeleton.add_edge(points[i], points[i + 1])
        return skeleton.copy(immutable=True)

    def skeleton_points(self, k=1):
        r"""
        Return the increasing list of indices of lattice points in
        k-skeleton of the polytope (k is 1 by default).

        EXAMPLES: We compute all skeleton points for the cube::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: c = o.polar()
            sage: c.skeleton_points()                                                   # needs palp sage.graphs
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        The default was 1-skeleton::

            sage: c.skeleton_points(k=1)                                                # needs palp sage.graphs
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        0-skeleton just lists all vertices::

            sage: c.skeleton_points(k=0)                                                # needs palp sage.graphs
            [0, 1, 2, 3, 4, 5, 6, 7]

        2-skeleton lists all points except for the origin (point #17)::

            sage: c.skeleton_points(k=2)                                                # needs palp sage.graphs
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
             18, 19, 20, 21, 22, 23, 24, 25, 26]

        3-skeleton includes all points::

            sage: c.skeleton_points(k=3)                                                # needs palp
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23, 24, 25, 26]

        It is OK to compute higher dimensional skeletons - you will get the
        list of all points::

            sage: c.skeleton_points(k=100)                                              # needs palp
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23, 24, 25, 26]
        """
        if k >= self.dim():
            return list(range(self.npoints()))
        skeleton = set([])
        for face in self.faces(dim=k):
            skeleton.update(face.ambient_point_indices())
        skeleton = sorted(skeleton)
        return skeleton

    def skeleton_show(self, normal=None):
        r"""Show the graph of one-skeleton of this polytope.
        Works only for polytopes in a 3-dimensional space.

        INPUT:

        - ``normal`` -- a 3-dimensional vector (can be given as a list), which
          should be perpendicular to the screen. If not given, will be selected
          randomly (new each time and it may be far from "nice").

        EXAMPLES: Show a pretty picture of the octahedron::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.skeleton_show([1,2,4])                                              # needs palp sage.plot

        Does not work for a diamond at the moment::

            sage: d = lattice_polytope.cross_polytope(2)
            sage: d.skeleton_show()
            Traceback (most recent call last):
            ...
            NotImplementedError: skeleton view is implemented only in 3-d space
        """
        if self.lattice_dim() != 3:
            raise NotImplementedError("skeleton view is implemented only in 3-d space")
        if normal is None:
            normal = [ZZ.random_element(20),ZZ.random_element(20),ZZ.random_element(20)]
        normal = matrix(QQ,3,1,list(normal))
        projectionm = normal.kernel().basis_matrix()
        positions = dict(enumerate([list(c) for c in (projectionm*self.points()).columns(copy=False)]))
        self.skeleton().show(pos=positions)

    def traverse_boundary(self):
        r"""
        Return a list of indices of vertices of a 2-dimensional polytope in their boundary order.

        Needed for plot3d function of polytopes.

        EXAMPLES::

            sage: p = lattice_polytope.cross_polytope(2).polar()
            sage: p.traverse_boundary()                                                 # needs sage.graphs
            [3, 0, 1, 2]
        """
        if self.dim() != 2:
            raise ValueError("Boundary can be traversed only for 2-polytopes!")
        zero_faces = set(self.faces(0))
        l = [self.faces(0)[0]]
        prev, next = sorted(zero_faces.intersection(l[0].adjacent()))
        l = [prev, l[0], next]
        while len(l) < self.nvertices():
            prev, next = zero_faces.intersection(l[-1].adjacent())
            if next == l[-2]:
                next = prev
            l.append(next)
        vertex_to_index = {v: i for i, v in enumerate(self.vertices())}
        return [vertex_to_index[v.vertex(0)] for v in l]

    def vertex(self, i):
        r"""
        Return the `i`-th vertex of this polytope, i.e. the `i`-th column of
        the matrix returned by ``vertices()``.

        EXAMPLES: Note that numeration starts with zero::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: o.vertex(3)
            M(-1, 0, 0)
        """
        return self._vertices[i]

    def vertex_facet_pairing_matrix(self):
        r"""
        Return the vertex facet pairing matrix `PM`.

        Return a matrix whose the `i, j^\text{th}` entry is the height
        of the `j^\text{th}` vertex over the `i^\text{th}` facet.
        The ordering of the vertices and facets is as in
        :meth:`vertices` and :meth:`facets`.

        EXAMPLES::

            sage: L = lattice_polytope.cross_polytope(3)
            sage: L.vertex_facet_pairing_matrix()
            [2 0 0 0 2 2]
            [2 2 0 0 0 2]
            [2 2 2 0 0 0]
            [2 0 2 0 2 0]
            [0 0 2 2 2 0]
            [0 0 0 2 2 2]
            [0 2 0 2 0 2]
            [0 2 2 2 0 0]
        """
        V = self.vertices()
        nv = self.nvertices()
        PM = matrix(ZZ, [n * V + vector(ZZ, [c] * nv)
            for n, c in zip(self.facet_normals(), self.facet_constants())])
        PM.set_immutable()
        return PM

    def vertices(self, *args, **kwds):
        r"""
        Return vertices of ``self``.

        INPUT:

        - any arguments given will be passed on to the returned object.

        OUTPUT: a :class:`point collection <PointCollection>`

        EXAMPLES:

        Vertices of the octahedron and its polar cube are in dual lattices::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: o.vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: cube = o.polar()
            sage: cube.vertices()
            N( 1, -1, -1),
            N( 1,  1, -1),
            N( 1,  1,  1),
            N( 1, -1,  1),
            N(-1, -1,  1),
            N(-1, -1, -1),
            N(-1,  1, -1),
            N(-1,  1,  1)
            in 3-d lattice N
        """
        if args or kwds:
            return self._vertices(*args, **kwds)
        else:
            return self._vertices


def is_NefPartition(x):
    r"""
    Check if ``x`` is a nef-partition.

    INPUT:

    - ``x`` -- anything

    OUTPUT:

    - ``True`` if ``x`` is a :class:`nef-partition <NefPartition>` and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.lattice_polytope import NefPartition
        sage: isinstance(1, NefPartition)
        False
        sage: o = lattice_polytope.cross_polytope(3)
        sage: np = o.nef_partitions()[0]; np                                            # needs palp
        Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
        sage: isinstance(np, NefPartition)                                              # needs palp
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38126,
                "The function is_NefPartition is deprecated; "
                "use 'isinstance(..., NefPartition)' instead.")
    return isinstance(x, NefPartition)


class NefPartition(SageObject, Hashable):
    r"""
    Create a nef-partition.

    INPUT:

    - ``data`` -- list of integers, the `i`-th element of this list must be
      the part of the `i`-th vertex of ``Delta_polar`` in this nef-partition;

    - ``Delta_polar`` -- a :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`;

    - ``check`` -- by default the input will be checked for correctness, i.e.
      that ``data`` indeed specify a nef-partition. If you are sure that the
      input is correct, you can speed up construction via ``check=False``
      option.

    OUTPUT: a nef-partition of ``Delta_polar``

    Let `M` and `N` be dual lattices. Let `\Delta \subset M_\RR` be a reflexive
    polytope with polar `\Delta^\circ \subset N_\RR`. Let `X_\Delta` be the
    toric variety associated to the normal fan of `\Delta`. A **nef-partition**
    is a decomposition of the vertex set `V` of `\Delta^\circ` into a disjoint
    union `V = V_0 \sqcup V_1 \sqcup \dots \sqcup V_{k-1}` such that divisors
    `E_i = \sum_{v\in V_i} D_v` are Cartier (here `D_v` are prime
    torus-invariant Weil divisors corresponding to vertices of `\Delta^\circ`).
    Equivalently, let `\nabla_i \subset N_\RR` be the convex hull of vertices
    from `V_i` and the origin. These polytopes form a nef-partition if their
    Minkowski sum `\nabla \subset N_\RR` is a reflexive polytope.

    The **dual nef-partition** is formed by polytopes `\Delta_i \subset M_\RR`
    of `E_i`, which give a decomposition of the vertex set of `\nabla^\circ
    \subset M_\RR` and their Minkowski sum is `\Delta`, i.e. the polar duality
    of reflexive polytopes switches convex hull and Minkowski sum for dual
    nef-partitions:

    .. MATH::

        \Delta^\circ
        &=
        \mathrm{Conv} \left(\nabla_0, \nabla_1, \dots, \nabla_{k-1}\right), \\
        \nabla^{\phantom{\circ}}
        &=
        \nabla_0 + \nabla_1 + \dots + \nabla_{k-1}, \\
        &
        \\
        \Delta^{\phantom{\circ}}
        &=
        \Delta_0 + \Delta_1 + \dots + \Delta_{k-1}, \\
        \nabla^\circ
        &=
        \mathrm{Conv} \left(\Delta_0, \Delta_1, \dots, \Delta_{k-1}\right).

    One can also interpret the duality of nef-partitions as the duality of the
    associated cones. Below `\overline{M} = M \times \ZZ^k` and
    `\overline{N} = N \times \ZZ^k` are dual lattices.

    The **Cayley polytope** `P \subset \overline{M}_\RR` of a nef-partition is
    given by `P = \mathrm{Conv}(\Delta_0 \times e_0, \Delta_1 \times e_1,
    \ldots, \Delta_{k-1} \times e_{k-1})`, where `\{e_i\}_{i=0}^{k-1}` is the
    standard basis of `\ZZ^k`. The **dual Cayley polytope**
    `P^* \subset \overline{N}_\RR` is the Cayley polytope of the dual
    nef-partition.

    The **Cayley cone** `C \subset \overline{M}_\RR` of a nef-partition is the
    cone spanned by its Cayley polytope. The **dual Cayley cone**
    `C^\vee \subset \overline{M}_\RR` is the usual dual cone of `C`. It turns
    out, that `C^\vee` is spanned by `P^*`.

    It is also possible to go back from the Cayley cone to the Cayley polytope,
    since `C` is a reflexive Gorenstein cone supported by `P`: primitive
    integral ray generators of `C` are contained in an affine hyperplane and
    coincide with vertices of `P`.

    See Section 4.3.1 in [CK1999]_ and references therein for further details, or
    [BN2008]_ for a purely combinatorial approach.

    EXAMPLES:

    It is very easy to create a nef-partition for the octahedron, since for
    this polytope any decomposition of vertices is a nef-partition. We create a
    3-part nef-partition with the 0th and 1st vertices belonging to the 0th
    part (recall that numeration in Sage starts with 0), the 2nd and 5th
    vertices belonging to the 1st part, and 3rd and 4th vertices belonging
    to the 2nd part::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: np = NefPartition([0,0,1,2,2,1], o)
        sage: np
        Nef-partition {0, 1} ⊔ {2, 5} ⊔ {3, 4}

    The octahedron plays the role of `\Delta^\circ` in the above description::

        sage: np.Delta_polar() is o
        True

    The dual nef-partition (corresponding to the "mirror complete
    intersection") gives decomposition of the vertex set of `\nabla^\circ`::

        sage: np.dual()
        Nef-partition {0, 1, 2} ⊔ {3, 4} ⊔ {5, 6, 7}
        sage: np.nabla_polar().vertices()
        N(-1, -1,  0),
        N(-1,  0,  0),
        N( 0, -1,  0),
        N( 0,  0, -1),
        N( 0,  0,  1),
        N( 1,  0,  0),
        N( 0,  1,  0),
        N( 1,  1,  0)
        in 3-d lattice N

    Of course, `\nabla^\circ` is `\Delta^\circ` from the point of view of the
    dual nef-partition::

        sage: np.dual().Delta_polar() is np.nabla_polar()
        True
        sage: np.Delta(1).vertices()
        N(0, 0, -1),
        N(0, 0,  1)
        in 3-d lattice N
        sage: np.dual().nabla(1).vertices()
        N(0, 0, -1),
        N(0, 0,  1)
        in 3-d lattice N

    Instead of constructing nef-partitions directly, you can request all 2-part
    nef-partitions of a given reflexive polytope (they will be computed using
    ``nef.x`` program from PALP)::

        sage: o.nef_partitions()                                                        # needs palp
        [Nef-partition {0, 1, 3} ⊔ {2, 4, 5},
         Nef-partition {0, 1, 3, 4} ⊔ {2, 5} (direct product),
         Nef-partition {0, 1, 2} ⊔ {3, 4, 5},
         Nef-partition {0, 1, 2, 3} ⊔ {4, 5},
         Nef-partition {0, 1, 2, 3, 4} ⊔ {5} (projection)]
    """

    def __init__(self, data, Delta_polar, check=True):
        r"""
        See :class:`NefPartition` for documentation.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = o.nef_partitions()[0]                                            # needs palp
            sage: TestSuite(np).run()                                                   # needs palp
        """
        if check and not Delta_polar.is_reflexive():
            raise ValueError("nef-partitions can be constructed for reflexive "
                             "polytopes ony!")
        self._vertex_to_part = tuple(int(el) for el in data)
        self._nparts = max(self._vertex_to_part) + 1
        self._Delta_polar = Delta_polar
        if check and not self.nabla().is_reflexive():
            raise ValueError("%s do not form a nef-partition!" % str(data))

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``True`` if ``other`` is a :class:`nef-partition <NefPartition>`
          equal to ``self``, ``False`` otherwise.

        .. NOTE::

            Two nef-partitions are equal if they correspond to equal polytopes
            and their parts are the same, including their order.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o)
            sage: np == np
            True
            sage: np == o.nef_partitions()[0]                                           # needs palp
            True
            sage: np == o.nef_partitions()[1]                                           # needs palp
            False
            sage: np2 = NefPartition(np._vertex_to_part, o)
            sage: np2 is np
            False
            sage: np2 == np
            True
            sage: np == 0
            False
        """
        return (isinstance(other, NefPartition)
                and self._Delta_polar == other._Delta_polar
                and self._vertex_to_part == other._vertex_to_part)

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT: integer

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o)
            sage: hash(np) == hash(np)
            True
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(self._vertex_to_part) + hash(self._Delta_polar)
            return self._hash

    def __ne__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        - ``False`` if ``other`` is a :class:`nef-partition <NefPartition>`
          equal to ``self``, ``True`` otherwise.

        .. NOTE::

            Two nef-partitions are equal if they correspond to equal polytopes
            and their parts are the same, including their order.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o)
            sage: np != np
            False
            sage: np != o.nef_partitions()[0]                                           # needs palp
            False
            sage: np != o.nef_partitions()[1]                                           # needs palp
            True
            sage: np2 = NefPartition(np._vertex_to_part, o)
            sage: np2 is np
            False
            sage: np2 != np
            False
            sage: np != 0
            True
        """
        return not (self == other)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o)
            sage: latex(np)  # indirect doctest
            \text{Nef-partition } \{0, 1, 3\} \sqcup \{2, 4, 5\}
        """
        result = r"\text{Nef-partition } "
        for i, part in enumerate(self.parts()):
            if i != 0:
                result += r" \sqcup "
            result += r"\{" + ", ".join("%d" % v for v in part) + r"\}"
        try:
            # We may or may not know the type of the partition
            if self._is_product:
                result += r" \text{ (direct product)}"
            if self._is_projection:
                result += r" \text{ (projection)}"
        except AttributeError:
            pass
        return result

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o)
            sage: repr(np)  # indirect doctest
            'Nef-partition {0, 1, 3} ⊔ {2, 4, 5}'
        """
        result = "Nef-partition "
        for i, part in enumerate(self.parts()):
            if i != 0:
                result += " ⊔ "
            result += "{" + ", ".join("%d" % v for v in part) + "}"
        try:
            # We may or may not know the type of the partition
            if self._is_product:
                result += " (direct product)"
            if self._is_projection:
                result += " (projection)"
        except AttributeError:
            pass
        return result

    def _sage_input_(self, sib, coerced):
        """
        Return Sage command to reconstruct ``self``.

        See :mod:`sage.misc.sage_input` for details.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: sage_input(np, verify=True)
            # Verified
            NefPartition([0, 0, 1, 0, 1, 1],
                         LatticePolytope(sage.geometry.point_collection.PointCollection((vector(ZZ, [1, 0, 0]),
                                                                                         vector(ZZ, [0, 1, 0]),
                                                                                         vector(ZZ, [0, 0, 1]),
                                                                                         vector(ZZ, [-1, 0, 0]),
                                                                                         vector(ZZ, [0, -1, 0]),
                                                                                         vector(ZZ, [0, 0, -1]))),
                                         compute_vertices=False))
        """
        vertex_to_part = [ZZ(i) for i in self._vertex_to_part]
        return sib.name('NefPartition')(vertex_to_part, sib(self.Delta_polar()))

    def Delta(self, i=None):
        r"""
        Return the polytope `\Delta` or `\Delta_i` corresponding to ``self``.

        INPUT:

        - ``i`` -- integer; if not given, `\Delta` will be returned

        OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.Delta().polar() is o
            True
            sage: np.Delta().vertices()
            N( 1, -1, -1),
            N( 1,  1, -1),
            N( 1,  1,  1),
            N( 1, -1,  1),
            N(-1, -1,  1),
            N(-1, -1, -1),
            N(-1,  1, -1),
            N(-1,  1,  1)
            in 3-d lattice N
            sage: np.Delta(0).vertices()
            N(-1, -1, 0),
            N(-1,  0, 0),
            N( 1,  0, 0),
            N( 1, -1, 0)
            in 3-d lattice N
        """
        if i is None:
            return self._Delta_polar.polar()
        else:
            return self.dual().nabla(i)

    def Delta_polar(self):
        r"""
        Return the polytope `\Delta^\circ` corresponding to ``self``.

        OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.Delta_polar() is o
            True
        """
        return self._Delta_polar

    def Deltas(self):
        r"""
        Return the polytopes `\Delta_i` corresponding to ``self``.

        OUTPUT: a tuple of :class:`lattice polytopes <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.Delta().vertices()
            N( 1, -1, -1),
            N( 1,  1, -1),
            N( 1,  1,  1),
            N( 1, -1,  1),
            N(-1, -1,  1),
            N(-1, -1, -1),
            N(-1,  1, -1),
            N(-1,  1,  1)
            in 3-d lattice N
            sage: [Delta_i.vertices() for Delta_i in np.Deltas()]
            [N(-1, -1, 0),
             N(-1,  0, 0),
             N( 1,  0, 0),
             N( 1, -1, 0)
             in 3-d lattice N,
             N(0, 0, -1),
             N(0, 1,  1),
             N(0, 0,  1),
             N(0, 1, -1)
             in 3-d lattice N]
            sage: np.nabla_polar().vertices()
            N(-1, -1,  0),
            N( 1, -1,  0),
            N( 1,  0,  0),
            N(-1,  0,  0),
            N( 0,  1, -1),
            N( 0,  1,  1),
            N( 0,  0,  1),
            N( 0,  0, -1)
            in 3-d lattice N
        """
        return self.dual().nablas()

    @cached_method
    def dual(self):
        r"""
        Return the dual nef-partition.

        OUTPUT: a :class:`nef-partition <NefPartition>`

        See the class documentation for the definition.

        ALGORITHM:

        See Proposition 3.19 in [BN2008]_.

        .. NOTE::

            Automatically constructed dual nef-partitions will be ordered, i.e.
            vertex partition of `\nabla` will look like
            `\{0, 1, 2\} \sqcup \{3, 4, 5, 6\} \sqcup \{7, 8\}`.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.dual()
            Nef-partition {0, 1, 2, 3} ⊔ {4, 5, 6, 7}
            sage: np.dual().Delta() is np.nabla()
            True
            sage: np.dual().nabla(0) is np.Delta(0)
            True
        """
        # Delta and nabla are interchanged compared to [BN2008]_.
        # The order of vertices of this nabla_polar will be adjusted.
        nabla_polar = LatticePolytope(
            reduce(minkowski_sum,
                   (nabla.vertices() for nabla in self.nablas())),
            lattice=self._Delta_polar.lattice()).polar()
        vertex_to_part = []
        nabla_polar_vertices = []
        for i in range(self._nparts):
            A = nabla_polar.vertices().matrix() * self.nabla(i).vertices()
            for j, row in enumerate(A):
                if min(row) == -1:
                    vertex_to_part.append(i)
                    nabla_polar_vertices.append(nabla_polar.vertex(j))
        # Make dual look "ordered", like {0,1,2} ⊔ {3,4,5,6} ⊔ {7,8}.
        nabla_polar = LatticePolytope(nabla_polar_vertices,
                                      compute_vertices=False)
        # If self is a valid nef-partition, the dual is as well.
        dual = NefPartition(vertex_to_part, nabla_polar, check=False)
        dual.dual.set_cache(self)
        return dual

    def hodge_numbers(self):
        r"""
        Return Hodge numbers corresponding to ``self``.

        OUTPUT: a tuple of integers (produced by ``nef.x`` program from PALP)

        EXAMPLES:

        Currently, you need to request Hodge numbers when you compute
        nef-partitions::

            sage: # long time, needs palp
            sage: p = lattice_polytope.cross_polytope(5)
            sage: np = p.nef_partitions()[0]                    # 4s on sage.math, 2011
            sage: np.hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
            sage: np = p.nef_partitions(hodge_numbers=True)[0]  # 13s on sage.math, 2011
            sage: np.hodge_numbers()
            (19, 19)
        """
        try:
            return self._hodge_numbers
        except AttributeError:
            self._Delta_polar._compute_hodge_numbers()
            return self._hodge_numbers

    def nabla(self, i=None):
        r"""
        Return the polytope `\nabla` or `\nabla_i` corresponding to ``self``.

        INPUT:

        - ``i`` -- integer; if not given, `\nabla` will be returned

        OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.Delta_polar().vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: np.nabla(0).vertices()
            M(-1, 0, 0),
            M( 1, 0, 0),
            M( 0, 1, 0)
            in 3-d lattice M
            sage: np.nabla().vertices()
            M(-1,  0,  1),
            M(-1,  0, -1),
            M( 1,  0,  1),
            M( 1,  0, -1),
            M( 0,  1,  1),
            M( 0,  1, -1),
            M( 1, -1,  0),
            M(-1, -1,  0)
            in 3-d lattice M
        """
        if i is None:
            return self.dual().Delta()
        else:
            return self.nablas()[i]

    def nabla_polar(self):
        r"""
        Return the polytope `\nabla^\circ` corresponding to ``self``.

        OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.nabla_polar().vertices()
            N(-1, -1,  0),
            N( 1, -1,  0),
            N( 1,  0,  0),
            N(-1,  0,  0),
            N( 0,  1, -1),
            N( 0,  1,  1),
            N( 0,  0,  1),
            N( 0,  0, -1)
            in 3-d lattice N
            sage: np.nabla_polar() is np.dual().Delta_polar()
            True
        """
        return self.nabla().polar()

    def nablas(self):
        r"""
        Return the polytopes `\nabla_i` corresponding to ``self``.

        OUTPUT: a tuple of :class:`lattice polytopes <LatticePolytopeClass>`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.Delta_polar().vertices()
            M( 1,  0,  0),
            M( 0,  1,  0),
            M( 0,  0,  1),
            M(-1,  0,  0),
            M( 0, -1,  0),
            M( 0,  0, -1)
            in 3-d lattice M
            sage: [nabla_i.vertices() for nabla_i in np.nablas()]
            [M(-1, 0, 0),
             M( 1, 0, 0),
             M( 0, 1, 0)
             in 3-d lattice M,
             M(0, -1,  0),
             M(0,  0, -1),
             M(0,  0,  1)
             in 3-d lattice M]
        """
        try:
            return self._nablas
        except AttributeError:
            Delta_polar = self._Delta_polar
            origin = [[0] * Delta_polar.dim()]
            self._nablas = tuple(LatticePolytope(
                                [Delta_polar.vertex(j) for j in part] + origin,
                                lattice=Delta_polar.lattice())
                                for part in self.parts())
            return self._nablas

    def nparts(self):
        r"""
        Return the number of parts in ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.nparts()
            2
        """
        return self._nparts

    def part(self, i, all_points=False):
        r"""
        Return the ``i``-th part of ``self``.

        INPUT:

        - ``i`` -- integer

        - ``all_points`` -- boolean (default: ``False``); whether to list all lattice points
          or just vertices

        OUTPUT:

        - a tuple of integers, indices of vertices (or all lattice points) of
          `\Delta^\circ` belonging to `V_i`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.part(0)
            (0, 1, 3)
            sage: np.part(0, all_points=True)                                           # needs palp
            (0, 1, 3)
            sage: np.dual().part(0)
            (0, 1, 2, 3)
            sage: np.dual().part(0, all_points=True)                                    # needs palp
            (0, 1, 2, 3, 8)
        """
        return self.parts(all_points)[i]

    @cached_method
    def parts(self, all_points=False):
        r"""
        Return all parts of ``self``.

        INPUT:

        - ``all_points`` -- boolean (default: ``False``); whether to list all lattice points
          or just vertices

        OUTPUT:

        - a tuple of tuples of integers. The `i`-th tuple contains indices of
          vertices (or all lattice points) of `\Delta^\circ` belonging to `V_i`

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.parts()
            ((0, 1, 3), (2, 4, 5))
            sage: np.parts(all_points=True)                                             # needs palp
            ((0, 1, 3), (2, 4, 5))
            sage: np.dual().parts()
            ((0, 1, 2, 3), (4, 5, 6, 7))
            sage: np.dual().parts(all_points=True)                                      # needs palp
            ((0, 1, 2, 3, 8), (4, 5, 6, 7, 10))
        """
        parts = [[] for _ in range(self._nparts)]
        if all_points:
            for point in range(self._Delta_polar.npoints()):
                if point != self._Delta_polar.origin():
                    parts[self.part_of_point(point)].append(point)
        else:
            for vertex, part in enumerate(self._vertex_to_part):
                parts[part].append(vertex)
        return tuple(tuple(part) for part in parts)

    def part_of(self, i):
        r"""
        Return the index of the part containing the ``i``-th vertex.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        - an integer `j` such that the ``i``-th vertex of `\Delta^\circ`
          belongs to `V_j`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = NefPartition([0, 0, 1, 0, 1, 1], o); np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: np.part_of(3)
            0
            sage: np.part_of(2)
            1
        """
        return self._vertex_to_part[i]

    @cached_method
    def part_of_point(self, i):
        r"""
        Return the index of the part containing the ``i``-th point.

        INPUT:

        - ``i`` -- integer

        OUTPUT:

        - an integer `j` such that the ``i``-th point of `\Delta^\circ`
          belongs to `\nabla_j`.

        .. NOTE::

            Since a nef-partition induces a partition on the set of boundary
            lattice points of `\Delta^\circ`, the value of `j` is well-defined
            for all `i` but the one that corresponds to the origin, in which
            case this method will raise a :exc:`ValueError` exception.
            (The origin always belongs to all `\nabla_j`.)

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES:

        We consider a relatively complicated reflexive polytope #2252 (easily
        accessible in Sage as ``ReflexivePolytope(3, 2252)``, we create it here
        explicitly to avoid loading the whole database)::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,1), (0,1,-1),
            ....:         (0,-1,1), (-1,1,0), (0,-1,-1), (-1,-1,0), (-1,-1,2)])
            sage: np = p.nef_partitions()[0]; np                                        # needs palp
            Nef-partition {1, 2, 5, 7, 8} ⊔ {0, 3, 4, 6}
            sage: p.nvertices()
            9
            sage: p.npoints()                                                           # needs palp
            15

        We see that the polytope has 6 more points in addition to vertices. One
        of them is the origin::

            sage: p.origin()                                                            # needs palp
            14
            sage: np.part_of_point(14)                                                  # needs palp
            Traceback (most recent call last):
            ...
            ValueError: the origin belongs to all parts!

        But the remaining 5 are partitioned by ``np``::

            sage: [n for n in range(p.npoints())                                        # needs palp
            ....:    if p.origin() != n and np.part_of_point(n) == 0]
            [1, 2, 5, 7, 8, 9, 11, 13]
            sage: [n for n in range(p.npoints())                                        # needs palp
            ....:    if p.origin() != n and np.part_of_point(n) == 1]
            [0, 3, 4, 6, 10, 12]
        """
        if i < self._Delta_polar.nvertices():
            return self.part_of(i)
        if i == self._Delta_polar.origin():
            raise ValueError("the origin belongs to all parts!")
        point = self._Delta_polar.point(i)
        for part, nabla in enumerate(self.nablas()):
            if point in nabla:
                return part


_palp_dimension = None


def _palp(command, polytopes, reduce_dimension=False):
    r"""
    Run ``command`` on vertices of given
    ``polytopes``.

    Returns the name of the file containing the output of
    ``command``. You should delete it after using.

    .. NOTE::

      PALP cannot be called for polytopes that do not span the ambient space.
      If you specify ``reduce_dimension=True`` argument, PALP will be
      called for vertices of this polytope in some basis of the affine space
      it spans.

    TESTS::

        sage: # needs palp
        sage: o = lattice_polytope.cross_polytope(3)
        sage: result_name = lattice_polytope._palp("poly.x -f", [o])
        sage: f = open(result_name)
        sage: f.readlines()
        ['M:7 6 N:27 8 Pic:17 Cor:0\n']
        sage: f.close()
        sage: os.remove(result_name)

        sage: p = LatticePolytope([(1,0,0), (0,1,0), (-1,0,0), (0,-1,0)])
        sage: lattice_polytope._palp("poly.x -f", [p])                                  # needs palp
        Traceback (most recent call last):
        ...
        ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!

        sage: # needs palp
        sage: result_name = lattice_polytope._palp("poly.x -f", [p],
        ....:                                      reduce_dimension=True)
        sage: f = open(result_name)
        sage: f.readlines()
        ['M:5 4 F:4\n']
        sage: f.close()
        sage: os.remove(result_name)
    """
    if _palp_dimension is not None:
        dot = command.find(".")
        command = command[:dot] + "-%dd" % _palp_dimension + command[dot:]
    executable, args = command.split(" ", 1)
    if executable.endswith('.x'):
        executable = executable[:-2]
    executable = PalpExecutable(executable).absolute_filename()
    command = " ".join([shlex.quote(executable), args])
    input_file_name = tmp_filename()
    input_file = open(input_file_name, "w")
    for p in polytopes:
        if p.dim() == 0:
            raise ValueError(("Cannot run \"%s\" for the zero-dimensional "
                + "polytope!\nPolytope: %s") % (command, p))
        if p.dim() < p.lattice_dim():
            if not reduce_dimension:
                raise ValueError(("Cannot run PALP for a %d-dimensional polytope " +
                "in a %d-dimensional space!") % (p.dim(), p.lattice_dim()))
            write_palp_matrix(p._pullback(p._vertices), input_file)
        else:
            p._vertices.write_for_palp(input_file)
    input_file.close()
    output_file_name = tmp_filename()
    c = "%s <%s >%s" % (command, input_file_name, output_file_name)
    p = Popen(c, shell=True, bufsize=2048,
              stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
    stderr = p.stderr
    err = stderr.read()
    if len(err):
        raise RuntimeError(("Error executing \"%s\" for a polytope sequence!"
            + "\nOutput:\n%s") % (command, err))
    os.remove(input_file_name)
    try:
        p.terminate()
    except OSError:
        pass
    return output_file_name


def _palp_canonical_order(V, PM_max, permutations):
    r"""
    Compute the PALP normal form of the vertices V
    using auxiliary data computed elsewhere.

    This is a helper function for
    :meth:`~sage.geometry.lattice_polytope.LatticePolytopeClass.normal_form`
    and should not be called directly.

    Given a matrix of vertices, the maximal vertex-facet pairing matrix
    and the permutations realizing this matrix, apply the last part of the
    PALP algorithm and return the normal form.

    INPUT:

    - ``V`` -- :class:`point collection <PointCollection>`. The vertices

    - ``PM_max`` -- the maximal vertex-facet pairing matrix

    - ``permutations`` -- the permutations of the vertices yielding
        ``PM_max``

    OUTPUT:

    The PALP normal form as a :class:`point collection <PointCollection>`
    and a permutation.

    TESTS::

        sage: L = lattice_polytope.cross_polytope(2)
        sage: V = L.vertices()
        sage: PM_max, permutations = L._palp_PM_max(check=True)                         # needs sage.groups
        sage: from sage.geometry.lattice_polytope import _palp_canonical_order
        sage: _palp_canonical_order(V, PM_max, permutations)                            # needs sage.groups
        (M( 1,  0),
         M( 0,  1),
         M( 0, -1),
         M(-1,  0)
         in 2-d lattice M, (1,3,2,4))
    """
    from .palp_normal_form import _palp_canonical_order

    Vmatrix = V.column_matrix()
    Vmodule = V.module()
    vertices, permutation = _palp_canonical_order(Vmatrix, PM_max, permutations)
    vertices = [Vmodule(v) for v in vertices]
    for v in vertices:
        v.set_immutable()
    return PointCollection(vertices, Vmodule), permutation


def _palp_convert_permutation(permutation):
    r"""
    Convert a permutation from PALPs notation to a Sage permutation.

    PALP specifies a permutation group element by its domain. Furthermore,
    it only supports permutations of up to 62 objects and labels these by
    `0 \dots 9`, `a \dots z`, and `A \dots Z`.

    INPUT:

    - ``permutation`` -- string specifying a PALP style permutation

    OUTPUT: a :class:`permutation group element <sage.groups.perm_gps.permgroup_element.PermutationGroupElement>`

    EXAMPLES::

        sage: from sage.geometry.lattice_polytope import _palp_convert_permutation
        sage: _palp_convert_permutation('1023')                                         # needs sage.groups
        (1,2)
        sage: _palp_convert_permutation('0123456789bac')                                # needs sage.groups
        (11,12)
    """
    def from_palp_index(i):
        if i.isdigit():
            i = int(i)
            i += 1
        else:
            o = ord(i)
            if o in range(97, 123):
                i = o - 86
            elif o in range(65, 91):
                i = o - 28
            else:
                raise ValueError('cannot convert PALP index '
                                 + i + ' to number')
        return i
    n = len(permutation)
    domain = [from_palp_index(i) for i in permutation]
    from sage.groups.perm_gps.permgroup_element import make_permgroup_element
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    S = SymmetricGroup(n)
    return make_permgroup_element(S, domain)


def _read_nef_x_partitions(data):
    r"""
    Read all nef-partitions for one polytope from a string or an open
    file.

    INPUT:

    - ``data`` -- should be an output of ``nef.x``

    Returns the sequence of nef-partitions. Each nef-partition is given
    as a sequence of integers.

    If there are no nef-partitions, returns the empty sequence. If the
    string is empty or EOF is reached, raises :exc:`ValueError`.

    TESTS::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: s = o.nef_x("-N -p")                                                      # needs palp
        sage: print(s)  # random                                                        # needs palp
        M:27 8 N:7 6  codim=2 #part=5
         P:0 V:2 4 5       0sec  0cpu
         P:2 V:3 4 5       0sec  0cpu
         P:3 V:4 5       0sec  0cpu
        np=3 d:1 p:1    0sec     0cpu
        sage: lattice_polytope._read_nef_x_partitions(s)                                # needs palp
        [[2, 4, 5], [3, 4, 5], [4, 5]]
    """
    if isinstance(data, str):
        f = StringIO(data)
        partitions = _read_nef_x_partitions(f)
        f.close()
        return partitions
    line = data.readline()
    if line == "":
        raise ValueError("Empty file!")
    partitions = []
    while line and line.find("np=") == -1:
        if line.find("V:") == -1:
            line = data.readline()
            continue
        start = line.find("V:") + 2
        end = line.find("  ", start)  # Find DOUBLE space
        partitions.append(Sequence(line[start:end].split(), int))
        line = data.readline()
    # Compare the number of found partitions with np in data.
    start = line.find("np=")
    if start != -1:
        # start += 3
        # end = line.find(" ", start)
        # np = int(line[start:end])
        # if False and np != len(partitions):
        #     raise ValueError("Found %d partitions, expected %d!" %
        #                          (len(partitions), np))
        pass
    else:
        raise ValueError("Wrong data format, cannot find \"np=\"!")
    return partitions


def _read_poly_x_incidences(data, dim):
    r"""
    Convert incidence data from binary numbers to sequences.

    INPUT:

    - ``data`` -- an opened file with incidence
      information. The first line will be skipped, each consecutive line
      contains incidence information for all faces of one dimension, the
      first word of each line is a comment and is dropped.

    - ``dim`` -- dimension of the polytope

    OUTPUT: a sequence F, such that F[d][i] is a sequence of vertices
    or facets corresponding to the `i`-th d-dimensional face

    TESTS::

        sage: # needs palp
        sage: p = lattice_polytope.cross_polytope(2)
        sage: result_name = lattice_polytope._palp("poly.x -fi", [p])
        sage: with open(result_name) as f:
        ....:     print(f.read())
        Incidences as binary numbers [F-vector=(4 4)]:
        v[d][i]: sum_j Incidence(i'th dim-d-face, j-th vertex) x 2^j
        v[0]: 1000 0001 0100 0010
        v[1]: 1001 1100 0011 0110
        f[d][i]: sum_j Incidence(i'th dim-d-face, j-th facet) x 2^j
        f[0]: 0011 0101 1010 1100
        f[1]: 0001 0010 0100 1000
        sage: f = open(result_name)
        sage: l = f.readline()
        sage: lattice_polytope._read_poly_x_incidences(f, 2)
        [[[3], [0], [2], [1]], [[0, 3], [2, 3], [0, 1], [1, 2]]]
        sage: f.close()
        sage: os.remove(result_name)
    """
    data.readline()
    lines = [data.readline().split() for i in range(dim)]
    if len(lines) != dim:
        raise ValueError("Not enough data!")
    n = len(lines[0][1])     # Number of vertices or facets
    result = []
    for line in lines:
        line.pop(0)
        subr = []
        for e in line:
            f = Sequence([j for j in range(n) if e[n-1-j] == '1'],
                         int, check=False)
            f.set_immutable()
            subr.append(f)
        result.append(subr)
    return result


def all_cached_data(polytopes):
    r"""
    Compute all cached data for all given ``polytopes`` and
    their polars.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data. None of the
    polytopes in the given sequence should be constructed as the polar
    polytope to another one.

    INPUT:

    - ``polytopes`` -- a sequence of lattice polytopes

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: lattice_polytope.all_cached_data([o])                                     # needs palp
    """
    all_polars(polytopes)
    all_points(polytopes)
    reflexive = [p for p in polytopes if p.is_reflexive()]
    all_nef_partitions(reflexive)
    polar = [p.polar() for p in reflexive]
    all_points(polar)
    all_nef_partitions(polar)


def all_nef_partitions(polytopes, keep_symmetric=False):
    r"""
    Compute nef-partitions for all given ``polytopes``.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    Note: member function ``is_reflexive`` will be called
    separately for each polytope. It is strictly recommended to call
    ``all_polars`` on the sequence of
    ``polytopes`` before using this function.

    INPUT:

    - ``polytopes`` -- a sequence of lattice polytopes

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: lattice_polytope.all_nef_partitions([o])                                  # needs palp
        sage: o.nef_partitions()                                                        # needs palp
        [Nef-partition {0, 1, 3} ⊔ {2, 4, 5},
         Nef-partition {0, 1, 3, 4} ⊔ {2, 5} (direct product),
         Nef-partition {0, 1, 2} ⊔ {3, 4, 5},
         Nef-partition {0, 1, 2, 3} ⊔ {4, 5},
         Nef-partition {0, 1, 2, 3, 4} ⊔ {5} (projection)]

    You cannot use this function for non-reflexive polytopes::

        sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,2),
        ....:                      (-1,0,0), (0,-1,0), (0,0,-1)])
        sage: lattice_polytope.all_nef_partitions([o, p])                               # needs palp
        Traceback (most recent call last):
        ...
        ValueError: nef-partitions can be computed for reflexive polytopes only
    """
    keys = "-N -V -D -P -p"
    if keep_symmetric:
        keys += " -s"
    result_name = _palp("nef.x -f " + keys, polytopes)
    result = open(result_name)
    for p in polytopes:
        if not p.is_reflexive():
            raise ValueError("nef-partitions can be computed for reflexive "
                             "polytopes only")
        p._read_nef_partitions(result)
        p._nef_partitions_s = keep_symmetric
    result.close()
    os.remove(result_name)


def all_points(polytopes):
    r"""
    Compute lattice points for all given ``polytopes``.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    INPUT:

    - ``polytopes`` -- a sequence of lattice polytopes

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: lattice_polytope.all_points([o])                                          # needs palp
        sage: o.points()                                                                # needs palp
        M( 1,  0,  0),
        M( 0,  1,  0),
        M( 0,  0,  1),
        M(-1,  0,  0),
        M( 0, -1,  0),
        M( 0,  0, -1),
        M( 0,  0,  0)
        in 3-d lattice M
    """
    result_name = _palp("poly.x -fp", polytopes, reduce_dimension=True)
    result = open(result_name)
    for p in polytopes:
        M = p.lattice()
        nv = p.nvertices()
        if p.dim() == p.lattice_dim():
            points = read_palp_point_collection(result, M)
            p._points = points if len(points) > nv else p.vertices()
        else:
            m = p._embed(read_palp_matrix(result))
            if m.ncols() == nv:
                p._points = p.vertices()
            else:
                points = list(p.vertices())
                for j in range(nv, m.ncols()):
                    current = M.element_class(
                        M, [m[i, j] for i in range(M.rank())])
                    current.set_immutable()
                    points.append(current)
                p._points = PointCollection(points, M)
    result.close()
    os.remove(result_name)


def all_polars(polytopes):
    r"""
    Compute polar polytopes for all reflexive and equations of facets
    for all non-reflexive ``polytopes``.

    ``all_facet_equations`` and ``all_polars`` are synonyms.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    INPUT:

    - ``polytopes`` -- a sequence of lattice polytopes

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: lattice_polytope.all_polars([o])                                          # needs palp
        sage: o.polar()                                                                 # needs palp
        3-d reflexive polytope in 3-d lattice N
    """
    result_name = _palp("poly.x -fe", polytopes)
    result = open(result_name)
    for p in polytopes:
        p._read_equations(result)
    result.close()
    os.remove(result_name)


# Synonym for the above function
all_facet_equations = all_polars


def convex_hull(points):
    r"""
    Compute the convex hull of the given points.

    .. NOTE::

       ``points`` might not span the space. Also, it fails for large
       numbers of vertices in dimensions 4 or greater

    INPUT:

    - ``points`` -- list that can be converted into
      vectors of the same dimension over `\ZZ`

    OUTPUT: list of vertices of the convex hull of the given points (as vectors)

    EXAMPLES: Let's compute the convex hull of several points on a line
    in the plane::

        sage: lattice_polytope.convex_hull([[1,2],[3,4],[5,6],[7,8]])
        [(1, 2), (7, 8)]
    """
    if len(points) == 0:
        return []
    vpoints = []
    for p in points:
        v = vector(ZZ, p)
        if v not in vpoints:
            vpoints.append(v)
    p0 = vpoints[0]
    vpoints = [p - p0 for p in vpoints]
    N = ZZ**p0.degree()
    H = N.submodule(vpoints)
    if H.rank() == 0:
        return [p0]
    elif H.rank() == N.rank():
        vpoints = list(LatticePolytope(vpoints, lattice=N).vertices())
    else:
        H_points = [H.coordinates(p) for p in vpoints]
        H_polytope = LatticePolytope(H_points)
        vpoints = (H_polytope.vertices() * H.basis_matrix()).rows(copy=False)
    vpoints = [p + p0 for p in vpoints]
    return vpoints


def cross_polytope(dim):
    r"""
    Return a cross-polytope of the given dimension.

    INPUT:

    - ``dim`` -- integer

    OUTPUT: a :class:`lattice polytope <LatticePolytopeClass>`

    EXAMPLES::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: o
        3-d reflexive polytope in 3-d lattice M
        sage: o.vertices()
        M( 1,  0,  0),
        M( 0,  1,  0),
        M( 0,  0,  1),
        M(-1,  0,  0),
        M( 0, -1,  0),
        M( 0,  0, -1)
        in 3-d lattice M
    """
    M = ZZ**dim
    vertices = list(M.gens())
    vertices += [-v for v in vertices]
    return LatticePolytope(vertices, compute_vertices=False)


def minkowski_sum(points1, points2):
    r"""
    Compute the Minkowski sum of two convex polytopes.

    .. NOTE::

       Polytopes might not be of maximal dimension.

    INPUT:

    - ``points1``, ``points2`` -- lists of objects that can be
      converted into vectors of the same dimension, treated as vertices
      of two polytopes.

    OUTPUT: list of vertices of the Minkowski sum, given as vectors

    EXAMPLES: Let's compute the Minkowski sum of two line segments::

        sage: lattice_polytope.minkowski_sum([[1,0],[-1,0]],[[0,1],[0,-1]])
        [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    """
    points1 = [vector(p) for p in points1]
    points2 = [vector(p) for p in points2]
    points = []
    for p1 in points1:
        for p2 in points2:
            points.append(p1+p2)
    return convex_hull(points)


def positive_integer_relations(points):
    r"""
    Return relations between given points.

    INPUT:

    - ``points`` -- lattice points given as columns of a matrix

    OUTPUT: matrix of relations between given points with nonnegative
    integer coefficients

    EXAMPLES: This is a 3-dimensional reflexive polytope::

        sage: p = LatticePolytope([(1,0,0), (0,1,0),
        ....:                      (-1,-1,0), (0,0,1), (-1,0,-1)])
        sage: p.points()                                                                # needs palp
        M( 1,  0,  0),
        M( 0,  1,  0),
        M(-1, -1,  0),
        M( 0,  0,  1),
        M(-1,  0, -1),
        M( 0,  0,  0)
        in 3-d lattice M

    We can compute linear relations between its points in the following
    way::

        sage: p.points().matrix().kernel().echelonized_basis_matrix()                   # needs palp
        [ 1  0  0  1  1  0]
        [ 0  1  1 -1 -1  0]
        [ 0  0  0  0  0  1]

    However, the above relations may contain negative and rational
    numbers. This function transforms them in such a way, that all
    coefficients are nonnegative integers::

        sage: points = p.points().column_matrix()
        sage: lattice_polytope.positive_integer_relations(points)                       # needs palp
        [1 0 0 1 1 0]
        [1 1 1 0 0 0]
        [0 0 0 0 0 1]

        sage: cm = ReflexivePolytope(2,1).vertices().column_matrix()
        sage: lattice_polytope.positive_integer_relations(cm)
        [2 1 1]
    """
    points = points.transpose().base_extend(QQ)
    relations = points.kernel().echelonized_basis_matrix()
    nonpivots = relations.nonpivots()
    nonpivot_relations = relations.matrix_from_columns(nonpivots)
    n_nonpivots = len(nonpivots)
    n = nonpivot_relations.nrows()
    a = matrix(QQ, n_nonpivots, n_nonpivots)
    for i in range(n_nonpivots):
        a[i, i] = -1
    a = nonpivot_relations.stack(a).transpose()
    new_relations = []
    for i in range(n_nonpivots):
        # Find a nonnegative linear combination of relations,
        # such that all components are nonnegative and the `i`-th one is 1
        MIP = MixedIntegerLinearProgram(maximization=False, base_ring=QQ)
        w = MIP.new_variable(integer=False, nonnegative=True)
        b = vector([0] * i + [1] + [0] * (n_nonpivots - i - 1))
        MIP.add_constraint(a * w == b)
        c = [0] * (n + i) + [1] + [0] * (n_nonpivots - i - 1)
        MIP.set_objective(sum(ci * w[i] for i, ci in enumerate(c)))
        MIP.solve()
        x = list(MIP.get_values(w).values())[:n]
        v = relations.linear_combination_of_rows(x)
        new_relations.append(v)

    relations = relations.stack(matrix(QQ, new_relations))
    # Use the new relation to remove negative entries in non-pivot columns
    for i in range(n_nonpivots):
        for j in range(n):
            coef = relations[j,nonpivots[i]]
            if coef < 0:
                relations.add_multiple_of_row(j, n + i, -coef)
    # Get a new basis
    relations = relations.matrix_from_rows(relations.transpose().pivots())
    # Switch to integers
    for i in range(n):
        relations.rescale_row(i, 1 / integral_length(relations[i]))
    return relations.change_ring(ZZ)


def read_all_polytopes(file_name):
    r"""
    Read all polytopes from the given file.

    INPUT:

    - ``file_name`` -- string with the name of a file with VERTICES of
      polytopes

    OUTPUT: a sequence of polytopes

    EXAMPLES:

    We use ``poly.x`` to compute two polar polytopes and read them::

        sage: # needs palp
        sage: d = lattice_polytope.cross_polytope(2)
        sage: o = lattice_polytope.cross_polytope(3)
        sage: result_name = lattice_polytope._palp("poly.x -fe", [d, o])
        sage: with open(result_name) as f:
        ....:     print(f.read())
        4 2  Vertices of P-dual <-> Equations of P
          -1   1
           1   1
          -1  -1
           1  -1
        8 3  Vertices of P-dual <-> Equations of P
          -1  -1   1
           1  -1   1
          -1   1   1
           1   1   1
          -1  -1  -1
           1  -1  -1
          -1   1  -1
           1   1  -1
        sage: lattice_polytope.read_all_polytopes(result_name)
        [2-d reflexive polytope #14 in 2-d lattice M,
         3-d reflexive polytope in 3-d lattice M]
        sage: os.remove(result_name)
    """
    polytopes = []
    with open(file_name) as f:
        pc = read_palp_point_collection(f)
        while pc is not None:
            polytopes.append(LatticePolytope(pc, compute_vertices=False))
            pc = read_palp_point_collection(f)
    return polytopes


def read_palp_matrix(data, permutation=False):
    r"""
    Read and return an integer matrix from a string or an opened file.

    First input line must start with two integers m and n, the number
    of rows and columns of the matrix. The rest of the first line is
    ignored. The next m lines must contain n numbers each.

    If m>n, returns the transposed matrix. If the string is empty or EOF
    is reached, returns the empty matrix, constructed by
    ``matrix()``.

    INPUT:

    - ``data`` -- either a string containing the filename or the file itself
      containing the output by PALP

    - ``permutation`` -- boolean (default: ``False``); if ``True``, try to retrieve
      the permutation output by PALP. This parameter makes sense only
      when PALP computed the normal form of a lattice polytope.

    OUTPUT: a matrix or a tuple of a matrix and a permutation

    EXAMPLES::

        sage: lattice_polytope.read_palp_matrix("2 3 comment \n 1 2 3 \n 4 5 6")
        [1 2 3]
        [4 5 6]
        sage: lattice_polytope.read_palp_matrix("3 2 Will be transposed \n 1 2 \n 3 4 \n 5 6")
        [1 3 5]
        [2 4 6]
    """
    if isinstance(data,str):
        f = StringIO(data)
        mat = read_palp_matrix(f, permutation=permutation)
        f.close()
        return mat
    # If data is not a string, try to treat it as a file.
    first_line = data.readline()
    if first_line == "":
        return matrix()
    first_line = first_line.split()
    m = int(first_line[0])
    n = int(first_line[1])
    seq = []
    for i in range(m):
        seq.extend(int(el) for el in data.readline().split())
    mat = matrix(ZZ,m,n,seq)
    if m > n:
        mat = mat.transpose()
    # In some cases there may be additional information to extract
    if permutation:
        last_piece = first_line[-1]
        last_piece = last_piece.split('=')
        if last_piece[0] != 'perm':
            raise ValueError('PALP did not return a permutation.')
        p = _palp_convert_permutation(last_piece[1])
        return (mat, p)
    else:
        return mat


def set_palp_dimension(d):
    r"""
    Set the dimension for PALP calls to ``d``.

    INPUT:

    - ``d`` -- integer from the list ``[4,5,6,11]`` or ``None``

    OUTPUT: none

    PALP has many hard-coded limits, which must be specified before
    compilation, one of them is dimension. Sage includes several versions with
    different dimension settings (which may also affect other limits and enable
    certain features of PALP). You can change the version which will be used by
    calling this function. Such a change is not done automatically for each
    polytope based on its dimension, since depending on what you are doing it
    may be necessary to use dimensions higher than that of the input polytope.

    EXAMPLES:

    Let's try to work with a 7-dimensional polytope::

        sage: p = lattice_polytope.cross_polytope(7)
        sage: p._palp("poly.x -fv")                                                     # needs palp
        Traceback (most recent call last):
        ...
        ValueError: Error executing 'poly.x -fv' for the given polytope!
        Output:
        Please increase POLY_Dmax to at least 7

    However, we can work with this polytope by changing PALP dimension to 11::

        sage: lattice_polytope.set_palp_dimension(11)
        sage: p._palp("poly.x -fv")                                                     # needs palp
        '7 14  Vertices of P...'

    Let's go back to default settings::

        sage: lattice_polytope.set_palp_dimension(None)
    """
    global _palp_dimension
    _palp_dimension = d


def skip_palp_matrix(data, n=1):
    r"""
    Skip matrix data in a file.

    INPUT:

    - ``data`` -- opened file with blocks of matrix data in
      the following format: A block consisting of m+1 lines has the
      number m as the first element of its first line.

    - ``n`` -- (default: 1) integer, specifies how many
      blocks should be skipped

    If EOF is reached during the process, raises :exc:`ValueError` exception.

    EXAMPLES: We create a file with vertices of the square and the cube,
    but read only the second set::

        sage: # needs palp
        sage: d = lattice_polytope.cross_polytope(2)
        sage: o = lattice_polytope.cross_polytope(3)
        sage: result_name = lattice_polytope._palp("poly.x -fe", [d, o])
        sage: with open(result_name) as f:
        ....:     print(f.read())
        4 2  Vertices of P-dual <-> Equations of P
          -1   1
           1   1
          -1  -1
           1  -1
        8 3  Vertices of P-dual <-> Equations of P
          -1  -1   1
           1  -1   1
          -1   1   1
           1   1   1
          -1  -1  -1
           1  -1  -1
          -1   1  -1
           1   1  -1
        sage: f = open(result_name)
        sage: lattice_polytope.skip_palp_matrix(f)
        sage: lattice_polytope.read_palp_matrix(f)
        [-1  1 -1  1 -1  1 -1  1]
        [-1 -1  1  1 -1 -1  1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        sage: f.close()
        sage: os.remove(result_name)
    """
    for i in range(n):
        line = data.readline()
        if line == "":
            raise ValueError("There are not enough data to skip!")
        for j in range(int(line.split()[0])):
            if data.readline() == "":
                raise ValueError("There are not enough data to skip!")


def write_palp_matrix(m, ofile=None, comment='', format=None):
    r"""
    Write ``m`` into ``ofile`` in PALP format.

    INPUT:

    - ``m`` -- a matrix over integers or a
      :class:`point collection <PointCollection>`

    - ``ofile`` -- a file opened for writing (default: ``stdout``)

    - ``comment`` -- string (default: empty); see output description

    - ``format`` -- a format string used to print matrix entries

    OUTPUT: nothing is returned, output written to ``ofile`` has the format

      * First line: number_of_rows number_of_columns comment
      * Next number_of_rows lines: rows of the matrix.

    EXAMPLES::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: lattice_polytope.write_palp_matrix(o.vertices(), comment="3D Octahedron")
        3 6 3D Octahedron
         1  0  0 -1  0  0
         0  1  0  0 -1  0
         0  0  1  0  0 -1
        sage: lattice_polytope.write_palp_matrix(o.vertices(), format='%4d')
        3 6
           1    0    0   -1    0    0
           0    1    0    0   -1    0
           0    0    1    0    0   -1
    """
    if isinstance(m, PointCollection):
        m = m.column_matrix()
    if format is None:
        n = max(len(str(m[i,j]))
                for i in range(m.nrows()) for j in range(m.ncols()))
        format = "%" + str(n) + "d"
    s = "%d %d %s\n" % (m.nrows(),m.ncols(),comment)
    if ofile is None:
        print(s, end=" ")
    else:
        ofile.write(s)
    for i in range(m.nrows()):
        s = " ".join(format % m[i,j] for j in range(m.ncols()))+"\n"
        if ofile is None:
            print(s, end=" ")
        else:
            ofile.write(s)
