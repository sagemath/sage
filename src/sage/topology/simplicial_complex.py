# sage.doctest: needs sage.graphs
r"""
Finite simplicial complexes

AUTHORS:

- John H. Palmieri (2009-04)

- D. Benjamin Antieau (2009-06): added is_connected, generated_subcomplex,
  remove_facet, and is_flag_complex methods;
  cached the output of the graph() method.

- Travis Scrimshaw (2012-08-17): Made :class:`SimplicialComplex` have an
  immutable option, and added ``__hash__()`` function which checks to make
  sure it is immutable. Made :meth:`SimplicialComplex.remove_face()` into a
  mutator. Deprecated the ``vertex_set`` parameter.

- Christian Stump (2011-06): implementation of is_cohen_macaulay

- Travis Scrimshaw (2013-02-16): Allowed :class:`SimplicialComplex` to make
  mutable copies.

- Simon King (2014-05-02): Let simplicial complexes be objects of the
  category of simplicial complexes.

- Jeremy Martin (2016-06-02): added cone_vertices, decone, is_balanced,
  is_partitionable, intersection methods

This module implements the basic structure of finite simplicial
complexes. Given a set `V` of "vertices", a simplicial complex on `V`
is a collection `K` of subsets of `V` satisfying the condition that if
`S` is one of the subsets in `K`, then so is every subset of `S`.  The
subsets `S` are called the 'simplices' of `K`.

.. NOTE::

   In Sage, the elements of the vertex set are determined
   automatically: `V` is defined to be the union of the sets in
   `K`. So in Sage's implementation of simplicial complexes, every
   vertex is included in some face.

A simplicial complex `K` can be viewed as a purely combinatorial
object, as described above, but it also gives rise to a topological
space `|K|` (its *geometric realization*) as follows: first, the
points of `V` should be in general position in euclidean space.  Next,
if `\{v\}` is in `K`, then the vertex `v` is in `|K|`.  If `\{v, w\}`
is in `K`, then the line segment from `v` to `w` is in `|K|`. If `\{u,
v, w\}` is in `K`, then the triangle with vertices `u`, `v`, and `w`
is in `|K|`.  In general, `|K|` is the union of the convex hulls of
simplices of `K`.  Frequently, one abuses notation and uses `K` to
denote both the simplicial complex and the associated topological
space.

.. image:: ../../media/simplices.png

For any simplicial complex `K` and any commutative ring `R` there is
an associated chain complex, with differential of degree `-1`.  The
`n`-th term is the free `R`-module with basis given by the
`n`-simplices of `K`.  The differential is determined by its value on
any simplex: on the `n`-simplex with vertices `(v_0, v_1, ..., v_n)`,
the differential is the alternating sum with `i`-th summand `(-1)^i`
multiplied by the `(n-1)`-simplex obtained by omitting vertex `v_i`.

In the implementation here, the vertex set must be finite. To define a
simplicial complex, specify its *facets*: the maximal subsets (with
respect to inclusion) of the vertex set belonging to `K`. Each facet
can be specified as a list, a tuple, or a set.

.. NOTE::

   This class derives from
   :class:`~sage.topology.cell_complex.GenericCellComplex`, and so
   inherits its methods.  Some of those methods are not listed here;
   see the :mod:`Generic Cell Complex <sage.topology.cell_complex>`
   page instead.

EXAMPLES::

    sage: SimplicialComplex([[1], [3, 7]])
    Simplicial complex with vertex set (1, 3, 7) and facets {(1,), (3, 7)}
    sage: SimplicialComplex()   # the empty simplicial complex
    Simplicial complex with vertex set () and facets {()}
    sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
    sage: X
    Simplicial complex with vertex set (0, 1, 2, 3) and
     facets {(0, 1), (0, 3), (1, 2), (2, 3)}

Sage can perform a number of operations on simplicial complexes, such
as the join and the product, and it can also compute homology::

    sage: S = SimplicialComplex([[0,1], [1,2], [0,2]]) # circle
    sage: T = S.product(S)  # torus
    sage: T
    Simplicial complex with 9 vertices and 18 facets
    sage: T.homology()   # this computes reduced homology                               # needs sage.modules
    {0: 0, 1: Z x Z, 2: Z}
    sage: T.euler_characteristic()
    0

Sage knows about some basic combinatorial data associated to a
simplicial complex::

    sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [0,3]])
    sage: X.f_vector()
    [1, 4, 4]
    sage: X.face_poset()
    Finite poset containing 8 elements
    sage: x0, x1, x2, x3 = X.stanley_reisner_ring().gens()                              # needs sage.libs.singular
    sage: x0*x2 == x1*x3 == 0                                                           # needs sage.libs.singular
    True
    sage: X.is_pure()
    True

Mutability (see :issue:`12587`)::

    sage: S = SimplicialComplex([[1,4], [2,4]])
    sage: S.add_face([1,3])
    sage: S.remove_face([1,3]); S
    Simplicial complex with vertex set (1, 2, 3, 4) and facets {(3,), (1, 4), (2, 4)}
    sage: hash(S)
    Traceback (most recent call last):
    ...
    ValueError: this simplicial complex must be immutable; call set_immutable()
    sage: S = SimplicialComplex([[1,4], [2,4]])
    sage: S.set_immutable()
    sage: S.add_face([1,3])
    Traceback (most recent call last):
    ...
    ValueError: this simplicial complex is not mutable
    sage: S.remove_face([1,3])
    Traceback (most recent call last):
    ...
    ValueError: this simplicial complex is not mutable
    sage: hash(S) == hash(S)
    True

    sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
    sage: hash(S2) == hash(S)
    True

We can also make mutable copies of an immutable simplicial complex
(see :issue:`14142`)::

    sage: S = SimplicialComplex([[1,4], [2,4]])
    sage: S.set_immutable()
    sage: T = copy(S)
    sage: T.is_mutable()
    True
    sage: S == T
    True
"""

# possible future directions for SimplicialComplex:
#
#  make compatible with GAP (see http://linalg.org/gap.html)
#  compare to and make compatible with polymake
#  see Macaulay: http://www.math.uiuc.edu/Macaulay2/doc/Macaulay2-1.1/share/doc/Macaulay2/SimplicialComplexes/html/___Simplicial__Complex.html; compare performance and make compatible
#  should + have any meaning?
#  cohomology: compute cup products (and Massey products?)

from copy import copy
from itertools import combinations, chain
from functools import total_ordering

from .cell_complex import GenericCellComplex
from sage.categories.fields import Fields
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.misc.superseded import deprecation
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygens
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.set import Set
from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject

lazy_import('sage.categories.simplicial_complexes', 'SimplicialComplexes')
lazy_import('sage.matrix.constructor', 'matrix')


def lattice_paths(t1, t2, length=None):
    r"""
    Given lists (or tuples or ...) ``t1`` and ``t2``, think of them as
    labelings for vertices: ``t1`` labeling points on the x-axis,
    ``t2`` labeling points on the y-axis, both increasing.  Return the
    list of rectilinear paths along the grid defined by these points
    in the plane, starting from ``(t1[0], t2[0])``, ending at
    ``(t1[last], t2[last])``, and at each grid point, going either
    right or up.  See the examples.

    INPUT:

    - ``t1`` -- list or other iterable; labeling for vertices
    - ``t2`` -- list or other iterable; labeling for vertices
    - ``length`` -- integer or ``None`` (default: ``None``); if not ``None``, then
      an integer, the length of the desired path

    OUTPUT: list of lists of vertices making up the paths

    This is used when triangulating the product of simplices.  The
    optional argument ``length`` is used for `\Delta`-complexes, to
    specify all simplices in a product: in the triangulation of a
    product of two simplices, there is a `d`-simplex for every path of
    length `d+1` in the lattice.  The path must start at the bottom
    left and end at the upper right, and it must use at least one
    point in each row and in each column, so if ``length`` is too
    small, there will be no paths.

    EXAMPLES::

        sage: from sage.topology.simplicial_complex import lattice_paths
        sage: lattice_paths([0,1,2], [0,1,2])
        [[(0, 0), (0, 1), (0, 2), (1, 2), (2, 2)],
         [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)],
         [(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)],
         [(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)],
         [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2)],
         [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2)]]
        sage: lattice_paths(('a', 'b', 'c'), (0, 3, 5))
        [[('a', 0), ('a', 3), ('a', 5), ('b', 5), ('c', 5)],
         [('a', 0), ('a', 3), ('b', 3), ('b', 5), ('c', 5)],
         [('a', 0), ('b', 0), ('b', 3), ('b', 5), ('c', 5)],
         [('a', 0), ('a', 3), ('b', 3), ('c', 3), ('c', 5)],
         [('a', 0), ('b', 0), ('b', 3), ('c', 3), ('c', 5)],
         [('a', 0), ('b', 0), ('c', 0), ('c', 3), ('c', 5)]]
        sage: lattice_paths(range(3), range(3), length=2)
        []
        sage: lattice_paths(range(3), range(3), length=3)
        [[(0, 0), (1, 1), (2, 2)]]
        sage: lattice_paths(range(3), range(3), length=4)
        [[(0, 0), (1, 1), (1, 2), (2, 2)],
         [(0, 0), (0, 1), (1, 2), (2, 2)],
         [(0, 0), (1, 1), (2, 1), (2, 2)],
         [(0, 0), (1, 0), (2, 1), (2, 2)],
         [(0, 0), (0, 1), (1, 1), (2, 2)],
         [(0, 0), (1, 0), (1, 1), (2, 2)]]
    """
    # Convert t1, t2 to tuples, in case they are (for example) Python 3 ranges.
    t1 = tuple(t1)
    t2 = tuple(t2)
    if length is None:
        # 0 x n (or k x 0) rectangle:
        if len(t1) == 0 or len(t2) == 0:
            return [[]]
        # 1 x n (or k x 1) rectangle:
        elif len(t1) == 1:
            return [[(t1[0], w) for w in t2]]
        elif len(t2) == 1:
            return [[(v, t2[0]) for v in t1]]
        else:
            # recursive: paths in rectangle with either one fewer row
            # or column, plus the upper right corner
            return ([path + [(t1[-1], t2[-1])] for path
                     in lattice_paths(t1[:-1], t2)] +
                    [path + [(t1[-1], t2[-1])] for path
                     in lattice_paths(t1, t2[:-1])])
    else:
        if length > len(t1) + len(t2) - 1:
            return []
        # as above, except make sure that lengths are correct.  if
        # not, return an empty list.
        #
        # 0 x n (or k x 0) rectangle:
        elif len(t1) == 0 or len(t2) == 0:
            if length == 0:
                return [[]]
            else:
                return []
        # 1 x n (or k x 1) rectangle:
        elif len(t1) == 1:
            if length == len(t2):
                return [[(t1[0], w) for w in t2]]
            else:
                return []
        elif len(t2) == 1:
            if length == len(t1):
                return [[(v, t2[0]) for v in t1]]
            else:
                return []
        else:
            # recursive: paths of length one fewer in rectangle with
            # either one fewer row, one fewer column, or one fewer of
            # each, and then plus the upper right corner
            return ([path + [(t1[-1], t2[-1])] for path
                     in lattice_paths(t1[:-1], t2, length=length-1)] +
                    [path + [(t1[-1], t2[-1])] for path
                     in lattice_paths(t1, t2[:-1], length=length-1)] +
                    [path + [(t1[-1], t2[-1])] for path
                     in lattice_paths(t1[:-1], t2[:-1], length=length-1)])


def rename_vertex(n, keep, left=True):
    """
    Rename a vertex: the vertices from the list ``keep`` get
    relabeled 0, 1, 2, ..., in order.  Any other vertex (e.g. 4) gets
    renamed to by prepending an 'L' or an 'R' (thus to either 'L4' or
    'R4'), depending on whether the argument left is ``True`` or ``False``.

    INPUT:

    - ``n`` -- a 'vertex'; either an integer or a string
    - ``keep`` -- list of three vertices
    - ``left`` -- boolean (default: ``True``); if ``True``, rename for use in
      left factor

    This is used by the :meth:`~SimplicialComplex.connected_sum` method for
    simplicial complexes.

    EXAMPLES::

        sage: from sage.topology.simplicial_complex import rename_vertex
        sage: rename_vertex(6, [5, 6, 7])
        1
        sage: rename_vertex(3, [5, 6, 7, 8, 9])
        'L3'
        sage: rename_vertex(3, [5, 6, 7], left=False)
        'R3'
    """
    lookup = {i: v for v, i in enumerate(keep)}
    try:
        return lookup[n]
    except KeyError:
        return ("L" + str(n)) if left else ("R" + str(n))


@total_ordering
class Simplex(SageObject):
    """
    Define a simplex.

    Topologically, a simplex is the convex hull of a collection of
    vertices in general position.  Combinatorially, it is defined just
    by specifying a set of vertices.  It is represented in Sage by the
    tuple of the vertices.

    INPUT:

    - ``X`` -- set of vertices (integer, list, or other iterable)

    OUTPUT: simplex with those vertices

    ``X`` may be a nonnegative integer `n`, in which case the
    simplicial complex will have `n+1` vertices `(0, 1, ..., n)`, or
    it may be anything which may be converted to a tuple, in which
    case the vertices will be that tuple.  In the second case, each
    vertex must be hashable, so it should be a number, a string, or a
    tuple, for instance, but not a list.

    .. WARNING::

       The vertices should be distinct, and no error checking is done
       to make sure this is the case.

    EXAMPLES::

        sage: Simplex(4)
        (0, 1, 2, 3, 4)
        sage: Simplex([3, 4, 1])
        (3, 4, 1)
        sage: X = Simplex((3, 'a', 'vertex')); X
        (3, 'a', 'vertex')
        sage: X == loads(dumps(X))
        True

    Vertices may be tuples but not lists::

        sage: Simplex([(1,2), (3,4)])
        ((1, 2), (3, 4))
        sage: Simplex([[1,2], [3,4]])
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'list'
    """

    def __init__(self, X):
        """
        Define a simplex.  See :class:`Simplex` for full documentation.

        EXAMPLES::

            sage: Simplex(2)
            (0, 1, 2)
            sage: Simplex(('a', 'b', 'c'))
            ('a', 'b', 'c')
            sage: Simplex(-1)
            ()
            sage: Simplex(-3)
            Traceback (most recent call last):
            ...
            ValueError: the n-simplex is only defined if n > -2
        """
        try:
            N = int(X) + 1
            if N < 0:
                raise ValueError('the n-simplex is only defined if n > -2')
            self.__tuple = tuple(range(N))
        except TypeError:
            self.__tuple = tuple(X)
        self.__set = frozenset(self.__tuple)

    def tuple(self):
        """
        The tuple attached to this simplex.

        EXAMPLES::

            sage: Simplex(3).tuple()
            (0, 1, 2, 3)

        Although simplices are printed as if they were tuples, they
        are not the same type::

            sage: type(Simplex(3).tuple())
            <... 'tuple'>
            sage: type(Simplex(3))
            <class 'sage.topology.simplicial_complex.Simplex'>
        """
        return self.__tuple

    def set(self):
        """
        The frozenset attached to this simplex.

        EXAMPLES::

            sage: Simplex(3).set()
            frozenset({0, 1, 2, 3})
        """
        return self.__set

    def is_face(self, other):
        """
        Return ``True`` iff this simplex is a face of other.

        EXAMPLES::

            sage: Simplex(3).is_face(Simplex(5))
            True
            sage: Simplex(5).is_face(Simplex(2))
            False
            sage: Simplex(['a', 'b', 'c']).is_face(Simplex(8))
            False
        """
        return self.__set.issubset(other.__set)

    def __contains__(self, x):
        """
        Return ``True`` iff ``x`` is a vertex of this simplex.

        EXAMPLES::

            sage: 3 in Simplex(5)
            True
            sage: 3 in Simplex(2)
            False
        """
        return x in self.__set

    def __getitem__(self, n):
        """
        Return the `n`-th vertex in this simplex.

        EXAMPLES::

            sage: Simplex(5)[2]
            2
            sage: Simplex(['a', 'b', 'c'])[1]
            'b'
        """
        return self.__tuple[n]

    def __iter__(self):
        """
        Iterator for the vertices of this simplex.

        EXAMPLES::

            sage: [v**2 for v in Simplex(3)]
            [0, 1, 4, 9]
        """
        return iter(self.__tuple)

    def __add__(self, other):
        """
        Simplex obtained by concatenating the underlying tuples of the
        two arguments.

        INPUT:

        - ``other`` -- another simplex

        EXAMPLES::

            sage: Simplex((1,2,3)) + Simplex((5,6))
            (1, 2, 3, 5, 6)
        """
        return Simplex(self.__tuple + other.__tuple)

    def face(self, n):
        """
        The `n`-th face of this simplex.

        INPUT:

        - ``n`` -- integer between 0 and the dimension of this simplex

        OUTPUT: the simplex obtained by removing the `n`-th vertex from this
        simplex

        EXAMPLES::

            sage: S = Simplex(4)
            sage: S.face(0)
            (1, 2, 3, 4)
            sage: S.face(3)
            (0, 1, 2, 4)
        """
        if n >= 0 and n <= self.dimension():
            return Simplex(self.__tuple[:n] + self.__tuple[n+1:])
        else:
            raise IndexError("{} does not have an n-th face for n={}".format(self, n))

    def faces(self):
        """
        The list of faces (of codimension 1) of this simplex.

        EXAMPLES::

            sage: S = Simplex(4)
            sage: S.faces()
            [(1, 2, 3, 4), (0, 2, 3, 4), (0, 1, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)]
            sage: len(Simplex(10).faces())
            11
        """
        return [self.face(i) for i in range(self.dimension() + 1)]

    def dimension(self):
        """
        The dimension of this simplex.

        The dimension of a simplex is the number of vertices minus 1.

        EXAMPLES::

            sage: Simplex(5).dimension() == 5
            True
            sage: Simplex(5).face(1).dimension()
            4
        """
        return len(self.__tuple) - 1

    def is_empty(self):
        """
        Return ``True`` iff this simplex is the empty simplex.

        EXAMPLES::

            sage: [Simplex(n).is_empty() for n in range(-1,4)]
            [True, False, False, False, False]
        """
        return self.dimension() < 0

    def join(self, right, rename_vertices=True):
        """
        The join of this simplex with another one.

        The join of two simplices `[v_0, ..., v_k]` and `[w_0, ...,
        w_n]` is the simplex `[v_0, ..., v_k, w_0, ..., w_n]`.

        INPUT:

        - ``right`` -- the other simplex (the right-hand factor)

        - ``rename_vertices`` -- boolean (default: ``True``); if this is ``True``,
          the vertices in the join will be renamed by this formula: vertex "v"
          in the left-hand factor --> vertex "Lv" in the join, vertex "w" in
          the right-hand factor --> vertex "Rw" in the join.  If this is
          ``False``, this tries to construct the join without renaming the
          vertices; this may cause problems if the two factors have any
          vertices with names in common.

        EXAMPLES::

            sage: Simplex(2).join(Simplex(3))
            ('L0', 'L1', 'L2', 'R0', 'R1', 'R2', 'R3')
            sage: Simplex(['a', 'b']).join(Simplex(['x', 'y', 'z']))
            ('La', 'Lb', 'Rx', 'Ry', 'Rz')
            sage: Simplex(['a', 'b']).join(Simplex(['x', 'y', 'z']), rename_vertices=False)
            ('a', 'b', 'x', 'y', 'z')
        """
        if rename_vertices:
            vertex_set = (["L" + str(v) for v in self]
                          + ["R" + str(w) for w in right])
        else:
            vertex_set = self.__tuple + right.__tuple
        return Simplex(vertex_set)

    def product(self, other, rename_vertices=True):
        r"""
        The product of this simplex with another one, as a list of simplices.

        INPUT:

        - ``other`` -- the other simplex

        - ``rename_vertices`` -- boolean (default: ``True``); if this is
          ``False``, then the vertices in the product are the set of ordered
          pairs `(v,w)` where `v` is a vertex in the left-hand factor
          (``self``) and `w` is a vertex in the right-hand factor (``other``).
          If this is ``True``, then the vertices are renamed as "LvRw" (e.g.,
          the vertex (1,2) would become "L1R2").  This is useful if you want to
          define the Stanley-Reisner ring of the complex: vertex names like
          (0,1) are not suitable for that, while vertex names like "L0R1" are.

        ALGORITHM: see Hatcher, p. 277-278 [Hat2002]_ (who in turn refers to
        Eilenberg-Steenrod, p. 68): given ``S = Simplex(m)`` and
        ``T = Simplex(n)``, then `S \times T` can be
        triangulated as follows: for each path `f` from `(0,0)` to
        `(m,n)` along the integer grid in the plane, going up or right
        at each lattice point, associate an `(m+n)`-simplex with
        vertices `v_0`, `v_1`, ..., where `v_k` is the `k`-th vertex
        in the path `f`.

        Note that there are `m+n` choose `n` such paths.  Note also
        that each vertex in the product is a pair of vertices `(v,w)`
        where `v` is a vertex in the left-hand factor and `w`
        is a vertex in the right-hand factor.

        .. NOTE::

           This produces a list of simplices -- not a :class:`Simplex`, not
           a :class:`SimplicialComplex`.

        EXAMPLES::

            sage: len(Simplex(2).product(Simplex(2)))
            6
            sage: Simplex(1).product(Simplex(1))
            [('L0R0', 'L0R1', 'L1R1'), ('L0R0', 'L1R0', 'L1R1')]
            sage: Simplex(1).product(Simplex(1), rename_vertices=False)
            [((0, 0), (0, 1), (1, 1)), ((0, 0), (1, 0), (1, 1))]
        """
        if not rename_vertices:
            return [Simplex(x) for x in lattice_paths(self.tuple(), other.tuple())]

        answer = []
        for x in lattice_paths(self.tuple(), other.tuple()):
            new = tuple(["L" + str(v) + "R" + str(w) for v, w in x])
            answer.append(Simplex(new))
        return answer

    def alexander_whitney(self, dim):
        r"""
        Subdivide this simplex into a pair of simplices.

        If this simplex has vertices `v_0`, `v_1`, ..., `v_n`, then
        subdivide it into simplices `(v_0, v_1, ..., v_{dim})` and
        `(v_{dim}, v_{dim + 1}, ..., v_n)`.

        INPUT:

        - ``dim`` -- integer between 0 and one more than the
          dimension of this simplex

        OUTPUT:

        - a list containing just the triple ``(1, left, right)``,
          where ``left`` and ``right`` are the two simplices described
          above.

        This method allows one to construct a coproduct from the
        `p+q`-chains to the tensor product of the `p`-chains and the
        `q`-chains. The number 1 (a Sage integer) is the coefficient
        of ``left tensor right`` in this coproduct. (The corresponding
        formula is more complicated for the cubes that make up a
        cubical complex, and the output format is intended to be
        consistent for both cubes and simplices.)

        Calling this method ``alexander_whitney`` is an abuse of
        notation, since the actual Alexander-Whitney map goes from
        `C(X \times Y) \to C(X) \otimes C(Y)`, where `C(-)` denotes
        the chain complex of singular chains, but this subdivision of
        simplices is at the heart of it.

        EXAMPLES::

            sage: s = Simplex((0,1,3,4))
            sage: s.alexander_whitney(0)
            [(1, (0,), (0, 1, 3, 4))]
            sage: s.alexander_whitney(2)
            [(1, (0, 1, 3), (3, 4))]
        """
        return [(ZZ.one(), Simplex(self.tuple()[:dim+1]),
                 Simplex(self.tuple()[dim:]))]

    def __eq__(self, other):
        """
        Return ``True`` iff this simplex is the same as ``other``: that
        is, if the vertices of the two are the same, even with a
        different ordering

        INPUT:

        - ``other`` -- the other simplex

        EXAMPLES::

            sage: Simplex([0,1,2]) == Simplex([0,2,1])
            True
            sage: Simplex([0,1,2]) == Simplex(['a','b','c'])
            False
            sage: Simplex([1]) < Simplex([2])
            True
            sage: Simplex([1]) > Simplex([2])
            False
        """
        if not isinstance(other, Simplex):
            return False
        return set(self) == set(other)

    def __ne__(self, other):
        """
        Return ``True`` iff this simplex is not equal to ``other``.

        INPUT:

        - ``other`` -- the other simplex

        EXAMPLES::

            sage: Simplex([0,1,2]) != Simplex([0,2,1])
            False
            sage: Simplex([0,1,2]) != Simplex(['a','b','c'])
            True
        """
        return not self == other

    def __lt__(self, other):
        """
        Return ``True`` iff the sorted tuple for this simplex is less than
        that for ``other``.

        INPUT:

        - ``other`` -- the other simplex

        EXAMPLES::

            sage: Simplex([1]) < Simplex([2])
            True
            sage: Simplex([2,3]) < Simplex([1])
            False
            sage: Simplex([0,1,2]) < Simplex([0,2,1])
            False

        Test ``@total_ordering`` by testing other comparisons::

            sage: Simplex([0,1,2]) <= Simplex([0,2,1])
            True
            sage: Simplex([1]) <= Simplex([2])
            True
            sage: Simplex([2]) <= Simplex([1])
            False
            sage: Simplex([0,1,2]) > Simplex([0,2,1])
            False
            sage: Simplex([1]) > Simplex([2])
            False
            sage: Simplex([2]) > Simplex([1])
            True
            sage: Simplex([0,1,2]) > Simplex([0,2,1])
            False
            sage: Simplex([0,1,2]) >= Simplex([0,2,1])
            True
            sage: Simplex([1]) >= Simplex([2])
            False
            sage: Simplex([2]) >= Simplex([1])
            True
        """
        if not isinstance(other, Simplex):
            return False
        try:
            return sorted(self) < sorted(other)
        except TypeError:
            return sorted(map(str, self)) < sorted(map(str, other))

    def __hash__(self):
        """
        Hash value for this simplex.  This computes the hash value of
        the Python frozenset of the underlying tuple, since this is
        what's important when testing equality.

        EXAMPLES::

            sage: Simplex([1,2,0]).__hash__() == Simplex(2).__hash__()
            True
            sage: Simplex([1,2,0,1,1,2]).__hash__() == Simplex(2).__hash__()
            True
        """
        return hash(self.__set)

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: S = Simplex(5)
            sage: S._repr_()
            '(0, 1, 2, 3, 4, 5)'
        """
        return repr(self.__tuple)

    def _latex_(self):
        r"""
        LaTeX representation.

        EXAMPLES::

            sage: Simplex(3)._latex_()
            \left(0,
            1,
            2,
            3\right)
        """
        return latex(self.__tuple)


class SimplicialComplex(Parent, GenericCellComplex):
    r"""
    Define a simplicial complex.

    INPUT:

    - ``maximal_faces`` -- set of maximal faces
    - ``from_characteristic_function`` -- see below
    - ``maximality_check`` -- boolean (default: ``True``); see below
    - ``sort_facets`` -- dictionary; see below
    - ``name_check`` -- boolean (default: ``False``); see below
    - ``is_mutable`` -- boolean (default: ``True``); set to ``False`` to make
      this immutable
    - ``category`` -- the category of the simplicial complex (default: finite
      simplicial complexes)

    OUTPUT: a simplicial complex

    ``maximal_faces`` should be a list or tuple or set (indeed,
    anything which may be converted to a set) whose elements are lists
    (or tuples, etc.) of vertices.  Maximal faces are also known as
    'facets'. ``maximal_faces`` can also be a list containing a single
    nonnegative integer `n`, in which case this constructs the
    simplicial complex with a single `n`-simplex as the only facet.

    Alternatively, the maximal faces can be defined from a monotone boolean
    function on the subsets of a set `X`. While defining ``maximal_faces=None``,
    you can thus set ``from_characteristic_function=(f,X)`` where ``X`` is the
    set of points and ``f`` a boolean monotone hereditary function that accepts
    a list of elements from ``X`` as input (see
    :func:`~sage.combinat.subsets_hereditary.subsets_with_hereditary_property`
    for more information).

    If ``maximality_check`` is ``True``, check that each maximal face is,
    in fact, maximal. In this case, when producing the internal
    representation of the simplicial complex, omit those that are not.
    It is highly recommended that this be ``True``; various methods for
    this class may fail if faces which are claimed to be maximal are
    in fact not.

    ``sort_facets``: if not set to ``None``, the default, this should
    be a dictionary, used for sorting the vertices in each facet. The
    keys must be the vertices for the simplicial complex, and the
    values should be distinct sortable objects, for example
    integers. This should not need to be specified except in very
    special circumstances; currently the only use in the Sage library
    is when defining the product of a simplicial complex with itself:
    in this case, the vertices in the product must be sorted
    compatibly with the vertices in each factor so that the diagonal
    map is properly defined.

    If ``name_check`` is ``True``, check the names of the vertices to see
    if they can be easily converted to generators of a polynomial ring
    -- use this if you plan to use the Stanley-Reisner ring for the
    simplicial complex.

    EXAMPLES::

        sage: SimplicialComplex([[1,2], [1,4]])
        Simplicial complex with vertex set (1, 2, 4) and facets {(1, 2), (1, 4)}
        sage: SimplicialComplex([[0,2], [0,3], [0]])
        Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2), (0, 3)}
        sage: SimplicialComplex([[0,2], [0,3], [0]], maximality_check=False)
        Simplicial complex with vertex set (0, 2, 3) and facets {(0,), (0, 2), (0, 3)}

    Finally, if the first argument is a simplicial complex, return
    that complex.  If it is an object with a built-in conversion to
    simplicial complexes (via a ``_simplicial_`` method), then the
    resulting simplicial complex is returned::

        sage: S = SimplicialComplex([[0,2], [0,3], [0,6]])
        sage: SimplicialComplex(S) == S
        True
        sage: Tc = cubical_complexes.Torus(); Tc
        Cubical complex with 16 vertices and 64 cubes
        sage: Ts = SimplicialComplex(Tc); Ts
        Simplicial complex with 16 vertices and 32 facets
        sage: Ts.homology()                                                             # needs sage.modules
        {0: 0, 1: Z x Z, 2: Z}

    In the situation where the first argument is a simplicial complex
    or another object with a built-in conversion, most of the other
    arguments are ignored. The only exception is ``is_mutable``::

        sage: S.is_mutable()
        True
        sage: SimplicialComplex(S, is_mutable=False).is_mutable()
        False

    From a characteristic monotone boolean function, e.g. the simplicial complex
    of all subsets `S\subseteq \{0,1,2,3,4\}` such that `sum(S)\leq 4`::

        sage: SimplicialComplex(from_characteristic_function=(lambda x:sum(x)<=4, range(5)))
        Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 4), (0, 1, 2), (0, 1, 3)}

    or e.g. the simplicial complex of all 168 hyperovals of the projective plane of order 4::

        sage: l = designs.ProjectiveGeometryDesign(2, 1, GF(4,name='a'))                # needs sage.rings.finite_rings
        sage: f = lambda S: not any(len(set(S).intersection(x))>2 for x in l)
        sage: SimplicialComplex(from_characteristic_function=(f, l.ground_set()))       # needs sage.rings.finite_rings, long time
        Simplicial complex with 21 vertices and 168 facets

    TESTS:

    Check that we can make mutable copies (see :issue:`14142`)::

        sage: S = SimplicialComplex([[0,2], [0,3]], is_mutable=False)
        sage: S.is_mutable()
        False
        sage: C = copy(S)
        sage: C.is_mutable()
        True
        sage: SimplicialComplex(S, is_mutable=True).is_mutable()
        True
        sage: SimplicialComplex(S, is_immutable=False).is_mutable()
        True

    .. WARNING::

        Simplicial complexes are not proper parents as they do
        not possess element classes. In particular, parents are assumed
        to be hashable (and hence immutable) by the coercion framework.
        However this is close enough to being a parent with elements
        being the faces of ``self`` that we currently allow this abuse.
    """

    def __init__(self,
                 maximal_faces=None,
                 from_characteristic_function=None,
                 maximality_check=True,
                 sort_facets=None,
                 name_check=False,
                 is_mutable=True,
                 is_immutable=False,
                 category=None):
        """
        Define a simplicial complex.  See ``SimplicialComplex`` for more
        documentation.

        EXAMPLES::

            sage: SimplicialComplex([[0,2], [0,3], [0]])
            Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2), (0, 3)}
            sage: SimplicialComplex((('a', 'b'), ['a', 'c'], ('b', 'c'))) == SimplicialComplex((('a', 'b'), ('b', 'c'), ('a', 'c')))
            True

        TESTS::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S == S2
            True
            sage: S3 = SimplicialComplex(maximal_faces=[[1,4], [2,4]])
            sage: S == S3
            True

        Test that we have fixed a problem revealed in :issue:`20718`;
        see also :issue:`20720`::

            sage: SimplicialComplex([2])
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}

            sage: S = SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')))
            sage: S == loads(dumps(S))
            True

            sage: TestSuite(S).run()
            sage: TestSuite(S3).run()

        Test ``sort_facets``::

            sage: S = SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')), sort_facets=True)
            Traceback (most recent call last):
            ...
            TypeError: sort_facets must be a dict
            sage: S = SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')), sort_facets={'a': 1, 6: 3, 'c': 2})
            Traceback (most recent call last):
            ...
            ValueError: the set of keys of sort_facets must equal the set of vertices
            sage: S = SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')), sort_facets={'a': 1, 'b': 3, 'c': 2})
            sage: S._vertex_to_index['b']
            3
        """
        if (maximal_faces is not None and
                from_characteristic_function is not None):
            raise ValueError("maximal_faces and from_characteristic_function cannot be both defined")
        category = SimplicialComplexes().Finite().or_subcategory(category)
        Parent.__init__(self, category=category)

        C = None
        vertices = ()
        if from_characteristic_function is not None:
            from sage.combinat.subsets_hereditary import subsets_with_hereditary_property
            f, X = from_characteristic_function
            maximal_faces = subsets_with_hereditary_property(f, X)

        if maximal_faces is None:
            maximal_faces = []
        elif isinstance(maximal_faces, SimplicialComplex):
            C = maximal_faces
        else:
            try:
                C = maximal_faces._simplicial_()
            except AttributeError:
                if not isinstance(maximal_faces, (list, tuple, Simplex)):
                    # Convert it into a list (in case it is an iterable)
                    maximal_faces = list(maximal_faces)
                if len(maximal_faces) == 1 and isinstance(maximal_faces[0], (int, Integer)):
                    # list containing a single nonnegative integer n;
                    # construct the simplicial complex with a single n-simplex as the only facet.
                    vertices = tuple(range(maximal_faces[0] + 1))
                    maximal_faces = [vertices]
                elif maximal_faces:
                    vertices = tuple(set(chain.from_iterable(maximal_faces)))
        if C is not None:
            self._facets = list(C.facets())
            self._faces = copy(C._faces)
            self._gen_dict = copy(C._gen_dict)
            self._complex = copy(C._complex)
            self.__contractible = copy(C.__contractible)
            self.__enlarged = copy(C.__enlarged)
            self._graph = copy(C._graph)
            self._vertex_to_index = copy(C._vertex_to_index)
            self._is_immutable = False
            if not is_mutable or is_immutable:
                self.set_immutable()
            self._bbn = C._bbn
            self._bbn_all_computed = C._bbn_all_computed
            return

        gen_dict = {}
        for v in vertices:
            if name_check:
                try:
                    if int(v) < 0:
                        raise ValueError("the vertex %s does not have an appropriate name" % v)
                except ValueError:  # v is not an integer
                    try:
                        normalize_names(1, v)
                    except ValueError:
                        raise ValueError("the vertex %s does not have an appropriate name" % v)
            # build dictionary of generator names
            try:
                gen_dict[v] = 'x%s' % int(v)
            except Exception:
                gen_dict[v] = v
        # build set of facets
        good_faces = []
        maximal_simplices = [Simplex(f) for f in maximal_faces]

        if maximality_check:  # Sorting is useful to filter maximal faces
            maximal_simplices.sort(key=lambda x: x.dimension(), reverse=True)
        # Translate vertices to numbers, for use in sorting
        # facets. Having a consistent ordering for the vertices in
        # each facet is necessary for homology computations.
        if sort_facets:
            if not isinstance(sort_facets, dict):
                raise TypeError("sort_facets must be a dict")
            if set(sort_facets.keys()) != set(vertices):
                raise ValueError("the set of keys of sort_facets must equal the set of vertices")
            vertex_to_index = sort_facets
        else:
            vertex_to_index = {v: i for i, v in enumerate(vertices)}

        for face in maximal_simplices:
            # check whether each given face is actually maximal
            if (maximality_check and
                    any(face.is_face(other) for other in good_faces)):
                continue
            # This sorting is crucial for homology computations:
            face = Simplex(sorted(face.tuple(), key=vertex_to_index.__getitem__))
            good_faces.append(face)

        # if no maximal faces, add the empty face as a facet
        if len(maximal_simplices) == 0:
            good_faces.append(Simplex(-1))
        # now record the attributes for self
        # self._vertex_to_index: dictionary to convert vertices to integers
        self._vertex_to_index = vertex_to_index
        # self._facets: unsorted list of facets
        self._facets = good_faces
        # self._faces: dictionary of dictionaries of faces.  The main
        # dictionary is keyed by subcomplexes, and each value is a
        # dictionary keyed by dimension.  This should be empty until
        # needed -- that is, until the faces method is called
        self._faces = {}
        # self._gen_dict: dictionary of names for the polynomial
        # generators of the Stanley-Reisner ring
        self._gen_dict = gen_dict
        # self._complex: dictionary indexed by dimension d, subcomplex,
        # etc.: differential from dim d to dim d-1 in the associated
        # chain complex.  thus to get the differential in the cochain
        # complex from dim d-1 to dim d, take the transpose of this
        # one.
        self._complex = {}
        # self.__contractible: if not None, a contractible subcomplex
        # of self, as found by the _contractible_subcomplex method.
        self.__contractible = None
        # self.__enlarged: dictionary of enlarged subcomplexes,
        # indexed by subcomplexes.  For use in the _enlarge_subcomplex
        # method.
        self.__enlarged = {}
        # initialize self._graph to None.
        self._graph = None

        # Handle mutability keywords
        self._is_immutable = False
        if not is_mutable or is_immutable:
            self.set_immutable()

        # self._bbn: a dictionary indexed by base_ring, whose value is a dictionary of
        # bigraded Betti numbers, indexed by tuples (-i, 2j).
        # For use in the bigraded_betti_numbers method.
        self._bbn = {}
        # self._bbn_all_computed: a set of base rings for which we called
        # bigraded_betti_numbers(base_ring=base_ring)
        self._bbn_all_computed = set()

    def __hash__(self):
        """
        Compute the hash value of ``self``.

        If this simplicial complex is immutable, it computes the hash value
        based upon the facets. Otherwise it raises a :class`ValueError`.

        EXAMPLES::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: hash(S)
            Traceback (most recent call last):
            ...
            ValueError: this simplicial complex must be immutable; call set_immutable()
            sage: S.set_immutable()
            sage: hash(S) == hash(S)
            True
            sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S == S2
            True
            sage: hash(S) == hash(S2)
            True
        """
        if not self._is_immutable:
            raise ValueError("this simplicial complex must be immutable; call set_immutable()")
        return hash(frozenset(self._facets))

    def __eq__(self, right):
        """
        Two simplicial complexes are equal iff their vertex sets are
        equal and their sets of facets are equal.

        EXAMPLES::

            sage: SimplicialComplex([[1,2], [2,3], [4]]) == SimplicialComplex([[4], [2,3], [3], [2,1]])
            True
            sage: X = SimplicialComplex()
            sage: X.add_face([1,3])
            sage: X == SimplicialComplex([[1,3]])
            True
        """
        return isinstance(right, SimplicialComplex) and set(self._facets) == set(right._facets)

    def __ne__(self, right):
        """
        Return ``True`` if ``self`` and ``right`` are not equal.

        EXAMPLES::

            sage: SimplicialComplex([[1,2], [2,3], [4]]) != SimplicialComplex([[4], [2,3], [3], [2,1]])
            False
            sage: X = SimplicialComplex()
            sage: X.add_face([1,3])
            sage: X != SimplicialComplex([[1,3]])
            False
        """
        return not self.__eq__(right)

    def __copy__(self):
        """
        Return a mutable copy of ``self``.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,2], [0,3]], is_mutable=False)
            sage: S.is_mutable()
            False
            sage: C = copy(S)
            sage: C.is_mutable()
            True
            sage: C == S
            True
            sage: S.is_mutable()
            False
            sage: T = copy(C)
            sage: T == C
            True
        """
        return SimplicialComplex(self, is_mutable=True)

    def vertices(self):
        """
        The vertex set, as a tuple, of this simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex([[i] for i in range(16)] + [[0,1], [1,2]])
            sage: S
            Simplicial complex with 16 vertices and 15 facets
            sage: sorted(S.vertices())
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        return tuple(self._vertex_to_index)

    def _an_element_(self):
        """
        The first facet of this complex.

        EXAMPLES::

            sage: SimplicialComplex()._an_element_()
            ()
            sage: simplicial_complexes.Sphere(3)._an_element_()
            (0, 1, 2, 3)
        """
        try:
            return sorted(self.facets())[0]
        except TypeError:
            return self.facets()[0]

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is a simplex which is contained in this complex.

        EXAMPLES::

            sage: K = SimplicialComplex([(0,1,2), (0,2,3)])
            sage: Simplex((0,2)) in K
            True
            sage: Simplex((1,3)) in K
            False
            sage: 0 in K  # not a simplex
            False
        """
        if not isinstance(x, Simplex):
            return False
        dim = x.dimension()
        return dim in self.faces() and x in self.faces()[dim]

    def __call__(self, simplex):
        """
        If ``simplex`` is a simplex in this complex, return it.

        Otherwise, this raises a :exc:`ValueError`.

        EXAMPLES::

            sage: K = SimplicialComplex([(0,1,2), (0,2,3)])
            sage: K(Simplex((1,2)))
            (1, 2)
            sage: K(Simplex((0,1,3)))
            Traceback (most recent call last):
            ...
            ValueError: the simplex is not in this complex
        """
        if simplex not in self:
            raise ValueError('the simplex is not in this complex')
        return simplex

    def maximal_faces(self):
        """
        The maximal faces (a.k.a. facets) of this simplicial complex.

        This just returns the set of facets used in defining the
        simplicial complex, so if the simplicial complex was defined
        with no maximality checking, none is done here, either.

        EXAMPLES::

            sage: Y = SimplicialComplex([[0,2], [1,4]])
            sage: sorted(Y.maximal_faces())
            [(0, 2), (1, 4)]

        ``facets`` is a synonym for ``maximal_faces``::

            sage: S = SimplicialComplex([[0,1], [0,1,2]])
            sage: S.facets()
            {(0, 1, 2)}
        """
        return Set(self._facets)

    facets = maximal_faces

    def faces(self, subcomplex=None):
        """
        The faces of this simplicial complex, in the form of a
        dictionary of sets keyed by dimension.  If the optional
        argument ``subcomplex`` is present, then return only the
        faces which are *not* in the subcomplex.

        INPUT:

        - ``subcomplex`` -- a subcomplex of this simplicial complex (default:
          ``None``); return faces which are not in this subcomplex

        EXAMPLES::

            sage: Y = SimplicialComplex([[1,2], [1,4]])
            sage: Y.faces()
            {-1: {()}, 0: {(1,), (2,), (4,)}, 1: {(1, 2), (1, 4)}}
            sage: L = SimplicialComplex([[1,2]])
            sage: Y.faces(subcomplex=L)
            {-1: set(), 0: {(4,)}, 1: {(1, 4)}}
        """
        # Make the subcomplex immutable if it is not
        if subcomplex is not None and not subcomplex._is_immutable:
            subcomplex = SimplicialComplex(subcomplex._facets, maximality_check=False,
                                           is_mutable=False)

        if subcomplex not in self._faces:
            # Faces is the dictionary of faces in self but not in
            # subcomplex, indexed by dimension
            Faces = {}
            # sub_facets is the dictionary of facets in the subcomplex
            sub_facets = {}
            dimension = max([face.dimension() for face in self._facets])
            for i in range(-1, dimension + 1):
                Faces[i] = set()
                sub_facets[i] = set()
            for f in self._facets:
                dim = f.dimension()
                Faces[dim].add(f)
            if subcomplex is not None:
                for g in subcomplex._facets:
                    dim = g.dimension()
                    Faces[dim].discard(g)
                    sub_facets[dim].add(g)
            # bad_faces is the set of faces in the subcomplex in the
            # current dimension
            bad_faces = sub_facets[dimension]
            for dim in range(dimension, -1, -1):
                # bad_bdries = boundaries of bad_faces: things to be
                # discarded in dim-1
                bad_bdries = sub_facets[dim-1]
                for f in bad_faces:
                    bad_bdries.update(f.faces())
                for f in Faces[dim]:
                    Faces[dim-1].update(set(f.faces()).difference(bad_bdries))
                bad_faces = bad_bdries
            self._faces[subcomplex] = Faces
        return self._faces[subcomplex]

    def face_iterator(self, increasing=True):
        """
        An iterator for the faces in this simplicial complex.

        INPUT:

        - ``increasing`` -- boolean (default: ``True``); if ``True``, return
          faces in increasing order of dimension, thus starting with
          the empty face. Otherwise it returns faces in decreasing order of
          dimension.

        .. NOTE::

            Among the faces of a fixed dimension, there is no sorting.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: sorted(S1.face_iterator())
            [(), (0,), (0, 1), (0, 2), (1,), (1, 2), (2,)]
        """
        Fs = self.faces()
        dim_index = range(-1, self.dimension() + 1)
        if not increasing:
            dim_index = reversed(dim_index)
        for i in dim_index:
            yield from Fs[i]

    cells = faces

    n_faces = GenericCellComplex.n_cells

    def is_pure(self):
        """
        Return ``True`` iff this simplicial complex is pure.

        A simplicial complex is pure if and only if all of its maximal faces
        have the same dimension.

        .. WARNING::

           This may give the wrong answer if the simplicial complex
           was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: U = SimplicialComplex([[1,2], [1, 3, 4]])
            sage: U.is_pure()
            False
            sage: X = SimplicialComplex([[0,1], [0,2], [1,2]])
            sage: X.is_pure()
            True

        Demonstration of the warning::

            sage: S = SimplicialComplex([[0,1], [0]], maximality_check=False)
            sage: S.is_pure()
            False
        """
        dims = [face.dimension() for face in self._facets]
        return max(dims) == min(dims)

    def h_vector(self):
        r"""
        The `h`-vector of this simplicial complex.

        If the complex has dimension `d` and `(f_{-1}, f_0, f_1, ...,
        f_d)` is its `f`-vector (with `f_{-1} = 1`, representing the
        empty simplex), then the `h`-vector `(h_0, h_1, ..., h_d,
        h_{d+1})` is defined by

        .. MATH::

           \sum_{i=0}^{d+1} h_i x^{d+1-i} = \sum_{i=0}^{d+1} f_{i-1} (x-1)^{d+1-i}.

        Alternatively,

        .. MATH::

           h_j = \sum_{i=-1}^{j-1} (-1)^{j-i-1} \binom{d-i}{j-i-1} f_i.

        EXAMPLES:

        The `f`- and `h`-vectors of the boundary of an octahedron are
        computed in :wikipedia:`Simplicial_complex`::

            sage: square = SimplicialComplex([[0,1], [1,2], [2,3], [0,3]])
            sage: S0 = SimplicialComplex([[0], [1]])
            sage: octa = square.join(S0) # boundary of an octahedron
            sage: octa.f_vector()
            [1, 6, 12, 8]
            sage: octa.h_vector()
            [1, 3, 3, 1]
        """
        from sage.arith.misc import binomial
        d = self.dimension()
        f = self.f_vector()  # indexed starting at 0, since it's a Python list
        h = []
        for j in range(d + 2):
            s = 0
            for i in range(-1, j):
                s += (-1)**(j-i-1) * binomial(d-i, j-i-1) * f[i+1]
            h.append(s)
        return h

    def g_vector(self):
        r"""
        The `g`-vector of this simplicial complex.

        If the `h`-vector of the complex is `(h_0, h_1, ..., h_d,
        h_{d+1})` -- see :meth:`h_vector` -- then its `g`-vector
        `(g_0, g_1, ..., g_{[(d+1)/2]})` is defined by `g_0 = 1` and
        `g_i = h_i - h_{i-1}` for `i > 0`.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: S3 = simplicial_complexes.Sphere(3).barycentric_subdivision()
            sage: S3.f_vector()
            [1, 30, 150, 240, 120]
            sage: S3.h_vector()
            [1, 26, 66, 26, 1]
            sage: S3.g_vector()
            [1, 25, 40]
        """
        d = self.dimension()
        h = self.h_vector()
        g = [1]
        for i in range(1, (d + 1) // 2 + 1):
            g.append(h[i] - h[i-1])
        return g

    def face(self, simplex, i):
        """
        The `i`-th face of ``simplex`` in this simplicial complex.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial complex
        - ``i`` -- integer

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1,4], [0,1,2]])
            sage: S.face(Simplex((0,2)), 0)
            (2,)

            sage: S.face(Simplex((0,3)), 0)
            Traceback (most recent call last):
            ...
            ValueError: this simplex is not in this simplicial complex
        """
        d = simplex.dimension()
        if d in self.faces() and simplex in self.faces()[d]:
            return simplex.face(i)
        else:
            raise ValueError('this simplex is not in this simplicial complex')

    def f_triangle(self):
        r"""
        Compute the `f`-triangle of ``self``.

        The `f`-triangle is given by `f_{i,j}` being the number of
        faces `F` of size `j` such that `i = \max_{G \subseteq F} |G|`.

        .. SEEALSO::

            Not to be confused with :meth:`F_triangle` .

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2,3], [3,4,5], [1,4], [1,5], [2,4], [2,5]])
            sage: X.f_triangle()   # this complex is not pure
            [[0],
             [0, 0],
             [0, 0, 4],
             [1, 5, 6, 2]]

        A complex is pure if and only if the last row is nonzero::

            sage: X = SimplicialComplex([[1,2,3], [3,4,5], [1,4,5]])
            sage: X.f_triangle()
            [[0], [0, 0], [0, 0, 0], [1, 5, 8, 3]]
        """
        ret = [[0]*(i+1) for i in range(self.dimension() + 2)]
        facets = [set(F) for F in self.facets()]
        faces = self.faces()
        for d in faces:
            for f in faces[d]:
                f = set(f)
                L = [len(F) for F in facets if f.issubset(F)]
                i = max(L)
                ret[i][len(f)] += 1
        return ret

    def h_triangle(self):
        r"""
        Compute the `h`-triangle of ``self``.

        The `h`-triangle of a simplicial complex `\Delta` is given by

        .. MATH::

            h_{i,j} = \sum_{k=0}^j (-1)^{j-k} \binom{i-k}{j-k} f_{i,k},

        where `f_{i,k}` is the `f`-triangle of `\Delta`.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2,3], [3,4,5], [1,4], [1,5], [2,4], [2,5]])
            sage: X.h_triangle()
            [[0],
             [0, 0],
             [0, 0, 4],
             [1, 2, -1, 0]]
        """
        from sage.arith.misc import binomial
        ret = [[0]*(i+1) for i in range(self.dimension() + 2)]
        f = self.f_triangle()
        for i, row in enumerate(ret):
            for j in range(i+1):
                row[j] = sum((-1)**(j-k) * binomial(i-k, j-k) * f[i][k]
                             for k in range(j+1))
        return ret

    def F_triangle(self, S):
        """
        Return the F-triangle of ``self`` with respect
        to one maximal simplex ``S``.

        This is the bivariate generating polynomial of all faces,
        according to the number of elements in ``S`` and outside ``S``.

        OUTPUT: an :class:`~sage.combinat.triangles_FHM.F_triangle`

        .. SEEALSO::

            Not to be confused with :meth:`f_triangle` .

        EXAMPLES::

            sage: cs = simplicial_complexes.Torus()
            sage: cs.F_triangle(cs.facets()[0])                                         # needs sage.combinat
            F: x^3 + 9*x^2*y + 3*x*y^2 + y^3 + 6*x^2 + 12*x*y
            + 3*y^2 + 4*x + 3*y + 1

        TESTS::

            sage: S = SimplicialComplex([])
            sage: S.F_triangle(S.facets()[0])                                           # needs sage.combinat
            F: 1
        """
        x, y = polygens(ZZ, 'x, y')
        from sage.combinat.triangles_FHM import F_triangle

        def nega(f):
            return sum(1 for v in f if v in S)

        def posi(f):
            return f.dimension() + 1 - nega(f)

        poly = sum(x**posi(fa) * y**nega(fa)
                   for fa in self.face_iterator())
        return F_triangle(poly)

    def flip_graph(self):
        """
        If ``self`` is pure, return the flip graph of ``self``,
        otherwise, return ``None``.

        The flip graph of a pure simplicial complex is the (undirected) graph
        with vertices being the facets, such that two facets are joined by
        an edge if they meet in a codimension `1` face.

        The flip graph is used to detect if ``self`` is a pseudomanifold.

        EXAMPLES::

            sage: S0 = simplicial_complexes.Sphere(0)
            sage: G = S0.flip_graph()
            sage: G.vertices(sort=True); G.edges(sort=True, labels=False)
            [(0,), (1,)]
            [((0,), (1,))]

            sage: G = (S0.wedge(S0)).flip_graph()
            sage: G.vertices(sort=True); G.edges(sort=True, labels=False)
            [(0,), ('L1',), ('R1',)]
            [((0,), ('L1',)), ((0,), ('R1',)), (('L1',), ('R1',))]

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: G = (S1.wedge(S1)).flip_graph()
            sage: len(G.vertices(sort=False))
            6
            sage: len(G.edges(sort=False))
            10

            sage: (S1.wedge(S2)).flip_graph() is None
            True

            sage: G = S2.flip_graph()
            sage: G.vertices(sort=True); G.edges(sort=True, labels=False)
            [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
            [((0, 1, 2), (0, 1, 3)),
             ((0, 1, 2), (0, 2, 3)),
             ((0, 1, 2), (1, 2, 3)),
             ((0, 1, 3), (0, 2, 3)),
             ((0, 1, 3), (1, 2, 3)),
             ((0, 2, 3), (1, 2, 3))]

            sage: T = simplicial_complexes.Torus()
            sage: G = T.suspension(4).flip_graph()
            sage: len(G.vertices(sort=False)); len(G.edges(sort=False, labels=False))
            46
            161
        """
        from collections import defaultdict
        from sage.graphs.graph import Graph

        if not self.is_pure():
            return None
        d = self.dimension()
        Fs = self.facets()
        flipG = Graph()
        flipG.add_vertices(Fs)
        edges = defaultdict(list)
        # go through all codim 1 faces to build the edge
        for F in Fs:
            try:
                F_tuple = sorted(F._Simplex__set)
            except TypeError:
                F_tuple = tuple(F._Simplex__set)
            for i in range(d+1):
                coF = tuple(F_tuple[:i]+F_tuple[i+1:])
                if coF in edges:
                    for G in edges[coF]:
                        flipG.add_edge((F, G))
                edges[coF].append(F)
        return flipG

    def is_pseudomanifold(self):
        """
        Return ``True`` if ``self`` is a pseudomanifold.

        A pseudomanifold is a simplicial complex with the following properties:

        - it is pure of some dimension `d` (all of its facets are `d`-dimensional)
        - every `(d-1)`-dimensional simplex is the face of exactly two facets
        - for every two facets `S` and `T`, there is a sequence of
          facets

          .. MATH::

            S = f_0, f_1, ..., f_n = T

          such that for each `i`, `f_i` and `f_{i-1}` intersect in a
          `(d-1)`-simplex.

        By convention, `S^0` is the only 0-dimensional pseudomanifold.

        EXAMPLES::

            sage: S0 = simplicial_complexes.Sphere(0)
            sage: S0.is_pseudomanifold()
            True
            sage: (S0.wedge(S0)).is_pseudomanifold()
            False
            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: (S1.wedge(S1)).is_pseudomanifold()
            False
            sage: (S1.wedge(S2)).is_pseudomanifold()
            False
            sage: S2.is_pseudomanifold()
            True
            sage: T = simplicial_complexes.Torus()
            sage: T.suspension(4).is_pseudomanifold()
            True
        """
        if not self.is_pure():
            return False
        d = self.dimension()
        if d == 0:
            return len(self.facets()) == 2
        F = self.facets()
        X = self.faces()[d-1]
        # is each (d-1)-simplex is the face of exactly two facets?
        for s in X:
            if len([a for a in [s.is_face(f) for f in F] if a]) != 2:
                return False
        # construct a graph with one vertex for each facet, one edge
        # when two facets intersect in a (d-1)-simplex, and see
        # whether that graph is connected.
        return self.flip_graph().is_connected()

    def product(self, right, rename_vertices=True, is_mutable=True):
        """
        The product of this simplicial complex with another one.

        INPUT:

        - ``right`` -- the other simplicial complex (the right-hand factor)

        - ``rename_vertices`` -- boolean (default: ``True``); if this is
          ``False``, then the vertices in the product are the set of ordered
          pairs `(v,w)` where `v` is a vertex in ``self`` and `w` is a vertex
          in ``right``. If this is ``True``, then the vertices are renamed as
          "LvRw" (e.g., the vertex (1,2) would become "L1R2"). This is useful
          if you want to define the Stanley-Reisner ring of the complex: vertex
          names like (0,1) are not suitable for that, while vertex names like
          "L0R1" are.

        - ``is_mutable`` -- boolean (default: ``True``); determines whether the
          output is mutable

        The vertices in the product will be the set of ordered pairs
        `(v,w)` where `v` is a vertex in ``self`` and `w` is a vertex in
        right.

        .. WARNING::

           If ``X`` and ``Y`` are simplicial complexes, then ``X*Y``
           returns their join, not their product.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1], [1,2], [0,2]]) # circle
            sage: K = SimplicialComplex([[0,1]])   # edge
            sage: Cyl = S.product(K)  # cylinder
            sage: sorted(Cyl.vertices())
            ['L0R0', 'L0R1', 'L1R0', 'L1R1', 'L2R0', 'L2R1']
            sage: Cyl2 = S.product(K, rename_vertices=False)
            sage: sorted(Cyl2.vertices())
            [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
            sage: T = S.product(S)  # torus
            sage: T
            Simplicial complex with 9 vertices and 18 facets
            sage: T.homology()                                                          # needs sage.modules
            {0: 0, 1: Z x Z, 2: Z}

        These can get large pretty quickly::

            sage: T = simplicial_complexes.Torus(); T
            Minimal triangulation of the torus
            sage: K = simplicial_complexes.KleinBottle(); K
            Minimal triangulation of the Klein bottle
            sage: T.product(K)      # long time: 5 or 6 seconds
            Simplicial complex with 56 vertices and 1344 facets
        """
        facets = []
        for f in self._facets:
            for g in right._facets:
                facets.extend(f.product(g, rename_vertices))
        if self != right:
            return SimplicialComplex(facets, is_mutable=is_mutable)
        else:
            # Need to sort the vertices compatibly with the sorting in
            # self, so that the diagonal map is defined properly.
            V = self._vertex_to_index
            L = len(V)
            d = {}
            for v in V.keys():
                for w in V.keys():
                    if rename_vertices:
                        d['L' + str(v) + 'R' + str(w)] = V[v] * L + V[w]
                    else:
                        d[(v, w)] = V[v] * L + V[w]
            return SimplicialComplex(facets, is_mutable=is_mutable, sort_facets=d)

    def join(self, right, rename_vertices=True, is_mutable=True):
        """
        The join of this simplicial complex with another one.

        The join of two simplicial complexes `S` and `T` is the
        simplicial complex `S*T` with simplices of the form `[v_0,
        ..., v_k, w_0, ..., w_n]` for all simplices `[v_0, ..., v_k]` in
        `S` and `[w_0, ..., w_n]` in `T`.

        INPUT:

        - ``right`` -- the other simplicial complex (the right-hand factor)

        - ``rename_vertices`` -- boolean (default: ``True``); if this is
          ``True``, the vertices in the join will be renamed by the formula:
          vertex "v" in the left-hand factor --> vertex "Lv" in the join,
          vertex "w" in the right-hand factor --> vertex "Rw" in the join.
          If this is ``False``, this tries to construct the join without
          renaming the vertices; this will cause problems if the two factors
          have any vertices with names in common.

        - ``is_mutable`` -- boolean (default: ``True``); determine whether the
          output is mutable

        EXAMPLES::

            sage: S = SimplicialComplex([[0], [1]])
            sage: T = SimplicialComplex([[2], [3]])
            sage: S.join(T)
            Simplicial complex with vertex set ('L0', 'L1', 'R2', 'R3') and 4 facets
            sage: S.join(T, rename_vertices=False)
            Simplicial complex with vertex set (0, 1, 2, 3)
            and facets {(0, 2), (0, 3), (1, 2), (1, 3)}

        The notation '*' may be used, as well::

            sage: S * S
            Simplicial complex with vertex set ('L0', 'L1', 'R0', 'R1') and 4 facets
            sage: S * S * S * S * S * S * S * S
            Simplicial complex with 16 vertices and 256 facets
        """
        facets = []
        for f in self._facets:
            for g in right._facets:
                facets.append(f.join(g, rename_vertices))
        return SimplicialComplex(facets, is_mutable=is_mutable)

    # Use * to mean 'join':
    __mul__ = join

    def cone(self, is_mutable=True):
        """
        The cone on this simplicial complex.

        INPUT:

        - ``is_mutable`` -- boolean (default: ``True``); determines whether
          the output is mutable

        The cone is the simplicial complex formed by adding a new
        vertex `C` and simplices of the form `[C, v_0, ..., v_k]` for
        every simplex `[v_0, ..., v_k]` in the original simplicial
        complex.  That is, the cone is the join of the original
        complex with a one-point simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex([[0], [1]])
            sage: CS = S.cone()
            sage: sorted(CS.vertices())
            ['L0', 'L1', 'R0']
            sage: len(CS.facets())
            2
            sage: CS.facets() == set([Simplex(['L0', 'R0']), Simplex(['L1', 'R0'])])
            True
        """
        return self.join(SimplicialComplex([["0"]], is_mutable=is_mutable),
                         rename_vertices=True)

    def suspension(self, n=1, is_mutable=True):
        r"""
        The suspension of this simplicial complex.

        INPUT:

        - ``n`` -- positive integer (default: 1); suspend this many times

        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        The suspension is the simplicial complex formed by adding two
        new vertices `S_0` and `S_1` and simplices of the form `[S_0,
        v_0, ..., v_k]` and `[S_1, v_0, ..., v_k]` for every simplex
        `[v_0, ..., v_k]` in the original simplicial complex.  That
        is, the suspension is the join of the original complex with a
        two-point simplicial complex.

        If the simplicial complex `M` happens to be a pseudomanifold
        (see :meth:`is_pseudomanifold`), then this instead constructs
        Datta's one-point suspension (see [Dat2007]_, p. 434):
        choose a vertex `u` in `M` and choose a new vertex
        `w` to add.  Denote the join of simplices by "`*`".  The
        facets in the one-point suspension are of the two forms

        - `u * \alpha` where `\alpha` is a facet of `M` not containing
          `u`

        - `w * \beta` where `\beta` is any facet of `M`.

        EXAMPLES::

            sage: S0 = SimplicialComplex([[0], [1]])
            sage: S0.suspension() == simplicial_complexes.Sphere(1)
            True
            sage: S3 = S0.suspension(3)  # the 3-sphere
            sage: S3.homology()                                                         # needs sage.modules
            {0: 0, 1: 0, 2: 0, 3: Z}

        For pseudomanifolds, the complex constructed here will be
        smaller than that obtained by taking the join with the
        0-sphere: the join adds two vertices, while this construction
        only adds one. ::

            sage: T = simplicial_complexes.Torus()
            sage: sorted(T.join(S0).vertices())      # 9 vertices
            ['L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'R0', 'R1']
            sage: T.suspension().vertices()  # 8 vertices
            (0, 1, 2, 3, 4, 5, 6, 7)
        """
        if n < 0:
            raise ValueError("n must be nonnegative")
        if n == 0:
            return self
        if n == 1:
            if self.is_pseudomanifold():
                # Use one-point compactification of Datta. The
                # construction is a bit slower, but the resulting
                # complex is smaller.
                V = self.vertices()
                u = V[0]
                w = 0
                while w in V:
                    w += 1
                w = Simplex([w])
                new_facets = []
                for f in self.facets():
                    if u not in f:
                        new_facets.append(f.join(Simplex([u]), rename_vertices=False))
                    new_facets.append(f.join(w, rename_vertices=False))
                return SimplicialComplex(new_facets)
            else:
                return self.join(SimplicialComplex([["0"], ["1"]], is_mutable=is_mutable),
                                 rename_vertices=True)
        return self.suspension(1, is_mutable).suspension(int(n-1), is_mutable)

    def disjoint_union(self, right, rename_vertices=None, is_mutable=True):
        """
        The disjoint union of this simplicial complex with another one.

        INPUT:

        - ``right`` -- the other simplicial complex (the right-hand factor)

        - ``rename_vertices`` -- boolean (default: ``True``); if this is
          ``True``, the vertices in the disjoint union will be renamed by the
          formula: vertex "v" in the left-hand factor --> vertex "Lv" in the
          disjoint union, vertex "w" in the right-hand factor --> vertex "Rw"
          in the disjoint union.  If this is false, this tries to construct the
          disjoint union without renaming the vertices; this will cause
          problems if the two factors have any vertices with names in common.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S1.disjoint_union(S2).homology()                                      # needs sage.modules
            {0: Z, 1: Z, 2: Z}
        """
        if rename_vertices is not None:
            from sage.misc.superseded import deprecation
            deprecation(35907, 'the "rename_vertices" argument is deprecated')

        facets = []
        for f in self._facets:
            facets.append(tuple(["L" + str(v) for v in f]))
        for f in right._facets:
            facets.append(tuple(["R" + str(v) for v in f]))
        return SimplicialComplex(facets, is_mutable=is_mutable)

    def wedge(self, right, rename_vertices=True, is_mutable=True):
        """
        The wedge (one-point union) of this simplicial complex with
        another one.

        INPUT:

        - ``right`` -- the other simplicial complex (the right-hand factor)

        - ``rename_vertices`` -- boolean (default: ``True``); if this is
          ``True``, the vertices in the wedge will be renamed by the formula:
          first vertex in each are glued together and called "0".  Otherwise,
          each vertex "v" in the left-hand factor --> vertex "Lv" in the wedge,
          vertex "w" in the right-hand factor --> vertex "Rw" in the wedge.  If
          this is ``False``, this tries to construct the wedge without renaming
          the vertices; this will cause problems if the two factors have any
          vertices with names in common.

        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        .. NOTE::

            This operation is not well-defined if ``self`` or
            ``other`` is not path-connected.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S1.wedge(S2).homology()                                               # needs sage.modules
            {0: 0, 1: Z, 2: Z}
        """
        left_vertices = list(self.vertices())
        left_0 = left_vertices.pop(0)
        right_vertices = list(right.vertices())
        right_0 = right_vertices.pop(0)
        left_dict = {left_0: 0}
        right_dict = {right_0: 0}
        if rename_vertices:
            facets = []
            for v in left_vertices:
                left_dict[v] = "L" + str(v)
            for v in right_vertices:
                right_dict[v] = "R" + str(v)

            for f in self._facets:
                facets.append(tuple([left_dict[v] for v in f]))
            for f in right._facets:
                facets.append(tuple([right_dict[v] for v in f]))
        else:
            facets = self._facets + right._facets
        return SimplicialComplex(facets, is_mutable=is_mutable)

    def chain_complex(self, subcomplex=None, augmented=False,
                      verbose=False, check=False, dimensions=None,
                      base_ring=ZZ, cochain=False):
        r"""
        The chain complex associated to this simplicial complex.

        INPUT:

        - ``dimensions`` -- if ``None``, compute the chain complex in all
          dimensions.  If a list or tuple of integers, compute the
          chain complex in those dimensions, setting the chain groups
          in all other dimensions to zero.
        - ``base_ring`` -- commutative ring (default: ``ZZ``)
        - ``subcomplex`` -- a subcomplex of this simplicial complex (default:
          empty); compute the chain complex relative to this subcomplex
        - ``augmented`` -- boolean (default: ``False``); if ``True``, return
          the augmented chain complex (that is, include a class in dimension
          `-1` corresponding to the empty cell). This is ignored if
          ``dimensions`` is specified
        - ``cochain`` -- boolean (default: ``False``); if ``True``, return the
          cochain complex (that is, the dual of the chain complex)
        - ``verbose`` -- boolean (default: ``False``); if ``True``, print some
          messages as the chain complex is computed
        - ``check`` -- boolean (default: ``False``); if ``True``, make sure
          that the chain complex is actually a chain complex: the differentials
          are composable and their product is zero

        .. NOTE::

           If subcomplex is nonempty, then the argument ``augmented``
           has no effect: the chain complex relative to a nonempty
           subcomplex is zero in dimension `-1`.

        The rows and columns of the boundary matrices are indexed by
        the lists given by the :meth:`_n_cells_sorted` method, which by
        default are sorted.

        EXAMPLES::

            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle.chain_complex()                                                # needs sage.modules
            Chain complex with at most 2 nonzero terms over Integer Ring
            sage: circle.chain_complex()._latex_()                                      # needs sage.modules
            '\\Bold{Z}^{3} \\xrightarrow{d_{1}} \\Bold{Z}^{3}'
            sage: circle.chain_complex(base_ring=QQ, augmented=True)                    # needs sage.modules
            Chain complex with at most 3 nonzero terms over Rational Field
        """
        # initialize subcomplex
        if subcomplex is None:
            subcomplex = SimplicialComplex(is_mutable=False)
        else:
            # subcomplex is not empty, so don't augment the chain complex
            augmented = False
            # Use an immutable copy of the subcomplex
            if subcomplex._is_immutable:
                subcomplex = SimplicialComplex(subcomplex._facets, maximality_check=False,
                                               is_mutable=False)
        # now construct the range of dimensions in which to compute
        if dimensions is None:
            dimensions = range(self.dimension() + 1)
            first = 0
        else:
            augmented = False
            first = dimensions[0]
        dimensions = list(dimensions)
        differentials = {}
        # in the chain complex, compute the first dimension by hand,
        # and don't cache it: it may be differ from situation to
        # situation because of boundary effects.
        current = None
        current_dim = None
        if augmented:  # then first == 0
            current = self._n_cells_sorted(0, subcomplex=subcomplex)
            current_dim = 0
            if cochain:
                differentials[-1] = matrix(base_ring, len(current), 1,
                                           [1]*len(current))
            else:
                differentials[0] = matrix(base_ring, 1, len(current),
                                          [1]*len(current))
        elif first == 0 and not augmented:
            current = self._n_cells_sorted(0, subcomplex=subcomplex)
            current_dim = 0
            if not cochain:
                differentials[0] = matrix(base_ring, 0, len(current))
        else:  # first > 0
            current = self._n_cells_sorted(first, subcomplex=subcomplex)
            current_dim = first
            if not cochain:
                differentials[first] = matrix(base_ring, 0, len(current))
        for n in dimensions[1:]:
            if verbose:
                print("  starting dimension %s" % n)
            if (n, subcomplex) in self._complex:
                if cochain:
                    differentials[n-1] = self._complex[(n, subcomplex)].transpose().change_ring(base_ring)
                    mat = differentials[n-1]
                else:
                    differentials[n] = self._complex[(n, subcomplex)].change_ring(base_ring)
                    mat = differentials[n]
                if verbose:
                    print("    boundary matrix (cached): it's {} by {}.".format(mat.nrows(), mat.ncols()))
            else:
                # 'current' is the list of faces in dimension n
                #
                # 'old' is a dictionary, with keys the faces in the
                # previous dimension (dim n-1 for the chain complex,
                # n+1 for the cochain complex), values the integers 0,
                # 1, 2, ... (the index of the face).  finding an entry
                # in a dictionary seems to be faster than finding the
                # index of an entry in a list.
                if current_dim == n-1:
                    old = dict(zip(current, range(len(current))))
                else:
                    set_of_faces = self._n_cells_sorted(n-1, subcomplex=subcomplex)
                    old = dict(zip(set_of_faces, range(len(set_of_faces))))
                current = self._n_cells_sorted(n, subcomplex=subcomplex)
                current_dim = n
                # construct matrix.  it is easiest to construct it as
                # a sparse matrix, specifying which entries are
                # nonzero via a dictionary.
                matrix_data = {}
                col = 0
                if len(old) and len(current):
                    for simplex in current:
                        for i in range(n + 1):
                            face_i = simplex.face(i)
                            try:
                                matrix_data[(old[face_i], col)] = (-1)**i
                            except KeyError:
                                pass
                        col += 1
                mat = matrix(ZZ, len(old), len(current), matrix_data)
                if cochain:
                    self._complex[(n, subcomplex)] = mat
                    differentials[n-1] = mat.transpose().change_ring(base_ring)
                else:
                    self._complex[(n, subcomplex)] = mat
                    differentials[n] = mat.change_ring(base_ring)
                if verbose:
                    print("    boundary matrix computed: it's {} by {}.".format(mat.nrows(), mat.ncols()))
        # now for the cochain complex, compute the last dimension by
        # hand, and don't cache it.
        if cochain:
            n = dimensions[-1] + 1
            if current_dim != n-1:
                current = self._n_cells_sorted(n-1, subcomplex=subcomplex)
            differentials[n-1] = matrix(base_ring, 0, len(current))
        # finally, return the chain complex
        from sage.homology.chain_complex import ChainComplex

        if cochain:
            return ChainComplex(data=differentials, degree=1,
                                base_ring=base_ring, check=check)
        else:
            return ChainComplex(data=differentials, degree=-1,
                                base_ring=base_ring, check=check)

    def _homology_(self, dim=None, base_ring=ZZ, subcomplex=None,
                   cohomology=False, enlarge=True, algorithm='pari',
                   verbose=False, reduced=True, generators=False):
        """
        The (reduced) homology of this simplicial complex.

        INPUT:

        - ``dim`` -- integer or list of integers or ``None`` (default:
          ``None``); if ``None``, then return the homology in every dimension.
          If ``dim`` is an integer or list, return the homology in the given
          dimensions.  (Actually, if ``dim`` is a list, return the homology in
          the range from ``min(dim)`` to ``max(dim)``.)

        - ``base_ring`` -- commutative ring (default: ``ZZ``); must be ``ZZ``
          or a field

        - ``subcomplex`` -- a subcomplex of this simplicial complex (default:
          ``None``); compute homology relative to this subcomplex

        - ``cohomology`` -- boolean (default: ``False``); if ``True``, compute
          cohomology rather than homology

        - ``enlarge`` -- boolean (default: ``True``); if ``True``, find a new
          subcomplex homotopy equivalent to, and probably larger than, the
          given one

        - ``algorithm`` -- string (default: ``'pari'``); the options are
          ``'auto'``, ``'dhsw'``, or ``'pari'``.  (``'no_chomp'`` is a synonym
          for ``'auto'``, maintained for backward compatibility.) If
          ``'auto'``, use the Dumas, Heckenbach, Saunders, and Welker
          elimination algorithm for large matrices, Pari for small ones. If
          ``'pari'``, then compute elementary divisors using Pari. If
          ``'dhsw'``, then use the DHSW algorithm to compute elementary
          divisors.  (As of this writing, ``'pari'`` is the fastest standard
          option.)

        - ``verbose`` -- boolean (default: ``False``); if ``True``, print some
          messages as the homology is computed

        - ``reduced`` -- boolean (default: ``Trues``); if ``True``, return the
          reduced homology

        - ``generators`` -- boolean (default: ``False``); if ``True``, return
          the homology groups and also generators for them

        ALGORITHM: if ``generators`` is ``True``, directly compute the
        chain complex, compute its homology along with its generators,
        and then convert the chain complex generators to chains in the
        simplicial complex.

        Otherwise: if ``subcomplex`` is ``None``, replace it with a
        facet -- a contractible subcomplex of the original complex.
        Then as long as ``enlarge`` is ``True``, no matter what
        ``subcomplex`` is, replace it with a subcomplex `L` which is
        homotopy equivalent and as large as possible.  Compute the
        homology of the original complex relative to `L`: if `L` is
        large, then the relative chain complex will be small enough to
        speed up computations considerably.

        EXAMPLES::

            sage: # needs sage.modules
            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle._homology_()
            {0: 0, 1: Z}
            sage: sphere = SimplicialComplex([[0,1,2,3]])
            sage: sphere.remove_face([0,1,2,3])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)}
            sage: sphere._homology_()
            {0: 0, 1: 0, 2: Z}
            sage: sphere._homology_(reduced=False)
            {0: Z, 1: 0, 2: Z}
            sage: sphere._homology_(base_ring=GF(2), reduced=False)
            {0: Vector space of dimension 1 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}

        We need an immutable complex to compute homology generators::

            sage: # needs sage.modules
            sage: sphere.set_immutable()
            sage: sphere._homology_(generators=True)
            {0: [], 1: [], 2: [(Z, (0, 1, 2) - (0, 1, 3) + (0, 2, 3) - (1, 2, 3))]}

        Another way to get a two-sphere: take a two-point space and take its
        three-fold join with itself::

            sage: S = SimplicialComplex([[0], [1]])
            sage: (S*S*S)._homology_(dim=2, cohomology=True)                            # needs sage.modules
            Z

        The same computation, done without finding a contractible subcomplex::

            sage: (S*S*S)._homology_(dim=2, cohomology=True, enlarge=False)             # needs sage.modules
            Z

        Relative homology::

            sage: T = SimplicialComplex([[0,1,2]])
            sage: U = SimplicialComplex([[0,1], [1,2], [0,2]])
            sage: T._homology_(subcomplex=U)                                            # needs sage.modules
            {0: 0, 1: 0, 2: Z}

        Generators::

            sage: simplicial_complexes.Torus().homology(generators=True)                # needs sage.modules
            {0: [],
             1: [(Z, (2, 4) - (2, 6) + (4, 6)), (Z, (1, 4) - (1, 6) + (4, 6))],
             2: [(Z, (0, 1, 2) - (0, 1, 5) + (0, 2, 6) - (0, 3, 4) + (0, 3, 5)
                      - (0, 4, 6) - (1, 2, 4) + (1, 3, 4) - (1, 3, 6) + (1, 5, 6)
                      - (2, 3, 5) + (2, 3, 6) + (2, 4, 5) - (4, 5, 6))]}
        """
        from sage.homology.homology_group import HomologyGroup

        if dim is not None:
            if isinstance(dim, (list, tuple, range)):
                low = min(dim) - 1
                high = max(dim) + 2
            else:
                low = dim - 1
                high = dim + 2
            dims = range(low, high)
        else:
            dims = None

        if generators:
            enlarge = False

        if verbose:
            print("starting calculation of the homology of this")
            print("%s-dimensional simplicial complex" % self.dimension())
        if subcomplex is None:
            if enlarge:
                if verbose:
                    print("Constructing contractible subcomplex...")
                L = self._contractible_subcomplex(verbose=verbose)
                if verbose:
                    print("Done finding contractible subcomplex.")
                    vec = [len(self.faces(subcomplex=L)[n-1]) for n in range(self.dimension()+2)]
                    print("The difference between the f-vectors is:")
                    print("  %s" % vec)
            else:
                L = SimplicialComplex([[self.vertices()[0]]])
        else:
            if enlarge:
                if verbose:
                    print("Enlarging subcomplex...")
                L = self._enlarge_subcomplex(subcomplex, verbose=verbose)
                if verbose:
                    print("Done enlarging subcomplex:")
            else:
                L = subcomplex
        L.set_immutable()

        if verbose:
            print("Computing the chain complex...")
        C = self.chain_complex(dimensions=dims, augmented=reduced,
                               cochain=cohomology, base_ring=base_ring,
                               subcomplex=L, verbose=verbose)
        if verbose:
            print(" Done computing the chain complex. ")
            print("Now computing homology...")
        answer = C.homology(base_ring=base_ring, verbose=verbose,
                            algorithm=algorithm, generators=generators)

        if generators:
            # Convert chain complex information to simplicial complex
            # information.
            for i in answer:
                H_with_gens = answer[i]
                if H_with_gens:
                    chains = self.n_chains(i, base_ring=base_ring)
                    new_H = []
                    for H, gen in H_with_gens:
                        v = gen.vector(i)
                        new_gen = chains.zero()
                        for coeff, chaine in zip(v, chains.gens()):
                            new_gen += coeff * chaine
                        new_H.append((H, new_gen))
                    answer[i] = new_H

        if dim is None:
            dim = range(self.dimension() + 1)
        zero = HomologyGroup(0, base_ring)
        if isinstance(dim, (list, tuple, range)):
            # Fix non-reduced answer.
            if subcomplex is None and not reduced and 0 in dim:
                try:
                    if base_ring in Fields():
                        rank = answer[0].dimension()
                    else:
                        rank = len(answer[0].invariants())
                except KeyError:
                    rank = 0
                answer[0] = HomologyGroup(rank + 1, base_ring)
            return dict([d, answer.get(d, zero)] for d in dim)
        return answer.get(dim, zero)

    # This is cached for speed reasons: it can be very slow to run
    # this function.
    @cached_method
    def algebraic_topological_model(self, base_ring=None):
        r"""
        Algebraic topological model for this simplicial complex with
        coefficients in ``base_ring``.

        The term "algebraic topological model" is defined by Pilarczyk
        and Réal [PR2015]_.

        INPUT:

        - ``base_ring`` -- coefficient ring (default: ``QQ``); must be a field

        Denote by `C` the chain complex associated to this simplicial
        complex. The algebraic topological model is a chain complex
        `M` with zero differential, with the same homology as `C`,
        along with chain maps `\pi: C \to M` and `\iota: M \to C`
        satisfying `\iota \pi = 1_M` and `\pi \iota` chain homotopic
        to `1_C`. The chain homotopy `\phi` must satisfy

        - `\phi \phi = 0`,
        - `\pi \phi = 0`,
        - `\phi \iota = 0`.

        Such a chain homotopy is called a *chain contraction*.

        OUTPUT: a pair consisting of

        - chain contraction ``phi`` associated to `C`, `M`, `\pi`, and
          `\iota`
        - the chain complex `M`

        Note that from the chain contraction ``phi``, one can recover the
        chain maps `\pi` and `\iota` via ``phi.pi()`` and
        ``phi.iota()``. Then one can recover `C` and `M` from, for
        example, ``phi.pi().domain()`` and ``phi.pi().codomain()``,
        respectively.

        EXAMPLES::

            sage: # needs sage.modules
            sage: RP2 = simplicial_complexes.RealProjectivePlane()
            sage: phi, M = RP2.algebraic_topological_model(GF(2))
            sage: M.homology()
            {0: Vector space of dimension 1 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}
            sage: T = simplicial_complexes.Torus()
            sage: phi, M = T.algebraic_topological_model(QQ)
            sage: M.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 2 over Rational Field,
             2: Vector space of dimension 1 over Rational Field}
        """
        from sage.homology.algebraic_topological_model import algebraic_topological_model
        if base_ring is None:
            base_ring = QQ
        return algebraic_topological_model(self, base_ring)

    def alexander_whitney(self, simplex, dim_left):
        r"""
        Subdivide this simplex into a pair of simplices.

        If this simplex has vertices `v_0`, `v_1`, ..., `v_n`, then
        subdivide it into simplices `(v_0, v_1, ..., v_{dim})` and
        `(v_{dim}, v_{dim + 1}, ..., v_n)`.

        See :meth:`Simplex.alexander_whitney` for more details. This
        method just calls that one.

        INPUT:

        - ``simplex`` -- a simplex in this complex
        - ``dim`` -- integer between 0 and one more than the
          dimension of this simplex

        OUTPUT: list containing just the triple ``(1, left,
        right)``, where ``left`` and ``right`` are the two simplices
        described above.

        EXAMPLES::

            sage: s = Simplex((0,1,3,4))
            sage: X = SimplicialComplex([s])
            sage: X.alexander_whitney(s, 0)
            [(1, (0,), (0, 1, 3, 4))]
            sage: X.alexander_whitney(s, 2)
            [(1, (0, 1, 3), (3, 4))]
        """
        return simplex.alexander_whitney(dim_left)

    def add_face(self, face):
        """
        Add a face to this simplicial complex.

        INPUT:

        - ``face`` -- a subset of the vertex set

        This *changes* the simplicial complex, adding a new face and all
        of its subfaces.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [0,2]])
            sage: X.add_face([0,1,2,]); X
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
            sage: Y = SimplicialComplex(); Y
            Simplicial complex with vertex set () and facets {()}
            sage: Y.add_face([0,1])
            sage: Y.add_face([1,2,3])
            sage: Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1), (1, 2, 3)}

        If you add a face which is already present, there is no effect::

            sage: Y.add_face([1,3]); Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1), (1, 2, 3)}

        TESTS:

        Check that the bug reported at :issue:`14354` has been fixed::

            sage: T = SimplicialComplex([range(1,5)]).n_skeleton(1)
            sage: T.homology()                                                          # needs sage.modules
            {0: 0, 1: Z x Z x Z}
            sage: T.add_face([1,2,3])
            sage: T.homology()                                                          # needs sage.modules
            {0: 0, 1: Z x Z, 2: 0}

        Check that the ``_faces`` cache is treated correctly
        (:issue:`20758`)::

            sage: T = SimplicialComplex([range(1,5)]).n_skeleton(1)
            sage: _ = T.faces()       # populate the _faces attribute
            sage: _ = T.homology()    # add more to _faces                              # needs sage.modules
            sage: T.add_face((1,2,3))
            sage: all(Simplex((1,2,3)) in T._faces[L][2] for L in T._faces)
            True

        Check that the ``__enlarged`` cache is treated correctly
        (:issue:`20758`)::

            sage: T = SimplicialComplex([range(1,5)]).n_skeleton(1)
            sage: T.homology()  # to populate the __enlarged attribute                  # needs sage.modules
            {0: 0, 1: Z x Z x Z}
            sage: T.add_face([1,2,3])
            sage: len(T._SimplicialComplex__enlarged) > 0                               # needs sage.modules
            True

        Check we've fixed the bug reported at :issue:`14578`::

            sage: t0 = SimplicialComplex()
            sage: t0.add_face(('a', 'b'))
            sage: t0.add_face(('c', 'd', 'e'))
            sage: t0.add_face(('e', 'f', 'c'))
            sage: t0.homology()                                                         # needs sage.modules
            {0: Z, 1: 0, 2: 0}

        Check that we've fixed the bug reported at :issue:`22880`::

            sage: X = SimplicialComplex([[0], [1]])
            sage: temp = X.faces(SimplicialComplex(()))
            sage: X.add_face([0,1])
        """
        if self._is_immutable:
            raise ValueError("this simplicial complex is not mutable")

        vertex_to_index = self._translation_to_numeric()

        # Update vertex_to_index by giving each new vertex a larger
        # entry than the existing ones.
        if vertex_to_index:
            idx = max(vertex_to_index.values()) + 1
        else:
            idx = 0
        new_vertices = []
        for v in face:
            if v not in self.vertices():
                new_vertices.append(v)
                vertex_to_index[v] = idx
                idx += 1

        new_face = Simplex(sorted(face, key=vertex_to_index.__getitem__))

        face_is_maximal = True
        for other in self._facets:
            if face_is_maximal:
                face_is_maximal = not new_face.is_face(other)
        if face_is_maximal:
            # remove any old facets which are no longer maximal
            Facets = list(self._facets)
            for old_face in self._facets:
                if old_face.is_face(new_face):
                    Facets.remove(old_face)
            # add new_face to facet list
            Facets.append(new_face)
            self._facets = Facets

            # Update the vertex set
            self._vertex_to_index = vertex_to_index

            # Update self._faces.
            all_new_faces = SimplicialComplex([new_face]).faces()
            for L in self._faces:
                L_complex = self._faces[L]
                for dim in range(new_face.dimension()+1):
                    if dim in L_complex:
                        if L is None:
                            new_faces = all_new_faces[dim]
                        else:
                            new_faces = all_new_faces[dim].difference(L.n_cells(dim))
                        L_complex[dim] = L_complex[dim].union(new_faces)
                    else:
                        L_complex[dim] = all_new_faces[dim]
            # update self._graph if necessary
            if self._graph is not None:
                d = new_face.dimension()+1
                for i in range(d):
                    for j in range(i + 1, d):
                        self._graph.add_edge(new_face[i], new_face[j])
            self._complex = {}
            self.__contractible = None
            self._bbn = {}
            self._bbn_all_computed = set()

    def remove_face(self, face, check=False):
        """
        Remove a face from this simplicial complex.

        INPUT:

        - ``face`` -- a face of the simplicial complex

        - ``check`` -- boolean (default: ``False``); if
          ``True``, raise an error if ``face`` is not a
          face of this simplicial complex

        This does not return anything; instead, it *changes* the
        simplicial complex.

        ALGORITHM:

        The facets of the new simplicial complex are
        the facets of the original complex not containing ``face``,
        together with those of ``link(face)*boundary(face)``.

        EXAMPLES::

            sage: S = range(1,5)
            sage: Z = SimplicialComplex([S]); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.remove_face([1,2])
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and
             facets {(1, 3, 4), (2, 3, 4)}

            sage: S = SimplicialComplex([[0,1,2],[2,3]])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(2, 3), (0, 1, 2)}
            sage: S.remove_face([0,1,2])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1), (0, 2), (1, 2), (2, 3)}

        TESTS:

        Check that the ``_faces`` cache is treated properly: see
        :issue:`20758`::

            sage: T = SimplicialComplex([range(1,5)]).n_skeleton(1)
            sage: _ = T.faces()     # populate the _faces attribute
            sage: _ = T.homology()  # add more to _faces                                # needs sage.modules
            sage: T.add_face((1,2,3))
            sage: T.remove_face((1,2,3))
            sage: len(T._faces)                                                         # needs sage.modules
            2
            sage: T.remove_face((1,2))
            sage: len(T._faces)
            1

        Check that the face to be removed can be given with a
        different vertex ordering::

            sage: S = SimplicialComplex([[1,2], [1,3]])
            sage: S.remove_face([3,1])
            sage: S
            Simplicial complex with vertex set (1, 2, 3) and facets {(3,), (1, 2)}
        """
        if self._is_immutable:
            raise ValueError("this simplicial complex is not mutable")

        getindex = self._translation_to_numeric().__getitem__
        simplex = Simplex(sorted(face, key=getindex))
        facets = self.facets()
        if all(not simplex.is_face(F) for F in facets):
            # face is not in self
            if check:
                raise ValueError('trying to remove a face which is not in the simplicial complex')
            return
        link = self.link(simplex)
        join_facets = []
        for f in simplex.faces():
            for g in link.facets():
                join_facets.append(f.join(g, rename_vertices=False))
        # join_facets is the list of facets in the join bdry(face) * link(face)
        remaining = join_facets + [elem for elem in facets if not simplex.is_face(elem)]

        # Check to see if there are any non-maximal faces
        # build set of facets
        self._facets = []
        for f in remaining:
            face2 = Simplex(f)
            face_is_maximal = True
            faces_to_be_removed = []
            for other in self._facets:
                if other.is_face(face2):
                    faces_to_be_removed.append(other)
                elif face_is_maximal:
                    face_is_maximal = not face2.is_face(other)
            for x in faces_to_be_removed:
                self._facets.remove(x)
            face2 = Simplex(sorted(face2.tuple()))
            if face_is_maximal:
                self._facets.append(face2)
        # if no maximal faces, add the empty face as a facet
        if len(remaining) == 0:
            self._facets.append(Simplex(-1))

        # Recreate the vertex set
        vertices = set(chain.from_iterable(self._facets))
        for v in self.vertices():
            if v not in vertices:
                del self._vertex_to_index[v]

        # Update self._faces.
        # Note: can't iterate over self._faces, because the dictionary
        # size may change during iteration.
        for L in list(self._faces):
            del self._faces[L]
            if L is None or Simplex(face) not in L:
                self.faces(L)
        # Update self._graph if necessary.
        if self._graph is not None:
            # Only if removing a 1 or 2 dim face will the graph be affected
            if len(face) == 1:
                self._graph.delete_vertex(face[0])
                self._graph.add_vertex(face[0])
            elif len(face) == 2:
                self._graph.delete_edge(face[0], face[1])
        self._complex = {}
        self.__contractible = None
        self.__enlarged = {}
        self._bbn = {}
        self._bbn_all_computed = set()

    def remove_faces(self, faces, check=False):
        """
        Remove a collection of faces from this simplicial complex.

        INPUT:

        - ``faces`` -- list (or any iterable) of faces of the simplicial
          complex

        - ``check`` -- boolean (default: ``False``); if ``True``, raise an
          error if any element of ``faces`` is not a face of this simplicial
          complex

        This does not return anything; instead, it *changes* the
        simplicial complex.

        ALGORITHM:

        Run ``self.remove_face(f)`` repeatedly, for ``f`` in ``faces``.

        EXAMPLES::

            sage: S = range(1,5)
            sage: Z = SimplicialComplex([S]); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.remove_faces([[1,2]])
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 3, 4), (2, 3, 4)}

            sage: Z = SimplicialComplex([S]); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.remove_faces([[1,2], [2,3]])
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(2, 4), (1, 3, 4)}

        TESTS:

        Check the ``check`` argument::

            sage: Z = SimplicialComplex([[1,2,3,4]])
            sage: Z.remove_faces([[1,2], [3,4]])
            sage: Z.remove_faces([[1,2]])
            sage: Z.remove_faces([[1,2]], check=True)
            Traceback (most recent call last):
            ...
            ValueError: trying to remove a face which is not in the simplicial complex
        """
        for f in faces:
            self.remove_face(f, check=check)

    def is_subcomplex(self, other):
        """
        Return ``True`` if this is a subcomplex of ``other``.

        INPUT:

        - ``other`` -- another simplicial complex

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1.is_subcomplex(S1)
            True
            sage: Empty = SimplicialComplex()
            sage: Empty.is_subcomplex(S1)
            True
            sage: S1.is_subcomplex(Empty)
            False

            sage: sorted(S1.facets())
            [(0, 1), (0, 2), (1, 2)]
            sage: T = S1.product(S1)
            sage: sorted(T.facets())[0] # typical facet in T
            ('L0R0', 'L0R1', 'L1R1')
            sage: S1.is_subcomplex(T)
            False
            sage: T._contractible_subcomplex().is_subcomplex(T)
            True
        """
        return all(f in other for f in self.facets())

    def connected_sum(self, other, is_mutable=True):
        """
        The connected sum of this simplicial complex with another one.

        INPUT:

        - ``other`` -- another simplicial complex
        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        OUTPUT: the connected sum ``self # other``

        .. WARNING::

           This does not check that ``self`` and ``other`` are manifolds,
           only that their facets all have the same dimension.  Since a
           (more or less) random facet is chosen from each complex and
           then glued together, this method may return random
           results if applied to non-manifolds, depending on which
           facet is chosen.

        Algorithm: a facet is chosen from each surface, and removed.
        The vertices of these two facets are relabeled to
        ``(0,1,...,dim)``.  Of the remaining vertices, the ones from
        the left-hand factor are renamed by prepending an "L", and
        similarly the remaining vertices in the right-hand factor are
        renamed by prepending an "R".

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1.connected_sum(S1.connected_sum(S1)).homology()                     # needs sage.modules
            {0: 0, 1: Z}
            sage: P = simplicial_complexes.RealProjectivePlane(); P
            Minimal triangulation of the real projective plane
            sage: P.connected_sum(P)    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets

        The notation '+' may be used for connected sum, also::

            sage: P + P    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets
            sage: (P + P).homology()[1]                                                 # needs sage.modules
            Z x C2
        """
        if not (self.is_pure() and other.is_pure() and
                self.dimension() == other.dimension()):
            raise ValueError("complexes are not pure of the same dimension")
        # first find a top-dimensional simplex to remove from each surface
        keep_left = self._facets[0]
        keep_right = other._facets[0]
        # construct the set of facets:
        left = set(self._facets).difference({keep_left})
        right = set(other._facets).difference({keep_right})
        facet_set = ([[rename_vertex(v, keep=list(keep_left))
                       for v in face] for face in left]
                     + [[rename_vertex(v, keep=list(keep_right), left=False)
                         for v in face] for face in right])
        # return the new surface
        return SimplicialComplex(facet_set, is_mutable=is_mutable)

    __add__ = connected_sum

    def link(self, simplex, is_mutable=True):
        r"""
        The link of a simplex in this simplicial complex.

        The link of a simplex `F` is the simplicial complex formed by
        all simplices `G` which are disjoint from `F` but for which `F
        \cup G` is a simplex.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial complex
        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2], [1,2,3]])
            sage: X.link(Simplex([0]))
            Simplicial complex with vertex set (1, 2) and facets {(1, 2)}
            sage: X.link([1,2])
            Simplicial complex with vertex set (0, 3) and facets {(0,), (3,)}
            sage: Y = SimplicialComplex([[0,1,2,3]])
            sage: Y.link([1])
            Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2, 3)}
        """
        faces = []
        s = Simplex(simplex)
        for f in self._facets:
            if s.is_face(f):
                faces.append(Simplex(f.set().difference(s.set())))
        return SimplicialComplex(faces, is_mutable=is_mutable)

    def star(self, simplex, is_mutable=True):
        """
        Return the star of a simplex in this simplicial complex.

        The star of ``simplex`` is the simplicial complex formed by
        all simplices which contain ``simplex``.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial complex
        - ``is_mutable`` -- boolean (default: ``True``); determines if the output
          is mutable

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2], [1,2,3]])
            sage: X.star(Simplex([0]))
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
            sage: X.star(Simplex([1]))
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (1, 2, 3)}
            sage: X.star(Simplex([1,2]))
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (1, 2, 3)}
            sage: X.star(Simplex([]))
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (1, 2, 3)}
        """
        faces = []
        s = Simplex(simplex)
        for f in self._facets:
            if s.is_face(f):
                faces.append(f)
        return SimplicialComplex(faces, is_mutable=is_mutable)

    def is_cohen_macaulay(self, base_ring=QQ, ncpus=0):
        r"""
        Return ``True`` if ``self`` is Cohen-Macaulay.

        A simplicial complex `\Delta` is Cohen-Macaulay over `R` iff
        `\tilde{H}_i(\mathrm{lk}_\Delta(F);R) = 0` for all
        `F \in \Delta` and `i < \dim\mathrm{lk}_\Delta(F)`.
        Here, `\Delta` is ``self`` and `R` is ``base_ring``, and
        `\mathrm{lk}` denotes the link operator on ``self``.

        INPUT:

        - ``base_ring`` -- (default: ``QQ``) the base ring

        - ``ncpus`` -- (default: 0) number of cpus used for the
          computation. If this is 0, determine the number of cpus
          automatically based on the hardware being used.

        For finite simplicial complexes, this is equivalent to the
        statement that the Stanley-Reisner ring of ``self`` is
        Cohen-Macaulay.

        EXAMPLES:

        Spheres are Cohen-Macaulay::

            sage: S = SimplicialComplex([[1,2],[2,3],[3,1]])
            sage: S.is_cohen_macaulay(ncpus=3)                                          # needs sage.modules
            True

        The following example is taken from Bruns, Herzog - Cohen-Macaulay
        rings, Figure 5.3::

            sage: S = SimplicialComplex([[1,2,3],[1,4,5]])
            sage: S.is_cohen_macaulay(ncpus=3)                                          # needs sage.modules
            False

        The choice of base ring can matter.  The real projective plane `\RR P^2`
        has `H_1(\RR P^2) = \ZZ/2`, hence is CM over `\QQ` but not over `\ZZ`. ::

            sage: X = simplicial_complexes.RealProjectivePlane()
            sage: X.is_cohen_macaulay()                                                 # needs sage.modules
            True
            sage: X.is_cohen_macaulay(ZZ)                                               # needs sage.modules
            False
        """
        from sage.parallel.decorate import parallel

        if not ncpus:
            from sage.parallel.ncpus import ncpus as get_ncpus
            ncpus = get_ncpus()

        facs = list(self.face_iterator())
        n = len(facs)
        facs_divided = [[] for i in range(ncpus)]
        for i in range(n):
            facs_divided[i % ncpus].append(facs[i])

        def all_homologies_vanish(F):
            S = self.link(F)
            H = S.homology(base_ring=base_ring)
            if base_ring in Fields():
                return all(H[j].dimension() == 0 for j in range(S.dimension()))
            else:
                return not any(H[j].invariants() for j in range(S.dimension()))

        @parallel(ncpus=ncpus)
        def all_homologies_in_list_vanish(Fs):
            return all(all_homologies_vanish(F) for F in Fs)

        return all(answer[1] for answer in all_homologies_in_list_vanish(facs_divided))

    def generated_subcomplex(self, sub_vertex_set, is_mutable=True):
        """
        Return the largest sub-simplicial complex of ``self`` containing
        exactly ``sub_vertex_set`` as vertices.

        INPUT:

        - ``sub_vertex_set`` -- the sub-vertex set
        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: S
            Minimal triangulation of the 2-sphere
            sage: S.generated_subcomplex([0,1,2])
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
        """
        if not set(self.vertices()).issuperset(sub_vertex_set):
            raise ValueError("input must be a subset of the vertex set")
        faces = []
        for i in range(self.dimension() + 1):
            for j in self.faces()[i]:
                if j.set().issubset(sub_vertex_set):
                    faces.append(j)
        return SimplicialComplex(faces, maximality_check=True,
                                 is_mutable=is_mutable)

    def is_shelling_order(self, shelling_order, certificate=False):
        r"""
        Return if the order of the facets given by ``shelling_order``
        is a shelling order for ``self``.

        A sequence of facets `(F_i)_{i=1}^N` of a simplicial
        complex of dimension `d` is a *shelling order* if for all
        `i = 2, 3, 4, \ldots`, the complex

        .. MATH::

            X_i = \left( \bigcup_{j=1}^{i-1} F_j \right) \cap F_i

        is pure and of dimension `\dim F_i - 1`.

        INPUT:

        - ``shelling_order`` -- an ordering of the facets of ``self``
        - ``certificate`` -- boolean (default: ``False``); if ``True`` then returns
          the index of the first facet that violate the condition

        .. SEEALSO::

            :meth:`is_shellable`

        EXAMPLES::

            sage: facets = [[1,2,5],[2,3,5],[3,4,5],[1,4,5]]
            sage: X = SimplicialComplex(facets)
            sage: X.is_shelling_order(facets)
            True

            sage: b = [[1,2,5], [3,4,5], [2,3,5], [1,4,5]]
            sage: X.is_shelling_order(b)
            False
            sage: X.is_shelling_order(b, True)
            (False, 1)

        A non-pure example::

            sage: facets = [[1,2,3], [3,4], [4,5], [5,6], [4,6]]
            sage: X = SimplicialComplex(facets)
            sage: X.is_shelling_order(facets)
            True

        REFERENCES:

        - [BW1996]_
        """
        # Quick check by Lemma 2.2 in [BW1996]
        if self.dimension() != len(list(shelling_order[0])) - 1:
            return False

        cur_complex = SimplicialComplex([])
        for i, F in enumerate(shelling_order):
            if i > 0:
                # The shelling condition is precisely that intersection is
                #    a pure complex of one dimension less and stop if this fails
                common = set(F).intersection(set(cur_complex.vertices()))
                intersection = cur_complex.generated_subcomplex(list(common))

                dim = len(list(F)) - 1
                if not intersection.is_pure() or dim - 1 != intersection.dimension():
                    if certificate:
                        return (False, i)
                    return False
            cur_complex.add_face(F)
        return True

    @cached_method
    def is_shellable(self, certificate=False):
        r"""
        Return if ``self`` is shellable.

        A simplicial complex is shellable if there exists a shelling
        order.

        .. NOTE::

            1. This method can check all orderings of the facets by brute
               force, hence can be very slow.

            2. This is shellability in the general (nonpure) sense of
               Bjorner and Wachs [BW1996]_. This method does not check purity.

        .. SEEALSO::

            :meth:`is_shelling_order`

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); if ``True`` then
          returns the shelling order (if it exists)

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2,5], [2,3,5], [3,4,5], [1,4,5]])
            sage: X.is_shellable()
            True
            sage: order = X.is_shellable(True); order
            ((1, 2, 5), (2, 3, 5), (1, 4, 5), (3, 4, 5))
            sage: X.is_shelling_order(order)
            True

            sage: X = SimplicialComplex([[1,2,3], [3,4,5]])
            sage: X.is_shellable()
            False

        Examples from Figure 1 in [BW1996]_::

            sage: X = SimplicialComplex([[1,2,3], [3,4], [4,5], [5,6], [4,6]])
            sage: X.is_shellable()
            True

            sage: X = SimplicialComplex([[1,2,3], [3,4], [4,5,6]])
            sage: X.is_shellable()
            False

        REFERENCES:

        - :wikipedia:`Shelling_(topology)`
        """
        if not certificate:
            return bool(self.is_shellable(certificate=True))

        if self.is_pure():
            if any(x < 0 for x in self.h_vector()):
                return False
        else:  # Non-pure complex
            if any(x < 0 for row in self.h_triangle() for x in row):
                return False

        facets = set(self.facets())
        cur_order = []
        # For consistency when using different Python versions, for example, sort 'faces'.
        it = [iter(sorted(facets, key=str))]
        cur_complex = SimplicialComplex([])
        while facets:
            try:
                F = next(it[-1])
            except StopIteration:
                # Backtrace
                if not cur_order:
                    return False
                it.pop()
                facets.add(cur_order.pop())
                cur_complex = SimplicialComplex(cur_order)
                continue

            # First facet must be top dimensional
            if not cur_order:
                if self.dimension() == F.dimension():
                    cur_complex.add_face(F)
                    cur_order.append(F)
                    facets.remove(F)
                    it.append(iter(set(facets)))
                continue

            # The shelling condition is precisely that intersection is
            #    a pure complex of one dimension less and stop if this fails
            common = set(F).intersection(set(cur_complex.vertices()))
            intersection = cur_complex.generated_subcomplex(list(common))

            if (not intersection.is_pure()
                    or F.dimension() - 1 != intersection.dimension()):
                continue
            cur_complex.add_face(F)
            cur_order.append(F)
            facets.remove(F)
            it.append(iter(set(facets)))  # Iterate over a copy of the current facets

        return tuple(cur_order)

    def restriction_sets(self, order):
        """
        Return the restriction sets of the facets according to ``order``.

        A restriction set of a shelling order is the sequence of
        smallest new faces that are created during the shelling order.

        .. SEEALSO::

            :meth:`is_shelling_order`

        EXAMPLES::

            sage: facets = [[1,2,5], [2,3,5], [3,4,5], [1,4,5]]
            sage: X = SimplicialComplex(facets)
            sage: X.restriction_sets(facets)
            [(), (3,), (4,), (1, 4)]

            sage: b = [[1,2,5], [3,4,5], [2,3,5], [1,4,5]]
            sage: X.restriction_sets(b)
            Traceback (most recent call last):
            ...
            ValueError: not a shelling order
        """
        # It starts with the first empty
        restrictions = [()]

        # Each time we hit a facet, the complement goes to the restriction
        cur_complex = SimplicialComplex([])
        for i, F in enumerate(order):
            if i > 0:
                # The shelling condition is precisely that intersection is
                #    a pure complex of one dimension less and stop if this fails
                common = set(F).intersection(set(cur_complex.vertices()))
                intersection = cur_complex.generated_subcomplex(list(common))

                if not intersection.is_pure() or self.dimension() - 1 > intersection.dimension():
                    raise ValueError("not a shelling order")
                faces = SimplicialComplex([F]).faces()
                for k, v in intersection.faces().items():
                    faces[k] = faces[k].difference(v)
                for k in sorted(faces.keys()):
                    if faces[k]:
                        restrictions.append(faces[k].pop())
                        break
            cur_complex.add_face(F)

        return restrictions

    def _complement(self, simplex):
        """
        Return the complement of a simplex in the vertex set of this
        simplicial complex.

        INPUT:

        - ``simplex`` -- a simplex (need not be in the simplicial complex)

        OUTPUT: its complement: the simplex formed by the vertices not
        contained in ``simplex``.

        Note that this only depends on the vertex set of the
        simplicial complex, not on its simplices.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2,3,4,5]])
            sage: X._complement([1,2,3])
            (0, 4, 5)
            sage: X._complement([0,1,3,4])
            (2, 5)
            sage: X._complement([0,4,1,3])
            (2, 5)
        """
        return Simplex(set(self.vertices()).difference(simplex))

    def _transpose_simplices(self, *simplices):
        """
        Given tuple ``L`` of simplices, returns new list, where each
        simplex is formed by taking a vertex from each simplex from
        ``L``.

        INPUT:

        - ``simplices`` -- a bunch of simplices

        If ``simplices`` consists of `(f_0, f_1, f_2, ...)`, then the
        output consists of all possible simplices of the form `(v_0,
        v_1, v_2, ...)`, where `v_i` is a vertex of `f_i`.  If a
        vertex appears more than once in such a simplex, remove all
        but one of its appearances.  If such a simplex contains others
        already produced, then ignore that larger simplex -- the
        output should be a list of minimal simplices constructed in
        this way.

        This is used in computing the minimal nonfaces and hence the
        Stanley-Reisner ring.

        Note that this only depends on the vertex set of the
        simplicial complex, not on its simplices.

        I don't know if there is a standard name for this, but it
        looked sort of like the transpose of a matrix; hence the name
        for this method.

        EXAMPLES::

            sage: X = SimplicialComplex()
            sage: X._transpose_simplices([1,2])
            [(1,), (2,)]
            sage: X._transpose_simplices([1,2], [3,4])
            [(1, 3), (1, 4), (2, 3), (2, 4)]

        In the following example, one can construct the simplices
        ``(1,2)`` and ``(1,3)``, but you can also construct ``(1,1) = (1,)``,
        which is a face of both of the others.  So the answer omits
        ``(1,2)`` and ``(1,3)``::

            sage: X._transpose_simplices([1,2], [1,3])
            [(1,), (2, 3)]
        """
        answer = []
        if len(simplices) == 1:
            answer = [Simplex((v,)) for v in simplices[0]]
        elif len(simplices) > 1:
            face = simplices[0]
            rest = simplices[1:]
            for v in face:
                for partial in self._transpose_simplices(*rest):
                    if v not in partial:
                        L = sorted([v] + list(partial))
                        simplex = Simplex(L)
                    else:
                        simplex = partial
                    add_simplex = True
                    simplices_to_delete = []
                    for already in answer:
                        if add_simplex:
                            if already.is_face(simplex):
                                add_simplex = False
                            if add_simplex and simplex.is_face(already):
                                simplices_to_delete.append(already)
                    if add_simplex:
                        answer.append(simplex)
                    for x in simplices_to_delete:
                        answer.remove(x)
        return answer

    def minimal_nonfaces(self):
        """
        Set consisting of the minimal subsets of the vertex set of
        this simplicial complex which do not form faces.

        Algorithm: Proceeds through the faces of the complex increasing the
        dimension, starting from dimension 0, and add the faces that are not
        contained in the complex and that are not already contained in a
        previously seen minimal non-face.

        This is used in computing the
        :meth:`Stanley-Reisner ring<stanley_reisner_ring>` and the
        :meth:`Alexander dual<alexander_dual>`.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,3],[1,2]])
            sage: X.minimal_nonfaces()
            {(2, 3)}
            sage: Y = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: sorted(Y.minimal_nonfaces())
            [(0, 2), (1, 3)]

        TESTS::

            sage: SC = SimplicialComplex([(0,1,2),(0,2,3),(2,3,4),(1,2,4), \
            ....:                         (1,4,5),(0,3,6),(3,6,7),(4,5,7)])

        This was taking a long time before :issue:`20078`::

            sage: sorted(SC.minimal_nonfaces())
            [(0, 4),
             (0, 5),
             (0, 7),
             (1, 3),
             (1, 6),
             (1, 7),
             (2, 5),
             (2, 6),
             (2, 7),
             (3, 4, 7),
             (3, 5),
             (4, 6),
             (5, 6)]
        """
        face_dict = self.faces()
        vertices = self.vertices()
        dimension = self.dimension()
        set_mnf = set()

        for dim in range(dimension + 1):
            face_sets = frozenset(f.set() for f in face_dict[dim])
            for candidate in combinations(vertices, dim + 1):
                set_candidate = frozenset(candidate)
                if set_candidate not in face_sets:
                    new = not any(set_candidate.issuperset(mnf) for mnf in set_mnf)
                    if new:
                        set_mnf.add(set_candidate)

        for candidate in combinations(vertices, dimension+2):  # Checks for minimal nonfaces in the remaining dimension
            set_candidate = frozenset(candidate)
            new = not any(set_candidate.issuperset(mnf) for mnf in set_mnf)
            if new:
                set_mnf.add(set_candidate)

        min_non_faces = Set([Simplex(mnf) for mnf in set_mnf])

        return min_non_faces

    def _stanley_reisner_base_ring(self, base_ring=ZZ):
        """
        The polynomial algebra of which the Stanley-Reisner ring is a
        quotient.

        INPUT:

        - ``base_ring`` -- a commutative ring (default: ``ZZ``)

        OUTPUT: a polynomial algebra with coefficients in base_ring,
        with one generator for each vertex in the simplicial complex.

        See the documentation for :meth:`stanley_reisner_ring` for a
        warning about the names of the vertices.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2], [0], [3]])
            sage: X._stanley_reisner_base_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring
            sage: Y = SimplicialComplex([['a', 'b', 'c']])
            sage: Y._stanley_reisner_base_ring(base_ring=QQ)
            Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        verts = self._gen_dict.values()
        try:
            verts = sorted(verts)
        except TypeError:
            verts = sorted(verts, key=str)
        return PolynomialRing(base_ring, verts)

    def stanley_reisner_ring(self, base_ring=ZZ):
        """
        The Stanley-Reisner ring of this simplicial complex.

        INPUT:

        - ``base_ring`` -- a commutative ring (default: ``ZZ``)

        OUTPUT: a quotient of a polynomial algebra with coefficients
        in ``base_ring``, with one generator for each vertex in the
        simplicial complex, by the ideal generated by the products
        of those vertices which do not form faces in it

        Thus the ideal is generated by the products corresponding to
        the minimal nonfaces of the simplicial complex.

        .. WARNING::

           This may be quite slow!

           Also, this may behave badly if the vertices have the
           'wrong' names. To avoid this, define the simplicial complex
           at the start with the flag ``name_check`` set to ``True``.

           More precisely, this is a quotient of a polynomial ring
           with one generator for each vertex.  If the name of a
           vertex is a nonnegative integer, then the corresponding
           polynomial generator is named ``'x'`` followed by that integer
           (e.g., ``'x2'``, ``'x3'``, ``'x5'``, ...).  Otherwise, the
           polynomial generators are given the same names as the vertices.
           Thus if the vertex set is ``(2, 'x2')``, there will be problems.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2], [0,2,3]])
            sage: X.stanley_reisner_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring
             by the ideal (x1*x3)
            sage: Y = SimplicialComplex([[0,1,2,3,4]]); Y
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 1, 2, 3, 4)}
            sage: Y.add_face([0,1,2,3,4])
            sage: Y.stanley_reisner_ring(base_ring=QQ)
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Rational Field
        """
        R = self._stanley_reisner_base_ring(base_ring)
        products = []
        for f in self.minimal_nonfaces():
            prod = 1
            for v in f:
                prod *= R(self._gen_dict[v])
            products.append(prod)
        return R.quotient(products)

    def alexander_dual(self, is_mutable=True):
        """
        The Alexander dual of this simplicial complex: according to
        the Macaulay2 documentation, this is the simplicial complex
        whose faces are the complements of its nonfaces.

        Thus find the minimal nonfaces and take their complements to
        find the facets in the Alexander dual.

        INPUT:

        - ``is_mutable`` -- boolean (default: ``True``); determine whether
          the output is mutable

        EXAMPLES::

            sage: Y = SimplicialComplex([[i] for i in range(5)]); Y
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and
             facets {(0,), (1,), (2,), (3,), (4,)}
            sage: Y.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and 10 facets
            sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: X.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2), (1, 3)}
        """
        nonfaces = self.minimal_nonfaces()
        return SimplicialComplex([self._complement(f) for f in nonfaces], is_mutable=is_mutable)

    def barycentric_subdivision(self):
        """
        The barycentric subdivision of this simplicial complex.

        See :wikipedia:`Barycentric_subdivision` for a
        definition.

        EXAMPLES::

            sage: triangle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: hexagon = triangle.barycentric_subdivision(); hexagon
            Simplicial complex with 6 vertices and 6 facets
            sage: hexagon.homology(1) == triangle.homology(1)                           # needs sage.modules
            True

        Barycentric subdivisions can get quite large, since each
        `n`-dimensional facet in the original complex produces
        `(n+1)!` facets in the subdivision::

            sage: S4 = simplicial_complexes.Sphere(4); S4
            Minimal triangulation of the 4-sphere
            sage: S4.barycentric_subdivision()
            Simplicial complex with 62 vertices and 720 facets
        """
        return self.face_poset().order_complex()

    def stellar_subdivision(self, simplex, inplace=False, is_mutable=True):
        """
        Return the stellar subdivision of a simplex in this simplicial complex.

        The stellar subdivision of a face is obtained by adding a new vertex to the
        simplicial complex ``self`` joined to the star of the face and then
        deleting the face ``simplex`` to the result.

        INPUT:

        - ``simplex`` -- a simplex face of ``self``
        - ``inplace`` -- boolean (default: ``False``); determines if the
          operation is done on ``self`` or on a copy
        - ``is_mutable`` -- boolean (default: ``True``); determines if the
          output is mutable

        OUTPUT:

        - A simplicial complex obtained by the stellar subdivision of the face
          ``simplex``

        EXAMPLES::

            sage: SC = SimplicialComplex([[0,1,2],[1,2,3]])
            sage: F1 = Simplex([1,2])
            sage: F2 = Simplex([1,3])
            sage: F3 = Simplex([1,2,3])
            sage: SC.stellar_subdivision(F1)
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and
             facets {(0, 1, 4), (0, 2, 4), (1, 3, 4), (2, 3, 4)}
            sage: SC.stellar_subdivision(F2)
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and
             facets {(0, 1, 2), (1, 2, 4), (2, 3, 4)}
            sage: SC.stellar_subdivision(F3)
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and
             facets {(0, 1, 2), (1, 2, 4), (1, 3, 4), (2, 3, 4)}
            sage: SC.stellar_subdivision(F3, inplace=True);SC
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and
             facets {(0, 1, 2), (1, 2, 4), (1, 3, 4), (2, 3, 4)}

        The simplex to subdivide should be a face of self::

            sage: SC = SimplicialComplex([[0,1,2],[1,2,3]])
            sage: F4 = Simplex([3,4])
            sage: SC.stellar_subdivision(F4)
            Traceback (most recent call last):
            ...
            ValueError: the face to subdivide is not a face of self

        One can not modify an immutable simplicial complex::

            sage: SC = SimplicialComplex([[0,1,2],[1,2,3]], is_mutable=False)
            sage: SC.stellar_subdivision(F1, inplace=True)
            Traceback (most recent call last):
            ...
            ValueError: this simplicial complex is not mutable
        """

        if inplace and self._is_immutable:
            raise ValueError("this simplicial complex is not mutable")

        if Simplex(simplex) not in self:
            raise ValueError("the face to subdivide is not a face of self")

        if inplace:
            working_complex = self
        else:
            working_complex = copy(self)

        vertices = working_complex.vertices()
        not_found = True
        vertex_label = 0
        while not_found:
            if vertex_label not in vertices:
                not_found = False
            else:
                vertex_label += 1
        new_vertex = SimplicialComplex([[vertex_label]])
        new_faces = new_vertex.join(working_complex.star(simplex), rename_vertices=False)
        for face in new_faces.facets():
            working_complex.add_face(face)

        working_complex.remove_face(simplex)

        if not is_mutable:
            working_complex.set_immutable()

        if not inplace:
            return working_complex

    def graph(self):
        """
        The 1-skeleton of this simplicial complex, as a graph.

        .. WARNING::

           This may give the wrong answer if the simplicial complex
           was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1,2,3]])
            sage: G = S.graph(); G
            Graph on 4 vertices
            sage: G.edges(sort=True)
            [(0, 1, None), (0, 2, None), (0, 3, None), (1, 2, None), (1, 3, None), (2, 3, None)]
        """
        if self._graph is None:
            from sage.graphs.graph import Graph

            edges = self.n_cells(1)
            vertices = [min(f) for f in self._facets if f.dimension() == 0]
            used_vertices = []  # vertices which are in an edge
            d = {}
            for e in edges:
                try:
                    v = min(e)
                    max_e = max(e)
                except TypeError:
                    v = min(e, key=str)
                    max_e = max(e, key=str)
                if v in d:
                    d[v].append(max_e)
                else:
                    d[v] = [max_e]
                used_vertices.extend(list(e))
            for v in vertices:
                if v not in used_vertices:
                    d[v] = []
            self._graph = Graph(d)
        return self._graph

    def delta_complex(self, sort_simplices=False):
        r"""
        Return ``self`` as a `\Delta`-complex.

        The `\Delta`-complex is essentially identical to the
        simplicial complex: it has same simplices with the same
        boundaries.

        INPUT:

        - ``sort_simplices`` -- boolean (default: ``False``); if ``True``, sort
          the list of simplices in each dimension

        EXAMPLES::

            sage: T = simplicial_complexes.Torus()
            sage: Td = T.delta_complex()
            sage: Td
            Delta complex with 7 vertices and 43 simplices
            sage: T.homology() == Td.homology()                                         # needs sage.modules
            True
        """
        from .delta_complex import DeltaComplex
        data = {}
        dim = self.dimension()
        n_cells = self._n_cells_sorted(dim)
        if sort_simplices:
            n_cells.sort()
        for n in range(dim, -1, -1):
            bdries = self._n_cells_sorted(n-1)
            if sort_simplices:
                bdries.sort()
            data[n] = []
            for f in n_cells:
                data[n].append([bdries.index(f.face(i)) for i in range(n+1)])
            n_cells = bdries
        return DeltaComplex(data)

    def is_flag_complex(self):
        """
        Return ``True`` if and only if ``self`` is a flag complex.

        A flag complex is a simplicial complex that is the largest simplicial
        complex on its 1-skeleton. Thus a flag complex is the clique complex
        of its graph.

        EXAMPLES::

            sage: h = Graph({0: [1,2,3,4], 1: [2,3,4], 2: [3]})
            sage: x = h.clique_complex(); x
            Simplicial complex with vertex set (0, 1, 2, 3, 4)
            and facets {(0, 1, 4), (0, 1, 2, 3)}
            sage: x.is_flag_complex()
            True

            sage: X = simplicial_complexes.ChessboardComplex(3,3)
            sage: X.is_flag_complex()
            True
        """
        return self == self.graph().clique_complex()

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this simplicial complex.

        The `n`-skeleton of a simplicial complex is obtained by discarding
        all of the simplices in dimensions larger than `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [1,2,3], [0,2,3]])
            sage: X.n_skeleton(1)
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)}
            sage: X.set_immutable()
            sage: X.n_skeleton(2)
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1), (0, 2, 3), (1, 2, 3)}
            sage: X.n_skeleton(4)
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1), (0, 2, 3), (1, 2, 3)}
        """
        if n >= self.dimension():
            return self
        # make sure it's a list (it will be a tuple if immutable)
        facets = [f for f in self._facets if f.dimension() < n]
        facets.extend(self.faces()[n])
        return SimplicialComplex(facets, is_immutable=self._is_immutable)

    def _contractible_subcomplex(self, verbose=False):
        """
        Find a contractible subcomplex `L` of this simplicial complex,
        preferably one which is as large as possible.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); if ``True``, print some
          messages as the simplicial complex is computed

        Motivation: if `K` is the original complex and if `L` is
        contractible, then the relative homology `H_*(K,L)` is
        isomorphic to the reduced homology of `K`.  If `L` is large,
        then the relative chain complex will be a good deal smaller
        than the augmented chain complex for `K`, and this leads to a
        speed improvement for computing the homology of `K`.

        This just passes an immutable subcomplex consisting of a facet to the
        method ``_enlarge_subcomplex``.

        .. NOTE::

           Thus when the simplicial complex is empty, so is the
           resulting 'contractible subcomplex', which is therefore not
           technically contractible.  In this case, that doesn't
           matter because the homology is computed correctly anyway.

        EXAMPLES::

            sage: sphere = SimplicialComplex([[0,1,2,3]])
            sage: sphere.remove_face([0,1,2,3])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)}
            sage: L = sphere._contractible_subcomplex(); L
            Simplicial complex with vertex set (0, 1, 2, 3) and
             facets {(0, 1, 2), (0, 1, 3), (0, 2, 3)}
            sage: L.homology()                                                          # needs sage.modules
            {0: 0, 1: 0, 2: 0}
        """
        facets = [sorted(self._facets, key=str)[0]]
        return self._enlarge_subcomplex(SimplicialComplex(facets, is_mutable=False), verbose=verbose)

    def _enlarge_subcomplex(self, subcomplex, verbose=False):
        """
        Given a subcomplex `S` of this simplicial complex `K`, find a
        subcomplex `L`, as large as possible, containing `S` which is
        homotopy equivalent to `S` (so that `H_{*}(K,S)` is isomorphic
        to `H_{*}(K,L)`).  This way, the chain complex for computing
        `H_{*}(K,L)` will be smaller than that for computing
        `H_{*}(K,S)`, so the computations should be faster.

        INPUT:

        - ``subcomplex`` -- a subcomplex of this simplicial complex
        - ``verbose`` -- boolean (default: ``False``); if ``True``, print some
          messages as the simplicial complex is computed

        OUTPUT: a complex `L` containing ``subcomplex`` and contained
        in ``self``, homotopy equivalent to ``subcomplex``

        Algorithm: start with the subcomplex `S` and loop through the
        facets of `K` which are not in `S`.  For each one, see whether
        its intersection with `S` is contractible, and if so, add it.
        This is recursive: testing for contractibility calls this
        routine again, via ``_contractible_subcomplex``.

        EXAMPLES::

            sage: T = simplicial_complexes.Torus(); T
            Minimal triangulation of the torus

        Inside the torus, define a subcomplex consisting of a loop::

            sage: S = SimplicialComplex([[0,1], [1,2], [0,2]], is_mutable=False)
            sage: S.homology()                                                          # needs sage.modules
            {0: 0, 1: Z}
            sage: L = T._enlarge_subcomplex(S)
            sage: L
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 8 facets
            sage: sorted(L.facets())
            [(0, 1), (0, 1, 5), (0, 2), (0, 2, 6), (0, 3, 4), (0, 3, 5), (0, 4, 6), (1, 2)]
            sage: L.homology()[1]                                                       # needs sage.modules
            Z
        """
        # Make the subcomplex immutable if not
        if subcomplex is not None and not subcomplex._is_immutable:
            subcomplex = SimplicialComplex(subcomplex._facets,
                                           maximality_check=False,
                                           is_mutable=False)

        if subcomplex in self.__enlarged:
            return self.__enlarged[subcomplex]
        faces = [x for x in list(self._facets) if x not in subcomplex._facets]
        # For consistency when using different Python versions, for example, sort 'faces'.
        faces = sorted(faces, key=str)
        done = False
        new_facets = sorted(subcomplex._facets, key=str)
        while not done:
            done = True
            remove_these = []
            if verbose:
                print(f"  looping through {len(faces)} facets")
            for f in faces:
                f_set = f.set()
                int_facets = {a.set().intersection(f_set) for a in new_facets}
                intersection = SimplicialComplex(int_facets)
                if not intersection._facets[0].is_empty():
                    if (len(intersection._facets) == 1 or
                            intersection == intersection._contractible_subcomplex()):
                        new_facets.append(f)
                        remove_these.append(f)
                        done = False
            if verbose and not done:
                print("    added %s facets" % len(remove_these))
            for f in remove_these:
                faces.remove(f)
        if verbose:
            print("  now constructing a simplicial complex with {} vertices and {} facets".format(len(self.vertices()), len(new_facets)))
        L = SimplicialComplex(new_facets, maximality_check=False,
                              is_immutable=self._is_immutable)
        self.__enlarged[subcomplex] = L
        # Use the same sorting on the vertices in L as in the ambient complex.
        L._vertex_to_index = self._vertex_to_index
        return L

    def _cubical_(self):
        r"""
        Cubical complex constructed from ``self``.

        ALGORITHM:

        The algorithm comes from a paper by Shtan'ko and Shtogrin, as
        reported by Bukhshtaber and Panov.  Let `I^m` denote the unit
        `m`-cube, viewed as a cubical complex.  Let `[m] = \{1, 2,
        ..., m\}`; then each face of `I^m` has the following form, for
        subsets `I \subset J \subset [m]`:

        .. MATH::

            F_{I \subset J} = \{ (y_1,...,y_m) \in I^m \,:\, y_i =0 \text{
            for } i \in I, y_j = 1 \text{ for } j \not \in J\}.

        If `K` is a simplicial complex on vertex set `[m]` and if `I
        \subset [m]`, write `I \in K` if `I` is a simplex of `K`.
        Then we associate to `K` the cubical subcomplex of `I^m` with
        faces

        .. MATH::

            \{F_{I \subset J} \,:\, J \in K, I \neq \emptyset \}

        The geometric realization of this cubical complex is
        homeomorphic to the geometric realization of the original
        simplicial complex.

        REFERENCES:

        - [BP2000]_
        - [SS1992]_

        EXAMPLES::

            sage: T = simplicial_complexes.Torus()
            sage: T.homology()                                                          # needs sage.modules
            {0: 0, 1: Z x Z, 2: Z}
            sage: Tc = T._cubical_()
            sage: Tc
            Cubical complex with 42 vertices and 168 cubes
            sage: Tc.homology()                                                         # needs sage.modules
            {0: 0, 1: Z x Z, 2: Z}
        """
        from .cubical_complex import CubicalComplex
        V = self.vertices()
        embed = len(V)
        # dictionary to translate vertices to the numbers 1, ..., embed
        vd = dict(zip(V, range(1, embed + 1)))
        cubes = []
        for JJ in self.facets():
            J = [vd[i] for i in JJ]
            for i in J:
                # loop over indices from 1 to embed.  if equal to i,
                # set to 0. if not in J, set to 1.  Otherwise, range
                # from 0 to 1
                cube = []
                for n in range(1, embed+1):
                    if n == i:
                        cube.append([0])
                    elif n not in J:
                        cube.append([1])
                    else:
                        cube.append([0, 1])
                cubes.append(cube)
        return CubicalComplex(cubes)

    def connected_component(self, simplex=None):
        """
        Return the connected component of this simplicial complex
        containing ``simplex``. If ``simplex`` is omitted, then return
        the connected component containing the zeroth vertex in the
        vertex list. (If the simplicial complex is empty, raise an
        error.)

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1 == S1.connected_component()
            True
            sage: X = S1.disjoint_union(S1)
            sage: X == X.connected_component()
            False
            sage: CL0 = X.connected_component(Simplex(['L0']))
            sage: CR0 = X.connected_component(Simplex(['R0']))
            sage: CL0 == CR0
            False

            sage: S0 = simplicial_complexes.Sphere(0)
            sage: S0.vertices()
            (0, 1)
            sage: S0.connected_component()
            Simplicial complex with vertex set (0,) and facets {(0,)}
            sage: S0.connected_component(Simplex((1,)))
            Simplicial complex with vertex set (1,) and facets {(1,)}

            sage: SimplicialComplex([[]]).connected_component()
            Traceback (most recent call last):
            ...
            ValueError: the empty simplicial complex has no connected components
        """
        if self.dimension() == -1:
            raise ValueError("the empty simplicial complex has no connected components")
        if simplex is None:
            v = self.vertices()[0]
        else:
            v = simplex[0]
        vertices = self.graph().connected_component_containing_vertex(v, sort=False)
        facets = [f for f in self.facets() if f.is_face(Simplex(vertices))]
        return SimplicialComplex(facets)

    def fundamental_group(self, base_point=None, simplify=True):
        r"""
        Return the fundamental group of this simplicial complex.

        INPUT:

        - ``base_point`` -- (default: ``None``) if this complex is
          not path-connected, then specify a vertex; the fundamental
          group is computed with that vertex as a base point. If the
          complex is path-connected, then you may specify a vertex or
          leave this as its default setting of ``None``. (If this
          complex is path-connected, then this argument is ignored.)

        - ``simplify`` -- boolean (default: ``True``); then return a
          presentation of the group in terms of generators and
          relations. If ``True``, the default, simplify as much as GAP is
          able to.

        Algorithm: we compute the edge-path group -- see
        :wikipedia:`Fundamental_group`. Choose a spanning tree for the
        1-skeleton, and then the group's generators are given by the
        edges in the 1-skeleton; there are two types of relations:
        `e=1` if `e` is in the spanning tree, and for every 2-simplex,
        if its edges are `e_0`, `e_1`, and `e_2`, then we impose the
        relation `e_0 e_1^{-1} e_2 = 1`.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S1.fundamental_group()                                                # needs sage.groups
            Finitely presented group < e |  >

        If we pass the argument ``simplify=False``, we get generators and
        relations in a form which is not usually very helpful. Here is the
        cyclic group of order 2, for instance::

            sage: RP2 = simplicial_complexes.RealProjectiveSpace(2)
            sage: C2 = RP2.fundamental_group(simplify=False); C2                        # needs sage.groups
            Finitely presented group < e0, e1, e2, e3, e4, e5, e6, e7, e8, e9 | e0, e3,
            e4, e7, e9, e5*e2^-1*e0, e7*e2^-1*e1, e8*e3^-1*e1, e8*e6^-1*e4, e9*e6^-1*e5 >
            sage: C2.simplified()                                                       # needs sage.groups
            Finitely presented group < e1 | e1^2 >

        This is the same answer given if the argument ``simplify`` is True
        (the default)::

            sage: RP2.fundamental_group()                                               # needs sage.groups
            Finitely presented group < e1 | e1^2 >

        You must specify a base point to compute the fundamental group
        of a non-connected complex::

            sage: # needs sage.groups
            sage: K = S1.disjoint_union(RP2)
            sage: K.fundamental_group()
            Traceback (most recent call last):
            ...
            ValueError: this complex is not connected, so you must specify a base point
            sage: K.fundamental_group(base_point='L0')
            Finitely presented group < e |  >
            sage: K.fundamental_group(base_point='R0').order()
            2

        Some other examples::

            sage: S1.wedge(S1).fundamental_group()                                      # needs sage.groups
            Finitely presented group < e0, e1 | >
            sage: simplicial_complexes.Torus().fundamental_group()                      # needs sage.groups
            Finitely presented group < e1, e4 | e4^-1*e1^-1*e4*e1 >

            sage: # needs sage.groups
            sage: G = simplicial_complexes.MooreSpace(5).fundamental_group()
            sage: G.ngens()
            1
            sage: x = G.gen(0)
            sage: [(x**n).is_one() for n in range(1,6)]
            [False, False, False, False, True]
        """
        if not self.is_connected():
            if base_point is None:
                raise ValueError("this complex is not connected, so you must specify a base point")
            return self.connected_component(Simplex([base_point])).fundamental_group(simplify=simplify)

        from sage.groups.free_group import FreeGroup
        from sage.libs.gap.libgap import libgap as gap
        G = self.graph()
        # Edges in the graph may be sorted differently than in the
        # simplicial complex, so convert the edges to frozensets so we
        # don't have to worry about it. Convert spanning_tree to a set
        # to make lookup faster.
        spanning_tree = {frozenset((u, v)) for u, v, _ in G.min_spanning_tree()}
        gens = [e for e in G.edge_iterator(labels=False)
                if frozenset(e) not in spanning_tree]
        if not gens:
            return gap.TrivialGroup()

        gens_dict = {frozenset(g): i for i, g in enumerate(gens)}
        FG = FreeGroup(len(gens), 'e')
        rels = []
        for f in self._n_cells_sorted(2):
            bdry = [tuple(e) for e in f.faces()]
            z = {}
            for i in range(3):
                x = frozenset(bdry[i])
                if x in spanning_tree:
                    z[i] = FG.one()
                else:
                    z[i] = FG.gen(gens_dict[x])
            rels.append(z[0]*z[1].inverse()*z[2])
        if simplify:
            return FG.quotient(rels).simplified()
        else:
            return FG.quotient(rels)

    def is_isomorphic(self, other, certificate=False):
        r"""
        Check whether two simplicial complexes are isomorphic.

        INPUT:

        - ``certificate`` -- if ``True``, then output is ``(a, b)``, where ``a``
          is a boolean and ``b`` is either a map or ``None``

        This is done by creating two graphs and checking whether they
        are isomorphic.

        EXAMPLES::

            sage: Z1 = SimplicialComplex([[0,1],[1,2],[2,3,4],[4,5]])
            sage: Z2 = SimplicialComplex([['a','b'],['b','c'],['c','d','e'],['e','f']])
            sage: Z3 = SimplicialComplex([[1,2,3]])
            sage: Z1.is_isomorphic(Z2)
            True
            sage: Z1.is_isomorphic(Z2, certificate=True)
            (True, {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f'})
            sage: Z3.is_isomorphic(Z2)
            False

        We check that :issue:`20751` is fixed::

            sage: C1 = SimplicialComplex([[1,2,3], [2,4], [3,5], [5,6]])
            sage: C2 = SimplicialComplex([['a','b','c'], ['b','d'], ['c','e'], ['e','f']])
            sage: C1.is_isomorphic(C2, certificate=True)
            (True, {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f'})
        """
        # Check easy invariants agree
        if (sorted(x.dimension() for x in self._facets)
            != sorted(x.dimension() for x in other._facets)
                or len(self.vertices()) != len(other.vertices())):
            return False

        from sage.graphs.graph import Graph

        g1 = Graph()
        g2 = Graph()
        # With Python 3, "is_isomorphic" for graphs works best if the
        # vertices and edges are sortable. So we translate them all to
        # ints and then if a certificate is needed, we translate
        # back at the end.
        self_to_int = {v: i for i, v in enumerate(list(self.vertices()) + list(self._facets))}
        other_to_int = {v: i for i, v in enumerate(list(other.vertices()) + list(other._facets))}
        g1.add_edges((self_to_int[v], self_to_int[f], "generic edge") for f in self._facets for v in f)
        g2.add_edges((other_to_int[v], other_to_int[f], "generic edge") for f in other._facets for v in f)
        fake = -1
        g1.add_edges((fake, self_to_int[v], "special_edge")
                     for v in self.vertices())
        g2.add_edges((fake, other_to_int[v], "special_edge")
                     for v in other.vertices())
        if not certificate:
            return g1.is_isomorphic(g2, edge_labels=True)
        isisom, tr = g1.is_isomorphic(g2, edge_labels=True, certificate=True)

        if isisom:
            for f in self.facets():
                tr.pop(self_to_int[f])
            tr.pop(fake)

        int_to_self = {idx: x for x, idx in self_to_int.items()}
        int_to_other = {idx: x for x, idx in other_to_int.items()}
        return isisom, {int_to_self[i]: int_to_other[tr[i]] for i in tr}

    def automorphism_group(self):
        r"""
        Return the automorphism group of the simplicial complex.

        This is done by creating a bipartite graph, whose vertices are
        vertices and facets of the simplicial complex, and computing
        its automorphism group.

        .. WARNING::

            Since :issue:`14319` the domain of the automorphism group is equal to
            the graph's vertex set, and the ``translation`` argument has become
            useless.

        EXAMPLES::

            sage: S = simplicial_complexes.Simplex(3)
            sage: S.automorphism_group().is_isomorphic(SymmetricGroup(4))               # needs sage.groups
            True

            sage: P = simplicial_complexes.RealProjectivePlane()
            sage: P.automorphism_group().is_isomorphic(AlternatingGroup(5))             # needs sage.groups
            True

            sage: Z = SimplicialComplex([['1','2'],['2','3','a']])
            sage: Z.automorphism_group().is_isomorphic(CyclicPermutationGroup(2))       # needs sage.groups
            True
            sage: group = Z.automorphism_group()                                        # needs sage.groups
            sage: sorted(group.domain())                                                # needs sage.groups
            ['1', '2', '3', 'a']

        Check that :issue:`17032` is fixed::

            sage: s = SimplicialComplex([[(0,1),(2,3)]])
            sage: s.automorphism_group().cardinality()                                  # needs sage.groups
            2
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.graphs.graph import Graph

        G = Graph()
        G.add_vertices(self.vertices())
        G.add_edges((f.tuple(), v) for f in self.facets() for v in f)
        group = G.automorphism_group(partition=[list(self.vertices()),
                                                [f.tuple()
                                                 for f in self.facets()]])

        gens = [[tuple(c) for c in g.cycle_tuples()
                 if c[0] in self.vertices()]
                for g in group.gens()]

        return PermutationGroup(gens=gens, domain=self.vertices())

    def fixed_complex(self, G):
        r"""
        Return the fixed simplicial complex `Fix(G)` for a subgroup `G`.

        INPUT:

        - ``G`` -- a subgroup of the automorphism group of the simplicial
          complex or a list of elements of the automorphism group

        OUTPUT:

        - a simplicial complex `Fix(G)`

        Vertices in `Fix(G)` are the orbits of `G` (acting on vertices
        of ``self``) that form a simplex in ``self``. More generally,
        simplices in `Fix(G)` correspond to simplices in ``self`` that
        are union of such orbits.

        A basic example::

            sage: S4 = simplicial_complexes.Sphere(4)
            sage: S3 = simplicial_complexes.Sphere(3)
            sage: fix = S4.fixed_complex([S4.automorphism_group()([(0,1)])]); fix       # needs sage.groups
            Simplicial complex with vertex set (0, 2, 3, 4, 5) and 5 facets
            sage: fix.is_isomorphic(S3)                                                 # needs sage.groups
            True

        Another simple example::

            sage: T = SimplicialComplex([[1,2,3],[2,3,4]])
            sage: G = T.automorphism_group()                                            # needs sage.groups
            sage: T.fixed_complex([G([(1,4)])])                                         # needs sage.groups
            Simplicial complex with vertex set (2, 3) and facets {(2, 3)}

        A more sophisticated example::

            sage: RP2 = simplicial_complexes.ProjectivePlane()
            sage: CP2 = simplicial_complexes.ComplexProjectivePlane()
            sage: G = CP2.automorphism_group()                                          # needs sage.groups
            sage: H = G.subgroup([G([(2,3),(5,6),(8,9)])])                              # needs sage.groups
            sage: CP2.fixed_complex(H).is_isomorphic(RP2)                               # needs sage.groups
            True
        """
        from sage.categories.groups import Groups
        if G in Groups():
            gens = G.gens()
        else:
            gens = G
            G = self.automorphism_group().subgroup(gens)

        invariant_f = [tuple(u) for u in self.face_iterator()
                       if all(sorted(sigma(j) for j in u) == sorted(u)
                              for sigma in gens)]
        new_verts = [min(o) for o in G.orbits() if o in invariant_f]
        return SimplicialComplex([[s for s in f if s in new_verts]
                                  for f in invariant_f])

    def _Hom_(self, other, category=None):
        """
        Return the set of simplicial maps between simplicial complexes
        ``self`` and ``other``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,T)  # indirect doctest
            sage: H
            Set of Morphisms from Minimal triangulation of the 1-sphere
             to Minimal triangulation of the 2-sphere
             in Category of finite simplicial complexes
            sage: f = {0:0,1:1,2:3}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism:
              From: Minimal triangulation of the 1-sphere
              To: Minimal triangulation of the 2-sphere
            Defn: 0 |--> 0
                  1 |--> 1
                  2 |--> 3

            sage: S._Hom_(T, Objects())
            Traceback (most recent call last):
            ...
            TypeError: Category of objects is not a subcategory of SimplicialComplexes()
            sage: type(Hom(S, T, Objects()))
            <class 'sage.categories.homset.Homset_with_category_with_equality_by_id'>
        """
        if not category.is_subcategory(SimplicialComplexes()):
            raise TypeError("{} is not a subcategory of SimplicialComplexes()".format(category))
        from sage.topology.simplicial_complex_homset import SimplicialComplexHomset
        return SimplicialComplexHomset(self, other)

    # @cached_method    when we switch to immutable SimplicialComplex
    def _is_numeric(self):
        """
        Test whether all vertices are labeled by integers.

        OUTPUT: boolean; whether all vertices are labeled by (not necessarily
        consecutive) integers

        EXAMPLES::

            sage: s = SimplicialComplex()
            sage: s._is_numeric()
            True
            sage: s.add_face(['a', 'b', 123])
            sage: s._is_numeric()
            False
        """
        return all(isinstance(v, (int, Integer))
                   for v in self.vertices())

    # @cached_method    when we switch to immutable SimplicialComplex
    def _translation_to_numeric(self):
        """
        Return a dictionary enumerating the vertices.

        See also :meth:`_translation_from_numeric`, which returns the
        inverse map.

        OUTPUT:

        A dictionary. The keys are the vertices, and the associated
        values are integers from 0 to number of vertices - 1.

        EXAMPLES::

            sage: s = SimplicialComplex()
            sage: s._translation_to_numeric()
            {}
            sage: s.add_face(['a', 'b', 123])
            sage: s._translation_to_numeric()   # random output
            {'a': 1, 123: 0, 'b': 2}
            sage: set(s._translation_to_numeric().keys()) == set(['a', 'b', 123])
            True
            sage: sorted(s._translation_to_numeric().values())
            [0, 1, 2]
        """
        return self._vertex_to_index

    # @cached_method    when we switch to immutable SimplicialComplex
    def _translation_from_numeric(self):
        """
        Return a dictionary mapping vertex indices to vertices.

        See also :meth:`_translation_to_numeric`, which returns the
        inverse map.

        OUTPUT:

        A dictionary. The keys are integers from 0 to the number of
        vertices - 1. The associated values are the vertices.

        EXAMPLES::

            sage: s = SimplicialComplex()
            sage: s._translation_from_numeric()
            {}
            sage: s.add_face(['a', 'b', 123])
            sage: s._translation_from_numeric()   # random output
            {0: 123, 1: 'a', 2: 'b'}
            sage: sorted(s._translation_from_numeric().keys())
            [0, 1, 2]
            sage: set(s._translation_from_numeric().values()) == set(['a', 'b', 123])
            True
        """
        d = self._vertex_to_index
        return {idx: v for v, idx in d.items()}

    # this function overrides the standard one for GenericCellComplex,
    # because it lists the maximal faces, not the total number of faces.
    def _repr_(self):
        """
        Print representation.

        If there are only a few vertices or faces, they are listed. If
        there are lots, the number is given.

        Facets are sorted in increasing order of dimension, and within
        each dimension, they are sorted using the underlying tuple.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [1,2]])
            sage: X._repr_()
            'Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1), (1, 2)}'
            sage: SimplicialComplex([[i for i in range(16)]])
            Simplicial complex with 16 vertices and 1 facets
        """
        vertex_limit = 45
        facet_limit = 55
        vertices = self.vertices()
        try:
            vertices = sorted(vertices)
        except TypeError:
            vertices = sorted(vertices, key=str)
        try:
            facets = sorted(self._facets, key=lambda f: (f.dimension(), f.tuple()))
        except TypeError:
            # Sorting failed.
            facets = self._facets

        vertex_string = "with vertex set {}".format(tuple(vertices))
        if len(vertex_string) > vertex_limit:
            vertex_string = "with %s vertices" % len(vertices)
        facet_string = 'facets {' + repr(facets)[1:-1] + '}'
        if len(facet_string) > facet_limit:
            facet_string = "%s facets" % len(facets)
        return "Simplicial complex " + vertex_string + " and " + facet_string

    def set_immutable(self):
        """
        Make this simplicial complex immutable.

        EXAMPLES::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: S.is_mutable()
            True
            sage: S.set_immutable()
            sage: S.is_mutable()
            False
        """
        self._is_immutable = True
        self._facets = tuple(self._facets)

    def is_mutable(self):
        """
        Return ``True`` if mutable.

        EXAMPLES::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: S.is_mutable()
            True
            sage: S.set_immutable()
            sage: S.is_mutable()
            False
            sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S2.is_mutable()
            False
            sage: S3 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S3.is_mutable()
            False
        """
        return not self._is_immutable

    def is_immutable(self):
        """
        Return ``True`` if immutable.

        EXAMPLES::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: S.is_immutable()
            False
            sage: S.set_immutable()
            sage: S.is_immutable()
            True
        """
        return self._is_immutable

    def cone_vertices(self):
        r"""
        Return the list of cone vertices of ``self``.

        A vertex is a cone vertex if and only if it appears in every facet.

        EXAMPLES::

            sage: SimplicialComplex([[1,2,3]]).cone_vertices()
            [1, 2, 3]
            sage: SimplicialComplex([[1,2,3], [1,3,4], [1,5,6]]).cone_vertices()
            [1]
            sage: SimplicialComplex([[1,2,3], [1,3,4], [2,5,6]]).cone_vertices()
            []
        """
        F = self.facets()
        C = set(self.vertices())
        for f in F:
            C = C.intersection(list(f))
            if not C:
                break
        return sorted(C)

    def decone(self):
        r"""
        Return the subcomplex of ``self`` induced by the non-cone vertices.

        EXAMPLES::

            sage: SimplicialComplex([[1,2,3]]).decone()
            Simplicial complex with vertex set () and facets {()}
            sage: SimplicialComplex([[1,2,3], [1,3,4], [1,5,6]]).decone()
            Simplicial complex with vertex set (2, 3, 4, 5, 6)
             and facets {(2, 3), (3, 4), (5, 6)}
            sage: X = SimplicialComplex([[1,2,3], [1,3,4], [2,5,6]])
            sage: X.decone() == X
            True
        """
        V = set(self.vertices()).difference(self.cone_vertices())
        return self.generated_subcomplex(V)

    def is_balanced(self, check_purity=False, certificate=False):
        r"""
        Determine whether ``self`` is balanced.

        A simplicial complex `X` of dimension `d-1` is balanced if and
        only if its vertices can be colored with `d` colors such that
        every face contains at most one vertex of each color.  An
        equivalent condition is that the 1-skeleton of `X` is
        `d`-colorable.  In some contexts, it is also required that `X`
        be pure (i.e., that all maximal faces of `X` have the same
        dimension).

        INPUT:

        - ``check_purity`` -- boolean (default: ``False``); if this is ``True``,
          require that ``self`` be pure as well as balanced

        - ``certificate`` -- boolean (default: ``False``); if this is ``True`` and
          ``self`` is balanced, then return a `d`-coloring of the 1-skeleton

        EXAMPLES:

        A 1-dim simplicial complex is balanced iff it is bipartite::

            sage: X = SimplicialComplex([[1,2], [1,4], [3,4], [2,5]])
            sage: X.is_balanced()
            True
            sage: sorted(X.is_balanced(certificate=True))
            [[1, 3, 5], [2, 4]]
            sage: X = SimplicialComplex([[1,2], [1,4], [3,4], [2,4]])
            sage: X.is_balanced()
            False

        Any barycentric division is balanced::

            sage: X = SimplicialComplex([[1,2,3], [1,2,4], [2,3,4]])
            sage: X.is_balanced()
            False
            sage: X.barycentric_subdivision().is_balanced()
            True

        A non-pure balanced complex::

            sage: X = SimplicialComplex([[1,2,3], [3,4]])
            sage: X.is_balanced(check_purity=True)
            False
            sage: sorted(X.is_balanced(certificate=True))
            [[1, 4], [2], [3]]
        """
        d = 1 + self.dimension()
        if check_purity and not self.is_pure():
            return False
        Skel = self.graph()
        if certificate:
            C = Skel.coloring()
            C = C if len(C) == d else False
            return C
        else:
            return Skel.chromatic_number() == d

    def is_partitionable(self, certificate=False,
                         *, solver=None, integrality_tolerance=1e-3):
        r"""
        Determine whether ``self`` is partitionable.

        A partitioning of a simplicial complex `X` is a decomposition
        of its face poset into disjoint Boolean intervals `[R,F]`,
        where `F` ranges over all facets of `X`.

        The method sets up an integer program with:

        - a variable `y_i` for each pair `(R,F)`, where `F` is a facet of `X`
          and `R` is a subface of `F`

        - a constraint `y_i+y_j \leq 1` for each pair `(R_i,F_i)`, `(R_j,F_j)`
          whose Boolean intervals intersect nontrivially (equivalent to
          `(R_i\subseteq F_j and R_j\subseteq F_i))`

        - objective function equal to the sum of all `y_i`

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); if ``True``,
          and ``self`` is partitionable, then return a list of pairs `(R,F)`
          that form a partitioning.

        - ``solver`` -- (default: ``None``) specifies a Mixed Integer Linear Programming
          (MILP) solver to be used. If set to ``None``, the default one is used. For
          more information on MILP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
          inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

        EXAMPLES:

        Simplices are trivially partitionable::

            sage: X = SimplicialComplex([[1,2,3,4]])
            sage: X.is_partitionable()                                                  # needs sage.numerical.mip
            True
            sage: X.is_partitionable(certificate=True)                                  # needs sage.numerical.mip
            [((), (1, 2, 3, 4), 4)]

        Shellable complexes are partitionable::

            sage: # needs sage.numerical.mip, long time
            sage: X = SimplicialComplex([[1,3,5], [1,3,6], [1,4,5], [1,4,6],
            ....:                        [2,3,5], [2,3,6], [2,4,5]])
            sage: X.is_partitionable()
            True
            sage: P = X.is_partitionable(certificate=True)
            sage: def n_intervals_containing(f):
            ....:     return len([RF for RF in P
            ....:                    if RF[0].is_face(f) and f.is_face(RF[1])])
            sage: all(n_intervals_containing(f) == 1
            ....:     for k in X.faces().keys() for f in X.faces()[k])
            True

        A non-shellable, non-Cohen-Macaulay, partitionable example, constructed by Björner::

            sage: X = SimplicialComplex([[1,2,3], [1,2,4], [1,3,4], [2,3,4], [1,5,6]])
            sage: X.is_partitionable()                                                  # needs sage.numerical.mip
            True

        The bowtie complex is not partitionable::

            sage: X = SimplicialComplex([[1,2,3], [1,4,5]])
            sage: X.is_partitionable()                                                  # needs sage.numerical.mip
            False
        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        RFPairs = [(Simplex(r), f, f.dimension() - len(r) + 1)
                   for f in self.facets() for r in Set(f).subsets()]
        n = len(RFPairs)
        IP = MixedIntegerLinearProgram(solver=solver)
        y = IP.new_variable(binary=True)
        for i0, pair0 in enumerate(RFPairs):
            for i1, pair1 in enumerate(RFPairs):
                if (i0 < i1 and pair0[0].is_face(pair1[1]) and
                        pair1[0].is_face(pair0[1])):
                    IP.add_constraint(y[i0] + y[i1] <= 1)
        IP.set_objective(sum(2**RFPairs[i][2] * y[i] for i in range(n)))
        sol = round(IP.solve())
        if sol < sum(self.f_vector()):
            return False
        elif not certificate:
            return True
        else:
            x = IP.get_values(y, convert=bool, tolerance=integrality_tolerance)
            return [RFPairs[i] for i in range(n) if x[i]]

    def intersection(self, other):
        r"""
        Calculate the intersection of two simplicial complexes.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2,3], [1,2,4]])
            sage: Y = SimplicialComplex([[1,2,3], [1,4,5]])
            sage: Z = SimplicialComplex([[1,2,3], [1,4], [2,4]])
            sage: sorted(X.intersection(Y).facets())
            [(1, 2, 3), (1, 4)]
            sage: X.intersection(X) == X
            True
            sage: X.intersection(Z) == X
            False
            sage: X.intersection(Z) == Z
            True
        """
        F = []
        for k in range(1 + min(self.dimension(), other.dimension())):
            F = F + [s for s in self.faces()[k] if s in other.faces()[k]]
        return SimplicialComplex(F)

    def bigraded_betti_numbers(self, base_ring=ZZ, verbose=False):
        r"""
        Return a dictionary of the bigraded Betti numbers of ``self``,
        with keys `(-a, 2b)`.

        INPUT:

        - ``base_ring`` -- (default: ``ZZ``) the base ring used
          when computing homology
        - ``verbose`` -- boolean (default: ``False``); if ``True``, print
          messages during the computation, which indicate in which
          subcomplexes non-trivial homologies appear

        .. NOTE::

            If ``verbose`` is ``True``, then caching is avoided.

        .. SEEALSO::

            See :meth:`bigraded_betti_number` for more information.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1],[1,2],[1,3],[2,3]])
            sage: Y = SimplicialComplex([[1,2,3],[1,2,4],[3,5],[4,5]])
            sage: sorted(X.bigraded_betti_numbers(base_ring=QQ).items(), reverse=True)
            [((0, 0), 1), ((-1, 6), 1), ((-1, 4), 2), ((-2, 8), 1), ((-2, 6), 1)]
            sage: sorted(Y.bigraded_betti_numbers(verbose=True).items(), reverse=True)
            (-1, 4): Non-trivial homology Z in dimension 0 of the full
            subcomplex generated by a set of vertices (1, 5)
            (-1, 4): Non-trivial homology Z in dimension 0 of the full
            subcomplex generated by a set of vertices (2, 5)
            (-1, 4): Non-trivial homology Z in dimension 0 of the full
            subcomplex generated by a set of vertices (3, 4)
            (-2, 6): Non-trivial homology Z in dimension 0 of the full
            subcomplex generated by a set of vertices (1, 2, 5)
            (-2, 8): Non-trivial homology Z in dimension 1 of the full
            subcomplex generated by a set of vertices (1, 3, 4, 5)
            (-2, 8): Non-trivial homology Z in dimension 1 of the full
            subcomplex generated by a set of vertices (2, 3, 4, 5)
            (-3, 10): Non-trivial homology Z in dimension 1 of the full
            subcomplex generated by a set of vertices (1, 2, 3, 4, 5)
            [((0, 0), 1), ((-1, 4), 3), ((-2, 8), 2), ((-2, 6), 1), ((-3, 10), 1)]

        If we wish to view them in a form of a table, it is
        simple enough to create a function as such::

            sage: def print_table(bbns):
            ....:     max_a = max(-p[0] for p in bbns)
            ....:     max_b = max(p[1] for p in bbns)
            ....:     bbn_table = [[bbns.get((-a,b), 0) for a in range(max_a+1)]
            ....:                                       for b in range(max_b+1)]
            ....:     width = len(str(max(bbns.values()))) + 1
            ....:     print(' '*width, end=' ')
            ....:     for i in range(max_a+1):
            ....:         print(f'{-i:{width}}', end=' ')
            ....:     print()
            ....:     for j in range(len(bbn_table)):
            ....:         print(f'{j:{width}}', end=' ')
            ....:         for r in bbn_table[j]:
            ....:             print(f'{r:{width}}', end=' ')
            ....:         print()
            sage: print_table(X.bigraded_betti_numbers())
                0 -1 -2
             0  1  0  0
             1  0  0  0
             2  0  0  0
             3  0  0  0
             4  0  2  0
             5  0  0  0
             6  0  1  1
             7  0  0  0
             8  0  0  1
        """
        if base_ring in self._bbn_all_computed and not verbose:
            return self._bbn[base_ring]

        from sage.homology.homology_group import HomologyGroup
        L = self.vertices()
        n = len(L)
        B = {}
        H0 = HomologyGroup(0, base_ring)

        B[(0, 0)] = ZZ.one()

        for j in range(n+1):
            for x in combinations(L, j):
                S = self.generated_subcomplex(x)
                H = S.homology(base_ring=base_ring)
                for k in range(j):
                    if j-k-1 in H and H[j-k-1] != H0:
                        ind = (-k, 2*j)
                        if ind not in B:
                            B[ind] = ZZ.zero()
                        B[ind] += len(H[j-k-1].gens())
                        if verbose:
                            print("{}: Non-trivial homology {} in dimension {} of the full subcomplex generated by a set of vertices {}".format(ind, H[j-k-1], j-k-1, x))

        self._bbn[base_ring] = B
        self._bbn_all_computed.add(base_ring)

        return B

    def bigraded_betti_number(self, a, b, base_ring=ZZ, verbose=False):
        r"""
        Return the bigraded Betti number indexed in the form `(-a, 2b)`.

        Bigraded Betti number with indices `(-a, 2b)` is defined as a
        sum of ranks of `(b-a-1)`-th (co)homologies of full subcomplexes
        with exactly `b` vertices.

        INPUT:

        - ``base_ring`` -- (default: ``ZZ``) the base ring used
          when computing homology
        - ``verbose`` -- boolean (default: ``False``); if ``True``, print
          messages during the computation, which indicate in which
          subcomplexes non-trivial homologies appear

        .. NOTE::

            If ``verbose`` is ``True``, then caching is avoided.

        EXAMPLES::

            sage: # needs sage.modules
            sage: X = SimplicialComplex([[0,1],[1,2],[2,0],[1,3]])
            sage: X.bigraded_betti_number(-1, 4, base_ring=QQ)
            2
            sage: X.bigraded_betti_number(-1, 8)
            0
            sage: X.bigraded_betti_number(-2, 5)
            0
            sage: X.bigraded_betti_number(0, 0)
            1
            sage: sorted(X.bigraded_betti_numbers().items(), reverse=True)
            [((0, 0), 1), ((-1, 6), 1), ((-1, 4), 2), ((-2, 8), 1), ((-2, 6), 1)]
            sage: X.bigraded_betti_number(-1, 4, base_ring=QQ)
            2
            sage: X.bigraded_betti_number(-1, 8)
            0
            sage: Y = SimplicialComplex([[1,2,3],[1,2,4],[3,5],[4,5]])
            sage: Y.bigraded_betti_number(-1, 4, verbose=True)
            Non-trivial homology Z in dimension 0 of the full subcomplex
            generated by a set of vertices (1, 5)
            Non-trivial homology Z in dimension 0 of the full subcomplex
            generated by a set of vertices (2, 5)
            Non-trivial homology Z in dimension 0 of the full subcomplex
            generated by a set of vertices (3, 4)
            3
        """
        if b % 2:
            return ZZ.zero()
        if a == 0 == b:
            return ZZ.one()
        if base_ring in self._bbn and not verbose:
            if base_ring in self._bbn_all_computed:
                return self._bbn[base_ring].get((a, b), ZZ.zero())
            elif (a, b) in self._bbn[base_ring]:
                return self._bbn[base_ring][a, b]

        from sage.homology.homology_group import HomologyGroup

        b //= 2
        L = self.vertices()
        H0 = HomologyGroup(0, base_ring)

        B = 0

        for x in combinations(L, b):
            S = self.generated_subcomplex(x)
            H = S.homology(base_ring=base_ring)
            if b+a-1 in H and H[b+a-1] != H0:
                B += len(H[b+a-1].gens())
                if verbose:
                    print("Non-trivial homology {} in dimension {} of the full subcomplex generated by a set of vertices {}".format(H[b+a-1], b+a-1, x))

        B = ZZ(B)

        if base_ring in self._bbn:
            self._bbn[base_ring][(a, 2*b)] = B
        else:
            self._bbn[base_ring] = {(a, 2*b): B}

        return B

    def is_golod(self) -> bool:
        r"""
        Return whether ``self`` is Golod.

        A simplicial complex is Golod if multiplication and all higher
        Massey operations in the associated Tor-algebra are trivial. This
        is done by checking the bigraded Betti numbers.

        EXAMPLES::

            sage: # needs sage.modules
            sage: X = SimplicialComplex([[0,1],[1,2],[2,3],[3,0]])
            sage: Y = SimplicialComplex([[0,1,2],[0,2],[0,4]])
            sage: X.is_golod()
            False
            sage: Y.is_golod()
            True
        """
        H = [a+b for a, b in self.bigraded_betti_numbers()]
        if 0 in H:
            H.remove(0)

        return not any(i+j in H for ii, i in enumerate(H) for j in H[ii:])

    def is_minimally_non_golod(self) -> bool:
        r"""
        Return whether ``self`` is minimally non-Golod.

        If a simplicial complex itself is not Golod, but deleting any vertex
        gives us a full subcomplex that is Golod, then we say that a simplicial
        complex is minimally non-Golod.

        .. SEEALSO::

            See :meth:`is_golod` for more information.

        EXAMPLES::

            sage: # needs sage.modules
            sage: X = SimplicialComplex([[0,1],[1,2],[2,3],[3,0]])
            sage: Y = SimplicialComplex([[1,2,3],[1,2,4],[3,5],[4,5]])
            sage: X.is_golod()
            False
            sage: X.is_minimally_non_golod()
            True
            sage: Y.is_golod()
            False
            sage: Y.is_minimally_non_golod()
            False
        """
        def test(v):
            X = copy(self)
            X.remove_face([v])
            return X.is_golod()

        return not self.is_golod() and all(test(v) for v in self.vertices())

    def moment_angle_complex(self):
        """
        Return the moment-angle complex of ``self``.

        A moment-angle complex is a topological space created
        from this simplicial complex, which holds a lot of
        information about the simplicial complex itself.

        .. SEEALSO::

            See :mod:`sage.topology.moment_angle_complex` for
            more information on moment-angle complexes.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2,3], [1,4], [3,2,4]])
            sage: X.moment_angle_complex()
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3, 4) and facets {(1, 4), (2, 3, 4), (0, 1, 2, 3)}
            sage: K = simplicial_complexes.KleinBottle()
            sage: K.moment_angle_complex()
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3, 4, 5, 6, 7) and 16 facets

        We can also create it explicitly::

            sage: Z = MomentAngleComplex(K); Z
            Moment-angle complex of Simplicial complex with vertex set
            (0, 1, 2, 3, 4, 5, 6, 7) and 16 facets
        """
        from .moment_angle_complex import MomentAngleComplex
        return MomentAngleComplex(self)


# Miscellaneous utility functions.

# The following two functions can be used to generate the facets for
# the corresponding examples in sage.homology.examples. These take a
# few seconds to run, so the actual examples have the facets
# hard-coded. Thus the following functions are not currently used in
# the Sage library.

def facets_for_RP4():
    """
    Return the list of facets for a minimal triangulation of 4-dimensional
    real projective space.

    We use vertices numbered 1 through 16, define two facets, and define
    a certain subgroup `G` of the symmetric group `S_{16}`. Then the set
    of all facets is the `G`-orbit of the two given facets.

    See the description in Example 3.12 in Datta [Dat2007]_.

    EXAMPLES::

        sage: from sage.topology.simplicial_complex import facets_for_RP4
        sage: A = facets_for_RP4()   # long time (1 or 2 seconds)
        sage: SimplicialComplex(A) == simplicial_complexes.RealProjectiveSpace(4) # long time
        True
    """
    # Define the group:
    from sage.groups.perm_gps.permgroup import PermutationGroup
    g1 = '(2,7)(4,10)(5,6)(11,12)'
    g2 = '(1, 2, 3, 4, 5, 10)(6, 8, 9)(11, 12, 13, 14, 15, 16)'
    G = PermutationGroup([g1, g2])
    # Define the two simplices:
    t1 = (1, 2, 4, 5, 11)
    t2 = (1, 2, 4, 11, 13)
    # Apply the group elements to the simplices:
    facets = []
    for g in G:
        d = g.dict()
        for t in [t1, t2]:
            new = tuple([d[j] for j in t])
            if new not in facets:
                facets.append(new)
    return facets


def facets_for_K3():
    """
    Return the facets for a minimal triangulation of the K3 surface.

    This is a pure simplicial complex of dimension 4 with 16
    vertices and 288 facets. The facets are obtained by constructing a
    few facets and a permutation group `G`, and then computing the
    `G`-orbit of those facets.

    See Casella and Kühnel in [CK2001]_ and Spreer and Kühnel [SK2011]_;
    the construction here uses the labeling from Spreer and Kühnel.

    EXAMPLES::

        sage: from sage.topology.simplicial_complex import facets_for_K3
        sage: A = facets_for_K3()   # long time (a few seconds)
        sage: SimplicialComplex(A) == simplicial_complexes.K3Surface()  # long time
        True
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup
    G = PermutationGroup([[(1, 3, 8, 4, 9, 16, 15, 2, 14, 12, 6, 7, 13, 5, 10)],
                          [(1, 11, 16), (2, 10, 14), (3, 12, 13),
                           (4, 9, 15), (5, 7, 8)]])
    return ([tuple([g(i) for i in (1, 2, 3, 8, 12)]) for g in G]
            + [tuple([g(i) for i in (1, 2, 5, 8, 14)]) for g in G])
