r"""
Finite simplicial complexes

AUTHORS:

- John H. Palmieri (2009-04)

- D. Benjamin Antieau (2009-06) - added is_connected, generated_subcomplex,
  remove_facet, and is_flag_complex methods;
  cached the output of the graph() method.

- Travis Scrimshaw (2012-08-17): Made :class:`SimplicialComplex` have an
  immutable option, and added ``__hash__()`` function which checks to make
  sure it is immutable. Made :meth:`SimplicialComplex.remove_face()` into a
  mutator. Deprecated the ``vertex_set`` parameter.

- Christian Stump (2011-06) - implementation of is_cohen_macaulay

- Travis Scrimshaw (2013-02-16): Allowed :class:`SimplicialComplex` to make
  mutable copies.

This module implements the basic structure of finite simplicial
complexes. Given a set `V` of "vertices", a simplicial complex on `V`
is a collection `K` of subsets of `V` satisfying the condition that if
`S` is one of the subsets in `K`, then so is every subset of `S`.  The
subsets `S` are called the 'simplices' of `K`.

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
`n^{th}` term is the free `R`-module with basis given by the
`n`-simplices of `K`.  The differential is determined by its value on
any simplex: on the `n`-simplex with vertices `(v_0, v_1, ..., v_n)`,
the differential is the alternating sum with `i^{th}` summand `(-1)^i`
multiplied by the `(n-1)`-simplex obtained by omitting vertex `v_i`.

In the implementation here, the vertex set must be finite. To define a
simplicial complex, specify its vertex set: this should be a list,
tuple, or set, or it can be a non-negative integer `n`, in which case
the vertex set is `(0, ..., n)`.  Also specify the facets: the maximal
faces.

.. NOTE::

   The elements of the vertex set are not automatically contained in
   the simplicial complex: each one is only included if and only if it
   is a vertex of at least one of the specified facets.

.. NOTE::

   This class derives from
   :class:`~sage.homology.cell_complex.GenericCellComplex`, and so
   inherits its methods.  Some of those methods are not listed here;
   see the :mod:`Generic Cell Complex <sage.homology.cell_complex>`
   page instead.

EXAMPLES::

    sage: SimplicialComplex([[1], [3, 7]])
    Simplicial complex with vertex set (1, 3, 7) and facets {(3, 7), (1,)}
    sage: SimplicialComplex()   # the empty simplicial complex
    Simplicial complex with vertex set () and facets {()}
    sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
    sage: X
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2), (2, 3), (0, 3), (0, 1)}
    sage: X.stanley_reisner_ring()
    Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)
    sage: X.is_pure()
    True

Sage can perform a number of operations on simplicial complexes, such
as the join and the product, and it can also compute homology::

    sage: S = SimplicialComplex([[0,1], [1,2], [0,2]]) # circle
    sage: T = S.product(S)  # torus
    sage: T
    Simplicial complex with 9 vertices and 18 facets
    sage: T.homology()   # this computes reduced homology
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
    sage: X.stanley_reisner_ring()
    Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)

Mutability (see :trac:`12587`)::

    sage: S = SimplicialComplex([[1,4], [2,4]])
    sage: S.add_face([1,3])
    sage: S.remove_face([1,3]); S
    Simplicial complex with vertex set (1, 2, 3, 4) and facets {(2, 4), (1, 4), (3,)}
    sage: hash(S)
    Traceback (most recent call last):
    ...
    ValueError: This simplicial complex must be immutable. Call set_immutable().
    sage: S = SimplicialComplex([[1,4], [2,4]])
    sage: S.set_immutable()
    sage: S.add_face([1,3])
    Traceback (most recent call last):
    ...
    ValueError: This simplicial complex is not mutable
    sage: S.remove_face([1,3])
    Traceback (most recent call last):
    ...
    ValueError: This simplicial complex is not mutable
    sage: hash(S) == hash(S)
    True

    sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
    sage: hash(S2) == hash(S)
    True

We can also make mutable copies of an immutable simplicial complex
(see :trac:`14142`)::

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
from sage.homology.cell_complex import GenericCellComplex
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.integer_ring import ZZ
from sage.structure.parent_gens import normalize_names
from sage.misc.latex import latex
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex
from sage.graphs.graph import Graph

def lattice_paths(t1, t2, length=None):
    """
    Given lists (or tuples or ...) ``t1`` and ``t2``, think of them as
    labelings for vertices: ``t1`` labeling points on the x-axis,
    ``t2`` labeling points on the y-axis, both increasing.  Return the
    list of rectilinear paths along the grid defined by these points
    in the plane, starting from ``(t1[0], t2[0])``, ending at
    ``(t1[last], t2[last])``, and at each grid point, going either
    right or up.  See the examples.

    :param t1: labeling for vertices
    :param t2: labeling for vertices
    :param length: if not ``None``, then an integer, the length of the desired
        path.
    :type length: integer or ``None``; optional, default ``None``
    :type t1: tuple, list, other iterable
    :type t2: tuple, list, other iterable
    :return: list of lists of vertices making up the paths as described above
    :rtype: list of lists

    This is used when triangulating the product of simplices.  The
    optional argument ``length`` is used for `\Delta`-complexes, to
    specify all simplices in a product: in the triangulation of a
    product of two simplices, there is a `d`-simplex for every path of
    length `d+1` in the lattice.  The path must start at the bottom
    left and end at the upper right, and it must use at least one
    point in each row and in each column, so if ``length`` is too
    small, there will be no paths.

    EXAMPLES::

        sage: from sage.homology.simplicial_complex import lattice_paths
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

def rename_vertex(n, keep, left = True):
    """
    Rename a vertex: the vertices from the list ``keep`` get
    relabeled 0, 1, 2, ..., in order.  Any other vertex (e.g. 4) gets
    renamed to by prepending an 'L' or an 'R' (thus to either 'L4' or
    'R4'), depending on whether the argument left is ``True`` or ``False``.

    :param n: a 'vertex': either an integer or a string
    :param keep: a list of three vertices
    :param left: if ``True``, rename for use in left factor
    :type left: boolean; optional, default ``True``

    This is used by the :meth:`~SimplicialComplex.connected_sum` method for
    simplicial complexes.

    EXAMPLES::

        sage: from sage.homology.simplicial_complex import rename_vertex
        sage: rename_vertex(6, [5, 6, 7])
        1
        sage: rename_vertex(3, [5, 6, 7, 8, 9])
        'L3'
        sage: rename_vertex(3, [5, 6, 7], left=False)
        'R3'
    """
    lookup = dict(zip(keep, range(len(keep))))
    try:
        return lookup[n]
    except KeyError:
        if left:
            return "L" + str(n)
        else:
            return "R" + str(n)

class Simplex(SageObject):
    """
    Define a simplex.

    Topologically, a simplex is the convex hull of a collection of
    vertices in general position.  Combinatorially, it is defined just
    by specifying a set of vertices.  It is represented in Sage by the
    tuple of the vertices.

    :param X: set of vertices
    :type X: integer or list, tuple, or other iterable
    :return: simplex with those vertices

    ``X`` may be a non-negative integer `n`, in which case the
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
        """
        try:
            N = int(X) + 1
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
            <type 'tuple'>
            sage: type(Simplex(3))
            <class 'sage.homology.simplicial_complex.Simplex'>
        """
        return self.__tuple

    def set(self):
        """
        The frozenset attached to this simplex.

        EXAMPLES::

            sage: Simplex(3).set()
            frozenset([0, 1, 2, 3])
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
        return self.__set.__contains__(x)

    def __getitem__(self, n):
        """
        Return the `n`-th vertex in this simplex.

        EXAMPLES::

            sage: Simplex(5)[2]
            2
            sage: Simplex(['a', 'b', 'c'])[1]
            'b'
        """
        return self.__tuple.__getitem__(n)

    def __iter__(self):
        """
        Iterator for the vertices of this simplex.

        EXAMPLES::

            sage: [v**2 for v in Simplex(3)]
            [0, 1, 4, 9]
        """
        return self.__tuple.__iter__()

    def __add__(self, other):
        """
        Simplex obtained by concatenating the underlying tuples of the
        two arguments.

        :param other: another simplex

        EXAMPLES::

            sage: Simplex((1,2,3)) + Simplex((5,6))
            (1, 2, 3, 5, 6)
        """
        return Simplex(self.__tuple.__add__(other.__tuple))

    def face(self, n):
        """
        The `n`-th face of this simplex.

        :param n: an integer between 0 and the dimension of this simplex
        :type n: integer
        :return: the simplex obtained by removing the `n`-th vertex from this
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
            raise IndexError, "%s does not have an nth face for n=%s." % (self, n)

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
        return [self.face(i) for i in range(self.dimension()+1)]

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

        :param right: the other simplex (the right-hand factor)

        :param rename_vertices: If this is ``True``, the vertices in the
            join will be renamed by this formula: vertex "v" in the
            left-hand factor --> vertex "Lv" in the join, vertex "w"
            in the right-hand factor --> vertex "Rw" in the join.  If
            this is false, this tries to construct the join without
            renaming the vertices; this may cause problems if the two
            factors have any vertices with names in common.

        :type rename_vertices: boolean; optional, default ``True``

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

        :param other: the other simplex

        :param rename_vertices: If this is ``False``, then the vertices in
            the product are the set of ordered pairs `(v,w)` where `v`
            is a vertex in the left-hand factor (``self``) and `w` is
            a vertex in the right-hand factor (``other``). If this is
            ``True``, then the vertices are renamed as "LvRw" (e.g., the
            vertex (1,2) would become "L1R2").  This is useful if you
            want to define the Stanley-Reisner ring of the complex:
            vertex names like (0,1) are not suitable for that, while
            vertex names like "L0R1" are.

        :type rename_vertices: boolean; optional, default ``True``

        Algorithm: see Hatcher, p. 277-278 [Hat]_ (who in turn refers to
        Eilenberg-Steenrod, p. 68): given ``S = Simplex(m)`` and
        ``T = Simplex(n)``, then `S \times T` can be
        triangulated as follows: for each path `f` from `(0,0)` to
        `(m,n)` along the integer grid in the plane, going up or right
        at each lattice point, associate an `(m+n)`-simplex with
        vertices `v_0`, `v_1`, ..., where `v_k` is the `k^{th}` vertex
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
            new = tuple(["L" + str(v) + "R" + str(w) for (v,w) in x])
            answer.append(Simplex(new))
        return answer

    def __cmp__(self, other):
        """
        Return ``True`` iff this simplex is the same as ``other``: that
        is, if the vertices of the two are the same, even with a
        different ordering

        :param other: the other simplex

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
            return -1
        if self.__set == other.__set:
            return 0
        return cmp(sorted(tuple(self.__set)), sorted(tuple(other.__set)))

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
        return self.__tuple.__repr__()

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

class SimplicialComplex(GenericCellComplex):
    r"""
    Define a simplicial complex.

    :param maximal_faces: set of maximal faces
    :param maximality_check: see below
    :type maximality_check: boolean; optional, default ``True``
    :param sort_facets: see below
    :type sort_facets: boolean; optional, default ``True``
    :param name_check: see below
    :type name_check: boolean; optional, default ``False``
    :param is_mutable: Set to ``False`` to make this immutable
    :type is_mutable: boolean; optional, default ``True``
    :return: a simplicial complex

    ``maximal_faces`` should be a list or tuple or set (indeed,
    anything which may be converted to a set) whose elements are lists
    (or tuples, etc.) of vertices.  Maximal faces are also known as
    'facets'.

    If ``maximality_check`` is ``True``, check that each maximal face is,
    in fact, maximal. In this case, when producing the internal
    representation of the simplicial complex, omit those that are not.
    It is highly recommended that this be ``True``; various methods for
    this class may fail if faces which are claimed to be maximal are
    in fact not.

    If ``sort_facets`` is ``True``, sort the vertices in each facet.  If
    the vertices in different facets are not ordered compatibly (e.g.,
    if you have facets ``(1, 3, 5)`` and ``(5, 3, 8)``), then homology
    calculations may have unpredictable results.

    If ``name_check`` is ``True``, check the names of the vertices to see
    if they can be easily converted to generators of a polynomial ring
    -- use this if you plan to use the Stanley-Reisner ring for the
    simplicial complex.

    .. WARNING::

        Earlier versions of Sage supported a ``vertex_set`` argument
        to specify the vertices. This is now deprecated -- see
        :trac:`12587` -- the set of vertices is determined from the
        maximal faces.

    EXAMPLES::

        sage: SimplicialComplex([[1,2], [1,4]])
        Simplicial complex with vertex set (1, 2, 4) and facets {(1, 2), (1, 4)}
        sage: SimplicialComplex([[0,2], [0,3], [0]])
        Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2), (0, 3)}
        sage: SimplicialComplex([[0,2], [0,3], [0]], maximality_check=False)
        Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2), (0, 3), (0,)}
        sage: S = SimplicialComplex((('a', 'b'), ['a', 'c'], ('b', 'c')))
        sage: S
        Simplicial complex with vertex set ('a', 'b', 'c') and facets {('b', 'c'), ('a', 'c'), ('a', 'b')}

    Finally, if there is only one argument and it is a
    simplicial complex, return that complex.  If it is an object with
    a built-in conversion to simplicial complexes (via a
    ``_simplicial_`` method), then the resulting simplicial complex is
    returned::

        sage: S = SimplicialComplex([[0,2], [0,3], [0,6]])
        sage: SimplicialComplex(S) == S
        True
        sage: Tc = cubical_complexes.Torus(); Tc
        Cubical complex with 16 vertices and 64 cubes
        sage: Ts = SimplicialComplex(Tc); Ts
        Simplicial complex with 16 vertices and 32 facets
        sage: Ts.homology()
        {0: 0, 1: Z x Z, 2: Z}

    TESTS:

    Check that we can make mutable copies (see :trac:`14142`)::

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
        """

    def __init__(self, vertex_set=None, maximal_faces=None, **kwds):
        """
        Define a simplicial complex.  See ``SimplicialComplex`` for more
        documentation.

        .. WARNING::

            We are deprecating the option ``vertex_set`` in :trac:`12587`.

        EXAMPLES::

            sage: SimplicialComplex([[0,2], [0,3], [0]])
            Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2), (0, 3)}
            sage: SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')))
            Simplicial complex with vertex set ('a', 'b', 'c') and facets {('b', 'c'), ('a', 'c'), ('a', 'b')}

        TESTS::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S == S2
            True
            sage: S3 = SimplicialComplex(maximal_faces=[[1,4], [2,4]])
            sage: S == S3
            True

            sage: S = SimplicialComplex((('a', 'b'), ('a', 'c'), ('b', 'c')))
            sage: S == loads(dumps(S))
            True

            sage: Y = SimplicialComplex([1,2,3,4], [[1,2], [2,3], [3,4]])
            doctest:1: DeprecationWarning: vertex_set is deprecated.
            See http://trac.sagemath.org/12587 for details.
            sage: Y = SimplicialComplex([1,2,3,4], [[1,2], [2,3], [3,4]], vertex_check=False)
            doctest:1: DeprecationWarning: vertex_check is deprecated.
            See http://trac.sagemath.org/12587 for details.
        """
        from sage.misc.misc import union
        # process kwds
        sort_facets = kwds.get('sort_facets', True)
        maximality_check = kwds.get('maximality_check', True)
        name_check = kwds.get('name_check', False)
        # done with kwds except mutability

        # For deprecation #12587
        if maximal_faces is None:
            maximal_faces = vertex_set
        elif vertex_set is not None:
            # We've passed in both vertex_set and maximal_faces
            from sage.misc.superseded import deprecation
            deprecation(12587, "vertex_set is deprecated.")

        if 'vertex_check' in kwds:
            from sage.misc.superseded import deprecation
            deprecation(12587, "vertex_check is deprecated.")

        C = None
        if maximal_faces is None:
            vertex_set = []
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
                if len(maximal_faces) != 0:
                    vertex_set = reduce(union, maximal_faces)
        if C is not None:
            self._vertex_set = copy(C.vertices())
            self._facets = list(C.facets())
            self._faces = copy(C._faces)
            self._gen_dict = copy(C._gen_dict)
            self._complex = copy(C._complex)
            self.__contractible = copy(C.__contractible)
            self.__enlarged = copy(C.__enlarged)
            self._graph = copy(C._graph)
            self._is_mutable = True
            return

        if sort_facets:
            try:  # vertex_set is an iterable
                vertices = Simplex(sorted(vertex_set))
            except TypeError:  # vertex_set is an integer
                vertices = Simplex(vertex_set)
        else:
            vertices = Simplex(vertex_set)
        gen_dict = {}
        for v in vertices:
            if name_check:
                try:
                    if int(v) < 0:
                        raise ValueError("The vertex %s does not have an appropriate name."%v)
                except ValueError:  # v is not an integer
                    try:
                        normalize_names(1, v)
                    except ValueError:
                        raise ValueError("The vertex %s does not have an appropriate name."%v)
            # build dictionary of generator names
            try:
                gen_dict[v] = 'x%s'%int(v)
            except StandardError:
                gen_dict[v] = v
        # build set of facets
        good_faces = []
        maximal_simplices = [Simplex(f) for f in maximal_faces]
        for face in maximal_simplices:
            # check whether each given face is actually maximal
            face_is_maximal = True
            if maximality_check:
                faces_to_be_removed = []
                for other in good_faces:
                    if other.is_face(face):
                        faces_to_be_removed.append(other)
                    elif face_is_maximal:
                        face_is_maximal = not face.is_face(other)
                for x in faces_to_be_removed:
                    good_faces.remove(x)
            if sort_facets:
                face = Simplex(sorted(face.tuple()))
            if face_is_maximal:
                good_faces += [face]
        # if no maximal faces, add the empty face as a facet
        if len(maximal_simplices) == 0:
            good_faces.append(Simplex(-1))
        # now record the attributes for self
        # self._vertex_set: the Simplex formed by the vertices
        self._vertex_set = vertices
        # self._facets: list of facets
        self._facets = good_faces
        # self._sorted: True if the vertex set should be sorted. This
        # gets used by the add_faces method.
        self._sorted = sort_facets
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
        self._is_mutable = True
        if not kwds.get('is_mutable', True) or kwds.get('is_immutable', False):
            self.set_immutable()

    def __hash__(self):
        """
        Compute the hash value of ``self``.

        If this simplicial complex is immutable, it computes the hash value
        based upon the facets. Otherwise it raises a ``ValueError``.

        EXAMPLES::

            sage: S = SimplicialComplex([[1,4], [2,4]])
            sage: hash(S)
            Traceback (most recent call last):
            ...
            ValueError: This simplicial complex must be immutable. Call set_immutable().
            sage: S.set_immutable()
            sage: hash(S) == hash(S)
            True
            sage: S2 = SimplicialComplex([[1,4], [2,4]], is_mutable=False)
            sage: S == S2
            True
            sage: hash(S) == hash(S2)
            True
        """
        if self._is_mutable:
            raise ValueError("This simplicial complex must be immutable. Call set_immutable().")
        return hash(self._facets)

    def __cmp__(self,right):
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
        if set(self._facets) == set(right._facets):
            return 0
        else:
            return -1

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
        The vertex set of this simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex([[i] for i in range(16)] + [[0,1], [1,2]])
            sage: S
            Simplicial complex with 16 vertices and 15 facets
            sage: S.vertices()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

        Note that this actually returns a simplex::

            sage: type(S.vertices())
            <class 'sage.homology.simplicial_complex.Simplex'>
        """
        return self._vertex_set

    def maximal_faces(self):
        """
        The maximal faces (a.k.a. facets) of this simplicial complex.

        This just returns the set of facets used in defining the
        simplicial complex, so if the simplicial complex was defined
        with no maximality checking, none is done here, either.

        EXAMPLES::

            sage: Y = SimplicialComplex([[0,2], [1,4]])
            sage: Y.maximal_faces()
            {(1, 4), (0, 2)}

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

        :param subcomplex: a subcomplex of this simplicial complex.
            Return faces which are not in this subcomplex.

        :type subcomplex: optional, default ``None``

        EXAMPLES::

            sage: Y = SimplicialComplex([[1,2], [1,4]])
            sage: Y.faces()
            {0: set([(4,), (2,), (1,)]), 1: set([(1, 2), (1, 4)]), -1: set([()])}
            sage: L = SimplicialComplex([[1,2]])
            sage: Y.faces(subcomplex=L)
            {0: set([(4,)]), 1: set([(1, 4)]), -1: set([])}
        """
        # Make the subcomplex immutable if it is not
        if subcomplex is not None and subcomplex._is_mutable:
            subcomplex = SimplicialComplex(subcomplex._facets, maximality_check=False,
                                           sort_facets=False, is_mutable=False)

        if subcomplex not in self._faces:
            # Faces is the dictionary of faces in self but not in
            # subcomplex, indexed by dimension
            Faces = {}
            # sub_facets is the dictionary of facets in the subcomplex
            sub_facets = {}
            dimension = max([face.dimension() for face in self._facets])
            for i in range(-1,dimension+1):
                Faces[i] = set([])
                sub_facets[i] = set([])
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

        INPUTS:

        - ``increasing`` -- (optional, default ``True``) if ``True``, return
          faces in increasing order of dimension, thus starting with
          the empty face. Otherwise it returns faces in decreasing order of
          dimension.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: [f for f in S1.face_iterator()]
            [(), (2,), (0,), (1,), (1, 2), (0, 2), (0, 1)]
        """
        Fs = self.faces()
        dim_index = xrange(-1,self.dimension()+1)
        if not increasing:
            dim_index = reversed(dim_index)
        for i in dim_index:
            for F in Fs[i]:
                yield F

    cells = faces

    def n_faces(self, n, subcomplex=None):
        """
        The set of simplices of dimension ``n`` of this simplicial complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional faces which are *not* in the
        subcomplex.

        :param n: non-negative integer
        :param subcomplex: a subcomplex of this simplicial complex.
           Return ``n``-dimensional faces which are not in this
           subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES::

            sage: S = Set(range(1,5))
            sage: Z = SimplicialComplex(S.subsets())
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.n_faces(2)
            set([(1, 3, 4), (1, 2, 3), (2, 3, 4), (1, 2, 4)])
            sage: K = SimplicialComplex([[1,2,3], [2,3,4]])
            sage: Z.n_faces(2, subcomplex=K)
            set([(1, 3, 4), (1, 2, 4)])
        """
        if n in self.faces(subcomplex):
            return self.faces(subcomplex)[n]
        else:
            return set([])

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
        empy simplex), then the `h`-vector `(h_0, h_1, ..., h_d,
        h_{d+1})` is defined by

        .. MATH::

           \sum_{i=0}^{d+1} h_i x^{d+1-i} = \sum_{i=0}^{d+1} f_{i-1} (x-1)^{d+1-i}.

        Alternatively,

        .. MATH::

           h_j = \sum_{i=-1}^{j-1} (-1)^{j-i-1} \binom{d-i}{j-i-1} f_i.

        EXAMPLES:

        The `f`- and `h`-vectors of the boundary of an octahedron are
        computed in Wikipedia's page on simplicial complexes,
        http://en.wikipedia.org/wiki/Simplicial_complex::

            sage: square = SimplicialComplex([[0,1], [1,2], [2,3], [0,3]])
            sage: S0 = SimplicialComplex([[0], [1]])
            sage: octa = square.join(S0) # boundary of an octahedron
            sage: octa.f_vector()
            [1, 6, 12, 8]
            sage: octa.h_vector()
            [1, 3, 3, 1]
        """
        from sage.rings.arith import binomial
        d = self.dimension()
        f = self.f_vector()  # indexed starting at 0, since it's a Python list
        h = []
        for j in range(0, d+2):
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

            sage: S3 = simplicial_complexes.Sphere(3).barycentric_subdivision()
            sage: S3.f_vector()
            [1, 30, 150, 240, 120]
            sage: S3.h_vector()
            [1, 26, 66, 26, 1]
            sage: S3.g_vector()
            [1, 25, 40]
        """
        from sage.functions.other import floor
        d = self.dimension()
        h = self.h_vector()
        g = [1]
        for i in range(1, floor((d+1)/2) + 1):
            g.append(h[i] - h[i-1])
        return g

    def is_pseudomanifold(self):
        """
        Return True if self is a pseudomanifold.

        A pseudomanifold is a simplicial complex with the following properties:

        - it is pure of some dimension `d` (all of its facets are `d`-dimensional)
        - every `(d-1)`-dimensional simplex is the face of exactly two facets
        - for every two facets `S` and `T`, there is a sequence of
          facets

          .. math::

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
        X = self.n_faces(d-1)
        # is each (d-1)-simplex is the face of exactly two facets?
        for s in X:
            if len([a for a in [s.is_face(f) for f in F] if a]) != 2:
                return False
        # construct a graph with one vertex for each facet, one edge
        # when two facets intersect in a (d-1)-simplex, and see
        # whether that graph is connected.
        V = [f.set() for f in self.facets()]
        E = (lambda a,b: len(a.intersection(b)) == d)
        g = Graph([V,E])
        return g.is_connected()

    def product(self, right, rename_vertices=True, is_mutable=True):
        """
        The product of this simplicial complex with another one.

        :param right: the other simplicial complex (the right-hand
           factor)

        :param rename_vertices: If this is False, then the vertices in
           the product are the set of ordered pairs `(v,w)` where `v`
           is a vertex in ``self`` and `w` is a vertex in
           ``right``. If this is ``True``, then the vertices are renamed
           as "LvRw" (e.g., the vertex (1,2) would become "L1R2").
           This is useful if you want to define the Stanley-Reisner
           ring of the complex: vertex names like (0,1) are not
           suitable for that, while vertex names like "L0R1" are.

        :type rename_vertices: boolean; optional, default ``True``

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        The vertices in the product will be the set of ordered pairs
        `(v,w)` where `v` is a vertex in self and `w` is a vertex in
        right.

        .. WARNING::

           If ``X`` and ``Y`` are simplicial complexes, then ``X*Y``
           returns their join, not their product.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1], [1,2], [0,2]]) # circle
            sage: K = SimplicialComplex([[0,1]])   # edge
            sage: S.product(K).vertices()  # cylinder
            ('L0R0', 'L0R1', 'L1R0', 'L1R1', 'L2R0', 'L2R1')
            sage: S.product(K, rename_vertices=False).vertices()
            ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
            sage: T = S.product(S)  # torus
            sage: T
            Simplicial complex with 9 vertices and 18 facets
            sage: T.homology()
            {0: 0, 1: Z x Z, 2: Z}

        These can get large pretty quickly::

            sage: T = simplicial_complexes.Torus(); T
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 14 facets
            sage: K = simplicial_complexes.KleinBottle(); K
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 16 facets
            sage: T.product(K)      # long time: 5 or 6 seconds
            Simplicial complex with 56 vertices and 1344 facets
        """
        facets = []
        for f in self._facets:
            for g in right._facets:
                facets.extend(f.product(g, rename_vertices))
        return SimplicialComplex(facets, is_mutable=is_mutable)

    def join(self, right, rename_vertices=True, is_mutable=True):
        """
        The join of this simplicial complex with another one.

        The join of two simplicial complexes `S` and `T` is the
        simplicial complex `S*T` with simplices of the form `[v_0,
        ..., v_k, w_0, ..., w_n]` for all simplices `[v_0, ..., v_k]` in
        `S` and `[w_0, ..., w_n]` in `T`.

        :param right: the other simplicial complex (the right-hand factor)

        :param rename_vertices: If this is True, the vertices in the
           join will be renamed by the formula: vertex "v" in the
           left-hand factor --> vertex "Lv" in the join, vertex "w" in
           the right-hand factor --> vertex "Rw" in the join.  If this
           is false, this tries to construct the join without renaming
           the vertices; this will cause problems if the two factors
           have any vertices with names in common.

        :type rename_vertices: boolean; optional, default ``True``

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        EXAMPLES::

            sage: S = SimplicialComplex([[0], [1]])
            sage: T = SimplicialComplex([[2], [3]])
            sage: S.join(T)
            Simplicial complex with vertex set ('L0', 'L1', 'R2', 'R3') and 4 facets
            sage: S.join(T, rename_vertices=False)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 3), (1, 2), (0, 2), (0, 3)}

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

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        The cone is the simplicial complex formed by adding a new
        vertex `C` and simplices of the form `[C, v_0, ..., v_k]` for
        every simplex `[v_0, ..., v_k]` in the original simplicial
        complex.  That is, the cone is the join of the original
        complex with a one-point simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex([[0], [1]])
            sage: S.cone()
            Simplicial complex with vertex set ('L0', 'L1', 'R0') and facets {('L0', 'R0'), ('L1', 'R0')}
        """
        return self.join(SimplicialComplex([["0"]], is_mutable=is_mutable),
                         rename_vertices = True)

    def suspension(self, n=1, is_mutable=True):
        r"""
        The suspension of this simplicial complex.

        :param n: positive integer -- suspend this many times.

        :type n: optional, default 1

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        The suspension is the simplicial complex formed by adding two
        new vertices `S_0` and `S_1` and simplices of the form `[S_0,
        v_0, ..., v_k]` and `[S_1, v_0, ..., v_k]` for every simplex
        `[v_0, ..., v_k]` in the original simplicial complex.  That
        is, the suspension is the join of the original complex with a
        two-point simplicial complex.

        If the simplicial complex `M` happens to be a pseudomanifold
        (see :meth:`is_pseudomanifold`), then this instead constructs
        Datta's one-point suspension (see p. 434 in the cited
        article): choose a vertex `u` in `M` and choose a new vertex
        `w` to add.  Denote the join of simplices by "`*`".  The
        facets in the one-point suspension are of the two forms

        - `u * \alpha` where `\alpha` is a facet of `M` not containing
          `u`

        - `w * \beta` where `\beta` is any facet of `M`.

        REFERENCES:

        - Basudeb Datta, "Minimal triangulations of manifolds",
          J. Indian Inst. Sci. 87 (2007), no. 4, 429-449.

        EXAMPLES::

            sage: S0 = SimplicialComplex([[0], [1]])
            sage: S0.suspension() == simplicial_complexes.Sphere(1)
            True
            sage: S3 = S0.suspension(3)  # the 3-sphere
            sage: S3.homology()
            {0: 0, 1: 0, 2: 0, 3: Z}

        For pseudomanifolds, the complex constructed here will be
        smaller than that obtained by taking the join with the
        0-sphere: the join adds two vertices, while this construction
        only adds one. ::

            sage: T = simplicial_complexes.Torus()
            sage: T.join(S0).vertices()      # 9 vertices
            ('L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'R0', 'R1')
            sage: T.suspension().vertices()  # 8 vertices
            (0, 1, 2, 3, 4, 5, 6, 7)
        """
        if n<0:
            raise ValueError, "n must be non-negative."
        if n==0:
            return self
        if n==1:
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
                                 rename_vertices = True)
        return self.suspension(1, is_mutable).suspension(int(n-1), is_mutable)

    def disjoint_union(self, right, rename_vertices=True, is_mutable=True):
        """
        The disjoint union of this simplicial complex with another one.

        :param right: the other simplicial complex (the right-hand factor)

        :param rename_vertices: If this is True, the vertices in the
           disjoint union will be renamed by the formula: vertex "v"
           in the left-hand factor --> vertex "Lv" in the disjoint
           union, vertex "w" in the right-hand factor --> vertex "Rw"
           in the disjoint union.  If this is false, this tries to
           construct the disjoint union without renaming the vertices;
           this will cause problems if the two factors have any
           vertices with names in common.

        :type rename_vertices: boolean; optional, default True

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S1.disjoint_union(S2).homology()
            {0: Z, 1: Z, 2: Z}
        """
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

        :param right: the other simplicial complex (the right-hand factor)

        :param rename_vertices: If this is ``True``, the vertices in the
           wedge will be renamed by the formula: first vertex in each
           are glued together and called "0".  Otherwise, each vertex
           "v" in the left-hand factor --> vertex "Lv" in the wedge,
           vertex "w" in the right-hand factor --> vertex "Rw" in the
           wedge.  If this is ``False``, this tries to construct the wedge
           without renaming the vertices; this will cause problems if
           the two factors have any vertices with names in common.

        :type rename_vertices: boolean; optional, default ``True``

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        .. NOTE::

            This operation is not well-defined if ``self`` or
            ``other`` is not path-connected.

        EXAMPLES::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S1.wedge(S2).homology()
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

    def chain_complex(self, **kwds):
        """
        The chain complex associated to this simplicial complex.

        :param dimensions: if ``None``, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.
        :param base_ring: commutative ring
        :type base_ring: optional, default ``ZZ``
        :param subcomplex: a subcomplex of this simplicial complex.
           Compute the chain complex relative to this subcomplex.
        :type subcomplex: optional, default empty
        :param augmented: If ``True``, return the augmented chain complex
           (that is, include a class in dimension `-1` corresponding
           to the empty cell).  This is ignored if ``dimensions`` is
           specified.
        :type augmented: boolean; optional, default ``False``
        :param cochain: If ``True``, return the cochain complex (that is,
           the dual of the chain complex).
        :type cochain: boolean; optional, default ``False``
        :param verbose: If ``True``, print some messages as the chain
           complex is computed.
        :type verbose: boolean; optional, default ``False``
        :param check_diffs: If ``True``, make sure that the chain complex
           is actually a chain complex: the differentials are
           composable and their product is zero.
        :type check_diffs: boolean; optional, default ``False``

        .. NOTE::

           If subcomplex is nonempty, then the argument ``augmented``
           has no effect: the chain complex relative to a nonempty
           subcomplex is zero in dimension `-1`.

        EXAMPLES::

            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle.chain_complex()
            Chain complex with at most 2 nonzero terms over Integer Ring
            sage: circle.chain_complex()._latex_()
            '\\Bold{Z}^{3} \\xrightarrow{d_{1}} \\Bold{Z}^{3}'
            sage: circle.chain_complex(base_ring=QQ, augmented=True)
            Chain complex with at most 3 nonzero terms over Rational Field
        """
        augmented = kwds.get('augmented', False)
        cochain = kwds.get('cochain', False)
        verbose = kwds.get('verbose', False)
        check_diffs = kwds.get('check_diffs', False)
        base_ring = kwds.get('base_ring', ZZ)
        dimensions = kwds.get('dimensions', None)
        subcomplex = kwds.get('subcomplex', None)

        # initialize subcomplex
        if subcomplex is None:
            subcomplex = SimplicialComplex(is_mutable=False)
        else:
            # subcomplex is not empty, so don't augment the chain complex
            augmented = False
            # Use an immutable copy of the subcomplex
            if not subcomplex._is_mutable:
                subcomplex = SimplicialComplex(subcomplex._facets, maximality_check=False,
                                               sort_facets=False, is_mutable=False)
        # now construct the range of dimensions in which to compute
        if dimensions is None:
            dimensions = range(0, self.dimension()+1)
            first = 0
        else:
            augmented = False
            first = dimensions[0]
        differentials = {}
        # in the chain complex, compute the first dimension by hand,
        # and don't cache it: it may be differ from situation to
        # situation because of boundary effects.
        current = None
        current_dim = None
        if augmented:  # then first == 0
            current = list(self.n_faces(0, subcomplex=subcomplex))
            current_dim = 0
            if cochain:
                differentials[-1] = matrix(base_ring, len(current), 1,
                                           [1]*len(current))
            else:
                differentials[0] = matrix(base_ring, 1, len(current),
                                          [1]*len(current))
        elif first == 0 and not augmented:
            current = list(self.n_faces(0, subcomplex=subcomplex))
            current_dim = 0
            if not cochain:
                differentials[0] = matrix(base_ring, 0, len(current))
        else:  # first > 0
            current = list(self.n_faces(first, subcomplex=subcomplex))
            current_dim = first
            if not cochain:
                differentials[first] = matrix(base_ring, 0, len(current))
        for n in dimensions[1:]:
            if verbose:
                print "  starting dimension %s" % n
            if (n, subcomplex) in self._complex:
                if cochain:
                    differentials[n-1] = self._complex[(n, subcomplex)].transpose().change_ring(base_ring)
                    mat = differentials[n-1]
                else:
                    differentials[n] = self._complex[(n, subcomplex)].change_ring(base_ring)
                    mat = differentials[n]
                if verbose:
                    print "    boundary matrix (cached): it's %s by %s." % (mat.nrows(), mat.ncols())
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
                    set_of_faces = list(self.n_faces(n-1, subcomplex=subcomplex))
                    old = dict(zip(set_of_faces, range(len(set_of_faces))))
                current = list(self.n_faces(n, subcomplex=subcomplex))
                current_dim = n
                # construct matrix.  it is easiest to construct it as
                # a sparse matrix, specifying which entries are
                # nonzero via a dictionary.
                matrix_data = {}
                col = 0
                if len(old) and len(current):
                    for simplex in current:
                        for i in range(n+1):
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
                    print "    boundary matrix computed: it's %s by %s." % (mat.nrows(), mat.ncols())
        # now for the cochain complex, compute the last dimension by
        # hand, and don't cache it.
        if cochain:
            n = dimensions[-1] + 1
            if current_dim != n-1:
                current = list(self.n_faces(n-1, subcomplex=subcomplex))
            differentials[n-1] = matrix(base_ring, 0, len(current))
        # finally, return the chain complex
        if cochain:
            return ChainComplex(data=differentials, degree=1, **kwds)
        else:
            return ChainComplex(data=differentials, degree=-1, **kwds)

    def _homology_(self, dim=None, **kwds):
        """
        The reduced homology of this simplicial complex.

        :param dim: If ``None``, then return the homology in every
           dimension.  If ``dim`` is an integer or list, return the
           homology in the given dimensions.  (Actually, if ``dim`` is
           a list, return the homology in the range from ``min(dim)``
           to ``max(dim)``.)

        :type dim: integer or list of integers or ``None``; optional,
                   default ``None``

        :param base_ring: commutative ring. Must be ``ZZ`` or a field.

        :type base_ring: optional, default ``ZZ``

        :param subcomplex: a subcomplex of this simplicial complex.
           Compute homology relative to this subcomplex.

        :type subcomplex: optional, default ``None``

        :param cohomology: If ``True``, compute cohomology rather than
           homology.

        :type cohomology: boolean; optional, default ``False``

        :param enlarge: If ``True``, find a new subcomplex homotopy
           equivalent to, and probably larger than, the given one.

        :type enlarge: boolean; optional, default ``True``

        :param algorithm: The options are ``'auto'``, ``'dhsw'``,
           ``'pari'`` or  ``'no_chomp'``.  If ``'auto'``, first try CHomP,
           then use the Dumas, Heckenbach, Saunders, and Welker elimination
           algorithm for large matrices, Pari for small ones.  If
           ``'no_chomp'``, then don't try CHomP, but behave the same
           otherwise.  If ``'pari'``, then compute elementary divisors
           using Pari.  If ``'dhsw'``, then use the DHSW algorithm to
           compute elementary divisors.  (As of this writing, CHomP is
           by far the fastest option, followed by the ``'auto'`` or
           ``'no_chomp'`` setting of using DHSW for large matrices and
           Pari for small ones.)

        :type algorithm: string; optional, default ``'auto'``

        :param verbose: If ``True``, print some messages as the homology
           is computed.

        :type verbose: boolean; optional, default ``False``

        Algorithm: if ``subcomplex`` is ``None``, replace it with a facet
        -- a contractible subcomplex of the original complex.  Then no
        matter what ``subcomplex`` is, replace it with a subcomplex
        `L` which is homotopy equivalent and as large as possible.
        Compute the homology of the original complex relative to `L`:
        if `L` is large, then the relative chain complex will be small
        enough to speed up computations considerably.

        EXAMPLES::

            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle._homology_()
            {0: 0, 1: Z}
            sage: sphere = SimplicialComplex([[0,1,2,3]])
            sage: sphere.remove_face([0,1,2,3])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: sphere._homology_()
            {0: 0, 1: 0, 2: Z}

        Another way to get a two-sphere: take a two-point space and take its
        three-fold join with itself::

            sage: S = SimplicialComplex([[0], [1]])
            sage: (S*S*S)._homology_(dim=2, cohomology=True)
            Z

        Relative homology::

            sage: T = SimplicialComplex([[0,1,2]])
            sage: U = SimplicialComplex([[0,1], [1,2], [0,2]])
            sage: T._homology_(subcomplex=U)
            {0: 0, 1: 0, 2: Z}
        """
        from sage.modules.all import VectorSpace
        from sage.homology.homology_group import HomologyGroup

        base_ring = kwds.get('base_ring', ZZ)
        cohomology = kwds.get('cohomology', False)
        enlarge = kwds.get('enlarge', True)
        verbose = kwds.get('verbose', False)
        subcomplex = kwds.get('subcomplex', None)

        if dim is not None:
            if isinstance(dim, (list, tuple)):
                low = min(dim) - 1
                high = max(dim) + 2
            else:
                low = dim - 1
                high = dim + 2
            dims = range(low, high)
        else:
            dims = None

        if verbose:
            print "starting calculation of the homology of this",
            print "%s-dimensional simplicial complex" % self.dimension()
        if subcomplex is None:
            if enlarge:
                if verbose:
                    print "Constructing contractible subcomplex..."
                L = self._contractible_subcomplex(verbose=verbose)
                if verbose:
                    print "Done finding contractible subcomplex."
                    vec = [len(self.n_faces(n-1, subcomplex=L)) for n in range(self.dimension()+2)]
                    print "The difference between the f-vectors is:"
                    print "  %s" % vec
            else:
                L = SimplicialComplex([[self.vertices().tuple()[0]]])
        else:
            if enlarge:
                if verbose:
                    print "Enlarging subcomplex..."
                L = self._enlarge_subcomplex(subcomplex, verbose=verbose)
                if verbose:
                    print "Done enlarging subcomplex:"
            else:
                L = subcomplex
        L.set_immutable()

        if verbose:
            print "Computing the chain complex..."
        kwds['subcomplex']=L
        C = self.chain_complex(dimensions=dims, augmented=True,
                               cochain=cohomology, **kwds)
        if verbose:
            print " Done computing the chain complex. "
            print "Now computing homology..."
        if 'subcomplex' in kwds:
            del kwds['subcomplex']
        answer = C.homology(**kwds)

        if dim is None:
            dim = range(self.dimension()+1)
        zero = HomologyGroup(0, base_ring)
        if isinstance(dim, (list, tuple)):
            return dict([d, answer.get(d, zero)] for d in dim)
        return answer.get(dim, zero)

    def add_face(self, face):
        """
        Add a face to this simplicial complex

        :param face: a subset of the vertex set

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
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2, 3), (0, 1)}

        If you add a face which is already present, there is no effect::

            sage: Y.add_face([1,3]); Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2, 3), (0, 1)}

        Check that the bug reported at :trac:`14354` has been fixed::

            sage: T = SimplicialComplex([range(1,5)]).n_skeleton(1)
            sage: T.homology()
            {0: 0, 1: Z x Z x Z}
            sage: T.add_face([1,2,3])
            sage: T.homology()
            {0: 0, 1: Z x Z, 2: 0}

        Check we've fixed the bug reported at :trac:`14578`::

            sage: t0 = SimplicialComplex()
            sage: t0.add_face(('a', 'b'))
            sage: t0.add_face(('c', 'd', 'e'))
            sage: t0.add_face(('e', 'f', 'c'))
            sage: t0.homology()
            {0: Z, 1: 0, 2: 0}
        """
        if not self._is_mutable:
            raise ValueError("This simplicial complex is not mutable")

        if self._sorted:
            new_face = Simplex(sorted(face))
        else:
            new_face = Simplex(face)

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
            from sage.misc.misc import union

            if self._sorted:
                self._vertex_set = Simplex(sorted(reduce(union, [self._vertex_set, new_face])))
            else:
                self._vertex_set = Simplex(reduce(union, [self._vertex_set, new_face]))

            # update self._faces if necessary
            if None in self._faces:
                all_new_faces = SimplicialComplex([new_face]).faces()
                for dim in range(0, new_face.dimension()+1):
                    if dim in self._faces[None]:
                        self._faces[None][dim] = self._faces[None][dim].union(all_new_faces[dim])
                    else:
                        self._faces[None][dim] = all_new_faces[dim]
            # update self._graph if necessary
            if self._graph is not None:
                d = new_face.dimension()+1
                for i in range(d):
                    for j in range(i+1,d):
                        self._graph.add_edge(new_face[i],new_face[j])
            self._complex = {}
            self.__contractible = None
            self.__enlarged = {}

    def remove_face(self, face):
        """
        Remove a face from this simplicial complex and return the
        resulting simplicial complex.

        :param face: a face of the simplicial complex

        This *changes* the simplicial complex.

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
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 3, 4), (2, 3, 4)}

            sage: S = SimplicialComplex([[0,1,2],[2,3]])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (2, 3)}
            sage: S.remove_face([0,1,2])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2), (2, 3), (0, 2), (0, 1)}
        """
        if not self._is_mutable:
            raise ValueError("This simplicial complex is not mutable")

        simplex = Simplex(face)
        facets = self.facets()
        if all([not simplex.is_face(F) for F in facets]):
            # face is not in self: nothing to remove
            return self
        link = self.link(simplex)
        join_facets = []
        for f in simplex.faces():
            for g in link.facets():
                join_facets.append(f.join(g, rename_vertices=False))
        # join_facets is the list of facets in the join bdry(face) * link(face)
        remaining = join_facets + [elem for elem in facets if not simplex.is_face(elem)]

        # Check to see if there are any non-maximial faces
        # build set of facets
        self._facets = []
        for f in remaining:
            face = Simplex(f)
            face_is_maximal = True
            faces_to_be_removed = []
            for other in self._facets:
                if other.is_face(face):
                    faces_to_be_removed.append(other)
                elif face_is_maximal:
                    face_is_maximal = not face.is_face(other)
            for x in faces_to_be_removed:
                self._facets.remove(x)
            face = Simplex(sorted(face.tuple()))
            if face_is_maximal:
                self._facets.append(face)
        # if no maximal faces, add the empty face as a facet
        if len(remaining) == 0:
            self._facets.append(Simplex(-1))

        # Recreate the vertex set
        from sage.misc.misc import union
        if self._sorted:
            self._vertex_set = Simplex(sorted(reduce(union, self._facets)))
        else:
            self._vertex_set = Simplex(reduce(union, self._facets))

        # Update self._faces and self._graph if necessary
        if None in self._faces:
            self._faces = {}
            self.faces()
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

    def connected_sum(self, other, is_mutable=True):
        """
        The connected sum of this simplicial complex with another one.

        :param other: another simplicial complex
        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``
        :return: the connected sum ``self # other``

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
            sage: S1.connected_sum(S1.connected_sum(S1)).homology()
            {0: 0, 1: Z}
            sage: P = simplicial_complexes.RealProjectivePlane(); P
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
            sage: P.connected_sum(P)    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets

        The notation '+' may be used for connected sum, also::

            sage: P + P    # the Klein bottle
            Simplicial complex with 9 vertices and 18 facets
            sage: (P + P).homology()[1]
            Z x C2
        """
        if not (self.is_pure() and other.is_pure() and
                self.dimension() == other.dimension()):
            raise ValueError, "Complexes are not pure of the same dimension."
        # first find a top-dimensional simplex to remove from each surface
        keep_left = self._facets[0]
        keep_right = other._facets[0]
        # construct the set of vertices:
        left = set(self.vertices()).difference(set(keep_left))
        right = set(other.vertices()).difference(set(keep_right))
        # construct the set of facets:
        left = set(self._facets).difference(set([keep_left]))
        right = set(other._facets).difference(set([keep_right]))
        facet_set = ([[rename_vertex(v, keep=list(keep_left))
                       for v in face] for face in left]
                     + [[rename_vertex(v, keep=list(keep_right), left=False)
                         for v in face] for face in right])
        # return the new surface
        return SimplicialComplex(facet_set, is_mutable=is_mutable)

    __add__ = connected_sum

    def link(self, simplex, is_mutable=True):
        """
        The link of a simplex in this simplicial complex.

        The link of a simplex `F` is the simplicial complex formed by
        all simplices `G` which are disjoint from `F` but for which `F
        \cup G` is a simplex.

        :param simplex: a simplex in this simplicial complex.
        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1,2], [1,2,3]])
            sage: X.link(Simplex([0]))
            Simplicial complex with vertex set (1, 2) and facets {(1, 2)}
            sage: X.link([1,2])
            Simplicial complex with vertex set (0, 3) and facets {(3,), (0,)}
            sage: Y = SimplicialComplex([[0,1,2,3]])
            sage: Y.link([1])
            Simplicial complex with vertex set (0, 2, 3) and facets {(0, 2, 3)}
        """
        faces = []
        s = Simplex(simplex)
        for f in self._facets:
            if s.is_face(f):
                faces.append(Simplex(list(f.set().difference(s.set()))))
        return SimplicialComplex(faces, is_mutable=is_mutable)

    def is_cohen_macaulay(self, ncpus=0):
        r"""
        Returns True if ``self`` is Cohen-Macaulay, i.e., if
        `\tilde{H}_i(\operatorname{lk}_\Delta(F);\ZZ) = 0` for all
        `F \in \Delta` and `i < \operatorname{dim}\operatorname{lk}_\Delta(F)`.
        Here, `\Delta` is ``self``, and `\operatorname{lk}` denotes the
        link operator on ``self``.

        INPUT:

        - ``ncpus`` -- (default: 0) number of cpus used for the
          computation. If this is 0, determine the number of cpus
          automatically based on the hardware being used.

        For finite simplicial complexes, this is equivalent to the
        statement that the Stanley-Reisner ring of ``self`` is
        Cohen-Macaulay.

        EXAMPLES:

        Spheres are Cohen-Macaulay::

            sage: S = SimplicialComplex([[1,2],[2,3],[3,1]])
            sage: S.is_cohen_macaulay(ncpus=3)
            True

        The following example is taken from Bruns, Herzog - Cohen-Macaulay
        rings, Figure 5.3::

            sage: S = SimplicialComplex([[1,2,3],[1,4,5]])
            sage: S.is_cohen_macaulay(ncpus=3)
            ...
            False
        """
        from sage.parallel.decorate import parallel
        from sage.rings.rational_field import QQ

        if ncpus == 0:
            import os
            try:
                ncpus = int(os.environ['SAGE_NUM_THREADS'])
            except KeyError:
                ncpus = 1

        facs = [ x for x in self.face_iterator() ]
        n = len(facs)
        facs_divided = [ [] for i in range(ncpus) ]
        for i in range(n):
            facs_divided[i%ncpus].append(facs[i])

        def all_homologies_vanish(F):
            S = self.link(F)
            H = S.homology(base_ring=QQ)
            return all( H[j].dimension() == 0 for j in xrange(S.dimension()) )

        @parallel(ncpus=ncpus)
        def all_homologies_in_list_vanish(Fs):
            return all( all_homologies_vanish(F) for F in Fs )

        return all( answer[1] for answer in all_homologies_in_list_vanish(facs_divided) )

    def effective_vertices(self):
        """
        The set of vertices belonging to some face. Returns the list of
        vertices.

        .. WARNING::

            This method is deprecated. See :trac:`12587`.

        EXAMPLES::

            sage: S = SimplicialComplex([[0,1,2,3],[6,7]])
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3, 6, 7) and facets {(6, 7), (0, 1, 2, 3)}
            sage: S.effective_vertices()
            doctest:1: DeprecationWarning: effective_vertices is deprecated. Use vertices instead
            See http://trac.sagemath.org/12587 for details.
            (0, 1, 2, 3, 6, 7)
        """
        from sage.misc.superseded import deprecation
        deprecation(12587, "effective_vertices is deprecated. Use vertices instead")
        return self._vertex_set

    def generated_subcomplex(self,sub_vertex_set, is_mutable=True):
        """
        Returns the largest sub-simplicial complex of ``self`` containing
        exactly ``sub_vertex_set`` as vertices.

        :param sub_vertex_set: The sub-vertex set.
        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(2)
            sage: S
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: S.generated_subcomplex([0,1,2])
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}

        """
        if not self.vertices().set().issuperset(sub_vertex_set):
            raise ValueError, "input must be a subset of the vertex set."
        faces = []
        for i in range(self.dimension()+1):
            for j in self.faces()[i]:
                if j.set().issubset(sub_vertex_set):
                    faces.append(j)
        return SimplicialComplex(faces, maximality_check=True,
                                 is_mutable=is_mutable)

    def _complement(self, simplex):
        """
        Return the complement of a simplex in the vertex set of this
        simplicial complex.

        :param simplex: a simplex (need not be in the simplicial complex)

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

        :param simplices: a bunch of simplices

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
                        L = [v] + list(partial)
                        L.sort()
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

        Algorithm: first take the complement (within the vertex set)
        of each facet, obtaining a set `(f_1, f_2, ...)` of simplices.
        Now form the set of all simplices of the form `(v_1, v_2,
        ...)` where vertex `v_i` is in face `f_i`.  This set will
        contain the minimal nonfaces and may contain some non-minimal
        nonfaces also, so loop through the set to find the minimal
        ones.  (The last two steps are taken care of by the
        ``_transpose_simplices`` routine.)

        This is used in computing the
        :meth:`Stanley-Reisner ring<stanley_reisner_ring>` and the
        :meth:`Alexander dual<alexander_dual>`.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,3],[1,2]])
            sage: X.minimal_nonfaces()
            {(2, 3)}
            sage: Y = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: Y.minimal_nonfaces()
            {(1, 3), (0, 2)}
        """
        complements = [self._complement(facet) for facet in self._facets]
        return Set(self._transpose_simplices(*complements))

    def _stanley_reisner_base_ring(self, base_ring=ZZ):
        """
        The polynomial algebra of which the Stanley-Reisner ring is a
        quotient.

        :param base_ring: a commutative ring
        :type base_ring: optional, default ``ZZ``
        :return: a polynomial algebra with coefficients in base_ring,
          with one generator for each vertex in the simplicial complex.

        See the documentation for :meth:`stanley_reisner_ring` for a
        warning about the names of the vertices.

        EXAMPLES::

            sage: X = SimplicialComplex([[1,2], [0], [3]])
            sage: X._stanley_reisner_base_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring
            sage: Y = SimplicialComplex([['a', 'b', 'c']])
            sage: Y._stanley_reisner_base_ring(base_ring=QQ)
            Multivariate Polynomial Ring in a, c, b over Rational Field
        """
        return PolynomialRing(base_ring, self._gen_dict.values())

    def stanley_reisner_ring(self, base_ring=ZZ):
        """
        The Stanley-Reisner ring of this simplicial complex.

        :param base_ring: a commutative ring
        :type base_ring: optional, default ``ZZ``
        :return: a quotient of a polynomial algebra with coefficients
           in ``base_ring``, with one generator for each vertex in the
           simplicial complex, by the ideal generated by the products
           of those vertices which do not form faces in it.

        Thus the ideal is generated by the products corresponding to
        the minimal nonfaces of the simplicial complex.

        .. WARNING::

           This may be quite slow!

           Also, this may behave badly if the vertices have the
           'wrong' names. To avoid this, define the simplicial complex
           at the start with the flag ``name_check`` set to ``True``.

           More precisely, this is a quotient of a polynomial ring
           with one generator for each vertex.  If the name of a
           vertex is a non-negative integer, then the corresponding
           polynomial generator is named ``'x'`` followed by that integer
           (e.g., ``'x2'``, ``'x3'``, ``'x5'``, ...).  Otherwise, the
           polynomial generators are given the same names as the vertices.
           Thus if the vertex set is ``(2, 'x2')``, there will be problems.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [0,3]])
            sage: X.stanley_reisner_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)
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

        :param is_mutable: Determines if the output is mutable
        :type is_mutable: boolean; optional, default ``True``

        EXAMPLES::

            sage: Y = SimplicialComplex([[i] for i in range(5)]); Y
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(4,), (2,), (3,), (0,), (1,)}
            sage: Y.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and 10 facets
            sage: X = SimplicialComplex([[0,1], [1,2], [2,3], [3,0]])
            sage: X.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 3), (0, 2)}
        """
        nonfaces = self.minimal_nonfaces()
        return SimplicialComplex([self._complement(f) for f in nonfaces], is_mutable=is_mutable)

    def barycentric_subdivision(self):
        """
        The barycentric subdivision of this simplicial complex.

        See http://en.wikipedia.org/wiki/Barycentric_subdivision for a
        definition.

        EXAMPLES::

            sage: triangle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: hexagon = triangle.barycentric_subdivision()
            sage: hexagon
            Simplicial complex with 6 vertices and 6 facets
            sage: hexagon.homology(1) == triangle.homology(1)
            True

        Barycentric subdivisions can get quite large, since each
        `n`-dimensional facet in the original complex produces
        `(n+1)!` facets in the subdivision::

            sage: S4 = simplicial_complexes.Sphere(4)
            sage: S4
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 6 facets
            sage: S4.barycentric_subdivision()
            Simplicial complex with 62 vertices and 720 facets
        """
        return self.face_poset().order_complex()

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
            sage: G.edges()
            [(0, 1, None), (0, 2, None), (0, 3, None), (1, 2, None), (1, 3, None), (2, 3, None)]
        """
        if self._graph is None:
            edges = self.n_faces(1)
            vertices = map(min, filter(lambda f: f.dimension() == 0, self._facets))
            used_vertices = []  # vertices which are in an edge
            d = {}
            for e in edges:
                v = min(e)
                if v in d:
                    d[v].append(max(e))
                else:
                    d[v] = [max(e)]
                used_vertices.extend(list(e))
            for v in vertices:
                if v not in used_vertices:
                    d[v] = []
            self._graph = Graph(d)
        return self._graph

    def delta_complex(self, sort_simplices=False):
        r"""
        Returns ``self`` as a `\Delta`-complex.  The `\Delta`-complex
        is essentially identical to the simplicial complex: it has
        same simplices with the same boundaries.

        :param sort_simplices: if ``True``, sort the list of simplices in
          each dimension
        :type sort_simplices: boolean; optional, default ``False``

        EXAMPLES::

            sage: T = simplicial_complexes.Torus()
            sage: Td = T.delta_complex()
            sage: Td
            Delta complex with 7 vertices and 43 simplices
            sage: T.homology() == Td.homology()
            True
        """
        from delta_complex import DeltaComplex
        data = {}
        dim = self.dimension()
        n_cells = self.n_cells(dim)
        if sort_simplices:
            n_cells.sort()
        for n in range(dim, -1, -1):
            bdries = self.n_cells(n-1)
            if sort_simplices:
                bdries.sort()
            data[n] = []
            for f in n_cells:
                data[n].append([bdries.index(f.face(i)) for i in range(n+1)])
            n_cells = bdries
        return DeltaComplex(data)

    def is_flag_complex(self):
        """
        Returns ``True`` if and only if ``self`` is a flag complex.

        A flag complex is a simplicial complex that is the largest simplicial
        complex on its 1-skeleton. Thus a flag complex is the clique complex
        of its graph.

        EXAMPLES::

            sage: h = Graph({0:[1,2,3,4],1:[2,3,4],2:[3]})
            sage: x = h.clique_complex()
            sage: x
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 1, 4), (0, 1, 2, 3)}
            sage: x.is_flag_complex()
            True

            sage: X = simplicial_complexes.ChessboardComplex(3,3)
            sage: X.is_flag_complex()
            True
        """
        return self == self.graph().clique_complex()

    def is_connected(self):
        """
        Returns ``True`` if and only if ``self`` is connected.

        .. WARNING::

           This may give the wrong answer if the simplicial complex
           was constructed with ``maximality_check`` set to ``False``.

        EXAMPLES::

            sage: V = SimplicialComplex([[0,1,2],[3]])
            sage: V
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (3,)}
            sage: V.is_connected()
            False

            sage: X = SimplicialComplex([[0,1,2]])
            sage: X.is_connected()
            True

            sage: U = simplicial_complexes.ChessboardComplex(3,3)
            sage: U.is_connected()
            True

            sage: W = simplicial_complexes.Sphere(3)
            sage: W.is_connected()
            True

            sage: S = SimplicialComplex([[0,1],[2,3]])
            sage: S.is_connected()
            False
        """
        return self.graph().is_connected()

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this simplicial complex.

        The `n`-skeleton of a simplicial complex is obtained by discarding
        all of the simplices in dimensions larger than `n`.

        :param n: non-negative integer

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [1,2,3], [0,2,3]])
            sage: X.n_skeleton(1)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(2, 3), (0, 2), (1, 3), (1, 2), (0, 3), (0, 1)}
            sage: X.set_immutable()
            sage: X.n_skeleton(2)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (1, 2, 3), (0, 1)}
        """
        # make sure it's a list (it will be a tuple if immutable)
        facets = list(filter(lambda f: f.dimension()<n, self._facets))
        facets.extend(self.n_faces(n))
        return SimplicialComplex(facets, is_mutable=self._is_mutable)

    def _contractible_subcomplex(self, verbose=False):
        """
        Find a contractible subcomplex `L` of this simplicial complex,
        preferably one which is as large as possible.

        :param verbose: If ``True``, print some messages as the simplicial
           complex is computed.
        :type verbose: boolean; optional, default ``False``

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
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: L = sphere._contractible_subcomplex(); L
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (1, 2, 3), (0, 1, 3)}
            sage: L.homology()
            {0: 0, 1: 0, 2: 0}
        """
        facets = [self._facets[0]]
        return self._enlarge_subcomplex(SimplicialComplex(facets, is_mutable=False), verbose=verbose)

    def _enlarge_subcomplex(self, subcomplex, verbose=False):
        """
        Given a subcomplex `S` of this simplicial complex `K`, find a
        subcomplex `L`, as large as possible, containing `S` which is
        homotopy equivalent to `S` (so that `H_{*}(K,S)` is isomorphic
        to `H_{*}(K,L)`).  This way, the chain complex for computing
        `H_{*}(K,L)` will be smaller than that for computing
        `H_{*}(K,S)`, so the computations should be faster.

        :param subcomplex: a subcomplex of this simplicial complex
        :param verbose: If ``True``, print some messages as the simplicial
           complex is computed.
        :type verbose: boolean; optional, default ``False``
        :return: a complex `L` containing ``subcomplex`` and contained
           in ``self``, homotopy equivalent to ``subcomplex``.

        Algorithm: start with the subcomplex `S` and loop through the
        facets of `K` which are not in `S`.  For each one, see whether
        its intersection with `S` is contractible, and if so, add it.
        This is recursive: testing for contractibility calls this
        routine again, via ``_contractible_subcomplex``.

        EXAMPLES::

            sage: T = simplicial_complexes.Torus(); T
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 14 facets

        Inside the torus, define a subcomplex consisting of a loop::

            sage: S = SimplicialComplex([[0,1], [1,2], [0,2]], is_mutable=False)
            sage: S.homology()
            {0: 0, 1: Z}
            sage: L = T._enlarge_subcomplex(S)
            sage: L
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 8 facets
            sage: L.facets()
            {(0, 1, 5), (1, 3, 6), (1, 2), (1, 2, 4), (1, 3, 4), (0, 2), (1, 5, 6), (0, 1)}
            sage: L.homology()[1]
            Z
        """
        # Make the subcomplex immutable if not
        if subcomplex is not None and subcomplex._is_mutable:
            subcomplex = SimplicialComplex(subcomplex._facets,
                                           maximality_check=False,
                                           sort_facets=False,
                                           is_mutable=False)

        if subcomplex in self.__enlarged:
            return self.__enlarged[subcomplex]
        faces = filter(lambda x: x not in subcomplex._facets, list(self._facets))
        done = False
        new_facets = list(subcomplex._facets)
        while not done:
            done = True
            remove_these = []
            if verbose:
                print "  looping through %s facets" % len(faces)
            for f in faces:
                f_set = f.set()
                int_facets = set( a.set().intersection(f_set) for a in new_facets )
                intersection = SimplicialComplex(int_facets)
                if not intersection._facets[0].is_empty():
                    if (len(intersection._facets) == 1 or
                        intersection == intersection._contractible_subcomplex()):
                        new_facets.append(f)
                        remove_these.append(f)
                        done = False
            if verbose and not done:
                print "    added %s facets" % len(remove_these)
            for f in remove_these:
                faces.remove(f)
        if verbose:
            print "  now constructing a simplicial complex with %s vertices and %s facets" % (self.vertices().dimension()+1, len(new_facets))
        L = SimplicialComplex(new_facets, maximality_check=False,
                              sort_facets=False, is_mutable=self._is_mutable)
        self.__enlarged[subcomplex] = L
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

        .. [BP2000] V. M. Bukhshtaber and T. E. Panov, "Moment-angle complexes
           and combinatorics of simplicial manifolds," *Uspekhi
           Mat. Nauk* 55 (2000), 171--172.

        .. [SS1992] M. A. Shtan'ko and and M. I. Shtogrin, "Embedding cubic
           manifolds and complexes into a cubic lattice", *Uspekhi
           Mat. Nauk* 47 (1992), 219-220.

        EXAMPLES::

            sage: T = simplicial_complexes.Torus()
            sage: T.homology()
            {0: 0, 1: Z x Z, 2: Z}
            sage: Tc = T._cubical_()
            sage: Tc
            Cubical complex with 42 vertices and 168 cubes
            sage: Tc.homology()
            {0: 0, 1: Z x Z, 2: Z}
        """
        from sage.homology.cubical_complex import CubicalComplex
        V = self.vertices()
        embed = V.dimension() + 1
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
                        cube.append([0,])
                    elif n not in J:
                        cube.append([1,])
                    else:
                        cube.append([0,1])
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
            sage: v0 = X.vertices()[0]
            sage: v1 = X.vertices()[-1]
            sage: X.connected_component(Simplex([v0])) == X.connected_component(Simplex([v1]))
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
            ValueError: the empty simplicial complex has no connected components.
        """
        if self.dimension() == -1:
            raise ValueError("the empty simplicial complex has no connected components.")
        if simplex is None:
            v = self.vertices()[0]
        else:
            v = simplex[0]
        vertices = self.graph().connected_component_containing_vertex(v)
        facets = [f for f in self.facets() if f.is_face(Simplex(vertices))]
        return SimplicialComplex(facets)

    def fundamental_group(self, base_point=None, simplify=True):
        r"""
        Return the fundamental group of this simplicial complex.

        INPUT:

        - ``base_point`` (optional, default None) -- if this complex is
          not path-connected, then specify a vertex; the fundamental
          group is computed with that vertex as a base point. If the
          complex is path-connected, then you may specify a vertex or
          leave this as its default setting of ``None``. (If this
          complex is path-connected, then this argument is ignored.)

        - ``simplify`` (bool, optional True) -- if False, then return a
          presentation of the group in terms of generators and
          relations. If True, the default, simplify as much as GAP is
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
            sage: S1.fundamental_group()
            Finitely presented group < e |  >

        If we pass the argument ``simplify=False``, we get generators and
        relations in a form which is not usually very helpful. Here is the
        cyclic group of order 2, for instance::

            sage: RP2 = simplicial_complexes.RealProjectiveSpace(2)
            sage: C2 = RP2.fundamental_group(simplify=False)
            sage: C2
            Finitely presented group < e0, e1, e2, e3, e4, e5, e6, e7, e8, e9 | e6, e5, e3, e9, e4*e7^-1*e6, e9*e7^-1*e0, e0*e1^-1*e2, e5*e1^-1*e8, e4*e3^-1*e8, e2 >
            sage: C2.simplified()
            Finitely presented group < e0 | e0^2 >

        This is the same answer given if the argument ``simplify`` is True
        (the default)::

            sage: RP2.fundamental_group()
            Finitely presented group < e0 | e0^2 >

        You must specify a base point to compute the fundamental group
        of a non-connected complex::

            sage: K = S1.disjoint_union(RP2)
            sage: K.fundamental_group()
            Traceback (most recent call last):
            ...
            ValueError: this complex is not connected, so you must specify a base point.
            sage: v0 = list(K.vertices())[0]
            sage: K.fundamental_group(base_point=v0)
            Finitely presented group < e |  >
            sage: v1 = list(K.vertices())[-1]
            sage: K.fundamental_group(base_point=v1)
            Finitely presented group < e0 | e0^2 >

        Some other examples::

            sage: S1.wedge(S1).fundamental_group()
            Finitely presented group < e0, e1 | >
            sage: simplicial_complexes.Torus().fundamental_group()
            Finitely presented group < e0, e3 | e0*e3^-1*e0^-1*e3 >
            sage: simplicial_complexes.MooreSpace(5).fundamental_group()
            Finitely presented group < e1 | e1^5 >
        """
        if not self.is_connected():
            if base_point is None:
                raise ValueError("this complex is not connected, so you must specify a base point.")
            return self.connected_component(Simplex([base_point])).fundamental_group(simplify=simplify)

        from sage.groups.free_group import FreeGroup
        from sage.interfaces.gap import gap
        spanning_tree = [e[:2] for e in self.graph().min_spanning_tree()]
        gens = [tuple(e) for e in self.n_cells(1) if tuple(e) not in spanning_tree]

        if len(gens) == 0:
            return gap.TrivialGroup()

        gens_dict = dict(zip(gens, range(len(gens))))
        FG = FreeGroup(len(gens), 'e')
        rels = []
        for f in self.n_cells(2):
            bdry = [tuple(e) for e in f.faces()]
            z = dict()
            for i in range(3):
                if bdry[i] in spanning_tree:
                    z[i] = FG.one()
                else:
                    z[i] = FG.gen(gens_dict[bdry[i]])
            rels.append(z[0]*z[1].inverse()*z[2])
        if simplify:
            return FG.quotient(rels).simplified()
        else:
            return FG.quotient(rels)

    def category(self):
        """
        Return the category to which this simplicial complex belongs: the
        category of all simplicial complexes.

        EXAMPLES::

            sage: SimplicialComplex([[0,1], [1,2,3,4,5]]).category()
            Category of simplicial complexes
        """
        import sage.categories.all
        return sage.categories.all.SimplicialComplexes()

    def is_isomorphic(self,other, certify = False):
        r"""
        Checks whether two simplicial complexes are isomorphic

        INPUT:

        - ``certify`` - if ``True``, then output is ``(a,b)``, where ``a``
          is a boolean and ``b`` is either a map or ``None``.

        This is done by creating two graphs and checking whether they
        are isomorphic.

        EXAMPLES::

            sage: Z1 = SimplicialComplex([[0,1],[1,2],[2,3,4],[4,5]])
            sage: Z2 = SimplicialComplex([['a','b'],['b','c'],['c','d','e'],['e','f']])
            sage: Z3 = SimplicialComplex([[1,2,3]])
            sage: Z1.is_isomorphic(Z2)
            True
            sage: Z1.is_isomorphic(Z2, certify=True)
            (True, {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f'})
            sage: Z3.is_isomorphic(Z2)
            False
        """
        g1 = Graph()
        g2 = Graph()
        g1.add_edges((v,f) for f in self.facets() for v in f)
        g2.add_edges((v,f) for f in other.facets() for v in f)
        g1.add_edges(("fake_vertex",v,"special_edge") for v in self.vertices())
        g2.add_edges(("fake_vertex",v,"special_edge") for v in other.vertices())
        if not certify:
            return g1.is_isomorphic(g2)
        isisom, tr = g1.is_isomorphic(g2, certify = True)

        if isisom:
            for f in self.facets():
                tr.pop(f)
            tr.pop("fake_vertex")

        return isisom,tr

    def automorphism_group(self):
        r"""
        Returns the automorphism group of the simplicial complex

        This is done by creating a bipartite graph, whose vertices are
        vertices and facets of the simplicial complex, and computing
        its automorphism group.

        .. WARNING::

            Since :trac:`14319` the domain of the automorphism group is equal to
            the graph's vertex set, and the ``translation`` argument has become
            useless.

        EXAMPLES::

            sage: S = simplicial_complexes.Simplex(3)
            sage: S.automorphism_group().is_isomorphic(SymmetricGroup(4))
            True

            sage: P = simplicial_complexes.RealProjectivePlane()
            sage: P.automorphism_group().is_isomorphic(AlternatingGroup(5))
            True

            sage: Z = SimplicialComplex([['1','2'],['2','3','a']])
            sage: Z.automorphism_group().is_isomorphic(CyclicPermutationGroup(2))
            True
            sage: group = Z.automorphism_group()
            sage: group.domain()
            {'1', '2', '3', 'a'}
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup

        G = Graph()
        G.add_vertices(self.vertices())
        G.add_edges((f.tuple(),v) for f in self.facets() for v in f)
        groupe = G.automorphism_group(partition=[list(self.vertices()),
                                                 [f.tuple() for f in self.facets()]])


        gens = [ [tuple(c) for c in g.cycle_tuples() if not isinstance(c[0],tuple)]
                 for g in groupe.gens()]

        permgroup = PermutationGroup(gens = gens, domain = self.vertices())

        return permgroup

    def _Hom_(self, other, category=None):
        """
        Return the set of simplicial maps between simplicial complexes
        ``self`` and ``other``.

        EXAMPLES::

            sage: S = simplicial_complexes.Sphere(1)
            sage: T = simplicial_complexes.Sphere(2)
            sage: H = Hom(S,T)  # indirect doctest
            sage: H
            Set of Morphisms from Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)} in Category of simplicial complexes
            sage: f = {0:0,1:1,2:3}
            sage: x = H(f)
            sage: x
            Simplicial complex morphism {0: 0, 1: 1, 2: 3} from Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 2), (0, 1)} to Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
        """
        from sage.homology.simplicial_complex_homset import SimplicialComplexHomset
        return SimplicialComplexHomset(self, other)

    # @cached_method    when we switch to immutable SimplicialComplex
    def _is_numeric(self):
        """
        Test whether all vertices are labeled by integers

        OUTPUT:

        Boolean. Whether all vertices are labeled by (not necessarily
        consecutive) integers.

        EXAMPLES::

            sage: s = SimplicialComplex()
            sage: s._is_numeric()
            True
            sage: s.add_face(['a', 'b', 123])
            sage: s._is_numeric()
            False
        """
        return all([isinstance(v, (int, Integer, long)) for v in self._vertex_set])

    # @cached_method    when we switch to immutable SimplicialComplex
    def _translation_to_numeric(self):
        """
        Return a dictionary enumerating the vertices

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
        return dict((vertex, i) for i, vertex in enumerate(self._vertex_set))

    # @cached_method    when we switch to immutable SimplicialComplex
    def _translation_from_numeric(self):
        """
        Return a dictionary mapping vertex indices to vertices

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
        return dict(enumerate(self._vertex_set))

    def _chomp_repr_(self):
        r"""
        String representation of ``self`` suitable for use by the CHomP
        program.  This lists each facet on its own line, and makes
        sure vertices are listed as numbers.

        EXAMPLES::

            sage: S = SimplicialComplex([(0,1,2), (2,3,5)])
            sage: print S._chomp_repr_()
            (2, 3, 5)
            (0, 1, 2)

        A simplicial complex whose vertices are tuples, not integers::

            sage: S = SimplicialComplex([[(0,1), (1,2), (3,4)]])
            sage: S._chomp_repr_()
            '(0, 1, 2)\n'
        """
        s = ""
        numeric = self._is_numeric()
        if not numeric:
            d = self._translation_to_numeric()
        for f in self.facets():
            if numeric:
                s += str(f)
            else:
                s += '(' + ', '.join(str(d[a]) for a in f) + ')'
            s += '\n'
        return s

    # this function overrides the standard one for GenericCellComplex,
    # because it lists the maximal faces, not the total number of faces.
    def _repr_(self):
        """
        Print representation.

        If there are only a few vertices or faces, they are listed. If
        there are lots, the number is given.

        EXAMPLES::

            sage: X = SimplicialComplex([[0,1], [1,2]])
            sage: X._repr_()
            'Simplicial complex with vertex set (0, 1, 2) and facets {(1, 2), (0, 1)}'
            sage: SimplicialComplex([[i for i in range(16)]])
            Simplicial complex with 16 vertices and 1 facets
        """
        vertex_limit = 45
        facet_limit = 55
        vertices = self.vertices()
        facets = Set(self._facets)
        vertex_string = "with vertex set %s" % vertices
        if len(vertex_string) > vertex_limit:
            vertex_string = "with %s vertices" % str(1+vertices.dimension())
        facet_string = "facets %s" % facets
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
        self._is_mutable = False
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
        return self._is_mutable

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
        return not self._is_mutable

