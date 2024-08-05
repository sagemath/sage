r"""
Quiver Representations

AUTHORS:

- Jim Stark (2012-03-04): Initial implementation of acyclic quivers without
  relations.
- Simon King (2013-05, 2014-02): Split code up. Allow cyclic quivers where
  possible.

A quiver is a directed graph used for representation theory. In our
representation theoretic code, it is assumed that

- the vertices of the quiver are labelled by integers, and

- each edge of the quiver is labelled with a nonempty string.  The label
  cannot begin with ``'e_'`` or contain ``'*'`` and distinct edges must have
  distinct labels.

As far as the :class:`~sage.graphs.digraph.DiGraph` class is concerned, a
path is a finite list of pairwise distinct vertices `v_1, ..., v_n` such
that there exists an edge from `v_i` to `v_{i + 1}`.  If there are
multiple edges between the same two vertices this does not contribute
additional paths as listed by the :class:`~sage.graphs.digraph.DiGraph`
class; for example only two paths are listed from 1 to 3 in ``Q``::

    sage: Q = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
    sage: Q.edges(sort=True)
    [(1, 2, 'a'), (1, 2, 'b'), (1, 3, 'c'), (2, 3, 'd')]
    sage: Q.all_paths(1, 3)
    [[1, 2, 3], [1, 3]]

The notion of a path in a quiver (in representation theory) is
fundamentally different in several aspects. First, paths are no longer
required to have distinct vertices, or even distinct edges; thus, "path"
in quiver theory is closer to the notion of "walk" in graph theory.
Furthermore, paths in quiver theory "know" their edges, so parallel edges
between the same two vertices of a Quiver make different paths. But
paths in quiver theory also "know" their vertices, so that a length-`0`
path from `a` to `a` is not the same as a length-`0` path from `b` to `b`
for `a \neq b`.
Formally, we say that a path is given by two vertices, ``start`` and
``end``, and a finite (possibly empty) list of edges `e_1, e_2, \ldots, e_n`
such that the initial vertex of `e_1` is ``start``, the final vertex of `e_i`
is the initial vertex of `e_{i + 1}`, and the final vertex of `e_n` is
``end``.  In the case where no edges are specified, we must have
``start = end`` and the path is called the trivial path at the given vertex.

Quiver paths in the sense stated above correspond to the elements of a
partial semigroup, with multiplication of paths given by concatenation. Hence,
rather than overloading the method name inherited from
:class:`~sage.graphs.digraph.DiGraph` or inventing a
new method name, we move this functionality to this so-called *path
semigroup*.  Note that with this definition there are three paths from 1 to 3
in our example::

    sage: Q.path_semigroup().all_paths(1, 3)
    [a*d, b*d, c]

The returned paths are of type :class:`~sage.quivers.paths.QuiverPath`, which
are elements in the path semigroup that is associated with the quiver (a
partial semigroup, which does not generally have a neutral element). You can
specify a :class:`~sage.quivers.paths.QuiverPath` by giving an edge or a
list of edges, passed as arguments to the path semigroup containing this path.
Here an edge is a tuple of the form ``(i, j, l)``, where ``i`` and ``j``
are vertices and ``l`` is the label of an edge from i to j::

    sage: p = Q.path_semigroup()([(1, 2, 'a'), (2, 3, 'd')])
    sage: p
    a*d

Trivial paths are indicated by passing a list containing the tuple ``(vertex, vertex)``::

    sage: Q.path_semigroup()([(3, 3)])
    e_3

Here is an alternative way to define a path::

    sage: PQ = Q.path_semigroup()
    sage: q = PQ(['a', 'd'])
    sage: p == q
    True

If the vertices along the path do not match, a value error is raised::

    sage: inv1 = PQ([(2, 3, 'd'), (1, 2, 'a')])
    Traceback (most recent call last):
    ...
    ValueError: edge d ends at 3, but edge a starts at 1
    sage: inv2 = PQ([(1, 2, 'a'), (1, 2, 'a')])
    Traceback (most recent call last):
    ...
    ValueError: edge a ends at 2, but edge a starts at 1
    sage: inv3 = PQ([(1, 2, 'x')])
    Traceback (most recent call last):
    ...
    ValueError: (1, 2, 'x') is not an edge

The ``*`` operator is concatenation of paths. If the two paths do not compose,
then the result is ``None`` (whence the "partial" in "partial semigroup").  ::

    sage: print(p*q)
    None

Let us now construct a larger quiver::

    sage: Qbig = DiGraph({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}, 3:{4:['e']}, 4:{5:['f']}, 5:{1:['g']} })
    sage: Pbig = Qbig.path_semigroup()

Since ``Q`` is a sub-digraph of ``Qbig``, we have a coercion of the associated
path semigroups::

    sage: Pbig.has_coerce_map_from(PQ)
    True

In particular, ``p`` is considered to be an element of ``Pbig``, and can be
composed with paths that were defined for the larger quiver::

    sage: p in Pbig
    True
    sage: p*Pbig([(3, 4, 'e')])
    a*d*e
    sage: Pbig([(4, 5, 'f'), (5, 1, 'g')])*p
    f*g*a*d

The length of a path is the number of edges in that path::

    sage: len(p)
    2
    sage: triv = PQ([(1, 1)])
    sage: len(triv)
    0

List index and slice notation can be used to access the edges in a path.
QuiverPaths can also be iterated over.  Trivial paths have no elements::

    sage: for x in p: print(x)
    (1, 2, 'a')
    (2, 3, 'd')
    sage: triv[:]
    e_1

There are methods giving the initial and terminal vertex of a path::

    sage: p.initial_vertex()
    1
    sage: p.terminal_vertex()
    3

:class:`~sage.quivers.paths.QuiverPath` form the basis of the quiver
algebra of a quiver.  Given a field `k` and a quiver `Q`, the quiver
algebra `kQ` is, as a vector space, the free `k`-vector space whose basis
is the set of all paths in `Q`. Multiplication is defined on this basis
and extended bilinearly.  The product of two basis elements is given by
path composition when it makes sense and is set to be zero otherwise.
Specifically, if the terminal vertex of the left path equals the initial
vertex of the right path, then their product is the concatenation of the
two paths, and otherwise their product is zero. In sage, quiver algebras
are handled by the :class:`~sage.quivers.algebra.QuiverAlgebra` class::

    sage: A = PQ.algebra(GF(7))
    sage: A
    Path algebra of Multi-digraph on 3 vertices over Finite Field of size 7

Quivers have a method that creates their algebra over a given field (or,
more generally, commutative ring).  Note that
:class:`~sage.quivers.algebra.QuiverAlgebras` are uniquely defined by
their quiver and field, and play nicely with coercions of the underlying
path semigroups::

    sage: A is PQ.algebra(GF(7))
    True
    sage: A is PQ.algebra(RR)
    False
    sage: Q1 = Q.copy()
    sage: Q1.add_vertex(4)
    sage: PQ1 = Q1.path_semigroup()
    sage: A is PQ1.algebra(GF(7))
    False
    sage: Pbig.algebra(GF(7)).has_coerce_map_from(A)
    True

The :class:`~sage.quivers.algebra.QuiverAlgebra` can create elements
from :class:`~sage.quivers.paths.QuiverPaths` or from elements of the
base ring::

    sage: A(5)
    5*e_1 + 5*e_2 + 5*e_3
    sage: r = PQ([(1, 2, 'b'), (2, 3, 'd')])
    sage: e2 = PQ([(2, 2)])
    sage: x = A(p) + A(e2)
    sage: x
    a*d + e_2
    sage: y = A(p) + A(r)
    sage: y
    b*d + a*d

:class:`~sage.quivers.algebra.QuiverAlgebras` are `\NN`-graded algebras.
The grading is given by assigning to each basis element the length of the
path corresponding to that basis element::

    sage: x.is_homogeneous()
    False
    sage: x.degree()
    Traceback (most recent call last):
    ...
    ValueError: element is not homogeneous
    sage: y.is_homogeneous()
    True
    sage: y.degree()
    2
    sage: A[1]
    Free module spanned by [a, b, c, d] over Finite Field of size 7
    sage: A[2]
    Free module spanned by [a*d, b*d] over Finite Field of size 7

The category of right modules over a given quiver algebra is equivalent to the
category of representations of that quiver.  A quiver representation is a
diagram in the category of vector spaces whose underlying graph is the quiver.
So to each vertex of the quiver we assign a vector space and to each edge of
the quiver a linear map between the vector spaces assigned to the start and end
vertices of that edge.  To create the zero representation we just specify the
base ring and the path semigroup::

    sage: Z = Q1.path_semigroup().representation(GF(5))
    sage: Z.is_zero()
    True

To each vertex of a quiver there is associated a simple module, an
indecomposable projective, and an indecomposable injective, and these can
be created from the quiver::

    sage: S = PQ.S(GF(3), 1)
    sage: I = PQ.I(QQ, 2)
    sage: P = PQ.P(GF(3), 1)

Radicals, socles, tops, and quotients can all be computed and we can test if
modules are simple or semisimple, get their dimension, and test for equality.
Like quivers, :class:`~sage.quivers.representation.QuiverRep` objects
are unique and therefore equal if and only if they are identical::

    sage: P.is_simple()
    False
    sage: P.dimension()
    6
    sage: R = P.radical()
    sage: P.socle()
    Representation with dimension vector (0, 0, 3)
    sage: (P/R).is_simple()
    True
    sage: P == R
    False
    sage: P.top() is P/R
    True

There are special methods to deal with modules that are given as right
ideals in the quiver algebra.  To create such a module pass the keyword
``option='paths'`` along with a path or list of paths that generate the
desired ideal::

    sage: M = PQ.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='paths')
    sage: M.dimension_vector()
    (1, 2, 3)

There are also special methods to deal with modules that are given as the
linear dual of a right ideal in the quiver algebra.  To create such a
module, pass the keyword ``option='dual paths'`` to the constructor along
with a path or list of paths.  The module returned is the dual of the
ideal created in the opposite quiver by the reverses of the given paths::

    sage: D = PQ.representation(QQ, [[(1, 1)], [(1, 2, 'a')]], option='dual paths')
    sage: D.dimension_vector()
    (2, 0, 0)

For modules that are not a standard module or an ideal of the quiver algebra
:class:`~sage.quivers.representation.QuiverRep` can take as input two
dictionaries.  The first associates to each vertex a vector space or an
integer (the desired dimension of the vector space), the second associates to
each edge a map or a matrix or something from which sage can construct a map::

    sage: PQ2 = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
    sage: M2 = PQ2.representation(QQ, {1: QQ^2, 2: QQ^1}, {(1, 2, 'a'): [1, 0], (1, 2, 'b'): [0, 1]})
    sage: M.get_space(2)
    Vector space of dimension 2 over Rational Field
    sage: M2.get_space(2)
    Vector space of dimension 1 over Rational Field
    sage: M.get_map((1, 2, 'a'))
    Vector space morphism represented by the matrix:
    [1 0]
    Domain: Vector space of dimension 1 over Rational Field
    Codomain: Vector space of dimension 2 over Rational Field

A homomorphism between two quiver representations is given by homomorphisms
between the spaces assigned to the vertices of those representations such that
those homomorphisms commute with the edge maps of the representations. The
homomorphisms are created in the usual Sage syntax, the defining data given by
a dictionary associating maps to vertices::

    sage: P2 = PQ2.P(QQ, 1)
    sage: f = P2.hom({1:[1, 1], 2:[[1], [1]]}, M2)

When the domain is given as a right ideal in the quiver algebra we can also
create a homomorphism by just giving a single element in the codomain.
The map is then induced by acting on that element::

    sage: x = P2.gens('x')[0]
    sage: x
    x_0
    sage: f == P2.hom(f(x), M2)
    True

As you can see, the above homomorphisms can be applied to elements.  Just
like elements, addition is defined via the ``+`` operator.  On elements scalar
multiplication is defined via the ``*`` operator but on homomorphisms ``*``
defines composition, so scalar multiplication is done using a method::

    sage: g = f + f
    sage: g == f.scalar_mult(2)
    True
    sage: g == 2*f       # This multiplies the map with the scalar 2
    True
    sage: g(x) == 2*f(x) # This applies the map, then multiplies by the scalar
    True

The ``direct_sum`` method for modules returns only the resulting module by
default.  But can also return the projection and inclusion homomorphisms into
the various factors::

    sage: N2, inclusions, projections = M2.direct_sum([P2], return_maps=True)
    sage: inclusions[0].domain() is M2
    True
    sage: projections[0].codomain() is M2
    True
    sage: (projections[0]*inclusions[0]).is_isomorphism()
    True

As you see above we can determine if a given map is an isomorphism.  Testing
for injectivity and surjectivity works as well::

    sage: f.is_injective()
    False
    sage: f.is_surjective()
    False

We can create all the standard modules associated to maps::

    sage: f.kernel()
    Representation with dimension vector (0, 1)
    sage: f.cokernel()
    Representation with dimension vector (1, 0)
    sage: im = f.image()
    sage: im
    Representation with dimension vector (1, 1)

These methods, as well as the ``submodule`` and ``quotient`` methods that are
defined for representations, return only the resulting representation.  To get
the inclusion map of a submodule or the factor homomorphism of a quotient use
``coerce_map_from``::

    sage: incl = M2.coerce_map_from(im)
    sage: incl.domain() is im
    True
    sage: incl.codomain() is M2
    True
    sage: incl.is_injective()
    True

Both :class:`~sage.quivers.representation.QuiverRep` objects and
:class:`~sage.quivers.homspace.QuiverRepHom` objects have ``linear_dual`` and
``algebraic_dual`` methods.  The ``linear_dual`` method applies the functor
`Hom_k(..., k)` where `k` is the base ring of the representation, and the
``algebraic_dual`` method applies the functor `Hom_Q(..., kQ)` where `kQ`
is the quiver algebra.  Both these functors yield left modules.  A left
module is equivalent to a right module over the opposite algebra, and the
opposite of a quiver algebra is the algebra of the opposite quiver, so both
these methods yield modules and representations of the opposite quiver::

    sage: f.linear_dual()
    Homomorphism of representations of Reverse of (): Multi-digraph on 2 vertices
    sage: D = M2.algebraic_dual()
    sage: D.quiver() is PQ2.reverse().quiver()
    True

.. TODO::

    Change the wording ``Reverse of ()`` into something more meaningful.

There is a method returning the projective cover of any module.  Note that
this method returns the homomorphism; to get the module take the domain of
the homomorphism::

    sage: cov = M2.projective_cover()
    sage: cov
    Homomorphism of representations of Multi-digraph on 2 vertices
    sage: cov.domain()
    Representation with dimension vector (2, 4)

As projective covers are computable, so are the transpose and Auslander-Reiten
translates of modules::

    sage: M2.transpose()
    Representation with dimension vector (4, 3)
    sage: PQ2.I(QQ, 1).AR_translate()
    Representation with dimension vector (3, 2)

We have already used the ``gens`` method above to get an element of a quiver
representation.  An element of a quiver representation is simply a choice of
element from each of the spaces assigned to the vertices of the quiver.
Addition, subtraction, and scalar multiplication are performed pointwise and
implemented by the usual operators::

    sage: M2.dimension_vector()
    (2, 1)
    sage: x, y, z = M2.gens('xyz')
    sage: 2*x + y != x + 2*y
    True

To create a specific element of a given representation we just specify the
representation and a dictionary associating to each vertex an element of the
space associated to that vertex in the representation::

    sage: w = M2({1:(1, -1), 2:(3,)})
    sage: w.get_element(1)
    (1, -1)

The right action of a quiver algebra on an element is implemented via
the ``*`` operator::

    sage: A2 = x.quiver().path_semigroup().algebra(QQ)
    sage: a = A2('a')
    sage: x*a == z
    True
"""

# ****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from sage.structure.factory import UniqueFactory
from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById


class QuiverRepFactory(UniqueFactory):
    r"""
    A quiver representation is a diagram in the category of vector spaces whose
    underlying graph is the quiver.  Giving a finite dimensional representation
    is equivalent to giving a finite dimensional right module for the path
    algebra of the quiver.

    INPUT:

    The first two arguments specify the base ring and the quiver,
    and they are always required:

    - ``k`` -- ring, the base ring of the representation

    - ``P`` -- the partial semigroup formed by the paths of the quiver of the
      representation

    Then to specify the spaces and maps associated to the quiver
    there are three possible options.  The first is the ``'values'`` option,
    where the next two arguments give the data to be assigned.  The following
    can either be the next two entries in the argument list or they can be
    passed by keyword.  If the argument list is long enough the keywords
    are ignored; the keywords are only checked in the event that the argument
    list does not have enough entries after ``P``.

    - ``spaces`` -- dictionary (default: empty) associating to each vertex a
      free module over the base ring `k`.  Not all vertices must be specified;
      unspecified vertices are automatically set to `k^0`.  Keys of the
      dictionary  that don't correspond to vertices are ignored.

    - ``maps`` -- dictionary (default: empty) associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges must be
      specified; unspecified edges are automatically set to the zero map. Keys
      of the dictionary that don't correspond to edges are ignored.

    The second option is the ``paths`` option which creates a module by
    generating a right ideal from a list of paths.  Thus the basis elements
    of this module correspond to paths of the quiver and the maps are given
    by right multiplication by the corresponding edge.  As above this can be
    passed either as the next entry in the argument list or as a keyword.
    The keyword is only checked if there is no entry in the argument list
    after ``Q``.

    - ``basis`` -- list; a nonempty list of paths in the quiver ``Q``.
      Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  There must be at least one valid path in the list
      or a :exc:`ValueError` is raised.  The closure of this list under right
      multiplication forms the basis of the resulting representation.

    The third option is the ``dual paths`` option which creates the dual of
    a left ideal in the quiver algebra.  Thus the basis elements of this
    module correspond to paths of the quiver and the maps are given by
    deleting the corresponding edge from the start of the path (the edge map
    is zero on a path if that edge is not the initial edge of the path).
    As above this can be passed either as the next entry in the argument
    list or as a keyword.

    - ``basis`` -- list; a nonempty list of paths in the quiver ``Q``.
      Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  There must be at least one valid path in the list
      or a :exc:`ValueError` is raised.  The closure of this list under left
      multiplication of edges forms the basis of  the resulting representation.

    Using the second and third options requires that the following keyword be
    passed to the constructor.  This must be passed as a keyword.

    - ``option`` -- string (default: ``None``); either ``'values'`` or
      ``'paths'`` or ``'dual paths'``. ``None`` is equivalent to ``'values'``

    OUTPUT: :class:`QuiverRep`

    EXAMPLES::

        sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()

    When the ``option`` keyword is not supplied the constructor uses the
    ``'values'`` option and expects the spaces and maps to be specified.
    If no maps or spaces are given the zero module is created::

        sage: M = Q1.representation(GF(5))
        sage: M.is_zero()
        True

    The simple modules, indecomposable projectives, and indecomposable
    injectives are examples of quiver representations::

        sage: S = Q1.S(GF(3), 1)
        sage: I = Q1.I(QQ, 2)
        sage: P = Q1.P(GF(3), 1)

    Various standard submodules can be computed, such as radicals and socles.
    We can also form quotients and test for certain attributes such as
    semisimplicity::

        sage: R = P.radical()
        sage: R.is_zero()
        False
        sage: (P/R).is_simple()
        True
        sage: P == R
        False

    With the option ``'paths'`` the input data should be a list of
    :class:`QuiverPaths` or things that :class:`QuiverPaths` can be
    constructed from.  The resulting module is the submodule generated by
    these paths in the quiver algebra, when considered as a right module
    over itself::

        sage: P1 = Q1.representation(QQ, [[(1, 1)]], option='paths')
        sage: P1.dimension()
        2

    In the following example, the 3rd and 4th paths are actually the same,
    so the duplicate is removed::

        sage: N = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 2, 'a')],
        ....:                           [(1, 2, 'a')]], option='paths')
        sage: N.dimension()
        3

    The dimension at each vertex equals the number of paths in the closed
    basis whose terminal point is that vertex::

        sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
        sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
        sage: M.dimension_vector()
        (0, 2, 2)
        sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
        sage: N.dimension_vector()
        (0, 1, 2)
    """
    def create_key(self, k, P, *args, **kwds):
        """
        Return a key for the specified module.

        The key is a tuple.  The first and second entries are the base ring
        ``k`` and the partial semigroup ``P`` formed by the paths of a quiver.
        The third entry is the ``option`` and the remaining entries depend on
        that option.  If the option is ``'values'`` and the quiver
        has `n` vertices then the next `n` entries are the vector spaces
        to be assigned to those vertices.  After that are the matrices of
        the maps assigned to edges, listed in the same order that
        ``Q.edges(sort=True)`` uses.  If the option is ``'paths'`` or ``'dual paths'``
        then the next entry is a tuple containing a sorted list of the
        paths that form a basis of the quiver.

        INPUT:

        See the class documentation.

        OUTPUT: tuple

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: from sage.quivers.representation import QuiverRep
            sage: QuiverRep.create_key(GF(5), P)
            (Finite Field of size 5,
             Partial semigroup formed by the directed paths of Multi-digraph on 2 vertices,
             'values',
             Vector space of dimension 0 over Finite Field of size 5,
             Vector space of dimension 0 over Finite Field of size 5,
             [])
        """
        key = [k, P]
        Q = P.quiver()
        if 'option' in kwds and (kwds['option'] == 'paths' or kwds['option'] == 'dual paths'):
            # Follow the 'paths' specification for the input
            key.append(kwds['option'])
            if args:
                basis = args[0]
            else:
                basis = kwds['basis']

            # Add as QuiverPaths to a set
            paths = set()
            for p in basis:
                paths.add(P(p))

            if kwds['option'] == 'paths':
                # Close the set under right mult by edges
                edges = list(P.arrows())
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            pe = p * e
                            if pe is not None and pe not in paths:
                                to_be_added.append(pe)

                    paths.update(to_be_added)
                    just_added = to_be_added

            if kwds['option'] == 'dual paths':
                # Close the set under left mult by edges
                edges = list(P.arrows())
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            ep = e * p
                            if ep is not None and ep not in paths:
                                to_be_added.append(ep)

                    paths.update(to_be_added)
                    just_added = to_be_added

            # Add to the key
            key.append(tuple(sorted(paths)))

        else:
            # Assume the input type the 'values' option
            key.append('values')
            if args:
                spaces = args[0]
            elif 'spaces' in kwds:
                spaces = kwds['spaces']
            else:
                spaces = {}
            if len(args) > 1:
                maps = args[1]
            elif 'maps' in kwds:
                maps = kwds['maps']
            else:
                maps = {}

            # If the vertex is not specified set it as a free module of rank 0, if
            # an integer is given set it as a free module of that rank, otherwise
            # assume the object is a module and assign it to the vertex.
            from sage.rings.finite_rings.integer_mod_ring import Integers
            verts = Q.vertices(sort=True)
            for x in verts:
                if x not in spaces:
                    key.append(k**0)
                elif spaces[x] in Integers():
                    key.append(k**spaces[x])
                else:
                    key.append(spaces[x])

            # The preferred method of specifying an edge is as a tuple
            # (i, t, l) where i is the initial vertex, t is the terminal
            # vertex, and l is the label.  This is the form in which
            # quiver.edges(sort=True) and other such functions give the edge.  But here
            # edges can be specified by giving only the two vertices or giving
            # only the edge label.
            #
            # Note that the first space is assigned to key[3] and the first
            # vertex is 1 so the space assigned to vertex v is key[2 + v]
            from sage.matrix.constructor import Matrix
            from sage.categories.morphism import Morphism
            for x in P._sorted_edges:
                if x in maps:
                    e = maps[x]
                elif (x[0], x[1]) in maps:
                    e = maps[(x[0], x[1])]
                elif x[2] in maps:
                    e = maps[x[2]]
                else:
                    e = Matrix(k, key[3 + verts.index(x[0])].dimension(), key[3 + verts.index(x[1])].dimension())

                # If a morphism is specified take it's matrix.  Create one if
                # needed.  Otherwise assume the Matrix function can convert the
                # object to a Matrix.
                if isinstance(e, Morphism):
                    if hasattr(e, 'matrix'):
                        key.append(e.matrix())
                    else:
                        gens_images = [key[3 + verts.index(x[1])].coordinate_vector(e(x))
                                       for x in key[3 + verts.index(x[0])].gens()]
                        key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                          key[3 + verts.index(x[1])].dimension(), gens_images))
                else:
                    key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                      key[3 + verts.index(x[1])].dimension(), e))

                # Make sure the matrix is immutable so it hashes
                key[-1].set_immutable()

        # Wrap as a tuple and return
        return tuple(key)

    def create_object(self, version, key, **extra_args):
        """
        Create a :class:`QuiverRep_generic` or
        :class:`QuiverRep_with_path_basis` object from the key.

        The key is a tuple.  The first and second entries are the base ring
        ``k`` and the quiver ``Q``.  The third entry is the
        ``'option'`` and the remaining entries depend on that option.
        If the option is ``'values'`` and the quiver has `n`
        vertices then the next `n` entries are the vector spaces to be
        assigned to those vertices.  After that are the matrices
        of the maps assigned to edges, listed in the same order that
        ``Q.edges(sort=True)`` uses.  If the option is ``'paths'`` or ``'dual paths'``
        then the next entry is a tuple containing a sorted list of the
        paths that form a basis of the quiver.

        INPUT:

        - ``version`` -- the version of sage, this is currently ignored
        - ``key`` -- tuple

        OUTPUT: :class:`QuiverRep_generic` or :class:`QuiverRep_with_path_basis`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: from sage.quivers.representation import QuiverRep
            sage: key = QuiverRep.create_key(GF(5), Q)
            sage: QuiverRep.create_object(0, key)
            Representation with dimension vector (0, 0)
        """
        if len(key) < 4:
            raise ValueError("invalid key used in QuiverRepFactory")

        # Get the quiver
        P = key[1]
        Q = P.quiver()

        if key[2] == 'values':
            # Get the spaces
            spaces = {}
            i = 3
            for v in Q:
                spaces[v] = key[i]
                i += 1

            # Get the maps
            maps = {}
            for e in P._sorted_edges:
                maps[e] = key[i]
                i += 1

            # Create and return the module
            return QuiverRep_generic(key[0], P, spaces, maps)

        elif key[2] == 'paths':
            # Create and return the module
            return QuiverRep_with_path_basis(key[0], P, key[3])

        elif key[2] == 'dual paths':
            # Create and return the module
            return QuiverRep_with_dual_path_basis(key[0], P, key[3])

        raise ValueError("invalid key used in QuiverRepFactory")


QuiverRep = QuiverRepFactory("sage.quivers.representation.QuiverRep")

#########################################################################
#  Elements


class QuiverRepElement(ModuleElement):
    r"""
    An element of a quiver representation is a choice of element from each
    of the spaces assigned to the vertices of the quiver.  Addition,
    subtraction, and scalar multiplication of these elements is done
    pointwise within these spaces.

    INPUT:

    - ``module`` -- :class:`QuiverRep` (default: ``None``); the module to
      which the element belongs

    - ``elements`` -- dictionary (default: empty) associating to each vertex a
      vector or an object from which sage can create a vector. Not all vertices
      must be specified, unspecified vertices will be assigned the zero vector
      of the space associated to that vertex in the given module.  Keys that do
      not correspond to a vertex are ignored.

    - ``name`` -- string (default: ``None``); the name of the element

    OUTPUT: :class:`QuiverRepElement`

    .. NOTE::

        The constructor needs to know the quiver in order to create an
        element of a representation over that quiver.  The default is to
        read this information from ``module`` as well as to fill in
        unspecified vectors with the zeros of the spaces in ``module``.
        If ``module`` is ``None`` then ``quiver`` *MUST* be a quiver and each
        vertex *MUST* be specified or an error will result.  If both
        ``module`` and ``quiver`` are given then ``quiver`` is ignored.

    EXAMPLES::

        sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
        sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
        sage: M = Q.representation(GF(3), spaces)
        sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
        sage: M(elems)
        Element of quiver representation
        sage: v = M(elems, 'v')
        sage: v
        v
        sage: (v + v + v).is_zero()
        True
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, parent, elements=None, name=None):
        """
        Initialize ``self``. Type ``QuiverRepElement?`` for more information.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems, 'v')
            sage: TestSuite(v).run()
        """
        # In the default call method, the default value of the first argument is zero
        if not elements:
            elements = {}
        # The data describing an element is held in the following private
        # variables:
        #
        # * _elems
        #      A dictionary that assigns to each vertex of the quiver a choice
        #      of element from the space assigned to that vertex in the parent
        #      representation.
        # * _quiver
        #      The quiver of the representation.

        super().__init__(parent)

        self._elems = {}
        self._quiver = parent._quiver
        for v in self._quiver:
            if v in elements:
                self._elems[v] = parent._spaces[v](elements[v])
            else:
                self._elems[v] = parent._spaces[v].zero()

        # Assign a name if supplied
        if name is not None:
            self.rename(name)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: Q.P(QQ, 3).an_element() # indirect doctest
            Element of quiver representation
        """
        return "Element of quiver representation"

    def _add_(left, right):
        """
        This overrides the ``+`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: (v + v + v).is_zero() # indirect doctest
            True
        """
        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] + right._elems[v]

        return left.parent()(elements)

    def _sub_(left, right):
        """
        This overrides the ``-`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: (v - v).is_zero() # indirect doctest
            True
        """
        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] - right._elems[v]

        return left.parent()(elements)

    def _neg_(self):
        """
        This overrides the unary ``-`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v == -v # indirect doctest
            False
        """
        elements = {}
        for v in self._quiver:
            elements[v] = -self._elems[v]

        return self.parent()(elements)

    def __mul__(self, other):
        """
        Implement ``*`` for right multiplication by quiver algebra elements.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: A = Q.algebra(QQ)
            sage: m = P.an_element()
            sage: a = A('a')
            sage: e1 = A([(1, 1)])
            sage: m.support()
            [1, 2]
            sage: (m*a).support()
            [2]
            sage: (m*e1).support()
            [1]
            sage: (m*(a + e1)).support()
            [1, 2]
        """
        # Make sure the input is an element of the quiver algebra and get the
        # coefficients of the monomials in it
        parent = self.parent()
        mons = parent._actor(other).monomial_coefficients()

        # this must not be the cached parent.zero(),
        # since otherwise it gets changed in place!!
        result = parent()

        for path in mons:
            # Multiply by the scalar
            x = mons[path] * self._elems[path.initial_vertex()]

            # If the edge isn't trivial apply the corresponding maps
            if path:
                for edge in path:
                    x = parent.get_map(edge)(x)

            # Sum the results of acting by each monomial
            result._elems[path.terminal_vertex()] += x

        return result

    @cached_method
    def __hash__(self):
        """
        The hash only depends on the elements assigned to each vertex of the
        underlying quiver.

        .. NOTE::

            The default hash inherited from
            :class:`~sage.structure.element.Element` would depend on
            the name of the element and would thus be mutable.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}}).path_semigroup()
            sage: P = Q.P(QQ, 2)
            sage: v = P()   # let's not muddy P.zero(), which is cached
            sage: h1 = hash(v)
            sage: v.rename('foobar')
            sage: h1 == hash(v)
            True
        """
        return hash(frozenset((v, tuple(self._elems[v]))
                              for v in self._quiver))

    def __eq__(self, other):
        """
        This overrides the ``==`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = M(elems)
            sage: v == w
            True
            sage: v += w
            sage: v == w
            False
        """
        # Return False if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return False

        # Return False if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return False

        return True

    def __ne__(self, other):
        """
        This overrides the ``!=`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = M(elems)
            sage: v != w
            False
            sage: v += w
            sage: v != w
            True
        """
        # Return True if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return True

        # Return True if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return True

        return False

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view and modify the representation data. #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Return the quiver of the representation.

        OUTPUT: :class:`DiGraph`, the quiver of the representation

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: v = P.an_element()
            sage: v.quiver() is Q.quiver()
            True
        """
        return self._quiver

    def get_element(self, vertex):
        """
        Return the element at the given vertex.

        INPUT:

        - ``vertex`` -- integer; a vertex of the quiver

        OUTPUT: the vector assigned to the given vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v.get_element(3)
            (2, 1)
        """

        return self._elems[vertex]

    def _set_element(self, vector, vertex):
        r"""
        Set the value of the element ``self`` at the given vertex
        ``vertex`` to ``vector`` (while the values at all other vertices
        remain unchanged).

        INPUT:

        - ``vector`` -- a vector or an object from which the space
          associated to the given vertex in the parent can create a vector

        - ``vertex`` -- integer; a vertex of the quiver

        .. WARNING::

            Only use this method if you know what you are doing. In particular,
            do not apply it to an element that is cached somewhere (such as
            :meth:`zero`).

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v._set_element((1, 1), 1)
            sage: v.get_element(1)
            (1, 1)
        """
        self._elems[vertex] = self.parent()._spaces[vertex](vector)

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the homomorphism.         #
    #                                                                         #
    ###########################################################################

    def is_zero(self):
        """
        Test whether ``self`` is zero.

        OUTPUT: boolean, ``True`` if the element is the zero element, ``False``
        otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.is_zero()
            False
            sage: w = M()
            sage: w.is_zero()
            True

        TESTS::

            sage: M = Q.P(QQ, 1)
            sage: M.zero().is_zero()
            True
        """
        for v in self._quiver:
            if not self._elems[v].is_zero():
                return False

        return True

    def support(self):
        """
        Return the support of ``self`` as a list.

        The support is the set of vertices to which a nonzero vector is
        associated.

        OUTPUT: list; the support

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 0), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.support()
            [1, 3]
        """

        sup = []
        for v in self._quiver:
            if not self._elems[v].is_zero():
                sup.append(v)

        return sup

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def copy(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = v.copy()
            sage: w._set_element((0, 0), 1)
            sage: w.get_element(1)
            (0, 0)
            sage: v.get_element(1)
            (1, 0)
        """
        name = self.get_custom_name()
        return self.parent()(self._elems.copy(), name)

####################################################################
# The representations


class QuiverRep_generic(WithEqualityById, Module):
    """
    A generic quiver representation.

    This class should not be called by the user.

    Call :class:`QuiverRep` with ``option='values'`` (which is the default)
    instead.

    INPUT:

    - ``k`` -- ring, the base ring of the representation

    - ``P`` -- the path semigroup of the quiver `Q` of the representation

    - ``spaces`` -- dictionary (default: empty) associating to each vertex a
      free module over the base ring `k`.  Not all vertices need to be
      specified, unspecified vertices are automatically set to `k^0`.  Keys of
      the dictionary that don't correspond to vertices are ignored.

    - ``maps`` -- dictionary (default: empty) associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges need to be
      specified, unspecified edges are automatically set to the zero map.  Keys
      of the dictionary that don't correspond to edges are ignored.

    OUTPUT: :class:`QuiverRep`

    EXAMPLES::

        sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
        sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
        sage: M = Q.representation(QQ, spaces, maps)

    ::

        sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
        sage: P = Q.P(GF(3), 1)
        sage: I = Q.I(QQ, 1)
        sage: P.an_element() in P
        True
        sage: I.an_element() in P
        False

    TESTS::

        sage: TestSuite(M).run()
        sage: TestSuite(P).run()
        sage: TestSuite(I).run()
        sage: TestSuite(Q.S(ZZ,2)).run()
    """
    Element = QuiverRepElement

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, P, spaces, maps):
        """
        Initialize ``self``. Type ``QuiverRep?`` for more information.

        .. NOTE::

            Use :meth:`~sage.quivers.quiver.Quiver.representation` to create a
            representation. The following is just a test, but it is not the
            recommended way of creating a representation.

        TESTS::

            sage: from sage.quivers.representation import QuiverRep_generic
            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: QuiverRep_generic(GF(5), Q, {1: GF(5)^2, 2: GF(5)^3}, {})
            Representation with dimension vector (2, 3)
        """
        # This class can handle representations over
        # an arbitrary base ring, not necessarily a field, so long as sage can
        # construct free modules over that ring.  The data of a representation is held
        # in the following private variables:
        #
        # * _quiver
        #      The quiver of the representation.
        # * _base
        #       The quiver algebra of _quiver over _base_ring
        # * _base_ring
        #      The base ring of the representation.
        # * _spaces
        #      A dictionary which associates to each vertex of the quiver a free
        #      module over the base ring.
        # * _maps
        #      A dictionary which associates to each edge of the quiver a homomorphism
        #      whose domain and codomain equal the initial and terminal vertex of the
        #      edge.
        Q = P.quiver()
        self._semigroup = P
        self._actor = P.algebra(k)
        self._quiver = Q
        self._spaces = {}
        self._maps = {}

        # If the vertex is not specified set it as a free module of rank 0, if
        # an integer is given set it as a free module of that rank, otherwise
        # assume the object is a module and assign it to the vertex.
        from sage.rings.finite_rings.integer_mod_ring import Integers
        for x in Q:
            if x not in spaces:
                self._spaces[x] = k**0
            elif spaces[x] in Integers():
                self._spaces[x] = k**spaces[x]
            else:
                self._spaces[x] = spaces[x]

        # The preferred method of specifying an edge is as a tuple (i, t, l)
        # where i is the initial vertex, t is the terminal vertex, and l is the
        # label.  This is the form in which quiver.edges(sort=True) and other such
        # functions give the edge.  But here edges can be specified by giving
        # only the two vertices or giving only the edge label.
        for x in P._sorted_edges:
            if x in maps:
                e = maps[x]
            elif (x[0], x[1]) in maps:
                e = maps[(x[0], x[1])]
            elif x[2] in maps:
                e = maps[x[2]]
            else:
                e = self._spaces[x[0]].Hom(self._spaces[x[1]]).zero()

            # If a morphism is specified use it, otherwise assume the hom
            # function can convert the object to a morphism.  Matrices and the
            # zero and one of the base ring are valid inputs (one is valid only
            # when the domain and codomain are equal).
            from sage.categories.morphism import Morphism
            if isinstance(e, Morphism):
                self._maps[x] = e
            else:
                self._maps[x] = self._spaces[x[0]].hom(e, self._spaces[x[1]])

        self._assert_valid_quiverrep()

        super().__init__(k)  # Or explicitly Module.__init__(self, k)?

    def _assert_valid_quiverrep(self):
        r"""
        Raise an error if the representation is not well defined.

        Specifically it checks the map assigned to each edge.  The domain
        and codomain must equal the modules assigned to the initial and
        terminal vertices of the edge.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.P(GF(3), 2) # indirect doctest

        Due to unique representation, we will cause bugs in later code if
        we modify ``M`` to be an invalid representation.  So we make sure
        to store the original values and replace them after testing::

            sage: sv = M._spaces[1]
            sage: M._spaces[1] = 0
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: domain of map at edge 'a' does not match
            sage: M._spaces[1] = sv
            sage: M = Q.P(GF(3), 2)
            sage: sv = M._maps[(1, 2, 'a')]
            sage: M._maps[(1, 2, 'a')] = (QQ^2).Hom(QQ^1).zero()
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: domain of map at edge 'a' does not match
            sage: M._maps[(1, 2, 'a')] = sv
        """

        for x in self._semigroup._sorted_edges:
            if self._maps[x].domain() != self._spaces[x[0]]:
                raise ValueError("domain of map at edge '{}' does not match".format(x[2]))
            if self._maps[x].codomain() != self._spaces[x[1]]:
                raise ValueError("codomain of map at edge '{}' does not match".format(x[2]))

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: Q.P(GF(3), 2) # indirect doctest
            Representation with dimension vector (0, 1)
        """
        return "Representation with dimension vector {}".format(self.dimension_vector())

    def __truediv__(self, sub):
        """
        Overload the ``/`` operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(GF(3), 1)
            sage: R = P.radical()
            sage: (P/R).is_simple()
            True
        """
        return self.quotient(sub)

    def _submodule(self, spaces=None):
        """
        Return the submodule specified by the data.

        This differs from ``self.submodule`` in that it assumes the data
        correctly specifies a submodule whereas ``self.submodule`` returns
        the smallest submodule containing the data.

        TESTS::

            sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v]) # indirect doctest
            Representation with dimension vector (1, 1, 1)
            sage: M.submodule(spaces={1: QQ^2}) # indirect doctest
            Representation with dimension vector (2, 0, 2)
            sage: M.submodule(spaces={2: QQ^3}) # indirect doctest
            Representation with dimension vector (0, 3, 1)
            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3}) # indirect doctest
            Representation with dimension vector (1, 3, 1)
            sage: M.submodule().is_zero() # indirect doctest
            True
            sage: M.submodule(M.gens()) is M # indirect doctest
            True
        """
        if spaces is None:
            spaces = {}
        # Add zero submodules
        for v in self._quiver:
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

        # Create edge homomorphisms restricted to the new domains and codomains
        maps = {}
        for e in self._semigroup._sorted_edges:
            maps[e] = self._maps[e].restrict_domain(spaces[e[0]]).restrict_codomain(spaces[e[1]])

        return self._semigroup.representation(self.base_ring(), spaces, maps)

    def _coerce_map_from_(self, domain):
        """
        Return either a QuiverRepHom from ``domain`` to ``self``, or
        ``False``.

        .. NOTE::

            This function simply tries to coerce a map at each vertex and
            then check if the result is a valid homomorphism.  If it is,
            then that homomorphism is returned.  If it is not or if no
            coercion was possible then it returns ``False``.

        INPUT:

        - ``domain`` -- a Sage object

        OUTPUT: :class:`QuiverRepHom` or boolean

        TESTS::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: S = M.radical()
            sage: M.coerce_map_from(S) # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices
            sage: (M/S).coerce_map_from(M) # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices

        Here sage coerces a map but the result is not a homomorphism::

            sage: S.coerce_map_from(M) # indirect doctest

        Here sage cannot coerce a map::

            sage: N = Q.P(QQ, 3)
            sage: N.coerce_map_from(M) # indirect doctest
        """

        # Domain must be a QuiverRep
        if not isinstance(domain, QuiverRep_generic):
            return False

        # Coerce a map at each vertex, return false if it fails
        maps = {}
        for v in self._quiver:
            maps[v] = self._spaces[v].coerce_map_from(domain._spaces[v])
            if maps[v] is None or maps[v] is False:
                return False

        # Create and return the hom, return False if it wasn't valid
        try:
            return domain.hom(maps, self)
        except ValueError:
            return False

    def _Hom_(self, codomain, category):
        """
        This function is used by the coercion model.

        INPUT:

        - ``codomain`` -- :class:`QuiverRepHom`

        - ``category`` -- this input is (currently) ignored

        EXAMPLES::

            sage: from sage.categories.homset import Hom
            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c', 'd']}}).path_semigroup()
            sage: P = Q.P(GF(3), 2)
            sage: S = P/P.radical()
            sage: Hom(P, S) # indirect doctest
            Dimension 1 QuiverHomSpace
        """
        from sage.quivers.homspace import QuiverHomSpace

        if not isinstance(codomain, QuiverRep_generic):
            raise TypeError("codomain must be of type QuiverRep_generic")

        return QuiverHomSpace(self, codomain)

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view the representation data.            #
    #                                                                         #
    ###########################################################################

    def get_space(self, vertex):
        """
        Return the module associated to the given vertex ``vertex``.

        INPUT:

        - ``vertex`` -- integer; a vertex of the quiver of the module

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}}).path_semigroup()
            sage: Q.P(QQ, 1).get_space(1)
            Vector space of dimension 1 over Rational Field
        """
        return self._spaces[vertex]

    def get_map(self, edge):
        """
        Return the map associated to the given edge ``edge``.

        INPUT:

        - ``edge`` -- tuple of the form
          (initial vertex, terminal vertex, label) specifying the edge
          whose map is returned

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.P(ZZ, 1).get_map((1, 2, 'a'))
            Free module morphism defined by the matrix
            [1 0]
            Domain: Ambient free module of rank 1 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """

        return self._maps[edge]

    def quiver(self):
        """
        Return the quiver of the representation.

        OUTPUT: :class:`DiGraph`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.representation(GF(5))
            sage: M.quiver() is Q.quiver()
            True
        """
        return self._quiver

    def actor(self):
        r"""
        Return the quiver path algebra acting on this representation.

        OUTPUT: a quiver path algebra

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.representation(GF(5))
            sage: M.base_ring()
            Finite Field of size 5
            sage: M.actor()
            Path algebra of Multi-digraph on 2 vertices over Finite Field of size 5
        """
        return self._actor

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the representation.       #
    #                                                                         #
    ###########################################################################

    def dimension(self, vertex=None):
        """
        Return the dimension of the space associated to the given vertex
        ``vertex``.

        INPUT:

        - ``vertex`` -- integer or ``None`` (default: ``None``), the given
          vertex

        OUTPUT:

        - integer, the dimension over the base ring of the space
          associated to the given vertex.  If ``vertex=None`` then the
          dimension over the base ring of the module is returned

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension(1)
            1
            sage: P.dimension(2)
            2

        The total dimension of the module is the sum of the dimensions
        at each vertex::

            sage: P.dimension()
            3
        """
        if vertex is None:
            # Sum the dimensions of each vertex
            dim = 0
            for x in self._quiver:
                dim += self._spaces[x].dimension()
            return dim
        else:
            # Return the dimension of just the one vertex
            return self._spaces[vertex].dimension()

    def dimension_vector(self):
        """
        Return the dimension vector of the representation.

        OUTPUT: tuple

        .. NOTE::

            The order of the entries in the tuple matches the order given
            by calling the ``vertices()`` method on the quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension_vector()
            (1, 2)

        Each coordinate of the dimension vector is the dimension of the space
        associated to that coordinate::

            sage: P.get_space(2).dimension()
            2
        """
        return tuple(self._spaces[x].dimension() for x in self._quiver)

    def is_zero(self):
        """
        Test whether the representation is zero.

        OUTPUT: boolean

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(ZZ)
            sage: N = Q.representation(ZZ, {1: 1})
            sage: M
            Representation with dimension vector (0, 0)
            sage: N
            Representation with dimension vector (1, 0)
            sage: M.is_zero()
            True
            sage: N.is_zero()
            False
        """
        return self.dimension() == 0

    def is_simple(self):
        """
        Test whether the representation is simple.

        OUTPUT: boolean

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.P(RR, 1).is_simple()
            False
            sage: Q.S(RR, 1).is_simple()
            True
        """
        # A module for an acyclic quiver is simple if and only if it has total
        # dimension 1.
        return self.dimension() == 1

    def is_semisimple(self):
        """
        Test whether the representation is semisimple.

        OUTPUT: boolean

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: (M/M.radical()).is_semisimple()
            True
        """
        # A quiver representation is semisimple if and only if the zero map is
        # assigned to each edge.
        for x in self._semigroup._sorted_edges:
            if not self._maps[x].is_zero():
                return False
        return True

    def an_element(self):
        """
        Return an element of ``self``.

        OUTPUT: :class:`QuiverRepElement`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.an_element()
            Element of quiver representation
        """
        # Here we just use the an_element function from each space.
        elements = {v: self._spaces[v].an_element() for v in self._quiver}
        return self(elements)

    def support(self):
        """
        Return the support of ``self`` as a list.

        OUTPUT:

        - list, the vertices of the representation that have nonzero
          spaces associated to them

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 3:{2:['b'], 4:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 3)
            sage: M
            Representation with dimension vector (0, 1, 1, 1)
            sage: M.support()
            [2, 3, 4]
        """
        sup = []
        for v in self._quiver:
            if self._spaces[v].dimension() != 0:
                sup.append(v)

        return sup

    def gens(self, names='v') -> tuple:
        """
        Return a tuple of generators of ``self`` as a `k`-module.

        INPUT:

        - ``names`` -- an iterable variable of length equal to the number
          of generators, or a string (default: ``'v'``); gives the names of
          the generators either by giving a name to each generator or by
          giving a name to which an index will be appended

        OUTPUT:

        - tuple of :class:`QuiverRepElement` objects, the linear generators
          of the module (over the base ring)

        .. NOTE::

            The generators are ordered first by vertex and then by the order
            given by the ``gens()`` method of the space associated to
            that vertex.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.gens()
            (v_0, v_1, v_2)

        If a string is given then it is used as the name of each generator,
        with the index of the generator appended in order to differentiate
        them::

            sage: M.gens('generator')
            (generator_0, generator_1, generator_2)

        If a list or other iterable variable is given then each generator
        is named using the appropriate entry.  The length of the variable
        must equal the number of generators (the dimension of the module)::

            sage: M.gens(['w', 'x', 'y', 'z'])
            Traceback (most recent call last):
            ...
            TypeError: can only concatenate list (not "str") to list
            sage: M.gens(['x', 'y', 'z'])
            (x, y, z)

        Strings are iterable, so if the length of the string is equal to the
        number of generators then the characters of the string will be used
        as the names::

            sage: M.gens('xyz')
            (x, y, z)
        """
        # Use names as a list if and only if it is the correct length
        uselist = (len(names) == self.dimension())
        i = 0

        # Create bases for each space and construct QuiverRepElements from
        # them.  Other functions depend on the ordering of generators produced
        # by this code.  Make sure you know which they are before you change
        # anything
        basis = []
        for v in self._quiver:
            for m in self._spaces[v].gens():
                if uselist:
                    basis.append(self({v: m}, names[i]))
                else:
                    basis.append(self({v: m}, names + "_" + str(i)))
                i += 1

        return tuple(basis)

    def coordinates(self, vector):
        """
        Return the coordinates when ``vector`` is expressed in terms of
        the gens.

        INPUT:

        - ``vector`` -- :class:`QuiverRepElement`

        OUTPUT:

        - list, the coefficients when the vector is expressed as a linear
          combination of the generators of the module

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: M.coordinates(x - y + z)
            [1, -1, 1]
            sage: M.coordinates(M.an_element())
            [1, 1, 0]
            sage: M.an_element() == x + y
            True
        """
        # Just use the coordinates functions on each space
        coords = []
        for v in self._quiver:
            coords += self._spaces[v].coordinates(vector._elems[v])

        return coords

    def linear_combination_of_basis(self, coordinates):
        """
        Return the linear combination of the basis for ``self`` given
        by ``coordinates``.

        INPUT:

        - ``coordinates`` -- list; a list whose length is the dimension of
          ``self``.  The `i`-th element of this list defines the
          coefficient of the `i`-th basis vector in the linear
          combination.

        OUTPUT: :class:`QuiverRepElement`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: e = x - y + z
            sage: M.coordinates(e)
            [1, -1, 1]
            sage: M.linear_combination_of_basis([1, -1, 1]) == e
            True
        """
        # Make sure the input is valid
        gens = self.gens()
        if len(gens) != len(coordinates):
            raise ValueError("the coordinates do not match the dimension of the module")

        result = self()  # this must not be self.zero(), which is cached
        for ci, gi in zip(coordinates, gens):
            result += ci * gi
        return result

    ###########################################################################
    #                                                                         #
    # CONSTRUCTION FUNCTIONS                                                  #
    #    These functions create and return submodules and homomorphisms.      #
    #                                                                         #
    ###########################################################################

    def submodule(self, elements=[], spaces=None):
        """
        Return the submodule generated by the data.

        INPUT:

        - ``elements`` -- a collection of QuiverRepElements (default:
          empty list), each should be an element of ``self``

        - ``spaces`` -- dictionary (default: empty); should contain entries of
          the form ``{v: S}`` where `v` is a vertex of the quiver and `S` is a
          subspace of the vector space associated to `v`

        OUTPUT:

        - :class:`QuiverRep`, the smallest subrepresentation of ``self``
          containing the given elements and the given subspaces

        .. NOTE::

            This function returns only a :class:`QuiverRep` object ``sub``.
            The inclusion map of ``sub`` into ``M = self`` can be obtained
            by calling ``M.coerce_map_from(sub)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v])
            Representation with dimension vector (1, 1, 1)

        The smallest submodule containing the vector space at vertex 1
        also contains the entire vector space associated to vertex 3
        because there is an isomorphism associated to the edge
        ``(1, 3, 'a')``::

            sage: M.submodule(spaces={1: QQ^2})
            Representation with dimension vector (2, 0, 2)

        The smallest submodule containing the vector space at vertex 2
        also contains the image of the rank 1 homomorphism associated to
        the edge ``(2, 3, 'b')``::

            sage: M.submodule(spaces={2: QQ^3})
            Representation with dimension vector (0, 3, 1)

        As ``v`` is not already contained in this submodule, adding it as
        a generator yields a larger submodule::

            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3})
            Representation with dimension vector (1, 3, 1)

        Giving no generating data yields the zero submodule::

            sage: M.submodule().is_zero()
            True

        If the given data generates all of M then the result is M::

            sage: M.submodule(M.gens()) is M
            True
        """
        if spaces is None:
            spaces = {}

        # For each vertex generate a submodule from the given data
        dim = old_dim = 0
        for v in self._quiver:
            # Start with the zero submodule if no space is specified
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

            # Sum this with the submodule generated by the given elements.
            # Note that we are only taking the part of the element at the
            # vertex v.  We can always multiply an element of a quiver
            # representation by a constant path so we don't need to worry about
            # subspaces being embedded diagonally across multiple vertices.
            spaces[v] += self._spaces[v].submodule([m._elems[v] for m in elements])
            dim += spaces[v].dimension()

        # Now to enlarge the subspace to a submodule we sum a subspace at a
        # vertex with the images of the subspaces at adjacent vertices.  The
        # dimension of the subspace will strictly increase until we generate a
        # submodule.  At that point the dimension stabilizes and we can exit
        # the loop.
        while old_dim != dim:
            old_dim, dim = dim, 0

            # First sum the subspaces
            for e in self._semigroup._sorted_edges:
                spaces[e[1]] += self._maps[e](spaces[e[0]])

            # Then get the resulting dimensions
            for v in self._quiver:
                dim += spaces[v].dimension()

        # Return self if the entire module was generated, otherwise return a
        # submodule
        if dim == self.dimension():
            return self
        else:
            return self._submodule(spaces)

    def quotient(self, sub, check=True):
        """
        Return the quotient of ``self`` by the submodule ``sub``.

        INPUT:

        - ``sub`` -- :class:`QuiverRep`; this must be a submodule of ``self``,
          meaning the space associated to each vertex `v` of ``sub`` is a
          subspace of the space associated to `v` in ``self`` and the map
          associated to each edge `e` of ``sub`` is the restriction of
          the map associated to `e` in ``self``

        - ``check`` -- boolean; if ``True`` then ``sub`` is checked to verify
          that it is indeed a submodule of ``self`` and an error is raised
          if it is not

        OUTPUT: :class:`QuiverRep`, the quotient module ``self / sub``

        .. NOTE::

            This function returns only a QuiverRep object ``quot``.  The
            projection map from ``self`` to ``quot`` can be obtained by
            calling ``quot.coerce_map_from(self)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.I(GF(3), 3)
            sage: N = Q.S(GF(3), 3)
            sage: M.quotient(N)
            Representation with dimension vector (2, 1, 0)
            sage: M.quotient(M.radical())
            Representation with dimension vector (2, 0, 0)
            sage: M.quotient(M)
            Representation with dimension vector (0, 0, 0)
        """

        # First form the quotient space at each vertex
        spaces = {}
        for v in self._quiver:
            spaces[v] = self._spaces[v].quotient(sub._spaces[v], check)

        # Check the maps of sub if desired
        if check:
            for e in self._semigroup._sorted_edges:
                for x in sub._spaces[e[0]].gens():
                    if sub._maps[e](x) != self._maps[e](x):
                        raise ValueError("the quotient method was not passed a submodule")

        # Then pass the edge maps to the quotient
        maps = {}
        for e in self._semigroup._sorted_edges:
            # Sage can automatically coerce an element of a module to an
            # element of a quotient of that module but not the other way
            # around.  So in order to pass a map to the quotient we need to
            # construct the quotient map for the domain so that we can take
            # inverse images to lift elements.  As sage can coerce to a quotient
            # this is easy, we just send generators to themselves and set the
            # domain to be the quotient.

            # TODO: This 'if' shouldn't need to be here, but sage crashes when
            # coercing elements into a quotient that is zero.  Once Github issue
            # 12413 gets fixed only the else should need to execute.
            # NOTE: This is no longer necessary. Keep around for speed or
            # remove? -- darij, 16 Feb 2014
            if spaces[e[1]].dimension() == 0:
                maps[e] = spaces[e[0]].Hom(spaces[e[1]]).zero()
            else:
                factor = self._spaces[e[0]].hom(self._spaces[e[0]].gens(), spaces[e[0]])
                # Now we create a homomorphism by specifying the images of
                # generators.  Each generator is lifted to the original domain and
                # mapped over using the original map.  The codomain is set as the
                # quotient so sage will take care of pushing the result to the
                # quotient in the codomain.
                maps[e] = spaces[e[0]].hom([self._maps[e](factor.lift(x))
                                            for x in spaces[e[0]].gens()], spaces[e[1]])

        return self._semigroup.representation(self.base_ring(), spaces, maps)

    def socle(self):
        """
        The socle of ``self``.

        OUTPUT: :class:`QuiverRep`; the socle

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.socle()
            Representation with dimension vector (0, 0, 2)
        """
        # The socle of a representation is the intersection of the kernels of
        # all the edge maps.  The empty intersection is defined to be the
        # entire space so this is what we start with.
        spaces = self._spaces.copy()
        for e in self._semigroup._sorted_edges:
            spaces[e[0]] = spaces[e[0]].intersection(self._maps[e].kernel())

        return self._submodule(spaces)

    def radical(self):
        """
        Return the Jacobson radical of ``self``.

        OUTPUT: :class:`QuiverRep`; the Jacobson radical

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.radical()
            Representation with dimension vector (0, 2, 2)
        """
        # The Jacobson radical of a representation is the sum of the images of
        # all of the edge maps.  The empty sum is defined to be zero so this is
        # what we start with.
        spaces = {v: self._spaces[v].zero_submodule() for v in self._quiver}
        for e in self._semigroup._sorted_edges:
            spaces[e[1]] += self._maps[e].image()

        return self._submodule(spaces)

    def top(self):
        """
        Return the top of ``self``.

        OUTPUT: :class:`QuiverRep`; the quotient of ``self`` by its radical

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.top()
            Representation with dimension vector (1, 0, 0)
            sage: M.top() == M/M.radical()
            True
        """
        return self.quotient(self.radical())

    def zero_submodule(self):
        """
        Return the zero submodule of ``self``.

        OUTPUT: :class:`QuiverRep`; the zero submodule of ``self``

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.zero_submodule()
            Representation with dimension vector (0, 0, 0)
            sage: M.zero_submodule().is_zero()
            True
        """
        # When no data is specified this constructor automatically returns the
        # zero submodule
        return self._submodule()

    def linear_dual(self):
        """
        Compute the linear dual `Hom_k(M, k)` of the module
        `M =` ``self`` over the base ring `k`.

        OUTPUT: :class:`QuiverRep`; the dual representation

        .. NOTE::

            If `e` is an edge of the quiver `Q` then we let
            `(fe)(m) = f(me)`.  This gives `Hom_k(M, k)` a module
            structure over the opposite quiver ``Q.reverse()``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.linear_dual()
            Representation with dimension vector (1, 2, 2)
            sage: M.linear_dual().quiver() is Q.reverse().quiver()
            True
        """
        # This module is formed by taking the transpose of the edge maps.
        spaces = self._spaces.copy()
        maps = {(e[1], e[0], e[2]):
                self._spaces[e[1]].hom(self._maps[e].matrix().transpose(),
                                       self._spaces[e[0]])
                for e in self._semigroup._sorted_edges}

        # Reverse the bases if present
        if hasattr(self, '_bases'):
            bases = {}
            basis = []
            for v in self._bases:
                bases[v] = [p.reversal() for p in self._bases[v]]
                basis.extend(bases[v])

        if isinstance(self, QuiverRep_with_path_basis):
            result = self._semigroup.reverse().representation(self.base_ring(), basis, option='dual paths')
            result._maps = maps
            result._bases = bases
        elif isinstance(self, QuiverRep_with_dual_path_basis):
            result = self._semigroup.reverse().representation(self.base_ring(), basis, option='paths')
            result._maps = maps
            result._bases = bases
        else:
            result = self._semigroup.reverse().representation(self.base_ring(), spaces, maps)

        return result

    def algebraic_dual(self, basis=False):
        """
        Compute the algebraic dual `Hom_Q(M, kQ)` of the module
        `M` = ``self``.

        INPUT:

        - ``basis`` -- boolean; if ``False``, then only the module is
          returned.  If ``True``, then a tuple is returned.  The first
          element is the :class:`QuiverRep` and the second element is a
          dictionary which associates to each vertex a list.  The elements
          of this list are the homomorphisms which correspond to the basis
          elements of that vertex in the module.

        OUTPUT: :class:`QuiverRep` or tuple

        .. NOTE::

            Here `kQ` is the path algebra considered as a right module
            over itself.  If `e` is an edge of the quiver `Q` then we let
            `(fe)(m) = ef(m)`.  This gives `Hom_Q(M, kQ)` a module
            structure over the opposite quiver ``Q.reverse()``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}}).path_semigroup()
            sage: Q.free_module(GF(7)).algebraic_dual().dimension_vector()
            (7, 2, 1)
        """
        from sage.quivers.homspace import QuiverHomSpace
        return QuiverHomSpace(self, self._semigroup.free_module(self.base_ring())).left_module(basis)

    def Hom(self, codomain):
        """
        Return the hom space from ``self`` to ``codomain``.

        For more information see the :class:`QuiverHomSpace` documentation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            Dimension 2 QuiverHomSpace
        """
        from sage.quivers.homspace import QuiverHomSpace
        return QuiverHomSpace(self, codomain)

    def direct_sum(self, modules, return_maps=False):
        """
        Return the direct sum of ``self`` with the given modules
        ``modules``.

        The modules must be modules over the same quiver and base ring.

        INPUT:

        - ``modules`` -- :class:`QuiverRep` or list of :class:`QuiverRep`'s

        - ``return_maps`` -- boolean (default: ``False``); if ``False``, then
          the output is a single QuiverRep object which is the direct sum
          of ``self`` with the given module or modules.  If ``True``, then
          the output is a list ``[sum, iota, pi]``.  The first entry
          ``sum`` is the direct sum of ``self`` with the given module or
          modules.  Both ``iota`` and ``pi`` are lists of QuiverRepHoms
          with one entry for each summand; ``iota[i]`` is the inclusion
          map and ``pi[i]`` is the projection map of the `i`-th summand.
          The summands are ordered as given with ``self`` being the zeroth
          summand.

        OUTPUT: :class:`QuiverRep` or tuple

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}}).path_semigroup()
            sage: P1 = Q.P(QQ, 1)
            sage: P2 = Q.P(QQ, 2)
            sage: S = P1.direct_sum(P2)
            sage: P1.dimension_vector()
            (1, 1, 1, 2)
            sage: P2.dimension_vector()
            (0, 1, 0, 1)
            sage: S.dimension_vector()
            (1, 2, 1, 3)
            sage: S, iota, pi = P1.direct_sum(P2, return_maps=True)
            sage: iota[0].domain() is P1
            True
            sage: iota[1].domain() is P2
            True
            sage: pi[0].codomain() is P1
            True
            sage: pi[1].codomain() is P2
            True
            sage: iota[0].codomain() is S
            True
            sage: iota[1].codomain() is S
            True
            sage: pi[0].domain() is S
            True
            sage: pi[1].domain() is S
            True
            sage: iota[0].get_matrix(4)
            [1 0 0]
            [0 1 0]
            sage: pi[0].get_matrix(4)
            [1 0]
            [0 1]
            [0 0]
            sage: P1prime = S/iota[1].image()
            sage: f = P1prime.coerce_map_from(S)
            sage: g = f*iota[0]
            sage: g.is_isomorphism()
            True
        """
        # Convert a single module into a length 1 list and check validity
        if isinstance(modules, QuiverRep_generic):
            mods = [self, modules]
        else:
            mods = [self] + list(modules)
        for i in range(1, len(mods)):
            if not isinstance(mods[i], QuiverRep_generic):
                raise ValueError("modules must be a QuiverRep or a list of QuiverReps")
            if self._quiver is not mods[i]._quiver:
                raise ValueError("cannot direct sum modules of different quivers")
            if self.base_ring() is not mods[i].base_ring():
                raise ValueError("cannot direct sum modules with different base rings")

        # Get the dimensions of all spaces at each vertex
        dims = dict((v, [x._spaces[v].dimension() for x in mods]) for v in self._quiver)

        # Create spaces of the correct dimensions
        spaces = dict((v, self.base_ring()**sum(dims[v])) for v in self._quiver)

        # Take block sums of matrices to form the maps
        from sage.matrix.constructor import block_diagonal_matrix
        maps = {}
        for e in self._semigroup._sorted_edges:
            maps[e] = block_diagonal_matrix([x._maps[e].matrix() for x in mods], subdivide=False)

        # Create the QuiverRep, return if the maps aren't wanted
        result = self._semigroup.representation(self.base_ring(), spaces, maps)
        if not return_maps:
            return result

        # Create the inclusions and projections
        from sage.matrix.constructor import block_matrix, Matrix
        from sage.rings.integer import Integer
        iota = []
        pi = []
        for i in range(len(mods)):
            incl_maps = {}
            proj_maps = {}
            for v in self._quiver:
                # Create the maps using block matrices
                pre_dims = sum(dims[v][:i])
                post_dims = sum(dims[v][i + 1:])
                incl_maps[v] = block_matrix(1, 3, [Matrix(dims[v][i], pre_dims),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(dims[v][i], post_dims)])
                proj_maps[v] = block_matrix(3, 1, [Matrix(pre_dims, dims[v][i]),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(post_dims, dims[v][i])])
                # These matrices are over the integers, and get coerced
                # into the appropriate base ring at a later stage.
                # Might make trouble if the integer 1 does not coerce to
                # the 1 of the base ring; is that a real issue?
                # -- darij, 19 Feb 2014

            # Create and save the QuiverRepHom

            iota.append(mods[i].hom(incl_maps, result))
            pi.append(result.hom(proj_maps, mods[i]))

        # Return all the data
        return [result, iota, pi]

    def projective_cover(self, return_maps=False):
        """
        Return the projective cover of ``self``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: S1 = Q.S(GF(3), 1)
            sage: f1 = S1.projective_cover()
            sage: f1.is_surjective()
            True
            sage: f1._domain
            Representation with dimension vector (1, 2, 4)
            sage: Q.P(GF(3), 1)
            Representation with dimension vector (1, 2, 4)
        """
        # Deal with the zero module
        if self.dimension() == 0:
            return self.coerce_map_from(self)

        # Get the projection onto the top and generators
        top = self.top()
        proj_to_top = top.coerce_map_from(self)
        gens = top.gens()

        # Lift each of the generators
        lifts = [proj_to_top.lift(x) for x in gens]

        # Get projective summands of the cover
        Ps = [self._semigroup.P(self.base_ring(), x.support()[0]) for x in gens]

        # Factor the maps through self
        maps = [Ps[i].hom(lifts[i], self) for i in range(len(gens))]

        # Sum them and return
        return maps[0].direct_sum(maps[1:], return_maps, 'codomain')

    def transpose(self):
        r"""
        Return the transpose of ``self``.

        The transpose, `\mbox{Tr} M`, of a module `M` is defined as
        follows.  Let `p: P_1 \to P_2` be the second map in a minimal
        projective presentation `P_1 \to P_2 \to M \to 0` of `M`.
        If `p^t` is the algebraic dual of `p` then define
        `\mbox{Tr} M = \mbox{coker} p^t`.

        OUTPUT: :class:`QuiverRep`

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(GF(3), {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: M.transpose()
            Representation with dimension vector (1, 1)
        """

        # Get the second map in the minimal projective resolution of self.
        # Note we create this map as the projective cover of a kernel so we
        # need to alter the codomain to be the projective and not the kernel.
        p1 = self.projective_cover()
        k = p1.kernel()
        p = p1.domain().coerce_map_from(k) * k.projective_cover()

        # Return the cokernel
        return p.algebraic_dual().cokernel()

    def AR_translate(self):
        """
        Return the Auslander-Reiten translate of ``self``.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(GF(3), {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: tauM = M.AR_translate()
            sage: tauM
            Representation with dimension vector (1, 1)
            sage: tauM.get_map((1, 2, 'a')).matrix()
            [1]
            sage: tauM.get_map((1, 2, 'b')).matrix()
            [0]

        The module ``M`` above is its own AR translate.  This is not
        always true::

            sage: Q2 = DiGraph({3:{1:['b']}, 5:{3:['a']}}).path_semigroup()
            sage: Q2.S(QQ, 5).AR_translate()
            Representation with dimension vector (0, 1, 0)
        """

        return self.transpose().linear_dual()

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def right_edge_action(self, element, path):
        r"""
        Return the result of ``element*path``.

        INPUT:

        - ``element`` -- :class:`QuiverRepElement`, an element of ``self``

        - ``path`` -- :class:`QuiverPath` or list of tuples

        OUTPUT:

        - :class:`QuiverRepElement`, the result of ``element*path`` when
          ``path`` is considered an element of the path algebra of the quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]
            sage: M.right_edge_action(v, [(1, 1)]).support()
            [1]
            sage: M.right_edge_action(v, [(1, 1)]).support()
            [1]
            sage: M.right_edge_action(v, [(1, 2, 'a')]).support()
            [2]
            sage: M.right_edge_action(v, 'a') == M.right_edge_action(v, [(1, 2, 'a')])
            True
        """
        # Convert to a QuiverPath
        qpath = self._semigroup(path)

        # Invalid paths are zero in the quiver algebra
        result = self()  # this must not be self.zero(), which is cached

        # Start with the element at the initial vertex
        x = element._elems[qpath.initial_vertex()]

        # Act by each edge
        for e in qpath:
            x = self.get_map(e)(x)

        # Assign the result to the terminal vertex and return
        result._elems[qpath.terminal_vertex()] = x
        return result


class QuiverRep_with_path_basis(QuiverRep_generic):
    r"""
    The basis of the module must be closed under right multiplication by
    an edge; that is, appending any edge to the end of any path in the
    basis must result in either an invalid path or a valid path also
    contained in the basis of the module.

    INPUT:

    - ``k`` -- ring, the base ring of the representation

    - ``P`` -- the path semigroup of the quiver `Q` of the representation

    - ``basis`` -- list (default: empty); should be a list of paths (also
      lists) in the quiver `Q`.  Entries that do not represent valid paths
      are ignored and duplicate paths are deleted.  The closure of this
      list under right multiplication forms the basis of the resulting
      representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are path multiplication.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the codomain.  This class derives
    # from the QuiverRep class and the following private methods and variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths (also lists)
    #       which correspond to the basis elements of the space assigned to that
    #       vertex.
    #
    # If the right closure of the basis is also closed under left mult by an
    # edge then the object will have the following private variable:
    #
    # * _left_action_mats
    #       A dictionary of dictionaries.  If e is an edge or trivial path of
    #       the quiver then _left_action_mats[e] is a dictionary that
    #       associates a matrix to each vertex.  The left action of e on a
    #       QuiverRepElement is given by multiplying, on the left, the vector
    #       associated to vertex v in the QuiverRepElement by the matrix
    #       _left_action_mats[e][v]

    def __init__(self, k, P, basis):
        """
        Initialize ``self``. Type ``QuiverRep_with_path_basis?`` for more
        information.

        TESTS::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P1 = Q1.representation(QQ, [[(1, 1)]], option='paths')
            sage: P1.dimension()
            2
            sage: kQ = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 2, 'a')], [(1, 2, 'a')]], option='paths')
            sage: kQ.dimension()
            3
            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.dimension_vector()
            (0, 2, 2)
            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.dimension_vector()
            (0, 1, 2)
            sage: TestSuite(M).run()
        """
        self._quiver = Q = P.quiver()

        # Add the paths to the basis dictionary.  The terminal vertex is the
        # key
        self._bases = {v: [] for v in Q}
        for path in basis:
            # if path:
            self._bases[path.terminal_vertex()].append(path)

        # Create the matrices of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in P._sorted_edges:
            arrow = P([e], check=False)
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(k, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(len(self._bases[e[0]])):
                # Add an entry to the matrix corresponding to where the new path is found
                j = self._bases[e[1]].index(self._bases[e[0]][i] * arrow)
                maps[e][i, j] = k.one()

        # Create the spaces and then the representation
        spaces = dict((v, len(self._bases[v])) for v in Q)
        super().__init__(k, P, spaces, maps)

        # Try and create the matrices for the left edge action of edges.  If it
        # fails just return, there's no edge action and the construction is
        # done
        action_mats = {}
        for e in P._sorted_edges:
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(k, l, l)

                for j in range(l):
                    if e[1] == self._bases[v][j].initial_vertex():
                        try:
                            action_mats[e][v][self._bases[v].index(P([e], check=False) * self._bases[v][j]), j] = k.one()
                        except ValueError:
                            # There is no left action
                            return

        # If haven't returned yet then there is a left action.  Create the
        # matrices for acting by trivial paths
        for vert in self._quiver:
            e = (vert, vert)
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(k, l, l)

                # Paths not beginning at vert are sent to zero, paths beginning
                # at vert are fixed
                for i in range(l):
                    if self._bases[v][i].initial_vertex() == vert:
                        action_mats[e][v][i, i] = k.one()

        # Define the method and save the matrices
        self.left_edge_action = self._left_edge_action
        self._left_action_mats = action_mats

    def _left_edge_action(self, edge, element):
        r"""
        Return the result of ``edge*element``.

        INPUT:

        - ``element`` -- :class:`QuiverRepElement`; an element of ``self``

        - ``edge`` -- an edge of the quiver (a tuple) or a list of edges in
          the quiver.  Such a list can be empty (in which case no action
          is performed) and can contain trivial paths (tuples of the form
          `(v, v)` where `v` is a vertex of the quiver)

        OUTPUT:

        - :class:`QuiverRepElement`, the result of ``edge*element`` when
          ``edge`` is considered an element of the path algebra of the quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.representation(QQ, [[(1, 1)], [(2, 2)], [(3, 3)]], option='paths')
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]

        The sum of all trivial paths is the identity element under
        this action::

            sage: x = M.left_edge_action((1, 1), v) # indirect doctest
            sage: y = M.left_edge_action((2, 2), v) # indirect doctest
            sage: z = M.left_edge_action((3, 3), v) # indirect doctest
            sage: x + y + z == v
            True

        Note that the action only depends on the element of the path algebra
        that the input specifies::

            sage: a = M.left_edge_action([(1, 1), (1, 2, 'a'), (2, 2)], v) # indirect doctest
            sage: b = M.left_edge_action((1, 2, 'a'), v) # indirect doctest
            sage: a == b
            True

        These two edges multiply to zero in the path algebra::

            sage: M.left_edge_action([(1, 2, 'a'), (1, 2, 'b')], v).is_zero() # indirect doctest
            True
        """

        # Deal with lists by calling this function recursively
        if isinstance(edge, list):
            if not edge:
                return element
            else:
                return self.left_edge_action(edge[:-1], self.left_edge_action(edge[-1], element))

        # Now we are just acting by a single edge
        elems = {v: self._left_action_mats[edge][v] * element._elems[v]
                 for v in self._quiver}
        return self(elems)

    def is_left_module(self):
        """
        Test whether the basis is closed under left multiplication.

        EXAMPLES::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P2 = Q1.representation(QQ, [[(2, 2)]], option='paths')
            sage: P2.is_left_module()
            False

        The supplied basis is not closed under left multiplication, but it's
        not closed under right multiplication either.  When the closure under
        right multiplication is taken the result is also closed under left
        multiplication and therefore produces a left module structure::

            sage: kQ = Q1.representation(QQ, [[(1, 1)], [(2, 2)]], option='paths')
            sage: kQ.is_left_module()
            True

        Taking the right closure of a left closed set produces another
        left closed set::

            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.is_left_module()
            True

        Note that the second path is length 2, so even though the edge (1, 2, 'a')
        appears in the input the path [(1, 2, 'a')] is not in the right closure::

            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.is_left_module()
            False
        """
        return hasattr(self, 'left_edge_action')


class QuiverRep_with_dual_path_basis(QuiverRep_generic):
    r"""
    The basis of the module must be closed under left deletion of an edge; that
    is, deleting any edge from the beginning of any path in the basis must
    result in a path also contained in the basis of the module.

    INPUT:

    - ``k`` -- ring; the base ring of the representation

    - ``P`` -- the path semigroup of the quiver `Q` of the representation

    - ``basis`` -- list (default: empty); should be a list of paths (also
      lists) in the quiver `Q`.  Entries that do not represent valid paths
      are ignored and duplicate paths are deleted.  The closure of this
      list under left deletion forms the basis of the resulting
      representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are edge deletion.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the domain.  This class derives
    # from the QuiverRep class and the following methods and private variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths which
    #       correspond to the basis elements of the space assigned to that
    #       vertex.

    def __init__(self, k, P, basis):
        """
        Initialize ``self``. Type ``QuiverRep_with_dual_path_basis?`` for
        more information.

        TESTS::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: I2 = Q1.representation(QQ, [[(2, 2)]], option='dual paths')
            sage: I2.dimension()
            2
            sage: kQdual = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 2, 'a')], [(1, 2, 'a')]], option='dual paths')
            sage: kQdual.dimension()
            3
            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(1, 2, 'a'), (2, 3, 'd')], [(1, 3, 'b')]], option='dual paths')
            sage: M.dimension_vector()
            (2, 0, 0)
            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='dual paths')
            sage: N.dimension_vector()
            (2, 1, 0)
            sage: TestSuite(N).run()
        """
        self._quiver = Q = P.quiver()

        # Add the paths to the basis dictionary.  The initial vertex is the
        # key
        self._bases = {v: [] for v in Q}
        for path in basis:
            # if path:
            self._bases[path.initial_vertex()].append(path)

        # Create the matrices of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in P._sorted_edges:
            arrow = P([e], check=False)
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(k, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(len(self._bases[e[0]])):
                # Add an entry to the matrix corresponding to where the new path is found
                if self._bases[e[0]][i] % arrow in self._bases[e[1]]:
                    j = self._bases[e[1]].index(self._bases[e[0]][i] % arrow)
                    maps[e][i, j] = k.one()

        # Create the spaces and then the representation
        spaces = {v: len(self._bases[v]) for v in Q}
        super().__init__(k, P, spaces, maps)
