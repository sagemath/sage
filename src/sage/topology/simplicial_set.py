# sage.doctest: needs sage.graphs
r"""
Simplicial sets

AUTHORS:

- John H. Palmieri (2016-07)

This module implements simplicial sets.

A *simplicial set* `X` is a collection of sets `X_n` indexed by the
nonnegative integers; the set `X_n` is called the set of
`n`-simplices. These sets are connected by maps

.. MATH::

    d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
    s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

satisfying the *simplicial identities*:

.. MATH::

    d_i d_j &= d_{j-1} d_i \ \  \text{if } i<j \\
    d_i s_j &= s_{j-1} d_i \ \  \text{if } i<j \\
    d_j s_j &= 1 = d_{j+1} s_j \\
    d_i s_j &= s_{j} d_{i-1} \ \  \text{if } i>j+1 \\
    s_i s_j &= s_{j+1} s_{i} \ \  \text{if } i<j+1

See :wikipedia:`Simplicial_set`, Peter May's seminal book [May1967]_, or
Greg Friedman's "Illustrated introduction" :arxiv:`0809.4221` for more
information.

Several simplicial sets are predefined, and users can construct others
either by hand (using :class:`SimplicialSet_finite`) or from existing
ones using pushouts, pullbacks, etc.

EXAMPLES:

Some of the predefined simplicial sets::

    sage: simplicial_sets.Torus()
    Torus
    sage: simplicial_sets.RealProjectiveSpace(7)                                        # needs sage.groups
    RP^7
    sage: S5 = simplicial_sets.Sphere(5); S5
    S^5
    sage: S5.nondegenerate_simplices()
    [v_0, sigma_5]

One class of infinite simplicial sets is available: classifying spaces
of groups, or more generally, nerves of finite monoids::

    sage: Sigma4 = groups.permutation.Symmetric(4)                                      # needs sage.groups
    sage: Sigma4.nerve()                                                                # needs sage.groups
    Nerve of Symmetric group of order 4! as a permutation group

The same simplicial set (albeit with a different name) can also be
constructed as ::

    sage: simplicial_sets.ClassifyingSpace(Sigma4)                                      # needs sage.groups
    Classifying space of Symmetric group of order 4! as a permutation group

Type ``simplicial_sets.`` and hit the :kbd:`Tab` key to get a full list
of the predefined simplicial sets.

You can construct new simplicial sets from old by taking quotients,
subsimplicial sets, disjoint unions, wedges (if they are pointed),
smash products (if they are pointed and finite), products, pushouts,
pullbacks, cones, and suspensions, most of which also have maps
associated with them. Wedges, for example::

    sage: T = simplicial_sets.Torus()
    sage: S3 = simplicial_sets.Sphere(3)
    sage: T.is_pointed() and S3.is_pointed()
    True
    sage: T.wedge(S3)
    Wedge: (Torus v S^3)
    sage: T.disjoint_union(S3) == T.coproduct(S3)
    False

    sage: W = T.wedge(S3)
    sage: W.inclusion_map(0).domain()
    Torus
    sage: W.projection_map(1).codomain()
    Quotient: (Wedge: (Torus v S^3)/Simplicial set with 6 non-degenerate simplices)

If the `1`-sphere were not already available via
``simplicial_sets.Sphere(1)``, you could construct it as follows::

    sage: pt = simplicial_sets.Simplex(0)
    sage: edge = pt.cone()
    sage: S1 = edge.quotient(edge.n_skeleton(0))
    sage: S1
    Quotient: (Cone of 0-simplex/Simplicial set with 2 non-degenerate simplices)

At this point, ``S1`` is pointed: every quotient is automatically
given a base point, namely the image of the subcomplex. So its
suspension is the reduced suspension, and therefore is small::

    sage: S5 = S1.suspension(4)
    sage: S5
    Sigma^4(Quotient: (Cone of 0-simplex/Simplicial set with 2 non-degenerate simplices))
    sage: S5.f_vector()
    [1, 0, 0, 0, 0, 1]

If we forget about the base point in ``S1``, we would get the
unreduced suspension instead::

    sage: Z1 = S1.unset_base_point()
    sage: Z1.suspension(4).f_vector()
    [2, 2, 2, 2, 1, 1]

The cone on a pointed simplicial set is the reduced cone. The
`n`-simplex in Sage is not pointed, but the simplicial set ``Point``
is. ::

    sage: simplicial_sets.Simplex(0).cone().f_vector()
    [2, 1]
    sage: simplicial_sets.Point().cone().f_vector()
    [1]

For most simplicial sets (the ``Point`` is the main exception), each
time it is constructed, it gives a distinct copy, and two distinct
simplicial sets are never equal::

    sage: T = simplicial_sets.Torus()
    sage: T == simplicial_sets.Torus()
    False
    sage: T == T
    True
    sage: simplicial_sets.Torus() == simplicial_sets.Torus()
    False
    sage: simplicial_sets.Point() == simplicial_sets.Point()
    True

You can construct subsimplicial sets by specifying a list of simplices,
and then you can define the quotient simplicial set::

    sage: X = simplicial_sets.Simplex(2)
    sage: e,f,g = X.n_cells(1)
    sage: Y = X.subsimplicial_set([e,f,g])
    sage: Z = X.quotient(Y)

Or equivalently::

    sage: Y = X.n_skeleton(1)
    sage: Z = X.quotient(Y)
    sage: Z
    Quotient: (2-simplex/Simplicial set with 6 non-degenerate simplices)

Note that subsimplicial sets and quotients come equipped with
inclusion and quotient morphisms::

    sage: inc = Y.inclusion_map()
    sage: inc.domain() == Y and inc.codomain() == X
    True
    sage: quo = Z.quotient_map()
    sage: quo.domain()
    2-simplex
    sage: quo.codomain() == Z
    True

You can compute homology groups and the fundamental group of
any simplicial set::

    sage: S1 = simplicial_sets.Sphere(1)
    sage: eight = S1.wedge(S1)
    sage: eight.fundamental_group()                                                     # needs sage.groups
    Finitely presented group < e0, e1 | >

    sage: # needs sage.groups
    sage: Sigma3 = groups.permutation.Symmetric(3)
    sage: BSigma3 = Sigma3.nerve()
    sage: pi = BSigma3.fundamental_group(); pi
    Finitely presented group < e1, e2 | e2^2, e1^3, (e2*e1)^2 >
    sage: pi.order()
    6
    sage: pi.is_abelian()
    False

    sage: RP6 = simplicial_sets.RealProjectiveSpace(6)                                  # needs sage.groups
    sage: RP6.homology(reduced=False, base_ring=GF(2))                                  # needs sage.groups sage.modules
    {0: Vector space of dimension 1 over Finite Field of size 2,
     1: Vector space of dimension 1 over Finite Field of size 2,
     2: Vector space of dimension 1 over Finite Field of size 2,
     3: Vector space of dimension 1 over Finite Field of size 2,
     4: Vector space of dimension 1 over Finite Field of size 2,
     5: Vector space of dimension 1 over Finite Field of size 2,
     6: Vector space of dimension 1 over Finite Field of size 2}
    sage: RP6.homology(reduced=False, base_ring=QQ)                                     # needs sage.groups sage.modules
    {0: Vector space of dimension 1 over Rational Field,
     1: Vector space of dimension 0 over Rational Field,
     2: Vector space of dimension 0 over Rational Field,
     3: Vector space of dimension 0 over Rational Field,
     4: Vector space of dimension 0 over Rational Field,
     5: Vector space of dimension 0 over Rational Field,
     6: Vector space of dimension 0 over Rational Field}

When infinite simplicial sets are involved, most computations are done
by taking an `n`-skeleton for an appropriate `n`, either implicitly or
explicitly::

    sage: # needs sage.groups
    sage: G = groups.misc.MultiplicativeAbelian([3])
    sage: B3 = simplicial_sets.ClassifyingSpace(G)
    sage: B3.disjoint_union(B3).n_skeleton(3)
    Disjoint union: (Simplicial set with 15 non-degenerate simplices
                      u Simplicial set with 15 non-degenerate simplices)
    sage: S1 = simplicial_sets.Sphere(1)
    sage: B3.product(S1).homology(range(4))                                             # needs sage.modules
    {0: 0, 1: Z x C3, 2: C3, 3: C3}

Without the ``range`` argument, this would raise an error, since
``B3`` is infinite::

    sage: B3.product(S1).homology()                                                     # needs sage.groups sage.modules
    Traceback (most recent call last):
    ...
    NotImplementedError: this simplicial set may be infinite,
    so specify dimensions when computing homology

It should be easy to construct many simplicial sets from the
predefined ones using pushouts, pullbacks, etc., but they can also be
constructed "by hand": first define some simplices, then define a
simplicial set by specifying their faces::

    sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
    sage: v = AbstractSimplex(0, name='v')
    sage: w = AbstractSimplex(0, name='w')
    sage: e = AbstractSimplex(1, name='e')
    sage: f = AbstractSimplex(1, name='f')
    sage: X = SimplicialSet({e: (v,w), f: (w,w)})

Now `e` is an edge from `v` to `w` and `f` is an edge starting and
ending at `w`. Therefore the first homology group of `X` should be a
copy of the integers::

    sage: X.homology(1)                                                                 # needs sage.modules
    Z
"""
# ****************************************************************************
#  Copyright (C) 2016 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  https://www.gnu.org/licenses/
#
# ****************************************************************************

import copy

from sage.misc.cachefunc import cached_method
from sage.misc.fast_methods import WithEqualityById
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject

from .cell_complex import GenericCellComplex
from .delta_complex import DeltaComplex
from .simplicial_complex import SimplicialComplex

lazy_import('sage.categories.simplicial_sets', 'SimplicialSets')
lazy_import('sage.matrix.constructor', 'matrix')


########################################################################
# The classes for simplices.

class AbstractSimplex_class(SageObject):
    """
    A simplex of dimension ``dim``.

    INPUT:

    - ``dim`` -- integer; the dimension
    - ``degeneracies`` -- (optional) iterable, the indices of the
      degeneracy maps
    - ``underlying`` -- (optional) a non-degenerate simplex
    - ``name`` -- (optional) string
    - ``latex_name`` -- (optional) string

    Users should not call this directly, but instead use
    :func:`AbstractSimplex`. See that function for more documentation.
    """

    def __init__(self, dim, degeneracies=(), underlying=None, name=None,
                 latex_name=None):
        """
        A simplex of dimension ``dim``.

        INPUT:

        - ``dim`` -- integer; the dimension
        - ``degeneracies`` -- (optional) iterable, the indices of the degeneracy maps
        - ``underlying`` -- (optional) a non-degenerate simplex
        - ``name`` -- (optional) string
        - ``latex_name`` -- (optional) string

        Users should not call this directly, but instead use
        :func:`AbstractSimplex`. See that function for more
        documentation.

        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (1,2))
            s_3 s_1 Delta^3
            sage: AbstractSimplex(3, None)
            Delta^3

            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(0, (0,), underlying=v)
            sage: e
            s_0 v
            sage: e.nondegenerate() is v
            True

            sage: AbstractSimplex(3.2, None)
            Traceback (most recent call last):
            ...
            ValueError: the dimension must be an integer
            sage: AbstractSimplex(-3, None)
            Traceback (most recent call last):
            ...
            ValueError: the dimension must be nonnegative

            sage: AbstractSimplex(0, (1,))
            Traceback (most recent call last):
            ...
            ValueError: invalid list of degeneracy maps on 0-simplex

        Distinct non-degenerate simplices should never be equal, even
        if they have the same starting data::

            sage: from sage.topology.simplicial_set import AbstractSimplex_class
            sage: AbstractSimplex_class(3) == AbstractSimplex_class(3)
            False
            sage: AbstractSimplex(3) == AbstractSimplex(3)
            False

        Hashing::

            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: hash(v) == hash(w)
            False
            sage: x = v.apply_degeneracies(2,1,0)
            sage: hash(x) == hash(v.apply_degeneracies(2,1,0))
            True

        Equality::

            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: v == w
            False
            sage: v.apply_degeneracies(2,1,0) is v.apply_degeneracies(2,1,0)
            False
            sage: v.apply_degeneracies(2,1,0) == v.apply_degeneracies(2,1,0)
            True
            sage: v == None
            False
        """
        try:
            Integer(dim)
        except TypeError:
            raise ValueError('the dimension must be an integer')
        if dim < 0:
            raise ValueError('the dimension must be nonnegative')
        self._dim = dim
        if degeneracies:
            self._degens = standardize_degeneracies(*degeneracies)
            for (d, s) in enumerate(reversed(self._degens)):
                if d + dim < s:
                    raise ValueError('invalid list of degeneracy maps '
                                     'on {}-simplex'.format(dim))
            if underlying is None:
                self._underlying = NonDegenerateSimplex(dim)
            else:
                self._underlying = underlying
        else:
            self._degens = ()
            if underlying is None:
                self._underlying = self
            else:
                self._underlying = underlying
        if name is not None:
            self.rename(name)
        self._latex_name = latex_name

    def __hash__(self):
        """
        If nondegenerate: return the id of this simplex.

        Otherwise, combine the id of its underlying nondegenerate
        simplex with the tuple of indeterminacies.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: hash(v) == hash(w)
            False
            sage: x = v.apply_degeneracies(2,1,0)
            sage: id(x) == id(v.apply_degeneracies(2,1,0))
            False
            sage: hash(x) == hash(v.apply_degeneracies(2,1,0))
            True
        """
        if self.is_nondegenerate():
            return id(self)
        return hash(self.nondegenerate()) ^ hash(self._degens)

    def __eq__(self, other):
        """
        Two nondegenerate simplices are equal if they are identical.
        Two degenerate simplices are equal if their underlying
        nondegenerate simplices are identical and their tuples of
        degeneracies are equal.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: v == w
            False
            sage: v.apply_degeneracies(2,1,0) is v.apply_degeneracies(2,1,0)
            False
            sage: v.apply_degeneracies(2,1,0) == v.apply_degeneracies(2,1,0)
            True

        TESTS::

            sage: v == None
            False
        """
        if not isinstance(other, AbstractSimplex_class):
            return False
        return (self._degens == other._degens
                and self.nondegenerate() is other.nondegenerate())

    def __ne__(self, other):
        """
        This returns the negation of `__eq__`.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: v != w
            True
            sage: x = v.apply_degeneracies(1, 0)
            sage: y = v.apply_degeneracies(1, 0)
            sage: x != y
            False
        """
        return not self == other

    def __lt__(self, other):
        """
        We implement sorting in the hopes that sorted lists of simplices,
        for example as defining data for a simplicial set, will be
        well-defined invariants.

        Sort by dimension first. If dimensions are equal, if only one
        has a custom name (as set by specifying ``name=NAME`` upon
        creation or by calling ``object.rename('NAME')``, put it
        first. If both have custom names, sort by the names. As a last
        resort, sort by their id.

        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)

        At this point, comparison between v and w is random, based on
        their location in memory. ::

            sage: v < w and w < v
            False
            sage: v < w or w < v
            True
            sage: (v < w and w > v) or (w < v and v > w)
            True

        Now we add names to force an ordering::

            sage: w.rename('w')
            sage: v < w
            False
            sage: v > w
            True
            sage: v >= w
            True
            sage: v.rename('v')
            sage: v < w
            True
            sage: v <= w
            True
            sage: v > w
            False

        Test other sorting. Dimensions::

            sage: AbstractSimplex(0) < AbstractSimplex(3)
            True
            sage: AbstractSimplex(0, ((0,0,0))) <= AbstractSimplex(2)
            False

        Degenerate comes after non-degenerate, and if both are
        degenerate, sort on the degeneracies::

            sage: AbstractSimplex(0, ((0,0))) <= AbstractSimplex(2)
            False
            sage: AbstractSimplex(1, ((0,))) < AbstractSimplex(1, ((1,)))
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False
            sage: w.rename('a')
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            True

        Testing `<=`, `>`, `>=`::

            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: v <= v
            True
            sage: w <= v
            False
            sage: v.apply_degeneracies(1,0) <= w.apply_degeneracies(1,0)
            True

            sage: v > v
            False
            sage: w > v
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False

            sage: v >= v
            True
            sage: w >= v
            True
            sage: v.apply_degeneracies(1,0) >= w.apply_degeneracies(1,0)
            False
        """
        if self.dimension() < other.dimension():
            return True
        if self.dimension() > other.dimension():
            return False
        if self.degeneracies() and not other.degeneracies():
            return False
        if other.degeneracies() and not self.degeneracies():
            return True
        if self.degeneracies() and other.degeneracies() and self.degeneracies() != other.degeneracies():
            return self.degeneracies() < other.degeneracies()
        if self.nondegenerate().get_custom_name() is not None:
            if other.nondegenerate().get_custom_name() is not None:
                return self.nondegenerate().get_custom_name() < other.nondegenerate().get_custom_name()
            return True

        if other.nondegenerate().get_custom_name() is not None:
            return False
        return id(self) < id(other)

    def __gt__(self, other):
        """
        See :meth:`__lt__` for more doctests.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: e = AbstractSimplex(1, (1,0), name='e')
            sage: f = AbstractSimplex(1, (2,1), name='f')
            sage: e > f
            False
        """
        return not (self < other or self == other)

    def __le__(self, other):
        """
        See :meth:`__lt__` for more doctests.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: e = AbstractSimplex(1, (1,0), name='e')
            sage: f = AbstractSimplex(1, (2,1), name='f')
            sage: e <= f
            True
        """
        return self < other or self == other

    def __ge__(self, other):
        """
        See :meth:`__lt__` for more doctests.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: e = AbstractSimplex(1, (1,0), name='e')
            sage: f = AbstractSimplex(1, (2,1), name='f')
            sage: e >= f
            False
        """
        return not self < other

    def nondegenerate(self):
        """
        The non-degenerate simplex underlying this one.

        Therefore return itself if this simplex is non-degenerate.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0, name='v')
            sage: sigma = v.apply_degeneracies(1, 0)
            sage: sigma.nondegenerate()
            v
            sage: tau = AbstractSimplex(1, (3,2,1))
            sage: x = tau.nondegenerate(); x
            Delta^1
            sage: x == tau.nondegenerate()
            True

            sage: AbstractSimplex(1, None)
            Delta^1
            sage: AbstractSimplex(1, None) == x
            False
            sage: AbstractSimplex(1, None) == tau.nondegenerate()
            False
        """
        return self._underlying

    def degeneracies(self):
        """
        Return the list of indices for the degeneracy maps for this
        simplex.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(4, (0,0,0)).degeneracies()
            [2, 1, 0]
            sage: AbstractSimplex(4, None).degeneracies()
            []
        """
        return list(self._degens)

    def is_degenerate(self):
        """
        Return ``True`` if this simplex is degenerate.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).is_degenerate()
            True
            sage: AbstractSimplex(3, None).is_degenerate()
            False
        """
        return bool(self.degeneracies())

    def is_nondegenerate(self):
        """
        Return ``True`` if this simplex is non-degenerate.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).is_nondegenerate()
            False
            sage: AbstractSimplex(3, None).is_nondegenerate()
            True
            sage: AbstractSimplex(5).is_nondegenerate()
            True
        """
        return not self.is_degenerate()

    def dimension(self):
        """
        The dimension of this simplex.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).dimension()
            5
            sage: AbstractSimplex(3, None).dimension()
            3
            sage: AbstractSimplex(7).dimension()
            7
        """
        return self._dim + len(self.degeneracies())

    def apply_degeneracies(self, *args):
        """
        Apply the degeneracies given by the arguments ``args`` to this simplex.

        INPUT:

        - ``args`` -- integer

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: e = v.apply_degeneracies(0)
            sage: e.nondegenerate() == v
            True
            sage: f = e.apply_degeneracies(0)
            sage: f
            s_1 s_0 Delta^0
            sage: f.degeneracies()
            [1, 0]
            sage: f.nondegenerate() == v
            True
            sage: v.apply_degeneracies(1, 0)
            s_1 s_0 Delta^0

        TESTS::

            sage: e.apply_degeneracies() == e
            True

        Do not pass an explicit list or tuple as the argument: call
        this with the syntax ``x.apply_degeneracies(1,0)``, not
        ``x.apply_degeneracies([1,0])``::

            sage: e.apply_degeneracies([1,0])
            Traceback (most recent call last):
            ...
            TypeError: degeneracies are indexed by nonnegative integers; do not use an explicit list or tuple
        """
        if not args:
            return self
        underlying = self.nondegenerate()
        return AbstractSimplex(underlying.dimension(),
                               degeneracies=list(args) + self.degeneracies(),
                               underlying=underlying)

    def __copy__(self):
        """
        Return a copy of this simplex.

        Forget the "underlying" non-degenerate simplex. If this
        simplex has a name, then its copy's name is obtained by adding
        a prime ``'`` at the end.

        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: copy(v) == v
            False
            sage: copy(v).nondegenerate() == v
            False
            sage: x = v.apply_degeneracies(1, 0)
            sage: y = copy(v).apply_degeneracies(1, 0)
            sage: z = copy(x)
            sage: x == y or x == z or y == z
            False
            sage: x.nondegenerate() == copy(v)
            False
            sage: y.nondegenerate() == v
            False

            sage: v.rename('v')
            sage: copy(v)
            v'
            sage: copy(copy(v))
            v''
        """
        # Don't preserve the underlying simplex when copying, just the
        # dimension, the degeneracies, and the name (with a prime
        # added).
        sigma = AbstractSimplex(self._dim, degeneracies=self.degeneracies())
        if self.get_custom_name() is not None:
            sigma.rename(self.get_custom_name() + "'")
        return sigma

    def __deepcopy__(self, memo):
        """
        Return a "deep" copy of this simplex.

        INPUT:

        - ``memo`` -- "memo" dictionary required by the ``copy.deepcopy`` method

        This returns the same object as the :meth:`__copy__` method
        and also updates ``memo``.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: import copy
            sage: v = AbstractSimplex(0)
            sage: copy.deepcopy(v) == v
            False

        TESTS:

        The purpose for this method is to be able to make distinct
        copies of simplicial sets::

            sage: # needs sage.groups
            sage: from sage.topology.simplicial_set import SimplicialSet
            sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
            sage: dict(copy.copy(RP3._data)) == dict(RP3._data)
            True
            sage: dict(copy.deepcopy(RP3._data)) == dict(RP3._data)
            False
            sage: SimplicialSet(RP3) == RP3
            False
            sage: copy.copy(RP3) == RP3
            False
        """
        underlying = self.nondegenerate()
        degens = self.degeneracies()
        try:
            return memo[underlying].apply_degeneracies(*degens)
        except KeyError:
            sigma = AbstractSimplex(underlying._dim)
            if underlying.get_custom_name() is not None:
                sigma.rename(underlying.get_custom_name() + "'")
            memo[underlying] = sigma
            return sigma.apply_degeneracies(*degens)

    def _repr_(self):
        """
        Print representation.

        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, None)
            Delta^3
            sage: AbstractSimplex(3, (0,))
            s_0 Delta^3
            sage: AbstractSimplex(3, (0, 0))
            s_1 s_0 Delta^3

        Test renaming::

            sage: v = AbstractSimplex(0)
            sage: v
            Delta^0
            sage: v.rename('v')
            sage: v
            v
            sage: v.apply_degeneracies(1, 0)
            s_1 s_0 v
        """
        if self.degeneracies():
            degens = ' '.join(f's_{i}' for i in self.degeneracies())
            return degens + ' {}'.format(self.nondegenerate())
        return 'Delta^{}'.format(self._dim)

    def _latex_(self):
        r"""
        LaTeX representation.

        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: latex(AbstractSimplex(18, None))
            \Delta^{18}
            sage: latex(AbstractSimplex(3, (0, 0,)))
            s_{1} s_{0} \Delta^{3}
            sage: latex(AbstractSimplex(3, (0, 0,), name='x'))
            x
            sage: latex(AbstractSimplex(3, name='x').apply_degeneracies(0, 0))
            s_{1} s_{0} x
            sage: latex(AbstractSimplex(3, (0, 0,), name='x', latex_name='y'))
            y
            sage: latex(AbstractSimplex(3, name='x', latex_name='y').apply_degeneracies(0, 0))
            s_{1} s_{0} y
        """
        if self._latex_name is not None:
            return self._latex_name
        if self.get_custom_name() is not None:
            return self.get_custom_name()
        if self.nondegenerate()._latex_name is not None:
            simplex = self.nondegenerate()._latex_name
        elif self.nondegenerate().get_custom_name() is not None:
            simplex = self.nondegenerate().get_custom_name()
        else:
            simplex = "\\Delta^{{{}}}".format(self._dim)
        if self.degeneracies():
            degens = ' '.join(f's_{{{i}}}' for i in self.degeneracies())
            return degens + ' ' + simplex
        return simplex


# If we inherit from AbstractSimplex_class first in the following,
# then we have to override __eq__ and __hash__. If we inherit from
# WithEqualityById first, then we have to override __lt__, __gt__,
# __ge__, __le__. Inheriting from AbstractSimplex_class first seems to
# be slightly faster.
class NonDegenerateSimplex(AbstractSimplex_class, WithEqualityById):
    def __init__(self, dim, name=None, latex_name=None):
        """
        A nondegenerate simplex.

        INPUT:

        - ``dim`` -- nonnegative integer; the dimension

        - ``name`` -- (optional) string; a name for this simplex

        - ``latex_name`` -- (optional) string; a name for this simplex to
          use in the LaTeX representation

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0, name='v')
            sage: v
            v
            sage: type(v)
            <class 'sage.topology.simplicial_set.NonDegenerateSimplex'>

        Distinct non-degenerate simplices should never be equal, even
        if they have the same starting data. ::

            sage: v == AbstractSimplex(0, name='v')
            False
            sage: AbstractSimplex(3) == AbstractSimplex(3)
            False

            sage: from sage.topology.simplicial_set import NonDegenerateSimplex
            sage: x = NonDegenerateSimplex(0, name='x')
            sage: x == NonDegenerateSimplex(0, name='x')
            False
        """
        AbstractSimplex_class.__init__(self, dim, name=name, latex_name=latex_name)

    __eq__ = WithEqualityById.__eq__
    __hash__ = WithEqualityById.__hash__


# The following function returns an instance of either
# AbstractSimplex_class or NonDegenerateSimplex.

def AbstractSimplex(dim, degeneracies=(), underlying=None,
                    name=None, latex_name=None):
    r"""
    An abstract simplex, a building block of a simplicial set.

    In a simplicial set, a simplex either is non-degenerate or is
    obtained by applying degeneracy maps to a non-degenerate simplex.

    INPUT:

    - ``dim`` -- nonnegative integer; the dimension of the
      underlying non-degenerate simplex

    - ``degeneracies`` -- (default: ``None``) list or tuple of
      nonnegative integers, the degeneracies to be applied

    - ``underlying`` -- (optional) a non-degenerate simplex to which
      the degeneracies are being applied

    - ``name`` -- (optional) string; a name for this simplex

    - ``latex_name`` -- (optional) string; a name for this simplex to
      use in the LaTeX representation

    So to define a simplex formed by applying the degeneracy maps `s_2
    s_1` to a 1-simplex, call ``AbstractSimplex(1, (2, 1))``.

    Specify ``underlying`` if you need to keep explicit track of the
    underlying non-degenerate simplex, for example when computing
    faces of another simplex. This is mainly for use by the method
    :meth:`AbstractSimplex_class.apply_degeneracies`.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import AbstractSimplex
        sage: AbstractSimplex(3, (3, 1))
        s_3 s_1 Delta^3
        sage: AbstractSimplex(3, None)
        Delta^3
        sage: AbstractSimplex(3)
        Delta^3

    Simplices may be named (or renamed), affecting how they are printed::

        sage: AbstractSimplex(0)
        Delta^0
        sage: v = AbstractSimplex(0, name='v')
        sage: v
        v
        sage: v.rename('w_0')
        sage: v
        w_0
        sage: latex(v)
        w_0
        sage: latex(AbstractSimplex(0, latex_name='\\sigma'))
        \sigma

    The simplicial identities are used to put the degeneracies in
    standard decreasing form::

        sage: x = AbstractSimplex(0, (0, 0, 0))
        sage: x
        s_2 s_1 s_0 Delta^0
        sage: x.degeneracies()
        [2, 1, 0]

    Use of the ``underlying`` argument::

        sage: v = AbstractSimplex(0, name='v')
        sage: e = AbstractSimplex(0, (0,), underlying=v)
        sage: e
        s_0 v
        sage: e.nondegenerate() is v
        True

        sage: e.dimension()
        1
        sage: e.is_degenerate()
        True

    Distinct non-degenerate simplices are never equal::

        sage: AbstractSimplex(0, None) == AbstractSimplex(0, None)
        False
        sage: AbstractSimplex(0, (2,1,0)) == AbstractSimplex(0, (2,1,0))
        False

        sage: e = AbstractSimplex(0, ((0,)))
        sage: f = AbstractSimplex(0, ((0,)))
        sage: e == f
        False
        sage: e.nondegenerate() == f.nondegenerate()
        False

    This means that if, when defining a simplicial set, you specify
    the faces of a 2-simplex as::

        (e, e, e)

    then the faces are the same degenerate vertex, but if you specify
    the faces as::

        (AbstractSimplex(0, ((0,))), AbstractSimplex(0, ((0,))), AbstractSimplex(0, ((0,))))

    then the faces are three different degenerate vertices.

    View a command like ``AbstractSimplex(0, (2,1,0))`` as first
    constructing ``AbstractSimplex(0)`` and then applying degeneracies
    to it, and you always get distinct simplices from different calls
    to ``AbstractSimplex(0)``. On the other hand, if you apply
    degeneracies to the same non-degenerate simplex, the resulting
    simplices are equal::

        sage: v = AbstractSimplex(0)
        sage: v.apply_degeneracies(1, 0) == v.apply_degeneracies(1, 0)
        True
        sage: AbstractSimplex(1, (0,), underlying=v) == AbstractSimplex(1, (0,), underlying=v)
        True
    """
    if degeneracies:
        if underlying is None:
            underlying = NonDegenerateSimplex(dim)
        return AbstractSimplex_class(dim, degeneracies=degeneracies,
                                     underlying=underlying,
                                     name=name,
                                     latex_name=latex_name)
    else:
        return NonDegenerateSimplex(dim, name=name,
                                    latex_name=latex_name)


########################################################################
# The main classes for simplicial sets.

class SimplicialSet_arbitrary(Parent):
    r"""
    A simplicial set.

    A simplicial set `X` is a collection of sets `X_n`, the
    *n-simplices*, indexed by the nonnegative integers, together with
    maps

    .. MATH::

        d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. MATH::

        d_i d_j &= d_{j-1} d_i \ \  \text{if } i<j \\
        d_i s_j &= s_{j-1} d_i \ \  \text{if } i<j \\
        d_j s_j &= 1 = d_{j+1} s_j \\
        d_i s_j &= s_{j} d_{i-1} \ \  \text{if } i>j+1 \\
        s_i s_j &= s_{j+1} s_{i} \ \  \text{if } i<j+1

    This class is not fully implemented and is not intended to be
    called directly by users. It is intended instead to be used by
    other classes which inherit from this one. See
    :class:`SimplicialSet_finite` and :class:`Nerve` for two
    examples. In particular, any such class must implement a method
    ``n_skeleton`` -- without this, most computations will be
    impossible. It must also implement an ``__init__`` method which
    should also set the category, so that methods defined at the
    category level, like ``is_pointed`` and ``is_finite``, work
    correctly.

    Note that the method :meth:`subsimplicial_set` calls
    :meth:`n_skeleton`, so to avoid circularity, the
    :meth:`n_skeleton` method should call
    :class:`.simplicial_set_constructions.SubSimplicialSet` directly,
    not :meth:`subsimplicial_set`.
    """

    # This is cached because it is used frequently in morphism
    # construction when verifying that the morphism commutes with the
    # face maps.
    @cached_method
    def faces(self, simplex):
        """
        Return the list of faces of ``simplex`` in this simplicial set.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial set, either
          degenerate or not

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: S2.faces(sigma)
            (s_0 v_0, s_0 v_0, s_0 v_0)
            sage: S2.faces(sigma.apply_degeneracies(0))
            [sigma_2, sigma_2, s_1 s_0 v_0, s_1 s_0 v_0]

            sage: # needs sage.groups
            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)
            sage: f2 = BC3.n_cells(1)[1]; f2
            f^2
            sage: BC3.faces(f2)
            (1, 1)

        TESTS::

            sage: v_0 = S2.n_cells(0)[0]
            sage: S2.faces(v_0) is None
            True

            sage: from sage.topology.simplicial_set import AbstractSimplex
            sage: w = AbstractSimplex(0)
            sage: S2.faces(w)
            Traceback (most recent call last):
            ...
            ValueError: this simplex is not in this simplicial set
        """
        dim = simplex.dimension()
        if simplex not in self:
            raise ValueError('this simplex is not in this simplicial set')
        if simplex.is_nondegenerate():
            if self.is_finite():
                return self.face_data()[simplex]
            else:
                return self.n_skeleton(dim).face_data()[simplex]
        underlying = simplex.nondegenerate()
        faces = []
        for J, t in [face_degeneracies(m, simplex.degeneracies())
                     for m in range(dim+1)]:
            if t is None:
                faces.append(underlying.apply_degeneracies(*J))
            else:
                faces.append(self.face(underlying, t).apply_degeneracies(*J))
        return faces

    def face(self, simplex, i):
        """
        Return the `i`-th face of ``simplex`` in this simplicial set.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial set
        - ``i`` -- integer

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: v_0 = S2.n_cells(0)[0]
            sage: S2.face(sigma, 0)
            s_0 v_0
            sage: S2.face(sigma, 0) == v_0.apply_degeneracies(0)
            True
            sage: S2.face(S2.face(sigma, 0), 0) == v_0
            True
        """
        if i < 0 or i > simplex.dimension():
            raise ValueError('cannot compute face {} of {}-dimensional '
                             'simplex'.format(i, simplex.dimension()))
        faces = self.faces(simplex)
        if faces is not None:
            return self.faces(simplex)[i]
        return None

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is a simplex which is contained in this complex.

        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S0.n_cells(0)[0]
            sage: v0 in S0
            True
            sage: v0 in S1
            False

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: K = SimplicialSet({e: (v, v)})  # the circle
            sage: v in K
            True
            sage: v0 in K
            False
            sage: S1.n_cells(1)[0] in K
            False

        TESTS:

        Make sure we answer gracefully for unexpected input::

            sage: 248 in K
            False
        """
        try:
            underlying = x.nondegenerate()
            return underlying in self.n_cells(underlying.dimension())
        except AttributeError:
            return False

    def alexander_whitney(self, simplex, dim_left):
        r"""
        Return the 'subdivision' of ``simplex`` in this simplicial set
        into a pair of simplices.

        The left factor should have dimension ``dim_left``, so the
        right factor should have dimension ``dim - dim_left``, if
        ``dim`` is the dimension of the starting simplex. The results
        are obtained by applying iterated face maps to
        ``simplex``. Writing `d` for ``dim`` and `j` for ``dim_left``:
        apply `d_{j+1} d_{j+2} ... d_{d}` to get the left factor,
        `d_0 ... d_0` to get the right factor.

        INPUT:

        - ``dim_left`` -- integer; the dimension of the left-hand factor

        OUTPUT: list containing the triple ``(c, left, right)``,
        where ``left`` and ``right`` are the two simplices described
        above. If either ``left`` or ``right`` is degenerate, ``c`` is
        0; otherwise, ``c`` is 1. This is so that, when used to
        compute cup products, it is easy to ignore terms which have
        degenerate factors.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: S2.alexander_whitney(sigma, 0)
            [(1, v_0, sigma_2)]
            sage: S2.alexander_whitney(sigma, 1)
            [(0, s_0 v_0, s_0 v_0)]
        """
        dim = simplex.dimension()
        if dim_left < 0 or dim_left > dim:
            raise ValueError('alexander_whitney is only valid if dim_left '
                             'is between 0 and the dimension of the simplex')
        left = simplex
        for i in range(dim, dim_left, -1):
            left = self.face(left, i)
        right = simplex
        for i in range(dim_left):
            right = self.face(right, 0)
        if left.is_degenerate() or right.is_degenerate():
            c = ZZ.zero()
        else:
            c = ZZ.one()
        return [(c, left, right)]

    def nondegenerate_simplices(self, max_dim=None):
        """
        Return the sorted list of non-degenerate simplices in this simplicial set.

        INPUT:

        - ``max_dim`` -- (default: ``None``) if specified,
          return the non-degenerate simplices of this dimension or
          smaller. This argument is required if this simplicial set is
          infinite.

        The sorting is in increasing order of dimension, and within
        each dimension, by the name (if present) of each simplex.

        .. NOTE::

            The sorting is done when the simplicial set is
            constructed, so changing the name of a simplex after
            construction will not affect the ordering.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.nondegenerate_simplices()
            [Delta^0, Delta^0]

        Name the vertices and reconstruct the simplicial set: they
        should be ordered alphabetically::

            sage: v.rename('v')
            sage: w.rename('w')
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.nondegenerate_simplices()
            [v, w]

        Rename but do not reconstruct the set; the ordering does not
        take the new names into account::

            sage: v.rename('z')
            sage: S0.nondegenerate_simplices() # old ordering is used
            [z, w]

            sage: X0 = SimplicialSet({v: None, w: None})
            sage: X0.nondegenerate_simplices() # new ordering is used
            [w, z]

        Test an infinite example::

            sage: # needs sage.groups
            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)
            sage: BC3.nondegenerate_simplices(2)
            [1, f, f^2, f * f, f * f^2, f^2 * f, f^2 * f^2]
            sage: BC3.nondegenerate_simplices()
            Traceback (most recent call last):
            ...
            NotImplementedError: this simplicial set may be infinite, so specify max_dim
        """
        if self.is_finite():
            if max_dim is None:
                return list(self._simplices)
            return [sigma for sigma in self._simplices if sigma.dimension() <= max_dim]
        if max_dim is None:
            raise NotImplementedError('this simplicial set may be '
                                      'infinite, so specify max_dim')
        return list(self.n_skeleton(max_dim)._simplices)

    def cells(self, subcomplex=None, max_dim=None):
        """
        Return a dictionary of all non-degenerate simplices.

        INPUT:

        - ``subcomplex`` -- (optional) a subsimplicial set of this
          simplicial set. If ``subcomplex`` is specified, then return the
          simplices in the quotient by the subcomplex.

        - ``max_dim`` -- (default: ``None``) if specified,
          return the non-degenerate simplices of this dimension or
          smaller. This argument is required if this simplicial set is
          infinite.

        Each key is a dimension, and the corresponding value is the
        list of simplices in that dimension.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.cells()
            {0: [Delta^0, Delta^0]}

            sage: v.rename('v')
            sage: w.rename('w')
            sage: S0.cells()
            {0: [v, w]}

            sage: e = AbstractSimplex(1, name='e')
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: S1.cells()
            {0: [v], 1: [e]}

            sage: S0.cells(S0.subsimplicial_set([v, w]))
            {0: [*]}

            sage: X = SimplicialSet({e: (v,w)})
            sage: X.cells(X.subsimplicial_set([v, w]))
            {0: [*], 1: [e]}

        Test an infinite example::

            sage: # needs sage.groups
            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)
            sage: BC3.cells(max_dim=2)
            {0: [1], 1: [f, f^2], 2: [f * f, f * f^2, f^2 * f, f^2 * f^2]}
            sage: BC3.cells()
            Traceback (most recent call last):
            ...
            NotImplementedError: this simplicial set may be infinite, so specify max_dim
        """
        if subcomplex is None:
            if self.is_finite():
                simplices = {}
                for sigma in self.nondegenerate_simplices():
                    if sigma.dimension() in simplices:
                        simplices[sigma.dimension()].append(sigma)
                    else:
                        simplices[sigma.dimension()] = [sigma]
                if max_dim is not None:
                    return {d: sorted(simplices[d]) for d in simplices
                            if d <= max_dim}
                return {d: sorted(simplices[d]) for d in simplices}
            # Infinite case:
            if max_dim is None:
                raise NotImplementedError('this simplicial set may be '
                                          'infinite, so specify max_dim')
            return self.n_skeleton(max_dim).cells()
        # subcomplex is not None:
        return self.quotient(subcomplex).cells(max_dim=max_dim)

    def n_cells(self, n, subcomplex=None):
        """
        Return the list of cells of dimension ``n`` of this cell complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional faces in the quotient by this
        subcomplex.

        INPUT:

        - ``n`` -- the dimension

        - ``subcomplex`` -- (default: ``None``) a subcomplex
          of this cell complex. Return the cells which are in the
          quotient by this subcomplex.

        EXAMPLES::

            sage: simplicial_sets.Sphere(3).n_cells(3)
            [sigma_3]
            sage: simplicial_sets.Sphere(3).n_cells(2)
            []
            sage: C2 = groups.misc.MultiplicativeAbelian([2])                           # needs sage.groups
            sage: BC2 = C2.nerve()                                                      # needs sage.groups
            sage: BC2.n_cells(3)                                                        # needs sage.groups
            [f * f * f]
        """
        cells = self.cells(subcomplex=subcomplex, max_dim=n)
        try:
            return list(cells[n])
        except KeyError:
            # Don't barf if someone asks for n_cells in a dimension
            # where there are none.
            return []

    def _an_element_(self):
        """
        Return an element: a vertex of this simplicial set.

        Return ``None`` if the simplicial set is empty.

        EXAMPLES::

            sage: S4 = simplicial_sets.Sphere(4)
            sage: S4._an_element_()
            v_0
            sage: S4._an_element_() in S4
            True
            sage: from sage.topology.simplicial_set_examples import Empty
            sage: Empty()._an_element_() is None
            True
        """
        vertices = self.n_cells(0)
        if vertices:
            return vertices[0]
        return None

    def all_n_simplices(self, n):
        """
        Return a list of all simplices, non-degenerate and degenerate, in dimension ``n``.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2, name='tau')
            sage: Y = SimplicialSet({tau: (degen, degen, degen), w: None})

        ``Y`` is the disjoint union of a 2-sphere, with vertex ``v``
        and non-degenerate 2-simplex ``tau``, and a point ``w``. ::

            sage: Y.all_n_simplices(0)
            [v, w]
            sage: Y.all_n_simplices(1)
            [s_0 v, s_0 w]
            sage: Y.all_n_simplices(2)
            [tau, s_1 s_0 v, s_1 s_0 w]

        An example involving an infinite simplicial set::

            sage: C3 = groups.misc.MultiplicativeAbelian([3])                           # needs sage.groups
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                            # needs sage.groups
            sage: BC3.all_n_simplices(2)                                                # needs sage.groups
            [f * f,
             f * f^2,
             f^2 * f,
             f^2 * f^2, s_0 f, s_0 f^2, s_1 f, s_1 f^2, s_1 s_0 1]
        """
        non_degen = list(self.nondegenerate_simplices(max_dim=n))
        ans = {_ for _ in non_degen if _.dimension() == n}
        for sigma in non_degen:
            d = sigma.dimension()
            ans.update([sigma.apply_degeneracies(*_)
                        for _ in all_degeneracies(d, n-d)])
        return sorted(ans)

    def _map_from_empty_set(self):
        """
        Return the unique map from the empty set to this simplicial set.

        This is used to in the method :meth:`disjoint_union` to
        construct disjoint unions as pushouts.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: T._map_from_empty_set()
            Simplicial set morphism:
              From: Empty simplicial set
              To:   Torus
              Defn: [] --> []
        """
        from sage.topology.simplicial_set_examples import Empty
        return Empty().Hom(self)({})

    def identity(self):
        """
        Return the identity map on this simplicial set.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: S3.identity()
            Simplicial set endomorphism of S^3
              Defn: Identity map

            sage: # needs sage.groups
            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)
            sage: one = BC3.identity()
            sage: [(sigma, one(sigma)) for sigma in BC3.n_cells(2)]
            [(f * f, f * f),
             (f * f^2, f * f^2),
             (f^2 * f, f^2 * f),
             (f^2 * f^2, f^2 * f^2)]
        """
        return self.Hom(self).identity()

    def constant_map(self, codomain=None, point=None):
        """
        Return a constant map with this simplicial set as its domain.

        INPUT:

        - ``codomain`` -- (default: ``None``) if ``None``, the
          codomain is the standard one-point space constructed by
          :func:`Point`. Otherwise, either the codomain must be a
          pointed simplicial set, in which case the map is constant at
          the base point, or ``point`` must be specified.
        - ``point`` -- (default: ``None``) if specified, it
          must be a 0-simplex in the codomain, and it will be the
          target of the constant map

        EXAMPLES::

            sage: S4 = simplicial_sets.Sphere(4)
            sage: S4.constant_map()
            Simplicial set morphism:
              From: S^4
              To:   Point
              Defn: Constant map at *
            sage: S0 = simplicial_sets.Sphere(0)
            sage: S4.constant_map(codomain=S0)
            Simplicial set morphism:
              From: S^4
              To:   S^0
              Defn: Constant map at v_0

            sage: Sigma3 = groups.permutation.Symmetric(3)                              # needs sage.groups
            sage: Sigma3.nerve().constant_map()                                         # needs sage.groups
            Simplicial set morphism:
              From: Nerve of Symmetric group of order 3! as a permutation group
              To:   Point
              Defn: Constant map at *

        TESTS::

            sage: S0 = S0.unset_base_point()
            sage: S4.constant_map(codomain=S0)
            Traceback (most recent call last):
            ...
            ValueError: codomain is not pointed, so specify a target for the constant map
        """
        from sage.topology.simplicial_set_examples import Point
        if codomain is None:
            codomain = Point()
        return self.Hom(codomain).constant_map(point)

    def is_reduced(self):
        """
        Return ``True`` if this simplicial set has only one vertex.

        EXAMPLES::

            sage: simplicial_sets.Sphere(0).is_reduced()
            False
            sage: simplicial_sets.Sphere(3).is_reduced()
            True
        """
        return len(self.n_cells(0)) == 1

    def graph(self):
        """
        Return the 1-skeleton of this simplicial set, as a graph.

        EXAMPLES::

            sage: Delta3 = simplicial_sets.Simplex(3)
            sage: G = Delta3.graph()
            sage: G.edges(sort=True)
            [((0,), (1,), (0, 1)),
             ((0,), (2,), (0, 2)),
             ((0,), (3,), (0, 3)),
             ((1,), (2,), (1, 2)),
             ((1,), (3,), (1, 3)),
             ((2,), (3,), (2, 3))]

            sage: T = simplicial_sets.Torus()
            sage: T.graph()
            Looped multi-graph on 1 vertex
            sage: len(T.graph().edges(sort=False))
            3

            sage: # needs pyparsing
            sage: CP3 = simplicial_sets.ComplexProjectiveSpace(3)
            sage: G = CP3.graph()
            sage: len(G.vertices(sort=False))
            1
            sage: len(G.edges(sort=False))
            0

            sage: Sigma3 = groups.permutation.Symmetric(3)                              # needs sage.groups
            sage: Sigma3.nerve().is_connected()                                         # needs sage.groups
            True
        """
        from sage.graphs.graph import Graph

        G = Graph(loops=True, multiedges=True)
        for e in self.n_cells(1):
            G.add_edge(self.face(e, 0), self.face(e, 1), e)
        for v in self.n_cells(0):
            G.add_vertex(v)
        return G

    def is_connected(self):
        """
        Return ``True`` if this simplicial set is connected.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: X = T.disjoint_union(K)
            sage: T.is_connected()
            True
            sage: K.is_connected()
            True
            sage: X.is_connected()
            False
            sage: simplicial_sets.Sphere(0).is_connected()
            False
        """
        return self.graph().is_connected()

    def subsimplicial_set(self, simplices):
        """
        Return the sub-simplicial set of this simplicial set
        determined by ``simplices``, a set of nondegenerate simplices.

        INPUT:

        - ``simplices`` -- set, list, or tuple of nondegenerate
          simplices in this simplicial set, or a simplicial
          complex -- see below.

        Each sub-simplicial set comes equipped with an inclusion map
        to its ambient space, and you can easily recover its ambient
        space.

        If ``simplices`` is a simplicial complex, then the original
        simplicial set should itself have been converted from a
        simplicial complex, and ``simplices`` should be a subcomplex
        of that.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')

            sage: X = SimplicialSet({e: (v, w), f: (w, v)})
            sage: Y = X.subsimplicial_set([e]); Y
            Simplicial set with 3 non-degenerate simplices
            sage: Y.nondegenerate_simplices()
            [v, w, e]

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: K = SimplicialSet(S3)
            sage: tau = K.n_cells(3)[0]
            sage: tau.dimension()
            3
            sage: K.subsimplicial_set([tau])
            Simplicial set with 15 non-degenerate simplices

        A subsimplicial set knows about its ambient space and the
        inclusion map into it::

            sage: # needs sage.groups
            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: M = RP4.n_skeleton(2); M
            Simplicial set with 3 non-degenerate simplices
            sage: M.ambient_space()
            RP^4
            sage: M.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 3 non-degenerate simplices
              To:   RP^4
              Defn: [1, f, f * f] --> [1, f, f * f]

        An infinite ambient simplicial set::

            sage: # needs sage.groups
            sage: G = groups.misc.MultiplicativeAbelian([2])
            sage: B = simplicial_sets.ClassifyingSpace(G)
            sage: BxB = B.product(B)
            sage: BxB.n_cells(2)[5:]
            [(s_0 f, s_1 f), (s_1 f, f * f), (s_1 f, s_0 f), (s_1 s_0 1, f * f)]
            sage: BxB.subsimplicial_set(BxB.n_cells(2)[5:])
            Simplicial set with 8 non-degenerate simplices

        TESTS:

        Make sure vertices are treated properly::

            sage: X.subsimplicial_set([v]).nondegenerate_simplices()
            [v]
            sage: X.subsimplicial_set([v, w]).nondegenerate_simplices()
            [v, w]
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.subsimplicial_set([w]).nondegenerate_simplices()
            [w]

        Raise an error if an element of ``simplices`` is not actually
        in the original simplicial set::

            sage: sigma = AbstractSimplex(2, name='sigma_2')
            sage: Z = X.subsimplicial_set([e, sigma])
            Traceback (most recent call last):
            ...
            ValueError: not all simplices are in the original simplicial set

        Simplicial complexes::

            sage: X = simplicial_complexes.ComplexProjectivePlane()
            sage: Y = X._contractible_subcomplex()
            sage: CP2 = SimplicialSet(X)
            sage: sub = CP2.subsimplicial_set(Y)
            sage: CP2.f_vector()
            [9, 36, 84, 90, 36]
            sage: K = CP2.quotient(sub)
            sage: K.f_vector()
            [1, 0, 16, 30, 16]
            sage: K.homology()                                                          # needs sage.modules
            {0: 0, 1: 0, 2: Z, 3: 0, 4: Z}

        Try to construct a subcomplex from a simplicial complex which
        is not actually contained in ``self``::

            sage: Z = SimplicialComplex([[0,1,2,3,4]])
            sage: CP2.subsimplicial_set(Z)
            Traceback (most recent call last):
            ...
            ValueError: not all simplices are in the original simplicial set
        """
        # If simplices is a simplicial complex, turn it into a list of
        # nondegenerate simplices.
        from .simplicial_set_constructions import SubSimplicialSet
        if isinstance(simplices, SimplicialComplex):
            new = []
            for f in simplices.facets():
                d = f.dimension()
                found = False
                for x in self.n_cells(d):
                    if str(x) == str(tuple(sorted(f, key=str))):
                        new.append(x)
                        found = True
                        break
                if not found:
                    raise ValueError('not all simplices are in the original simplicial set')
            simplices = new

        if not self.is_finite():
            max_dim = max(sigma.dimension() for sigma in simplices)
            data = self.n_skeleton(max_dim).face_data()
            nondegenerate_simplices = self.nondegenerate_simplices(max_dim)
        else:
            data = self.face_data()
            nondegenerate_simplices = self.nondegenerate_simplices()
        vertices = set()
        keep = set(simplices)
        old_keep = set()
        while keep != old_keep:
            old_keep = copy.copy(keep)
            for x in old_keep:
                underlying = x.nondegenerate()
                if underlying not in data.keys():
                    raise ValueError('not all simplices are in the original simplicial set')
                keep.add(underlying)
                if underlying in data and data[underlying]:
                    keep.update([f.nondegenerate() for f in data[underlying]])
                else:
                    # x is a vertex
                    assert underlying.dimension() == 0
                    vertices.add(underlying)
        missing = set(nondegenerate_simplices).difference(keep)
        for x in missing:
            if x in data:
                del data[x]
        for x in vertices:
            data[x] = None
        return SubSimplicialSet(data, self)

    def chain_complex(self, dimensions=None, base_ring=ZZ, augmented=False,
                      cochain=False, verbose=False, subcomplex=None,
                      check=False):
        r"""
        Return the normalized chain complex.

        INPUT:

        - ``dimensions`` -- if ``None``, compute the chain complex in all
          dimensions.  If a list or tuple of integers, compute the
          chain complex in those dimensions, setting the chain groups
          in all other dimensions to zero.

        - ``base_ring`` -- (default: `\ZZ`) commutative ring

        - ``augmented`` -- boolean (default: ``False``); if ``True``,
          return the augmented chain complex (that is, include a class
          in dimension `-1` corresponding to the empty cell)

        - ``cochain`` -- boolean (default: ``False``); if ``True``,
          return the cochain complex (that is, the dual of the chain
          complex)

        - ``verbose`` -- boolean (default: ``False``); ignored

        - ``subcomplex`` -- (default: ``None``) if present,
          compute the chain complex relative to this subcomplex

        - ``check`` -- boolean (default: ``False``); if ``True``, make
          sure that the chain complex is actually a chain complex:
          the differentials are composable and their product is zero

        .. NOTE::

            If this simplicial set is not finite, you must specify
            dimensions in which to compute its chain complex via the
            argument ``dimensions``.

        EXAMPLES::

            sage: simplicial_sets.Sphere(5).chain_complex()                             # needs sage.modules
            Chain complex with at most 3 nonzero terms over Integer Ring

            sage: C3 = groups.misc.MultiplicativeAbelian([3])                           # needs sage.groups
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                            # needs sage.groups
            sage: BC3.chain_complex(range(4), base_ring=GF(3))                          # needs sage.groups sage.modules
            Chain complex with at most 4 nonzero terms over Finite Field of size 3

        TESTS::

            sage: BC3.chain_complex()                                                   # needs sage.groups
            Traceback (most recent call last):
            ...
            NotImplementedError: this simplicial set may be infinite, so specify dimensions when computing its chain complex
        """
        kwds = {'base_ring': base_ring, 'augmented': augmented, 'cochain': cochain,
                'verbose': verbose, 'subcomplex': subcomplex, 'check': check}
        if not self.is_finite():
            if dimensions is None:
                raise NotImplementedError('this simplicial set may be infinite, '
                                          'so specify dimensions when computing '
                                          'its chain complex')
            else:
                max_dim = max(dimensions)
                return SimplicialSet_finite.chain_complex(self.n_skeleton(max_dim+1),
                                                          dimensions=dimensions,
                                                          **kwds)
        return SimplicialSet_finite.chain_complex(self, dimensions=dimensions,
                                                  **kwds)

    def homology(self, dim=None, **kwds):
        r"""
        Return the (reduced) homology of this simplicial set.

        INPUT:

        - ``dim`` -- (default: ``None``) if ``None``, then
          return the homology in every dimension.  If ``dim`` is an
          integer or list, return the homology in the given
          dimensions.  (Actually, if ``dim`` is a list, return the
          homology in the range from ``min(dim)`` to ``max(dim)``.)

        - ``base_ring`` -- (default: ``ZZ``) commutative
          ring; must be ``ZZ`` or a field

        Other arguments are also allowed: see the documentation for
        :meth:`.cell_complex.GenericCellComplex.homology`.

        .. NOTE::

            If this simplicial set is not finite, you must specify
            dimensions in which to compute homology via the argument
            ``dim``.

        EXAMPLES::

            sage: simplicial_sets.Sphere(5).homology()                                  # needs sage.modules
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z}

            sage: C3 = groups.misc.MultiplicativeAbelian([3])                           # needs sage.groups
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                            # needs sage.groups
            sage: BC3.homology(range(4), base_ring=GF(3))                               # needs sage.groups sage.modules
            {0: Vector space of dimension 0 over Finite Field of size 3,
             1: Vector space of dimension 1 over Finite Field of size 3,
             2: Vector space of dimension 1 over Finite Field of size 3,
             3: Vector space of dimension 1 over Finite Field of size 3}

            sage: # needs sage.groups
            sage: C2 = groups.misc.MultiplicativeAbelian([2])
            sage: BC2 = simplicial_sets.ClassifyingSpace(C2)
            sage: BK = BC2.product(BC2)
            sage: BK.homology(range(4))                                                 # needs sage.modules
            {0: 0, 1: C2 x C2, 2: C2, 3: C2 x C2 x C2}

        TESTS::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: S3.homology(0)                                                        # needs sage.modules
            0
            sage: S3.homology((0,))                                                     # needs sage.modules
            {0: 0}
            sage: S3.homology(0, reduced=False)                                         # needs sage.modules
            Z

            sage: BC3.homology()                                                        # needs sage.groups sage.modules
            Traceback (most recent call last):
            ...
            NotImplementedError: this simplicial set may be infinite, so specify dimensions when computing homology
        """
        if not self.is_finite():
            if dim is None:
                raise NotImplementedError('this simplicial set may be infinite, so '
                                          'specify dimensions when computing homology')
            else:
                if isinstance(dim, (list, tuple, range)):
                    dim = list(dim)
                    max_dim = max(dim)
                    space = self.n_skeleton(max_dim+1)
                    min_dim = min(dim)
                    H = GenericCellComplex.homology(space, **kwds)
                    return {n: H[n] for n in H if min_dim <= n <= max_dim}
                else:
                    max_dim = dim
            space = self.n_skeleton(max_dim+1)
        else:
            space = self
        return GenericCellComplex.homology(space, dim=dim, **kwds)

    def cohomology(self, dim=None, **kwds):
        r"""
        Return the cohomology of this simplicial set.

        INPUT:

        - ``dim`` -- (default: ``None``) if ``None``, then
          return the homology in every dimension.  If ``dim`` is an
          integer or list, return the homology in the given
          dimensions.  (Actually, if ``dim`` is a list, return the
          homology in the range from ``min(dim)`` to ``max(dim)``.)

        - ``base_ring`` -- (default: ``ZZ``) commutative
          ring; must be ``ZZ`` or a field

        Other arguments are also allowed, the same as for the
        :meth:`homology` method -- see
        :meth:`.cell_complex.GenericCellComplex.homology` for complete
        documentation -- except that :meth:`homology` accepts a
        ``cohomology`` key word, while this function does not:
        ``cohomology`` is automatically true here.  Indeed, this
        function just calls :meth:`homology` with argument
        ``cohomology=True``.

        .. NOTE::

            If this simplicial set is not finite, you must specify
            dimensions in which to compute homology via the argument
            ``dim``.

        EXAMPLES::

            sage: simplicial_sets.KleinBottle().homology(1)                             # needs sage.modules
            Z x C2
            sage: simplicial_sets.KleinBottle().cohomology(1)                           # needs sage.modules
            Z
            sage: simplicial_sets.KleinBottle().cohomology(2)                           # needs sage.modules
            C2

        TESTS::

            sage: C3 = groups.misc.MultiplicativeAbelian([3])                           # needs sage.groups
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                            # needs sage.groups
            sage: BC3.cohomology()                                                      # needs sage.groups
            Traceback (most recent call last):
            ...
            NotImplementedError: this simplicial set may be infinite,
            so specify dimensions when computing homology
        """
        return self.homology(dim=dim, cohomology=True, **kwds)

    def betti(self, dim=None, subcomplex=None):
        r"""
        The Betti numbers of this simplicial complex as a dictionary
        (or a single Betti number, if only one dimension is given):
        the `i`-th Betti number is the rank of the `i`-th homology group.

        INPUT:

        - ``dim`` -- (default: ``None``) if ``None``, then
          return the homology in every dimension.  If ``dim`` is an
          integer or list, return the homology in the given
          dimensions.  (Actually, if ``dim`` is a list, return the
          homology in the range from ``min(dim)`` to ``max(dim)``.)

        - ``subcomplex`` -- (default: ``None``) a subcomplex
          of this cell complex.  Compute the Betti numbers of the
          homology relative to this subcomplex.

        .. NOTE::

            If this simplicial set is not finite, you must specify
            dimensions in which to compute Betti numbers via the
            argument ``dim``.

        EXAMPLES:

        Build the two-sphere as a three-fold join of a
        two-point space with itself::

            sage: simplicial_sets.Sphere(5).betti()                                     # needs sage.modules
            {0: 1, 1: 0, 2: 0, 3: 0, 4: 0, 5: 1}

            sage: C3 = groups.misc.MultiplicativeAbelian([3])                           # needs sage.groups
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3)                            # needs sage.groups
            sage: BC3.betti(range(4))                                                   # needs sage.groups sage.modules
            {0: 1, 1: 0, 2: 0, 3: 0}
        """
        dic = {}
        H = self.homology(dim, base_ring=QQ, subcomplex=subcomplex)
        try:
            for n in H.keys():
                dic[n] = H[n].dimension()
                if n == 0:
                    dic[n] += 1
        except AttributeError:
            return H.dimension()
        else:
            return dic

    def n_chains(self, n, base_ring=ZZ, cochains=False):
        r"""
        Return the free module of (normalized) chains in degree ``n``
        over ``base_ring``.

        This is the free module on the nondegenerate simplices in the
        given dimension.

        INPUT:

        - ``n`` -- integer
        - ``base_ring`` -- ring (default: `\ZZ`)
        - ``cochains`` -- boolean (default: ``False``); if
          ``True``, return cochains instead

        The only difference between chains and cochains is notation:
        the generator corresponding to the dual of a simplex
        ``sigma`` is written as ``'\chi_sigma'`` in the group of
        cochains.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: C = S3.n_chains(3, cochains=True)                                     # needs sage.modules
            sage: list(C.basis())                                                       # needs sage.modules
            [\chi_sigma_3]

            sage: # needs sage.groups
            sage: Sigma3 = groups.permutation.Symmetric(3)
            sage: BSigma3 = simplicial_sets.ClassifyingSpace(Sigma3)
            sage: list(BSigma3.n_chains(1).basis())                                     # needs sage.modules
            [(1,2), (1,2,3), (1,3), (1,3,2), (2,3)]
            sage: list(BSigma3.n_chains(1, cochains=True).basis())                      # needs sage.modules
            [\chi_(1,2), \chi_(1,2,3), \chi_(1,3), \chi_(1,3,2), \chi_(2,3)]
        """
        if self.is_finite():
            return GenericCellComplex.n_chains(self, n=n,
                                               base_ring=base_ring,
                                               cochains=cochains)

        from sage.homology.chains import Chains, Cochains

        n_cells = tuple(self.n_cells(n))
        if cochains:
            return Cochains(self, n, n_cells, base_ring)
        else:
            return Chains(self, n, n_cells, base_ring)

    def quotient(self, subcomplex, vertex_name='*'):
        """
        Return the quotient of this simplicial set by ``subcomplex``.

        That is, ``subcomplex`` is replaced by a vertex.

        INPUT:

        - ``subcomplex`` -- subsimplicial set of this simplicial set,
          or a list, tuple, or set of simplices defining a
          subsimplicial set

        - ``vertex_name`` -- string (default: ``'*'``); name to be given to the
          new vertex

        In Sage, from a quotient simplicial set, you can recover the
        ambient space, the subcomplex, and (if the ambient space is
        finite) the quotient map.

        Base points: if the original simplicial set has a base point
        not contained in ``subcomplex`` and if the original simplicial
        set is finite, then use its image as the base point for the
        quotient. In all other cases, ``*`` is the base point.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, w), f: (v, w)})
            sage: Y = X.quotient([f])
            sage: Y.nondegenerate_simplices()
            [*, e]
            sage: Y.homology(1)                                                         # needs sage.modules
            Z

            sage: E = SimplicialSet({e: (v, w)})
            sage: Z = E.quotient([v, w])
            sage: Z.nondegenerate_simplices()
            [*, e]
            sage: Z.homology(1)                                                         # needs sage.modules
            Z

            sage: F = E.quotient([v])
            sage: F.nondegenerate_simplices()
            [*, w, e]
            sage: F.base_point()
            *

            sage: # needs sage.groups
            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.homology(base_ring=GF(2))                                       # needs sage.modules
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 0 over Finite Field of size 2,
             3: Vector space of dimension 1 over Finite Field of size 2,
             4: Vector space of dimension 1 over Finite Field of size 2,
             5: Vector space of dimension 1 over Finite Field of size 2}
            sage: RP5_2.ambient()
            RP^5
            sage: RP5_2.subcomplex()
            Simplicial set with 3 non-degenerate simplices
            sage: RP5_2.quotient_map()
            Simplicial set morphism:
              From: RP^5
              To:   Quotient: (RP^5/Simplicial set with 3 non-degenerate simplices)
              Defn: [1, f, f * f, f * f * f, f * f * f * f, f * f * f * f * f]
                    --> [*, s_0 *, s_1 s_0 *, f * f * f, f * f * f * f, f * f * f * f * f]

        Behavior of base points::

            sage: K = simplicial_sets.Simplex(3)
            sage: K.is_pointed()
            False
            sage: L = K.subsimplicial_set([K.n_cells(1)[-1]])
            sage: L.nondegenerate_simplices()
            [(2,), (3,), (2, 3)]
            sage: K.quotient([K.n_cells(1)[-1]]).base_point()
            *

            sage: K = K.set_base_point(K.n_cells(0)[0])
            sage: K.base_point()
            (0,)
            sage: L = K.subsimplicial_set([K.n_cells(1)[-1]])
            sage: L.nondegenerate_simplices()
            [(2,), (3,), (2, 3)]
            sage: K.quotient(L).base_point()
            (0,)

        TESTS::

            sage: pt = RP5.quotient(RP5.n_skeleton(5)); pt                              # needs sage.groups
            Quotient: (RP^5/RP^5)
            sage: len(pt.nondegenerate_simplices())                                     # needs sage.groups
            1
        """
        from .simplicial_set_constructions import SubSimplicialSet
        from .simplicial_set_constructions import QuotientOfSimplicialSet, \
            QuotientOfSimplicialSet_finite
        if not isinstance(subcomplex, SimplicialSet_finite):
            # If it's not a simplicial set, subcomplex should be a
            # list, tuple, or set of simplices, so form the actual
            # subcomplex:
            subcomplex = self.subsimplicial_set(subcomplex)
        else:
            # Test whether subcomplex is actually a subcomplex of
            # self.
            if (not isinstance(subcomplex, SubSimplicialSet)
                    and subcomplex.ambient_space() == self):
                raise ValueError('the "subcomplex" is not actually a subcomplex')
        if self.is_finite():
            return QuotientOfSimplicialSet_finite(subcomplex.inclusion_map(),
                                                  vertex_name=vertex_name)
        else:
            return QuotientOfSimplicialSet(subcomplex.inclusion_map(),
                                           vertex_name=vertex_name)

    def disjoint_union(self, *others):
        """
        Return the disjoint union of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        As long as the factors are all finite, the inclusion map from
        each factor is available. Any factors which are empty are
        ignored completely: they do not appear in the list of factors,
        etc.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({f: (v, w)})
            sage: Z = X.disjoint_union(Y)

        Since ``X`` and ``Y`` have simplices in common, Sage uses a
        copy of ``Y`` when constructing the disjoint union. Note the
        name conflict in the list of simplices: ``v`` appears twice::

            sage: Z = X.disjoint_union(Y)
            sage: Z.nondegenerate_simplices()
            [v, v, w, e, f]

        Factors and inclusion maps::

            sage: T = simplicial_sets.Torus()
            sage: S2 = simplicial_sets.Sphere(2)
            sage: A = T.disjoint_union(S2)
            sage: A.factors()
            (Torus, S^2)
            sage: i = A.inclusion_map(0)
            sage: i.domain()
            Torus
            sage: i.codomain()
            Disjoint union: (Torus u S^2)

        Empty factors are ignored::

            sage: from sage.topology.simplicial_set_examples import Empty
            sage: E = Empty()
            sage: K = S2.disjoint_union(S2, E, E, S2)
            sage: K == S2.disjoint_union(S2, S2)
            True
            sage: K.factors()
            (S^2, S^2, S^2)
        """
        from .simplicial_set_constructions import DisjointUnionOfSimplicialSets, \
            DisjointUnionOfSimplicialSets_finite
        if all(space.is_finite() for space in [self] + list(others)):
            return DisjointUnionOfSimplicialSets_finite((self,) + others)
        else:
            return DisjointUnionOfSimplicialSets((self,) + others)

    def coproduct(self, *others):
        """
        Return the coproduct of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        If these simplicial sets are pointed, return their wedge sum;
        if they are not, return their disjoint union. If some are
        pointed and some are not, raise an error: it is not clear in
        which category to work.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: D3 = simplicial_sets.Simplex(3)
            sage: Y = S2.unset_base_point()
            sage: Z = K.unset_base_point()

            sage: S2.coproduct(K).is_pointed()
            True
            sage: S2.coproduct(K)
            Wedge: (S^2 v Klein bottle)
            sage: D3.coproduct(Y, Z).is_pointed()
            False
            sage: D3.coproduct(Y, Z)
            Disjoint union: (3-simplex u Simplicial set with 2 non-degenerate simplices
                             u Simplicial set with 6 non-degenerate simplices)

        The coproduct comes equipped with an inclusion map from each
        summand, as long as the summands are all finite::

            sage: S2.coproduct(K).inclusion_map(0)
            Simplicial set morphism:
              From: S^2
              To:   Wedge: (S^2 v Klein bottle)
              Defn: [v_0, sigma_2] --> [*, sigma_2]
            sage: D3.coproduct(Y, Z).inclusion_map(2)
            Simplicial set morphism:
              From: Simplicial set with 6 non-degenerate simplices
              To:   Disjoint union: (3-simplex
                                     u Simplicial set with 2 non-degenerate simplices
                                     u Simplicial set with 6 non-degenerate simplices)
              Defn: [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}]
                    --> [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}]

        TESTS::

            sage: D3.coproduct(S2, Z)
            Traceback (most recent call last):
            ...
            ValueError: some, but not all, of the simplicial sets are pointed,
            so the categorical coproduct is not defined: the category is ambiguous
        """
        if self.is_pointed() and all(X.is_pointed() for X in others):
            return self.wedge(*others)
        if self.is_pointed() or any(X.is_pointed() for X in others):
            raise ValueError('some, but not all, of the simplicial sets are pointed, '
                             'so the categorical coproduct is not defined: the '
                             'category is ambiguous')
        return self.disjoint_union(*others)

    def product(self, *others):
        r"""
        Return the product of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        If `X` and `Y` are simplicial sets, then their product `X
        \times Y` is defined to be the simplicial set with
        `n`-simplices `X_n \times Y_n`. See
        :class:`.simplicial_set_constructions.ProductOfSimplicialSets`
        for more information.

        If a simplicial set is constructed as a product, the factors
        are recorded and are accessible via the method
        :meth:`.simplicial_set_constructions.Factors.factors`.
        If each factor is finite, then you can also construct the
        projection maps onto each factor, the wedge as a subcomplex,
        and the fat wedge as a subcomplex.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, w)})
            sage: square = X.product(X)

        ``square`` is now the standard triangulation of the square: 4
        vertices, 5 edges (the four on the border and the diagonal), 2
        triangles::

            sage: square.f_vector()
            [4, 5, 2]

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: T.homology(reduced=False)                                             # needs sage.modules
            {0: Z, 1: Z x Z, 2: Z}

        Since ``S1`` is pointed, so is ``T``::

            sage: S1.is_pointed()
            True
            sage: S1.base_point()
            v_0
            sage: T.is_pointed()
            True
            sage: T.base_point()
            (v_0, v_0)

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2xS3 = S2.product(S3)
            sage: S2xS3.homology(reduced=False)                                         # needs sage.modules
            {0: Z, 1: 0, 2: Z, 3: Z, 4: 0, 5: Z}

            sage: S2xS3.factors() == (S2, S3)
            True
            sage: S2xS3.factors() == (S3, S2)
            False

            sage: # needs sage.groups
            sage: G = groups.misc.MultiplicativeAbelian([2])
            sage: B = simplicial_sets.ClassifyingSpace(G)
            sage: B.rename('RP^oo')
            sage: X = B.product(B, S2); X
            RP^oo x RP^oo x S^2
            sage: X.factor(1)
            RP^oo
            sage: X.factors()
            (RP^oo, RP^oo, S^2)

        Projection maps and wedges::

            sage: S2xS3.projection_map(0)
            Simplicial set morphism:
              From: S^2 x S^3
              To:   S^2
              Defn: ...
            sage: S2xS3.wedge_as_subset().homology()                                    # needs sage.modules
            {0: 0, 1: 0, 2: Z, 3: Z}

        In the case of pointed simplicial sets, there is an inclusion
        of each factor into the product. These are not automatically
        defined in Sage, but they are easy to construct using identity
        maps and constant maps and the universal property of the
        product::

            sage: one = S2.identity()
            sage: const = S2.constant_map(codomain=S3)
            sage: S2xS3.universal_property(one, const)
            Simplicial set morphism:
              From: S^2
              To:   S^2 x S^3
              Defn: [v_0, sigma_2] --> [(v_0, v_0), (sigma_2, s_1 s_0 v_0)]
        """
        from .simplicial_set_constructions import ProductOfSimplicialSets, \
            ProductOfSimplicialSets_finite
        if self.is_finite() and all(X.is_finite() for X in others):
            return ProductOfSimplicialSets_finite((self,) + others)
        else:
            return ProductOfSimplicialSets((self,) + others)

    cartesian_product = product

    def pushout(self, *maps):
        r"""
        Return the pushout obtained from given ``maps``.

        INPUT:

        - ``maps`` -- several maps of simplicial sets, each of which
          has this simplicial set as its domain

        If only a single map `f: X \to Y` is given, then return
        `Y`. If more than one map is given, say `f_i: X \to Y_i` for
        `0 \leq i \leq m`, then return the pushout defined by those
        maps. If no maps are given, return the empty simplicial set.

        In addition to the defining maps `f_i` used to construct the
        pushout `P`, there are also maps `\bar{f}_i: Y_i \to P`, which
        we refer to as *structure maps*. The pushout also has a
        universal property: given maps `g_i: Y_i \to Z` such that `g_i
        f_i = g_j f_j` for all `i`, `j`, then there is a unique map
        `g: P \to Z` making the appropriate diagram commute: that is,
        `g \bar{f}_i = g_i` for all `i`.

        In Sage, a pushout is equipped with its defining maps, and as
        long as the simplicial sets involved are finite, you can also
        access the structure maps and the universal property.

        EXAMPLES:

        Construct the 4-sphere as a quotient of a 4-simplex::

            sage: K = simplicial_sets.Simplex(4)
            sage: L = K.n_skeleton(3)
            sage: S4 = L.pushout(L.constant_map(), L.inclusion_map()); S4
            Pushout of maps:
              Simplicial set morphism:
                From: Simplicial set with 30 non-degenerate simplices
                To:   Point
                Defn: Constant map at *
              Simplicial set morphism:
                From: Simplicial set with 30 non-degenerate simplices
                To:   4-simplex
                Defn: [(0,), (1,), (2,), (3,), (4,),
                       (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4),
                       (2, 3), (2, 4), (3, 4),
                       (0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4),
                       (0, 3, 4), (1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4),
                       (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4), (0, 2, 3, 4),
                       (1, 2, 3, 4)]
                       --> [(0,), (1,), (2,), (3,), (4,),
                            (0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4),
                            (2, 3), (2, 4), (3, 4),
                            (0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 2, 3), (0, 2, 4), (0, 3, 4),
                            (1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4),
                            (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4), (0, 2, 3, 4), (1, 2, 3, 4)]
            sage: len(S4.nondegenerate_simplices())
            2
            sage: S4.homology(4)                                                        # needs sage.modules
            Z

        The associated maps::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: W = S1.wedge(T)  # wedge, constructed as a pushout
            sage: W.defining_map(1)
            Simplicial set morphism:
              From: Point
              To:   S^1 x S^1
              Defn: Constant map at (v_0, v_0)
            sage: W.structure_map(0)
            Simplicial set morphism:
              From: S^1
              To:   Wedge: (S^1 v S^1 x S^1)
              Defn: [v_0, sigma_1] --> [*, sigma_1]

            sage: f = S1.Hom(T)({S1.n_cells(0)[0]: K.n_cells(0)[0],
            ....:                S1.n_cells(1)[0]: K.n_cells(1)[0]})

        The maps `f: S^1 \to T` and `1: T \to T` induce a map `S^1 \vee T \to T`::

            sage: g = W.universal_property(f, Hom(T,T).identity())
            sage: g.domain() == W
            True
            sage: g.codomain() == T
            True

        TESTS::

            sage: K = simplicial_sets.Simplex(5)
            sage: K.pushout()
            Empty simplicial set

            sage: S0 = simplicial_sets.Sphere(0)
            sage: pt_map = S0.base_point_map()
            sage: pt_map.domain().pushout(pt_map) == S0
            True

            sage: K.pushout(K.constant_map(), pt_map)
            Traceback (most recent call last):
            ...
            ValueError: the domains of the maps must be equal
        """
        from .simplicial_set_constructions import PushoutOfSimplicialSets, \
            PushoutOfSimplicialSets_finite
        if any(self != f.domain() for f in maps):
            raise ValueError('the domains of the maps must be equal')
        if not maps:
            return PushoutOfSimplicialSets_finite()
        if all(f.codomain().is_finite() for f in maps):
            return PushoutOfSimplicialSets_finite(maps)
        else:
            return PushoutOfSimplicialSets(maps)

    def pullback(self, *maps):
        r"""
        Return the pullback obtained from given ``maps``.

        INPUT:

        - ``maps`` -- several maps of simplicial sets, each of which
          has this simplicial set as its codomain

        If only a single map `f: X \to Y` is given, then return
        `X`. If more than one map is given, say `f_i: X_i \to Y` for
        `0 \leq i \leq m`, then return the pullback defined by those
        maps. If no maps are given, return the one-point simplicial
        set.

        In addition to the defining maps `f_i` used to construct the
        pullback `P`, there are also maps `\bar{f}_i: P \to X_i`,
        which we refer to as *structure maps* or *projection
        maps*. The pullback also has a universal property: given maps
        `g_i: Z \to X_i` such that `f_i g_i = f_j g_j` for all `i`,
        `j`, then there is a unique map `g: Z \to P` making the
        appropriate diagram commute: that is, `\bar{f}_i g = g_i` for
        all `i`. For example, given maps `f: X \to Y` and `g: X \to
        Z`, there is an induced map `g: X \to Y \times Z`.

        In Sage, a pullback is equipped with its defining maps, and as
        long as the simplicial sets involved are finite, you can also
        access the structure maps and the universal property.

        EXAMPLES:

        Construct a product as a pullback::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: pt = simplicial_sets.Point()
            sage: P = pt.pullback(S2.constant_map(), S2.constant_map())
            sage: P.homology(2)                                                         # needs sage.modules
            Z x Z

        If the pullback is defined via maps `f_i: X_i \to Y`, then
        there are structure maps `\bar{f}_i: Y_i \to P`. The structure
        maps are only available in Sage when all of the maps involved
        have finite domains. ::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: P = S2.pullback(one, one)
            sage: P.homology()                                                          # needs sage.modules
            {0: 0, 1: 0, 2: Z}

            sage: P.defining_map(0) == one
            True
            sage: P.structure_map(1)
            Simplicial set morphism:
              From: Pullback of maps:
              Simplicial set endomorphism of S^2
                Defn: Identity map
              Simplicial set endomorphism of S^2
                Defn: Identity map
              To:   S^2
              Defn: [(v_0, v_0), (sigma_2, sigma_2)] --> [v_0, sigma_2]
            sage: P.structure_map(0).domain() == P
            True
            sage: P.structure_map(0).codomain() == S2
            True

        The universal property::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]: K.n_cells(0)[0],
            ....:                S1.n_cells(1)[0]: K.n_cells(1)[0]})
            sage: D = S1.cone()          # the cone C(S^1)
            sage: g = D.map_from_base()  # map from S^1 to C(S^1)
            sage: P = T.product(D)
            sage: h = P.universal_property(f, g)
            sage: h.domain() == S1
            True
            sage: h.codomain() == P
            True

        TESTS::

            sage: pt.pullback(S2.constant_map(), S2.base_point_map())
            Traceback (most recent call last):
            ...
            ValueError: the codomains of the maps must be equal
        """
        from .simplicial_set_constructions import PullbackOfSimplicialSets, \
            PullbackOfSimplicialSets_finite
        if any(self != f.codomain() for f in maps):
            raise ValueError('the codomains of the maps must be equal')
        if not maps:
            return PullbackOfSimplicialSets_finite()
        if self.is_finite() and all(f.domain().is_finite() for f in maps):
            return PullbackOfSimplicialSets_finite(maps)
        else:
            return PullbackOfSimplicialSets(maps)

    # Ideally, this would be defined at the category level and only
    # for pointed simplicial sets, but the abstract_method "wedge" in
    # cell_complex.py would shadow that.
    def wedge(self, *others):
        r"""
        Return the wedge sum of this pointed simplicial set with ``others``.

        - ``others`` -- one or several simplicial sets

        This constructs the quotient of the disjoint union in which
        the base points of all of the simplicial sets have been
        identified. This is the coproduct in the category of pointed
        simplicial sets.

        This raises an error if any of the factors is not pointed.

        From the wedge, you can access the factors, and if the
        simplicial sets involved are all finite, you can also access
        the inclusion map of each factor into the wedge, as well as
        the projection map onto each factor.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: W = X.wedge(Y)
            sage: W.nondegenerate_simplices()
            [*, e, f]
            sage: W.homology()                                                          # needs sage.modules
            {0: 0, 1: Z x Z}
            sage: S2 = simplicial_sets.Sphere(2)
            sage: X.wedge(S2).homology(reduced=False)                                   # needs sage.modules
            {0: Z, 1: Z, 2: Z}
            sage: X.wedge(X).nondegenerate_simplices()
            [*, e, e]

            sage: S3 = simplicial_sets.Sphere(3)
            sage: W = S2.wedge(S3, S2)
            sage: W.inclusion_map(2)
            Simplicial set morphism:
              From: S^2
              To:   Wedge: (S^2 v S^3 v S^2)
              Defn: [v_0, sigma_2] --> [*, sigma_2]
            sage: W.projection_map(1)
            Simplicial set morphism:
              From: Wedge: (S^2 v S^3 v S^2)
              To:   Quotient: (Wedge: (S^2 v S^3 v S^2)/Simplicial set with 3 non-degenerate simplices)
              Defn: [*, sigma_2, sigma_2, sigma_3] --> [*, s_1 s_0 *, s_1 s_0 *, sigma_3]

        Note that the codomain of the projection map is not identical
        to the original ``S2``, but is instead a quotient of the wedge
        which is isomorphic to ``S2``::

            sage: S2.f_vector()
            [1, 0, 1]
            sage: W.projection_map(2).codomain().f_vector()
            [1, 0, 1]
            sage: (W.projection_map(2) * W.inclusion_map(2)).is_bijective()
            True

        TESTS::

            sage: Z = SimplicialSet({e: (v,w)})
            sage: X.wedge(Z)
            Traceback (most recent call last):
            ...
            ValueError: the simplicial sets must be pointed
        """
        from .simplicial_set_constructions import WedgeOfSimplicialSets, \
            WedgeOfSimplicialSets_finite
        if all(space.is_finite() for space in [self] + list(others)):
            return WedgeOfSimplicialSets_finite((self,) + others)
        else:
            return WedgeOfSimplicialSets((self,) + others)

    def cone(self):
        r"""
        Return the (reduced) cone on this simplicial set.

        If this simplicial set `X` is not pointed, construct the
        ordinary cone: add a point `v` (which will become the base
        point) and for each simplex `\sigma` in `X`, add both `\sigma`
        and a simplex made up of `v` and `\sigma` (topologically, form
        the join of `v` and `\sigma`).

        If this simplicial set is pointed, then construct the reduced
        cone: take the quotient of the unreduced cone by the 1-simplex
        connecting the old base point to the new one.

        In either case, as long as the simplicial set is finite, it
        comes equipped in Sage with a map from it into the cone.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone()  # unreduced cone, since X not pointed
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *

        `X` as a subset of the cone, and also the map from `X`, in the
        unreduced case::

            sage: CX.base_as_subset()
            Simplicial set with 2 non-degenerate simplices
            sage: CX.map_from_base()
            Simplicial set morphism:
            From: Simplicial set with 2 non-degenerate simplices
              To:   Cone of Simplicial set with 2 non-degenerate simplices
              Defn: [v, e] --> [v, e]

        In the reduced case, only the map from `X` is available::

            sage: X = X.set_base_point(v)
            sage: CX = X.cone()  # reduced cone
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]
            sage: CX.map_from_base()
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   Reduced cone of Simplicial set with 2 non-degenerate simplices
              Defn: [v, e] --> [*, e]
        """
        from .simplicial_set_constructions import \
            ConeOfSimplicialSet, ConeOfSimplicialSet_finite, \
            ReducedConeOfSimplicialSet, ReducedConeOfSimplicialSet_finite
        if self.is_pointed():
            if self.is_finite():
                return ReducedConeOfSimplicialSet_finite(self)
            else:
                return ReducedConeOfSimplicialSet(self)
        if self.is_finite():
            return ConeOfSimplicialSet_finite(self)
        else:
            return ConeOfSimplicialSet(self)

    def suspension(self, n=1):
        """
        Return the (reduced) `n`-th suspension of this simplicial set.

        INPUT:

        - ``n`` -- integer (default: 1); suspend this many times

        If this simplicial set `X` is not pointed, return the
        suspension: the quotient `CX/X`, where `CX` is the (ordinary,
        unreduced) cone on `X`. If `X` is pointed, then use the
        reduced cone instead, and so return the reduced suspension.

        EXAMPLES::

            sage: # needs sage.groups
            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: SigmaRP4 = RP4.suspension()
            sage: S1_smash_RP4 = S1.smash_product(RP4)
            sage: SigmaRP4.homology() == S1_smash_RP4.homology()
            True

        The version of the suspension obtained by the smash product is
        typically less efficient than the reduced suspension produced
        here::

            sage: SigmaRP4.f_vector()                                                   # needs sage.groups
            [1, 0, 1, 1, 1, 1]
            sage: S1_smash_RP4.f_vector()                                               # needs sage.groups
            [1, 1, 4, 6, 8, 5]

        TESTS::

            sage: RP4.suspension(-3)                                                    # needs sage.groups
            Traceback (most recent call last):
            ...
            ValueError: n must be nonnegative
        """
        from .simplicial_set_constructions import \
            SuspensionOfSimplicialSet, SuspensionOfSimplicialSet_finite
        if n < 0:
            raise ValueError('n must be nonnegative')
        if n == 0:
            return self
        if self.is_finite():
            Sigma = SuspensionOfSimplicialSet_finite(self)
        else:
            Sigma = SuspensionOfSimplicialSet(self)
        if n == 1:
            return Sigma
        return Sigma.suspension(n-1)

    def join(self, *others):
        """
        The join of this simplicial set with ``others``.

        Not implemented. See
        https://ncatlab.org/nlab/show/join+of+simplicial+sets for a
        few descriptions, for anyone interested in implementing
        this. See also P. J. Ehlers and Tim Porter, Joins for
        (Augmented) Simplicial Sets, Jour. Pure Applied Algebra, 145
        (2000) 37-44 :arxiv:`9904039`.

        - ``others`` -- one or several simplicial sets

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(2)
            sage: K.join(K)
            Traceback (most recent call last):
            ...
            NotImplementedError: joins are not implemented for simplicial sets
        """
        raise NotImplementedError('joins are not implemented for simplicial sets')

    def reduce(self):
        """
        Reduce this simplicial set.

        That is, take the quotient by a spanning tree of the
        1-skeleton, so that the resulting simplicial set has only one
        vertex. This only makes sense if the simplicial set is
        connected, so raise an error if not. If already reduced,
        return itself.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(2)
            sage: K.is_reduced()
            False
            sage: X = K.reduce()
            sage: X.is_reduced()
            True

        ``X`` is reduced, so calling ``reduce`` on it again
        returns ``X`` itself::

            sage: X is X.reduce()
            True
            sage: K is K.reduce()
            False

        Raise an error for disconnected simplicial sets::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: S0.reduce()
            Traceback (most recent call last):
            ...
            ValueError: this simplicial set is not connected
        """
        if self.is_reduced():
            return self
        if not self.is_connected():
            raise ValueError("this simplicial set is not connected")
        graph = self.graph()
        spanning_tree = [e[2] for e in graph.min_spanning_tree()]
        return self.quotient(spanning_tree)

    def _Hom_(self, other, category=None):
        """
        Return the set of simplicial maps between simplicial sets
        ``self`` and ``other``.

        INPUT:

        - ``other`` -- another simplicial set
        - ``category`` -- (optional) the category in which to compute
          the maps. By default this is ``SimplicialSets``, and it must
          be a subcategory of this or else an error is raised.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3._Hom_(S2)
            Set of Morphisms from S^3 to S^2 in Category of finite pointed simplicial sets
            sage: Hom(S3, S2)
            Set of Morphisms from S^3 to S^2 in Category of finite pointed simplicial sets
            sage: K4 = simplicial_sets.Simplex(4)
            sage: S3._Hom_(K4)
            Set of Morphisms from S^3 to 4-simplex in Category of finite simplicial sets
        """
        # Import this here to prevent circular imports.
        from sage.topology.simplicial_set_morphism import SimplicialSetHomset
        # Error-checking on the ``category`` argument is done when
        # calling Hom(X,Y), so no need to do it again here.
        if category is None:
            if self.is_finite() and other.is_finite():
                if self.is_pointed() and other.is_pointed():
                    category = SimplicialSets().Finite().Pointed()
                else:
                    category = SimplicialSets().Finite()
            else:
                if self.is_pointed() and other.is_pointed():
                    category = SimplicialSets().Pointed()
                else:
                    category = SimplicialSets()
        return SimplicialSetHomset(self, other, category=category)

    def rename_latex(self, s):
        """
        Rename or set the LaTeX name for this simplicial set.

        INPUT:

        - ``s`` -- string; the LaTeX representation. Or ``s`` can be
          ``None``, in which case the LaTeX name is unset.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None}, latex_name='*')
            sage: latex(X)
            *
            sage: X.rename_latex('x_0')
            sage: latex(X)
            x_0
        """
        self._latex_name = s

    def _latex_(self):
        r"""
        LaTeX representation.

        If ``latex_name`` is set when the simplicial set is defined,
        or if :meth:`rename_latex` is used to set the LaTeX name, use
        that. Otherwise, use its string representation.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None}, latex_name='*')
            sage: latex(X)
            *
            sage: X.rename_latex('y_0')
            sage: latex(X)
            y_0
            sage: X.rename_latex(None)
            sage: latex(X)
            Simplicial set with 1 non-degenerate simplex
            sage: X.rename('v')
            sage: latex(X)
            v
        """
        if hasattr(self, '_latex_name') and self._latex_name is not None:
            return self._latex_name
        return str(self)

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2)
            sage: SimplicialSet({tau: (degen, degen, degen), w: None})
            Simplicial set with 3 non-degenerate simplices
            sage: SimplicialSet({w: None})
            Simplicial set with 1 non-degenerate simplex

        Test names and renaming::

            sage: SimplicialSet({w: None}, name='pt')
            pt
            sage: K = SimplicialSet({w: None}, name='pt')
            sage: K.rename('point')
            sage: K
            point
        """
        num = len(self.nondegenerate_simplices())
        if num == 1:
            return "Simplicial set with 1 non-degenerate simplex"
        return "Simplicial set with {} non-degenerate simplices".format(num)


class SimplicialSet_finite(SimplicialSet_arbitrary, GenericCellComplex):
    r"""
    A finite simplicial set.

    A simplicial set `X` is a collection of sets `X_n`, the
    *n-simplices*, indexed by the nonnegative integers, together with
    face maps `d_i` and degeneracy maps `s_j`.  A simplex is
    *degenerate* if it is in the image of some `s_j`, and a simplicial
    set is *finite* if there are only finitely many non-degenerate
    simplices.

    INPUT:

    - ``data`` -- the data defining the simplicial set; see below for
      details

    - ``base_point`` -- (default: ``None``) 0-simplex in this
      simplicial set, its base point

    - ``name`` -- string (default: ``None``); the name of the
      simplicial set

    - ``check`` -- boolean (default: ``True``); if ``True``,
      check the simplicial identity on the face maps when defining the
      simplicial set

    - ``category`` -- (default: ``None``) the category in
      which to define this simplicial set. The default is either
      finite simplicial sets or finite pointed simplicial sets,
      depending on whether a base point is defined.

    - ``latex_name`` -- string (default: ``None``); the LaTeX
      representation of the simplicial set

    ``data`` should have one of the following forms: it could be a
    simplicial complex or `\Delta`-complex, in case it is converted to
    a simplicial set. Alternatively, it could be a dictionary. The
    keys are the nondegenerate simplices of the simplicial set, and
    the value corresponding to a simplex `\sigma` is a tuple listing
    the faces of `\sigma`. The 0-dimensional simplices may be omitted
    from ``data`` if they (or their degeneracies) are faces of other
    simplices; otherwise they must be included with value ``None``.

    See :mod:`.simplicial_set` and the methods for simplicial sets for
    more information and examples.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: u = AbstractSimplex(0, name='u')
        sage: v = AbstractSimplex(0, name='v')
        sage: w = AbstractSimplex(0, name='w')
        sage: e = AbstractSimplex(1, name='e')
        sage: f = AbstractSimplex(1, name='f')

    In the following simplicial set, ``u`` is an isolated vertex::

        sage: X = SimplicialSet({e: (v,w), f: (w,w), u: None})
        sage: X
        Simplicial set with 5 non-degenerate simplices
        sage: X.rename('X')
        sage: X
        X
        sage: X = SimplicialSet({e: (v,w), f: (w,w), u: None}, name='Y')
        sage: X
        Y
    """

    def __init__(self, data, base_point=None, name=None, check=True,
                 category=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: SimplicialSet({e: (v, v, v)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1
            sage: SimplicialSet({e: (v,)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1

        Base points::

            sage: SimplicialSet({e: (v,v)}, base_point=AbstractSimplex(0))
            Traceback (most recent call last):
            ...
            ValueError: the base point is not a simplex in this simplicial set
            sage: SimplicialSet({e: (v,v)}, base_point=e)
            Traceback (most recent call last):
            ...
            ValueError: the base "point" is not a zero-simplex

        Simplicial identity::

            sage: sigma = AbstractSimplex(2)
            sage: w = AbstractSimplex(0)
            sage: K = SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0))})
            sage: SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                        v.apply_degeneracies(0),
            ....:                        w.apply_degeneracies(0))})
            Traceback (most recent call last):
            ...
            ValueError: simplicial identity d_i d_j = d_{j-1} d_i fails in dimension 2

        Returning a copy of the original::

            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: SimplicialSet(S1) == S1
            False

        Test suites::

            sage: skip = ["_test_pickling", "_test_elements"]
            sage: TestSuite(S1).run(skip=skip)
            sage: TestSuite(simplicial_sets.Sphere(5)).run(skip=skip)
            sage: TestSuite(simplicial_sets.RealProjectiveSpace(6)).run(skip=skip)      # needs sage.groups
        """
        def face(sigma, i):
            """
            Return the i-th face of sigma, a simplex in this simplicial set.

            Once the simplicial set has been fully initialized, use
            the :meth:`face` method instead.
            """
            if sigma.is_nondegenerate():
                return data[sigma][i]
            else:
                underlying = sigma.nondegenerate()
                J, t = face_degeneracies(i, sigma.degeneracies())
                if t is None:
                    return underlying.apply_degeneracies(*J)
                else:
                    return data[underlying][t].apply_degeneracies(*J)

        if isinstance(data, GenericCellComplex):
            # Construct new data appropriately.
            if isinstance(data, SimplicialComplex):
                simplices = {}
                faces = {}
                for d in range(data.dimension()+1):
                    old_faces = faces
                    faces = {}
                    for idx, sigma in enumerate(data.n_cells(d)):
                        new_sigma = AbstractSimplex(d)
                        new_sigma.rename(str(tuple(sorted(sigma, key=str))))
                        if d > 0:
                            simplices[new_sigma] = [old_faces[_] for _ in sigma.faces()]
                        else:
                            simplices[new_sigma] = None
                        faces[sigma] = new_sigma
                data = simplices

            elif isinstance(data, DeltaComplex):
                simplices = {}
                current = []
                for d in range(data.dimension()+1):
                    faces = tuple(current)
                    current = []
                    for idx, sigma in enumerate(data.n_cells(d)):
                        new_sigma = AbstractSimplex(d)
                        # Name: Delta_{d,idx} where d is dimension,
                        # idx is its index in the list of d-simplices.
                        new_sigma.rename('Delta_{{{},{}}}'.format(d, idx))
                        if d > 0:
                            simplices[new_sigma] = [faces[_] for _ in sigma]
                        else:
                            simplices[new_sigma] = None
                        current.append(new_sigma)
                data = simplices
            elif isinstance(data, SimplicialSet_finite):
                data = dict(copy.deepcopy(data._data))
            else:
                raise NotImplementedError('I do not know how to convert this '
                                          'to a simplicial set')
        # Convert each value in data to a tuple, and then convert all
        # of data to a tuple, so that it is hashable.
        for x in data:
            if data[x]:
                if x.dimension() != len(data[x]) - 1:
                    raise ValueError('wrong number of faces for simplex '
                                     'in dimension {}'.format(x.dimension()))
                if not all(y.dimension() == x.dimension() - 1 for y in data[x]):
                    raise ValueError('faces of a {}-simplex have the wrong '
                                     'dimension'.format(x.dimension()))
                data[x] = tuple(data[x])

        # To obtain the non-degenerate simplices, look at both the
        # keys for data and also the underlying non-degenerate
        # simplices in its values.
        simplices = set(data.keys())
        for t in data.values():
            if t:
                simplices.update([_.nondegenerate() for _ in t])

        for x in simplices:
            if x not in data:
                # x had better be a vertex.
                assert x.dimension() == 0
                data[x] = None

        # Check the simplicial identity d_i d_j = d_{j-1} d_i.
        if check:
            for sigma in simplices:
                d = sigma.dimension()
                if d >= 2:
                    for j in range(d+1):
                        for i in range(j):
                            if face(face(sigma, j), i) != face(face(sigma, i), j-1):
                                raise ValueError('simplicial identity d_i d_j '
                                                 '= d_{{j-1}} d_i fails '
                                                 'in dimension {}'.format(d))

        # Now define the attributes for an instance of this class.
        # self._data: a tuple representing the defining data of the
        # simplicial set.
        self._data = tuple(data.items())
        # self._simplices: a sorted tuple of non-degenerate simplices.
        self._simplices = sorted(simplices)
        # self._basepoint: the base point, or None.
        if base_point is not None:
            if base_point not in simplices:
                raise ValueError('the base point is not a simplex in '
                                 'this simplicial set')
            if base_point.dimension() != 0:
                raise ValueError('the base "point" is not a zero-simplex')
            self._basepoint = base_point
        if category is None:
            if base_point is None:
                category = SimplicialSets().Finite()
            else:
                category = SimplicialSets().Finite().Pointed()
        Parent.__init__(self, category=category)
        if name:
            self.rename(name)
        self._latex_name = latex_name

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` are equal as simplicial sets.

        Two simplicial sets are equal if they have the same defining
        data.  This means that they have *the same* simplices in each
        dimension, not just that they have the same numbers of
        `n`-simplices for each `n` with corresponding face maps.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({e: (w, w)})
            sage: X == X
            True
            sage: X == SimplicialSet({e: (v, v)})
            True
            sage: X == Y
            False
        """
        if self.is_pointed():
            return (isinstance(other, SimplicialSet_finite)
                    and other.is_pointed()
                    and sorted(self._data) == sorted(other._data)
                    and self.base_point() == other.base_point())
        else:
            return (isinstance(other, SimplicialSet_finite)
                    and not other.is_pointed()
                    and sorted(self._data) == sorted(other._data))

    def __ne__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` are not equal as simplicial sets.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None})
            sage: Y = SimplicialSet({w: None})
            sage: X != X
            False
            sage: X != SimplicialSet({v: None})
            False
            sage: X != Y
            True
        """
        return not (self == other)

    # This is cached because it is used frequently in simplicial set
    # construction: the last two lines in the __init__ method access
    # dictionaries which use instances of SimplicialSet_finite as keys and so
    # computes their hash. If the tuple self._data is long, this can
    # take a long time.
    @cached_method
    def __hash__(self):
        """
        The hash is formed from that of the tuple ``self._data``.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None})
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2)
            sage: Y = SimplicialSet({tau: (degen, degen, degen)})

            sage: hash(X) # random
            17
            sage: hash(X) != hash(Y)
            True
        """
        if self.is_pointed():
            return hash(self._data) ^ hash(self.base_point())
        else:
            return hash(self._data)

    def __copy__(self):
        """
        Return a distinct copy of this simplicial set.

        The copy will not be equal to the original simplicial set.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: copy(T) == T
            False
            sage: T.n_cells(0)[0] == copy(T).n_cells(0)[0]
            False
            sage: T.homology() == copy(T).homology()                                    # needs sage.modules
            True
        """
        return SimplicialSet(dict(copy.deepcopy(self._data)))

    def face_data(self):
        """
        Return the face-map data -- a dictionary -- defining this simplicial set.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, w)})
            sage: X.face_data()[e]
            (v, w)

            sage: Y = SimplicialSet({v: None, w: None})
            sage: v in Y.face_data()
            True
            sage: Y.face_data()[v] is None
            True
        """
        return dict(self._data)

    def n_skeleton(self, n):
        """
        Return the `n`-skeleton of this simplicial set.

        That is, the subsimplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        INPUT:

        - ``n`` -- the dimension

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2, name='tau')
            sage: Y = SimplicialSet({tau: (degen, degen, degen), w: None})

        ``Y`` is the disjoint union of a 2-sphere, with vertex ``v``
        and non-degenerate 2-simplex ``tau``, and a point ``w``. ::

            sage: Y.nondegenerate_simplices()
            [v, w, tau]
            sage: Y.n_skeleton(1).nondegenerate_simplices()
            [v, w]
            sage: Y.n_skeleton(2).nondegenerate_simplices()
            [v, w, tau]
        """
        data = [x for x in self.nondegenerate_simplices()
                if x.dimension() <= n]
        return self.subsimplicial_set(data)

    def _facets_(self):
        r"""
        Return the list of facets of this simplicial set, where by
        "facet" we mean a non-degenerate simplex which is not a face
        of another non-degenerate simplex.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: T._facets_()
            [(s_0 sigma_1, s_1 sigma_1),
             (s_1 sigma_1, s_0 sigma_1)]
            sage: S5 = simplicial_sets.Sphere(5)
            sage: S5._facets_()
            [sigma_5]
            sage: simplicial_sets.Sphere(0)._facets_()
            [v_0, w_0]
        """
        faces = set()
        for dim in range(self.dimension(), 0, -1):
            for sigma in self.n_cells(dim):
                faces.update([tau.nondegenerate() for tau in self.faces(sigma)])
        return sorted(set(self.nondegenerate_simplices()).difference(faces))

    def f_vector(self):
        """
        Return the list of the number of non-degenerate simplices in each
        dimension.

        Unlike for some other cell complexes in Sage, this does not
        include the empty simplex in dimension `-1`; thus its `i`-th
        entry is the number of `i`-dimensional simplices.

        EXAMPLES::

            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.f_vector()
            [2]

            sage: e = AbstractSimplex(1)
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: S1.f_vector()
            [1, 1]
            sage: simplicial_sets.Sphere(3).f_vector()
            [1, 0, 0, 1]
        """
        return [len(self.n_cells(_)) for _ in range(self.dimension()+1)]

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic of this simplicial set: the
        alternating sum over `n \geq 0` of the number of
        nondegenerate `n`-simplices.

        EXAMPLES::

            sage: simplicial_sets.RealProjectiveSpace(4).euler_characteristic()         # needs sage.groups
            1
            sage: simplicial_sets.Sphere(6).euler_characteristic()
            2
            sage: simplicial_sets.KleinBottle().euler_characteristic()
            0
        """
        return sum([(-1)**n * num for (n, num) in enumerate(self.f_vector())])

    def chain_complex(self, dimensions=None, base_ring=ZZ, augmented=False,
                      cochain=False, verbose=False, subcomplex=None,
                      check=False):
        r"""
        Return the normalized chain complex.

        INPUT:

        - ``dimensions`` -- if ``None``, compute the chain complex in all
          dimensions.  If a list or tuple of integers, compute the
          chain complex in those dimensions, setting the chain groups
          in all other dimensions to zero.

        - ``base_ring`` -- commutative ring (default: `\ZZ`)

        - ``augmented`` -- boolean (default: ``False``); if ``True``,
          return the augmented chain complex (that is, include a class
          in dimension `-1` corresponding to the empty cell).

        - ``cochain`` -- boolean (default: ``False``); if ``True``,
          return the cochain complex (that is, the dual of the chain
          complex).

        - ``verbose`` -- boolean (default: ``False``); ignored

        - ``subcomplex`` -- (default: ``None``) if present,
          compute the chain complex relative to this subcomplex

        - ``check`` -- boolean (default: ``False``); if ``True``, make
          sure that the chain complex is actually a chain complex:
          the differentials are composable and their product is zero.

        The normalized chain complex of a simplicial set is isomorphic
        to the chain complex obtained by modding out by degenerate
        simplices, and the latter is what is actually constructed
        here.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.topology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: degen = v.apply_degeneracies(1, 0) # s_1 s_0 applied to v
            sage: sigma = AbstractSimplex(3)
            sage: S3 = SimplicialSet({sigma: (degen, degen, degen, degen)})  # the 3-sphere
            sage: S3.chain_complex().homology()
            {0: Z, 3: Z}
            sage: S3.chain_complex(augmented=True).homology()
            {-1: 0, 0: 0, 3: Z}
            sage: S3.chain_complex(dimensions=range(3), base_ring=QQ).homology()
            {0: Vector space of dimension 1 over Rational Field}

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)                          # needs sage.groups
            sage: RP2 = RP5.n_skeleton(2)                                               # needs sage.groups
            sage: RP5.chain_complex(subcomplex=RP2).homology()                          # needs sage.groups sage.modules
            {0: Z, 3: C2, 4: 0, 5: Z}

        TESTS:

        Convert some simplicial complexes and `\Delta`-complexes to
        simplicial sets, and compare homology calculations::

            sage: # needs sage.modules
            sage: T = simplicial_complexes.Torus()
            sage: T.homology() == SimplicialSet(T).homology()
            True
            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: RP2.homology() == SimplicialSet(RP2).homology()
            True
            sage: cohoRP2 = RP2.cohomology(base_ring=GF(2))
            sage: cohoRP2 == SimplicialSet(RP2).cohomology(base_ring=GF(2))
            True
        """
        from sage.homology.chain_complex import ChainComplex

        if dimensions is None:
            if not self.cells():  # Empty
                if cochain:
                    return ChainComplex({-1: matrix(base_ring, 0, 0)},
                                        degree_of_differential=1)
                return ChainComplex({0: matrix(base_ring, 0, 0)},
                                    degree_of_differential=-1)
            dimensions = list(range(self.dimension() + 1))
        else:
            if not isinstance(dimensions, (list, tuple, range)):
                dimensions = list(range(dimensions - 1, dimensions + 2))
            else:
                dimensions = [n for n in dimensions if n >= 0]
            if not dimensions:
                # Return the empty chain complex.
                if cochain:
                    return ChainComplex(base_ring=base_ring, degree=1)
                else:
                    return ChainComplex(base_ring=base_ring, degree=-1)

        differentials = {}
        # Convert the tuple self._data to a dictionary indexed by the
        # non-degenerate simplices.
        if subcomplex:
            X = self.quotient(subcomplex)
            face_data = X.face_data()
            nondegens = X.nondegenerate_simplices()
        else:
            face_data = self.face_data()
            nondegens = self.nondegenerate_simplices()
        # simplices: dictionary indexed by dimension, values the list
        # of non-degenerate simplices in that dimension.
        simplices = {}
        for sigma in nondegens:
            if sigma.dimension() in simplices:
                simplices[sigma.dimension()].append(sigma)
            else:
                simplices[sigma.dimension()] = [sigma]
        first = dimensions.pop(0)
        if first in simplices:
            rank = len(simplices[first])
            current = sorted(simplices[first])
        else:
            rank = 0
            current = []
        if augmented and first == 0:
            differentials[first-1] = matrix(base_ring, 0, 1)
            differentials[first] = matrix(base_ring, 1, rank,
                                          [1] * rank)
        else:
            differentials[first] = matrix(base_ring, 0, rank)

        for d in dimensions:
            old_rank = rank
            faces = {_[1]: _[0] for _ in enumerate(current)}
            if d in simplices:
                current = sorted(simplices[d])
                rank = len(current)
                # old_rank: number of simplices in dimension d-1.
                # faces: list of simplices in dimension d-1.
                # rank: number of simplices in dimension d.
                # current: list of simplices in dimension d.
                if not faces:
                    differentials[d] = matrix(base_ring, old_rank, rank)
                else:
                    matrix_data = {}
                    for col, sigma in enumerate(current):
                        sign = 1
                        for tau in face_data[sigma]:
                            if tau.is_nondegenerate():
                                row = faces[tau]
                                if (row, col) in matrix_data:
                                    matrix_data[(row, col)] += sign
                                else:
                                    matrix_data[(row, col)] = sign
                            sign *= -1

                    differentials[d] = matrix(base_ring, old_rank,
                                              rank, matrix_data)

            else:
                rank = 0
                current = []
                differentials[d] = matrix(base_ring, old_rank, rank)

        if cochain:
            new_diffs = {}
            for d in differentials:
                new_diffs[d-1] = differentials[d].transpose()
            return ChainComplex(new_diffs, degree_of_differential=1,
                                check=check)
        return ChainComplex(differentials, degree_of_differential=-1,
                            check=check)

    @cached_method
    def algebraic_topological_model(self, base_ring=None):
        r"""
        Return the algebraic topological model for this simplicial set
        with coefficients in ``base_ring``.

        The term "algebraic topological model" is defined by Pilarczyk
        and Réal [PR2015]_.

        INPUT:

        - ``base_ring`` -- coefficient ring (default: ``QQ``); must be a field

        Denote by `C` the chain complex associated to this simplicial
        set. The algebraic topological model is a chain complex `M`
        with zero differential, with the same homology as `C`, along
        with chain maps `\pi: C \to M` and `\iota: M \to C` satisfying
        `\iota \pi = 1_M` and `\pi \iota` chain homotopic to
        `1_C`. The chain homotopy `\phi` must satisfy

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

            sage: RP2 = simplicial_sets.RealProjectiveSpace(2)                          # needs sage.groups
            sage: phi, M = RP2.algebraic_topological_model(GF(2))                       # needs sage.groups
            sage: M.homology()                                                          # needs sage.groups sage.modules
            {0: Vector space of dimension 1 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}

            sage: T = simplicial_sets.Torus()
            sage: phi, M = T.algebraic_topological_model(QQ)                            # needs sage.modules
            sage: M.homology()                                                          # needs sage.modules
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 2 over Rational Field,
             2: Vector space of dimension 1 over Rational Field}
        """
        from sage.homology.algebraic_topological_model import algebraic_topological_model_delta_complex

        if base_ring is None:
            base_ring = QQ
        return algebraic_topological_model_delta_complex(self, base_ring)


# TODO: possibly turn SimplicialSet into a function, for example
# allowing for the construction of infinite simplicial sets.
SimplicialSet = SimplicialSet_finite


########################################################################
# Functions for manipulating face and degeneracy maps.

def standardize_degeneracies(*L):
    r"""
    Return list of indices of degeneracy maps in standard (decreasing)
    order.

    INPUT:

    - ``L`` -- list of integers representing a composition of
      degeneracies in a simplicial set

    OUTPUT:

    an equivalent list of degeneracies, standardized to be
    written in decreasing order, using the simplicial identity

    .. MATH::

        s_i s_j = s_{j+1} s_i \ \   \text{if }   i \leq j.

    For example, `s_0 s_2 = s_3 s_0` and `s_0 s_0 = s_1 s_0`.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import standardize_degeneracies
        sage: standardize_degeneracies(0, 0)
        (1, 0)
        sage: standardize_degeneracies(0, 0, 0, 0)
        (3, 2, 1, 0)
        sage: standardize_degeneracies(1, 2)
        (3, 1)

    TESTS::

        sage: standardize_degeneracies()
        ()
        sage: standardize_degeneracies(2, -1)
        Traceback (most recent call last):
        ...
        ValueError: degeneracies are indexed by nonnegative integers
        sage: standardize_degeneracies([2, 1])
        Traceback (most recent call last):
        ...
        TypeError: degeneracies are indexed by nonnegative integers; do not use an explicit list or tuple
    """
    J = list(L)
    for m in J:
        try:
            if Integer(m) < 0:
                raise ValueError('degeneracies are indexed by nonnegative integers')
        except TypeError:
            # Likely if called via standard_degeneracies([1,2,3])
            # rather than          standard_degeneracies(1,2,3).
            raise TypeError('degeneracies are indexed by nonnegative integers; do not use an explicit list or tuple')
    inadmissible = True
    while inadmissible:
        inadmissible = False
        for idx in range(len(J)-1):
            if J[idx] <= J[idx + 1]:
                inadmissible = True
                tmp = J[idx]
                J[idx] = J[idx + 1] + 1
                J[idx + 1] = tmp
    return tuple(J)


def all_degeneracies(n, l=1):
    r"""
    Return list of all composites of degeneracies (written in
    "admissible" form, i.e., as a strictly decreasing sequence) of
    length `l` on an `n`-simplex.

    INPUT:

    - ``n``, ``l`` -- integers

    On an `n`-simplex, one may apply the degeneracies `s_i` for `0
    \leq i \leq n`. Then on the resulting `n+1`-simplex, one may apply
    `s_i` for `0 \leq i \leq n+1`, and so on. But one also has to take
    into account the simplicial identity

    .. MATH::

        s_i s_j = s_{j+1} s_i \ \   \text{if }   i \leq j.

    There are `\binom{l+n}{n}` such composites: each non-degenerate
    `n`-simplex leads to `\binom{l+n}{n}` degenerate `l+n` simplices.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import all_degeneracies
        sage: all_degeneracies(0, 3)
        {(2, 1, 0)}
        sage: all_degeneracies(1, 1)
        {(0,), (1,)}
        sage: all_degeneracies(1, 3)
        {(2, 1, 0), (3, 1, 0), (3, 2, 0), (3, 2, 1)}
    """
    if l == 0:
        return set()
    if l == 1:
        return {(_,) for _ in range(n+1)}
    ans = set()
    for i in range(n+l):
        ans.update({tuple(standardize_degeneracies(*([i] + list(_))))
                        for _ in all_degeneracies(n, l-1)})
    return ans


def standardize_face_maps(*L):
    r"""
    Return list of indices of face maps in standard (non-increasing)
    order.

    INPUT:

    - ``L`` -- list of integers representing a composition of
      face maps in a simplicial set

    OUTPUT:

    an equivalent list of face maps, standardized to be
    written in non-increasing order, using the simplicial identity

    .. MATH::

        d_i d_j = d_{j-1} d_i \ \   \text{if }   i<j.

    For example, `d_0 d_2 = d_1 d_0` and `d_0 d_1 = d_0 d_0`.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import standardize_face_maps
        sage: standardize_face_maps(0, 1)
        (0, 0)
        sage: standardize_face_maps(0, 2)
        (1, 0)
        sage: standardize_face_maps(1, 3, 5)
        (3, 2, 1)
    """
    J = list(L)
    for m in J:
        if Integer(m) < 0:
            raise ValueError('faces are indexed by nonnegative integers')
    inadmissible = True
    while inadmissible:
        inadmissible = False
        for idx in range(len(J)-1):
            if J[idx] < J[idx + 1]:
                inadmissible = True
                tmp = J[idx]
                J[idx] = J[idx + 1] - 1
                J[idx + 1] = tmp
    return tuple(J)


def face_degeneracies(m, I):
    r"""
    Return the result of applying the face map `d_m` to the iterated
    degeneracy `s_I = s_{i_1} s_{i_2} ... s_{i_n}`.

    INPUT:

    - ``m`` -- integer
    - ``I`` -- tuple ``(i_1, i_2, ..., i_n)`` of integers; we assume
      that this sequence is strictly decreasing

    Using the simplicial identities (see :mod:`.simplicial_set`), we
    can rewrite

    .. MATH::

        d_m s_{i_1} s_{i_2} ... s_{i_n}

    in one of the forms

    .. MATH::

        s_{j_1} s_{j_2} ... s_{j_n} d_t, \quad
        s_{j_1} s_{j_2} ... s_{j_{n-1}}.

    OUTPUT: the pair ``(J, t)`` or ``(J, None)``; ``J`` is returned as
    a list

    EXAMPLES::

        sage: from sage.topology.simplicial_set import face_degeneracies
        sage: face_degeneracies(0, (1, 0))
        ([0], None)
        sage: face_degeneracies(1, (1, 0))
        ([0], None)
        sage: face_degeneracies(2, (1, 0))
        ([0], None)
        sage: face_degeneracies(3, (1, 0))
        ([1, 0], 1)
        sage: face_degeneracies(3, ())
        ([], 3)
    """
    if not I:
        return ([], m)
    J = []
    t = m
    for i in I:
        if t is None:
            J.append(i)
        elif t < i:
            J.append(i-1)
        elif t == i or t == i+1:
            t = None
        else:
            J.append(i)
            t -= 1
    return (J, t)


########################################################################

def shrink_simplicial_complex(K):
    """
    Convert the simplicial complex ``K`` to a "small" simplicial set.

    First convert ``K`` naively, then mod out by a large contractible
    subcomplex, as found by
    :meth:`.simplicial_complex.SimplicialComplex._contractible_subcomplex`.
    This will produce a simplicial set no larger than, and sometimes
    much smaller than, the initial simplicial complex.

    EXAMPLES::

        sage: from sage.topology.simplicial_set import shrink_simplicial_complex
        sage: K = simplicial_complexes.Simplex(3)
        sage: X = shrink_simplicial_complex(K)
        sage: X.f_vector()
        [1]

        sage: Y = simplicial_complexes.Sphere(2)
        sage: S2 = shrink_simplicial_complex(Y); S2
        Quotient: (Simplicial set with
                   14 non-degenerate simplices/Simplicial set with
                                               13 non-degenerate simplices)
        sage: S2.f_vector()
        [1, 0, 1]
        sage: S2.homology()                                                             # needs sage.modules
        {0: 0, 1: 0, 2: Z}

        sage: Z = simplicial_complexes.SurfaceOfGenus(3)
        sage: Z.f_vector()
        [1, 15, 57, 38]
        sage: Z.homology()                                                              # needs sage.modules
        {0: 0, 1: Z^6, 2: Z}
        sage: M = shrink_simplicial_complex(Z)
        sage: M.f_vector()  # random
        [1, 32, 27]
        sage: M.homology()                                                              # needs sage.modules
        {0: 0, 1: Z^6, 2: Z}
    """
    L = K._contractible_subcomplex()
    return SimplicialSet_finite(K).quotient(L)
