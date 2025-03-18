# sage.doctest: needs sage.graphs sage.modules
r"""
Cluster algebras

This file constructs cluster algebras using the Parent-Element framework.
The implementation mainly utilizes structural theorems from [FZ2007]_.

The key points being used here are these:

- cluster variables are parametrized by their g-vectors;

- g-vectors (together with c-vectors) provide a self-standing model for the
  combinatorics behind any cluster algebra;

- each cluster variable in any cluster algebra can be computed, by the
  separation of additions formula, from its g-vector and F-polynomial.

Accordingly this file provides three classes:

- :class:`ClusterAlgebra`

- :class:`ClusterAlgebraSeed`

- :class:`ClusterAlgebraElement`

:class:`ClusterAlgebra`, constructed as a subobject of
:class:`sage.rings.polynomial.laurent_polynomial_ring.LaurentPolynomialRing_generic`,
is the frontend of this implementation. It provides all the algebraic
features (like ring morphisms), it computes cluster variables, it is
responsible for controlling the exploration of the exchange graph and
serves as the repository for all the data recursively computed so far.
In particular, all g-vectors and all F-polynomials of known cluster
variables as well as a mutation path by which they can be obtained
are recorded. In the optic of efficiency, this implementation does not
store directly the exchange graph nor the exchange relations. Both of
these could be added to :class:`ClusterAlgebra` with minimal effort.

:class:`ClusterAlgebraSeed` provides the combinatorial backbone
for :class:`ClusterAlgebra`. It is an auxiliary class and therefore its
instances should **not** be directly created by the user. Rather it
should be accessed via :meth:`ClusterAlgebra.current_seed`
and :meth:`ClusterAlgebra.initial_seed`. The task of performing current
seed mutations is delegated to this class. Seeds are considered equal if
they have the same parent cluster algebra and they can be obtained from
each other by a permutation of their data (i.e. if they coincide as
unlabelled seeds).  Cluster algebras whose initial seeds are equal in the
above sense are not considered equal but are endowed with coercion maps
to each other.  More generally, a cluster algebra is endowed with coercion
maps from any cluster algebra which is obtained by freezing a collection
of initial cluster variables and/or permuting both cluster variables
and coefficients.

:class:`ClusterAlgebraElement` is a thin wrapper around
:class:`sage.rings.polynomial.laurent_polynomial.LaurentPolynomial`
providing all the functions specific to cluster variables.
Elements of a cluster algebra with principal coefficients have special methods
and these are grouped in the subclass :class:`PrincipalClusterAlgebraElement`.

One more remark about this implementation. Instances of
:class:`ClusterAlgebra` are built by identifying the initial cluster variables
with the generators of :meth:`ClusterAlgebra.ambient`. In particular, this
forces a specific embedding into the ambient field of rational expressions. In
view of this, although cluster algebras themselves are independent of the
choice of initial seed, :meth:`ClusterAlgebra.mutate_initial` is forced to
return a different instance of :class:`ClusterAlgebra`. At the moment there
is no coercion implemented among the two instances but this could in
principle be added to :meth:`ClusterAlgebra.mutate_initial`.

REFERENCES:

- [FZ2007]_
- [LLZ2014]_
- [NZ2012]_

AUTHORS:

- Dylan Rupel (2015-06-15): initial version

- Salvatore Stella (2015-06-15): initial version

EXAMPLES:

We begin by creating a simple cluster algebra and printing its
initial exchange matrix::

    sage: A = ClusterAlgebra(['A', 2]); A
    A Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
    sage: A.b_matrix()
    [ 0  1]
    [-1  0]

``A`` is of finite type so we can explore all its exchange graph::

    sage: A.explore_to_depth(infinity)

and get all its g-vectors, F-polynomials, and cluster variables::

    sage: sorted(A.g_vectors_so_far())
    [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0)]
    sage: sorted(A.F_polynomials_so_far(), key=str)
    [1, 1, u0 + 1, u0*u1 + u0 + 1, u1 + 1]
    sage: sorted(A.cluster_variables_so_far(), key=str)
    [(x0 + 1)/x1, (x0 + x1 + 1)/(x0*x1), (x1 + 1)/x0, x0, x1]

Simple operations among cluster variables behave as expected::

    sage: s = A.cluster_variable((0, -1)); s
    (x0 + 1)/x1
    sage: t = A.cluster_variable((-1, 1)); t
    (x1 + 1)/x0
    sage: t + s
    (x0^2 + x1^2 + x0 + x1)/(x0*x1)
    sage: _.parent() == A
    True
    sage: t - s
    (-x0^2 + x1^2 - x0 + x1)/(x0*x1)
    sage: _.parent() == A
    True
    sage: t*s
    (x0*x1 + x0 + x1 + 1)/(x0*x1)
    sage: _.parent() == A
    True
    sage: t/s
    (x1^2 + x1)/(x0^2 + x0)
    sage: _.parent() == A
    False

Division is not guaranteed to yield an element of ``A`` so it returns an
element of ``A.ambient().fraction_field()`` instead::

    sage: (t/s).parent() == A.ambient().fraction_field()
    True

We can compute denominator vectors of any element of ``A``::

    sage: (t*s).d_vector()
    (1, 1)

Since we are in rank 2 and we do not have coefficients we can compute the
greedy element associated to any denominator vector::

    sage: A.rank() == 2 and A.coefficients() == ()
    True
    sage: A.greedy_element((1, 1))
    (x0 + x1 + 1)/(x0*x1)
    sage: _ == t*s
    False

not surprising since there is no cluster in ``A`` containing
both ``t`` and ``s``::

    sage: seeds = A.seeds(mutating_F=false)
    sage: [ S for S in seeds if (0, -1) in S and (-1, 1) in S ]
    []

indeed::

    sage: A.greedy_element((1, 1)) == A.cluster_variable((-1, 0))
    True

Disabling F-polynomials in the computation just done was redundant because we
already explored the whole exchange graph before. Though in different
circumstances it could have saved us considerable time.

g-vectors and F-polynomials can be computed from elements of ``A`` only if
``A`` has principal coefficients at the initial seed::

    sage: (t*s).g_vector()
    Traceback (most recent call last):
    ...
    AttributeError: 'ClusterAlgebra_with_category.element_class' object has no attribute 'g_vector'...
    sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
    sage: A.explore_to_depth(infinity)
    sage: s = A.cluster_variable((0, -1)); s
    (x0*y1 + 1)/x1
    sage: t = A.cluster_variable((-1, 1)); t
    (x1 + y0)/x0
    sage: (t*s).g_vector()
    (-1, 0)
    sage: (t*s).F_polynomial()
    u0*u1 + u0 + u1 + 1
    sage: (t*s).is_homogeneous()
    True
    sage: (t+s).is_homogeneous()
    False
    sage: (t+s).homogeneous_components()
    {(-1, 1): (x1 + y0)/x0, (0, -1): (x0*y1 + 1)/x1}

Each cluster algebra is endowed with a reference to a current seed;
it could be useful to assign a name to it::

    sage: A = ClusterAlgebra(['F', 4])
    sage: len(A.g_vectors_so_far())
    4
    sage: A.current_seed()
    The initial seed of a Cluster Algebra with cluster variables x0, x1, x2, x3
     and no coefficients over Integer Ring
    sage: A.current_seed() == A.initial_seed()
    True
    sage: S = A.current_seed()
    sage: S.b_matrix()
    [ 0  1  0  0]
    [-1  0 -1  0]
    [ 0  2  0  1]
    [ 0  0 -1  0]
    sage: S.g_matrix()
    [1 0 0 0]
    [0 1 0 0]
    [0 0 1 0]
    [0 0 0 1]
    sage: S.cluster_variables()
    [x0, x1, x2, x3]

and use ``S`` to walk around the exchange graph of ``A``::

    sage: S.mutate(0); S
    The seed of a Cluster Algebra with cluster variables x0, x1, x2, x3
     and no coefficients over Integer Ring obtained from the initial
     by mutating in direction 0
    sage: S.b_matrix()
    [ 0 -1  0  0]
    [ 1  0 -1  0]
    [ 0  2  0  1]
    [ 0  0 -1  0]
    sage: S.g_matrix()
    [-1  0  0  0]
    [ 1  1  0  0]
    [ 0  0  1  0]
    [ 0  0  0  1]
    sage: S.cluster_variables()
    [(x1 + 1)/x0, x1, x2, x3]
    sage: S.mutate('sinks'); S
    The seed of a Cluster Algebra with cluster variables x0, x1, x2, x3
     and no coefficients over Integer Ring obtained from the initial
     by mutating along the sequence [0, 2]
    sage: S.mutate([2, 3, 2, 1, 0]); S
    The seed of a Cluster Algebra with cluster variables x0, x1, x2, x3
     and no coefficients over Integer Ring obtained from the initial
     by mutating along the sequence [0, 3, 2, 1, 0]
    sage: S.g_vectors()
    [(0, 1, -2, 0), (-1, 2, -2, 0), (0, 1, -1, 0), (0, 0, 0, -1)]
    sage: S.cluster_variable(3)
    (x2 + 1)/x3

Walking around by mutating ``S`` updates the information stored in ``A``::

    sage: len(A.g_vectors_so_far())
    10
    sage: A.current_seed().path_from_initial_seed()
    [0, 3, 2, 1, 0]
    sage: A.current_seed() == S
    True

Starting from ``A.initial_seed()`` still records data in ``A`` but does not
update ``A.current_seed()``::

    sage: S1 = A.initial_seed()
    sage: S1.mutate([2, 1, 3])
    sage: len(A.g_vectors_so_far())
    11
    sage: S1 == A.current_seed()
    False

Since :class:`ClusterAlgebra` inherits from :class:`UniqueRepresentation`,
computed data is shared across instances::

    sage: A1 = ClusterAlgebra(['F', 4])
    sage: A1 is A
    True
    sage: len(A1.g_vectors_so_far())
    11

It can be useful, at times to forget all computed data. Because of
:class:`UniqueRepresentation` this cannot be achieved by simply creating a
new instance; instead it has to be manually triggered by::

    sage: A.clear_computed_data()
    sage: len(A.g_vectors_so_far())
    4

Given a cluster algebra ``A`` we may be looking for a specific cluster
variable::

    sage: A = ClusterAlgebra(['E', 8, 1])
    sage: v = (-1, 1, -1, 1, -1, 1, 0, 0, 1)
    sage: A.find_g_vector(v, depth=2)
    sage: seq = A.find_g_vector(v); seq  # random
    [0, 1, 2, 4, 3]
    sage: v in A.initial_seed().mutate(seq, inplace=False).g_vectors()
    True

This also performs mutations of F-polynomials::

    sage: A.F_polynomial((-1, 1, -1, 1, -1, 1, 0, 0, 1))
    u0*u1*u2*u3*u4 + u0*u1*u2*u4 + u0*u2*u3*u4 + u0*u1*u2 + u0*u2*u4
     + u2*u3*u4 + u0*u2 + u0*u4 + u2*u4 + u0 + u2 + u4 + 1

which might not be a good idea in algebras that are too big. One workaround is
to first disable F-polynomials and then recompute only the desired mutations::

    sage: # long time
    sage: A.reset_exploring_iterator(mutating_F=False)
    sage: v = (-1, 1, -2, 2, -1, 1, -1, 1, 1)
    sage: seq = A.find_g_vector(v); seq  # random
    [1, 0, 2, 6, 5, 4, 3, 8, 1]
    sage: S = A.initial_seed().mutate(seq, inplace=False)
    sage: v in S.g_vectors()
    True
    sage: A.current_seed().mutate(seq)
    sage: A.F_polynomial((-1, 1, -2, 2, -1, 1, -1, 1, 1))
    u0*u1^2*u2^2*u3*u4*u5*u6*u8 +
    ...
    2*u2 + u4 + u6 + 1

We can manually freeze cluster variables and get coercions in between
the two algebras::

    sage: A = ClusterAlgebra(['F', 4]); A
    A Cluster Algebra with cluster variables x0, x1, x2, x3 and no coefficients
     over Integer Ring
    sage: A1 = ClusterAlgebra(A.b_matrix().matrix_from_columns([0, 1, 2]), coefficient_prefix='x'); A1
    A Cluster Algebra with cluster variables x0, x1, x2 and coefficient x3
     over Integer Ring
    sage: A.has_coerce_map_from(A1)
    True

and we also have an immersion of ``A.base()`` into ``A`` and of ``A``
into ``A.ambient()``::

    sage: A.has_coerce_map_from(A.base())
    True
    sage: A.ambient().has_coerce_map_from(A)
    True

but there is currently no coercion in between algebras obtained by
mutating at the initial seed::

    sage: A1 = A.mutate_initial(0); A1
    A Cluster Algebra with cluster variables x4, x1, x2, x3 and no coefficients
     over Integer Ring
    sage: A.b_matrix() == A1.b_matrix()
    False
    sage: [X.has_coerce_map_from(Y) for X, Y in [(A, A1), (A1, A)]]
    [False, False]
"""

# ****************************************************************************
#       Copyright (C) 2015 Dylan Rupel and Salvatore Stella
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from typing import Any

from sage.arith.misc import binomial
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.rings import Rings
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.combinat.permutation import Permutation
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.special import block_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import (LaurentPolynomialRing_generic,
                                                           LaurentPolynomialRing)
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.geometry.cone', 'Cone')
lazy_import('sage.geometry.fan', 'Fan')

##############################################################################
# Elements of a cluster algebra
##############################################################################


class ClusterAlgebraElement(ElementWrapper):
    """
    An element of a cluster algebra.
    """
    # AdditiveMagmas.Subobjects currently does not implements _add_
    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an element of ``self.parent()``

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.an_element() + A.an_element()
            2*x0
        """
        return self.parent().retract(self.lift() + other.lift())

    def _neg_(self):
        r"""
        Return  the negative of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: -A.an_element()
            -x0
        """
        return self.parent().retract(-self.lift())

    def _div_(self, other):
        r"""
        Return the quotient of ``self`` and ``other``.

        .. WARNING::

            The result of a division is not guaranteed to be inside
            :meth:`parent` therefore this method does not return an
            instance of :class:`ClusterAlgebraElement`.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: x = A.an_element()
            sage: x/x
            1
            sage: _.parent()
            Multivariate Laurent Polynomial Ring in x0, x1, x2, x3
             over Integer Ring
            sage: A.retract(x/x)
            1
            sage: _.parent()
            A Cluster Algebra with cluster variables x0, x1, x2, x3
             and no coefficients over Integer Ring
        """
        return self.lift() / other.lift()

    def d_vector(self) -> tuple:
        r"""
        Return the denominator vector of ``self`` as a tuple of integers.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4], principal_coefficients=True)
            sage: A.current_seed().mutate([0, 2, 1])
            sage: x = A.cluster_variable((-1, 2, -2, 2)) * A.cluster_variable((0, 0, 0, 1))**2
            sage: x.d_vector()
            (1, 1, 2, -2)
        """
        monomials = self.lift().monomial_coefficients()
        minimal = map(min, zip(*monomials))
        return tuple(-vector(minimal))[:self.parent().rank()]

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4], principal_coefficients=True)
            sage: A.current_seed().mutate([0, 2, 1])
            sage: A.cluster_variable((-1, 2, -2, 2))
            (x0*x2^2*y0*y1*y2^2 + x1^3*x3^2 + x1^2*x3^2*y0 + 2*x1^2*x3*y2 + 2*x1*x3*y0*y2 + x1*y2^2 + y0*y2^2)/(x0*x1*x2^2)
        """
        numer, denom = self.lift()._fraction_pair()
        return repr(numer / denom)


class PrincipalClusterAlgebraElement(ClusterAlgebraElement):
    """
    An element in a cluster algebra with principle coefficients.
    """
    def g_vector(self):
        r"""
        Return the g-vector of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: A.cluster_variable((1, 0)).g_vector() == (1, 0)
            True
            sage: sum(A.initial_cluster_variables()).g_vector()
            Traceback (most recent call last):
            ...
            ValueError: this element does not have a well defined g-vector
        """
        if not self.is_homogeneous():
            raise ValueError("this element does not have a well defined g-vector")
        return self._g_vector

    def F_polynomial(self):
        r"""
        Return the F-polynomial of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: S = A.initial_seed()
            sage: S.mutate([0, 1, 0])
            sage: S.cluster_variable(0).F_polynomial() == S.F_polynomial(0)
            True
            sage: sum(A.initial_cluster_variables()).F_polynomial()
            Traceback (most recent call last):
            ...
            ValueError: this element does not have a well defined g-vector
        """
        if not self.is_homogeneous():
            raise ValueError("this element does not have a well defined g-vector")
        subs_dict = {}
        A = self.parent()
        for x in A.initial_cluster_variables():
            subs_dict[x.lift()] = A._U(1)
        for i in range(A.rank()):
            subs_dict[A.coefficient(i).lift()] = A._U.gen(i)
        return A._U(self.lift().substitute(subs_dict))

    def is_homogeneous(self) -> bool:
        r"""
        Return ``True`` if ``self`` is a homogeneous element
        of ``self.parent()``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: A.cluster_variable((1, 0)).is_homogeneous()
            True
            sage: x = A.cluster_variable((1, 0)) + A.cluster_variable((0, 1))
            sage: x.is_homogeneous()
            False
        """
        return getattr(self, '_is_homogeneous', len(self.homogeneous_components()) == 1)

    def homogeneous_components(self) -> dict:
        r"""
        Return a dictionary of the homogeneous components of ``self``.

        OUTPUT:

        A dictionary whose keys are homogeneous degrees and whose values
        are the summands of ``self`` of the given degree.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: x = A.cluster_variable((1, 0)) + A.cluster_variable((0, 1))
            sage: x.homogeneous_components()
            {(0, 1): x1, (1, 0): x0}
        """
        deg_matrix = block_matrix([[identity_matrix(self.parent().rank()),
                                    -self.parent().b_matrix()]])
        components: dict[tuple, Any] = {}
        x = self.lift()
        monomials = x.monomials()
        for m in monomials:
            g_vect = tuple(deg_matrix * vector(m.exponents()[0]))
            if g_vect in components:
                components[g_vect] += self.parent().retract(x.monomial_coefficient(m) * m)
            else:
                components[g_vect] = self.parent().retract(x.monomial_coefficient(m) * m)
        for g_vect, compo in components.items():
            compo._is_homogeneous = True
            compo._g_vector = g_vect
        self._is_homogeneous = (len(components) == 1)
        if self._is_homogeneous:
            self._g_vector = next(iter(components))
        return components

    def theta_basis_decomposition(self):
        r"""
        Return the decomposition of ``self`` in the theta basis.

        OUTPUT:

        A dictionary whose keys are the g-vectors and whose values are the coefficients
        in the decomposition of ``self`` in the theta basis.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,-2],[3,0]]), principal_coefficients=True)
            sage: f = (A.theta_basis_element((1,0)) + A.theta_basis_element((0,1)))**2 + A.coefficient(1)* A.theta_basis_element((1,1))
            sage: decomposition = f.theta_basis_decomposition()
            sage: sum(decomposition[g] * A.theta_basis_element(g) for g in decomposition) == f
            True
            sage: f = A.theta_basis_element((4,-4))*A.theta_basis_element((1,-1))
            sage: decomposition =  f.theta_basis_decomposition()
            sage: sum(decomposition[g] * A.theta_basis_element(g) for g in decomposition) == f
            True
        """
        A = self.parent()
        B = A.b_matrix()
        U = A._U
        out = {}
        zero_A = A(0)
        zero_U = U(0)
        zero_t = (0,) * A.rank()

        components = self.homogeneous_components()

        for g_vect in components:
            f_poly = components[g_vect].F_polynomial()
            g_vect = vector(g_vect)
            while f_poly != zero_U:
                coeffs = f_poly.monomial_coefficients()
                y_exp = min(coeffs)
                coeff = coeffs[y_exp]
                g_theta = tuple(g_vect + B * vector(y_exp))
                out[g_theta] = out.get(g_theta, zero_A) + A({zero_t + tuple(y_exp): coeff})
                f_poly -= U({y_exp: coeff}) * A.theta_basis_F_polynomial(g_theta)

        return out


##############################################################################
# Seeds
##############################################################################

class ClusterAlgebraSeed(SageObject):
    """
    A seed in a Cluster Algebra.

    INPUT:

    - ``B`` -- a skew-symmetrizable integer matrix
    - ``C`` -- the matrix of c-vectors of ``self``
    - ``G`` -- the matrix of g-vectors of ``self``
    - ``parent`` -- :class:`ClusterAlgebra`; the algebra to which the
      seed belongs
    - ``path`` -- list (default: ``[]``); the mutation sequence from the
      initial seed of ``parent`` to ``self``

    .. WARNING::

        Seeds should **not** be created manually: no test is performed to
        assert that they are built from consistent data nor that they
        really are seeds of ``parent``. If you create seeds with
        inconsistent data all sort of things can go wrong, even
        :meth:`__eq__` is no longer guaranteed to give correct answers.
        Use at your own risk.
    """
    def __init__(self, B, C, G, parent, **kwargs):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: from sage.algebras.cluster_algebra import ClusterAlgebraSeed
            sage: ClusterAlgebraSeed(A.b_matrix(), identity_matrix(4), identity_matrix(4), A, path=[1, 2, 3])
            The seed of a Cluster Algebra with cluster variables x0, x1, x2, x3
             and no coefficients over Integer Ring obtained from the initial
             by mutating along the sequence [1, 2, 3]
        """
        self._B = copy(B)
        self._C = copy(C)
        self._G = copy(G)
        self._parent = parent
        self._path = kwargs.get('path', [])

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = copy(A.current_seed())
            sage: S == A.current_seed()
            True
            sage: S is not A.current_seed()
            True
        """
        other = type(self).__new__(type(self))
        other._B = copy(self._B)
        other._C = copy(self._C)
        other._G = copy(self._G)
        other._parent = self._parent
        other._path = copy(self._path)
        return other

    def __eq__(self, other):
        r"""
        Test equality of two seeds.

        INPUT:

        - ``other`` -- a :class:`ClusterAlgebraSeed`

        ALGORITHM:

        ``self`` and ``other`` are deemed to be equal if they have the same
        parent and their set of g-vectors coincide, i.e. this tests
        equality of unlabelled seeds.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: A.clear_computed_data()
            sage: S = copy(A.current_seed())
            sage: S.mutate([0, 2, 0])
            sage: S == A.current_seed()
            False
            sage: S.mutate(2)
            sage: S == A.current_seed()
            True

            sage: A = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: S = A.current_seed()
            sage: S.mutate(0)
            sage: S == A.current_seed()
            True
        """
        return (isinstance(other, ClusterAlgebraSeed) and
                self.parent() == other.parent() and
                frozenset(self.g_vectors()) == frozenset(other.g_vectors()))

    def __contains__(self, element) -> bool:
        r"""
        Test whether ``element`` belong to ``self``.

        INPUT:

        - ``element`` -- either a g-vector or an element of :meth:`parent`

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: (1, 0, 0) in S
            True
            sage: (1, 1, 0) in S
            False
            sage: A.cluster_variable((1, 0, 0)) in S
            True
        """
        if isinstance(element, ClusterAlgebraElement):
            cluster = self.cluster_variables()
        else:
            element = tuple(element)
            cluster = self.g_vectors()
        return element in cluster

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        ALGORITHM:

        For speed purposes the hash is computed on :meth:`self.g_vectors`.
        In particular it is guaranteed to be unique only within a given
        instance of :class:`ClusterAlgebra`. Moreover unlabelled seeds that
        have the same set of g-vectors have the same hash.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: T = S.mutate(1, inplace=False)
            sage: hash(S) == hash(S)
            True
            sage: hash(S) == hash(T)
            False
        """
        return hash(frozenset(self.g_vectors()))

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: A.clear_computed_data()
            sage: S = A.current_seed(); S
            The initial seed of a Cluster Algebra with cluster variables x0, x1, x2
             and no coefficients over Integer Ring
            sage: S.mutate(0); S
            The seed of a Cluster Algebra with cluster variables x0, x1, x2
             and no coefficients over Integer Ring obtained from the initial
             by mutating in direction 0
            sage: S.mutate(1); S
            The seed of a Cluster Algebra with cluster variables x0, x1, x2
             and no coefficients over Integer Ring obtained from the initial
             by mutating along the sequence [0, 1]
        """
        if not self._path:
            return "The initial seed of a %s" % str(self.parent())[2:]
        if len(self._path) == 1:
            return "The seed of a %s obtained from the initial by mutating in direction %s" % (str(self.parent())[2:], str(self._path[0]))
        return "The seed of a %s obtained from the initial by mutating along the sequence %s" % (str(self.parent())[2:], str(self._path))

    def parent(self):
        r"""
        Return the parent of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B', 3])
            sage: A.current_seed().parent() == A
            True
        """
        return self._parent

    def depth(self) -> int:
        r"""
        Return the length of a mutation sequence from the initial seed
         of :meth:`parent` to ``self``.

        .. WARNING::

            This is the length of the mutation sequence returned by
            :meth:`path_from_initial_seed`, which need not be the
            shortest possible.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: S1 = A.initial_seed()
            sage: S1.mutate([0, 1, 0, 1])
            sage: S1.depth()
            4
            sage: S2 = A.initial_seed()
            sage: S2.mutate(1)
            sage: S2.depth()
            1
            sage: S1 == S2
            True
        """
        return len(self._path)

    def path_from_initial_seed(self) -> list:
        r"""
        Return a mutation sequence from the initial seed of :meth:`parent`
        to ``self``.

        .. WARNING::

            This is the path used to compute ``self`` and it does not
            have to be the shortest possible.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: S1 = A.initial_seed()
            sage: S1.mutate([0, 1, 0, 1])
            sage: S1.path_from_initial_seed()
            [0, 1, 0, 1]
            sage: S2 = A.initial_seed()
            sage: S2.mutate(1)
            sage: S2.path_from_initial_seed()
            [1]
            sage: S1 == S2
            True
        """
        return copy(self._path)

    def b_matrix(self):
        r"""
        Return the exchange matrix of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
        """
        return copy(self._B)

    def c_matrix(self):
        r"""
        Return the matrix whose columns are the c-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.c_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return copy(self._C)

    def c_vector(self, j) -> tuple:
        r"""
        Return the ``j``-th c-vector of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the c-vector to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.c_vector(0)
            (1, 0, 0)
            sage: S.mutate(0)
            sage: S.c_vector(0)
            (-1, 0, 0)
            sage: S.c_vector(1)
            (1, 1, 0)
        """
        return tuple(self._C.column(j))

    def c_vectors(self) -> list[tuple]:
        r"""
        Return all the c-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.c_vectors()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return list(map(tuple, self._C.columns()))

    def g_matrix(self):
        r"""
        Return the matrix whose columns are the g-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.g_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return copy(self._G)

    def g_vector(self, j) -> tuple:
        r"""
        Return the ``j``-th g-vector of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the g-vector to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.g_vector(0)
            (1, 0, 0)
        """
        return tuple(self._G.column(j))

    def g_vectors(self) -> list[tuple]:
        r"""
        Return all the g-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.g_vectors()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return list(map(tuple, self._G.columns()))

    def F_polynomial(self, j):
        r"""
        Return the ``j``-th F-polynomial of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the F-polynomial to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.F_polynomial(0)
            1
        """
        return self.parent().F_polynomial(self.g_vector(j))

    def F_polynomials(self) -> list:
        r"""
        Return all the F-polynomials of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.F_polynomials()
            [1, 1, 1]
        """
        return [self.parent().F_polynomial(g) for g in self.g_vectors()]

    def cluster_variable(self, j):
        r"""
        Return the ``j``-th cluster variable of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the cluster variable to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.cluster_variable(0)
            x0
            sage: S.mutate(0)
            sage: S.cluster_variable(0)
            (x1 + 1)/x0
        """
        return self.parent().cluster_variable(self.g_vector(j))

    def cluster_variables(self) -> list:
        r"""
        Return all the cluster variables of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: S = A.initial_seed()
            sage: S.cluster_variables()
            [x0, x1, x2]
        """
        return [self.parent().cluster_variable(g) for g in self.g_vectors()]

    def mutate(self, direction, **kwargs):
        r"""
        Mutate ``self``.

        INPUT:

        - ``direction`` -- in which direction(s) to mutate, it can be:

          * an integer in ``range(self.rank())`` to mutate in one direction only
          * an iterable of such integers to mutate along a sequence
          * a string "sinks" or "sources" to mutate at all sinks or sources simultaneously

        - ``inplace`` -- boolean (default: ``True``); whether to mutate in place
          or to return a new object

        - ``mutating_F`` -- boolean (default: ``True``); whether to compute
          F-polynomials while mutating

        .. NOTE::

            While knowing F-polynomials is essential to computing
            cluster variables, the process of mutating them is quite slow.
            If you care only about combinatorial data like g-vectors and
            c-vectors, setting ``mutating_F=False`` yields significant
            benefits in terms of speed.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: S = A.initial_seed()
            sage: S.mutate(0); S
            The seed of a Cluster Algebra with cluster variables x0, x1
             and no coefficients over Integer Ring obtained from the initial
             by mutating in direction 0
            sage: S.mutate(5)
            Traceback (most recent call last):
            ...
            ValueError: cannot mutate in direction 5
        """
        n = self.parent().rank()

        # do we want to change self?
        inplace = kwargs.pop('inplace', True)
        if inplace:
            to_mutate = self
        else:
            to_mutate = copy(self)

        # construct mutation sequence
        # if you change this be considerate and change also :class:`ClusterAlgebra`.mutate_initial
        if direction == "sinks":
            B = self.b_matrix()
            seq = [i for i in range(n) if all(x <= 0 for x in B.column(i))]
        elif direction == "sources":
            B = self.b_matrix()
            seq = [i for i in range(n) if all(x >= 0 for x in B.column(i))]
        else:
            try:
                seq = iter(direction)
            except TypeError:
                seq = iter((direction, ))

        # are we mutating F-polynomials?
        mutating_F = kwargs.pop('mutating_F', True)

        for k in seq:
            if k not in range(n):
                raise ValueError(f'cannot mutate in direction {k}')

            # store new mutation path
            if to_mutate._path and to_mutate._path[-1] == k:
                to_mutate._path.pop()
            else:
                to_mutate._path.append(k)

            # find sign of k-th c-vector
            if any(x > 0 for x in to_mutate._C.column(k)):
                eps = +1
            else:
                eps = -1

            # store the g-vector to be mutated in case we are mutating F-polynomials also
            old_g_vector = to_mutate.g_vector(k)

            # compute new G-matrix
            J = identity_matrix(n)
            for j in range(n):
                J[j, k] += max(0, -eps * to_mutate._B[j, k])
            J[k, k] = -1
            to_mutate._G = to_mutate._G * J

            # path to new g-vector (we store the shortest encountered so far)
            g_vector = to_mutate.g_vector(k)
            if g_vector not in to_mutate.parent()._path_dict or len(to_mutate.parent()._path_dict[g_vector]) > len(to_mutate._path):
                to_mutate.parent()._path_dict[g_vector] = copy(to_mutate._path)

            # compute F-polynomials
            if mutating_F and g_vector not in to_mutate.parent()._F_poly_dict:
                to_mutate.parent()._F_poly_dict[g_vector] = to_mutate._mutated_F(k, old_g_vector)

            # compute new C-matrix
            J = identity_matrix(n)
            for j in range(n):
                J[k, j] += max(0, eps * to_mutate._B[k, j])
            J[k, k] = -1
            to_mutate._C = to_mutate._C * J

            # compute new B-matrix
            to_mutate._B.mutate(k)

        # return if we need to
        if not inplace:
            return to_mutate

    def _mutated_F(self, k, old_g_vector):
        r"""
        Compute new F-polynomial obtained by mutating in direction ``k``.

        INPUT:

        - ``k`` -- an integer in ``range(self.parent().rank())``;
          the direction in which we are mutating

        - ``old_g_vector`` -- tuple; the `k`-th g-vector of ``self``
          before mutating

        .. NOTE::

            This function is the bottleneck of :meth:`mutate`. The problem is
            that operations on polynomials are slow. One can get a significant
            speed boost by disabling this method calling :meth:`mutate` with
            ``mutating_F=False``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: S = A.initial_seed()
            sage: S.mutate(0)
            sage: S._mutated_F(0, (1, 0))
            u0 + 1

        Check that :issue:`28176` is fixed::

            sage: A = ClusterAlgebra(matrix([[0,2],[-2,0]]))
            sage: S = A.initial_seed()
            sage: S.mutate([1, 0, 1])
            sage: parent(S._mutated_F(1, (0, -1)))
            Multivariate Polynomial Ring in u0, u1 over Rational Field
        """
        alg = self.parent()
        pos = alg._U(1)
        neg = alg._U(1)
        for j in range(alg.rank()):
            if self._C[j, k] > 0:
                pos *= alg._U.gen(j) ** self._C[j, k]
            else:
                neg *= alg._U.gen(j) ** (-self._C[j, k])
            if self._B[j, k] > 0:
                pos *= self.F_polynomial(j) ** self._B[j, k]
            elif self._B[j, k] < 0:
                neg *= self.F_polynomial(j) ** (-self._B[j, k])
        return (pos + neg) // alg.F_polynomial(old_g_vector)

##############################################################################
# Cluster algebras
##############################################################################


class ClusterAlgebra(Parent, UniqueRepresentation):
    r"""
    A Cluster Algebra.

    INPUT:

    - ``data`` -- some data defining a cluster algebra; it can be anything
      that can be parsed by :class:`ClusterQuiver`

    - ``scalars`` -- a ring (default: `\ZZ`); the scalars over
      which the cluster algebra is defined

    - ``cluster_variable_prefix`` -- string (default: ``'x'``); it needs to be
      a valid variable name

    - ``cluster_variable_names`` -- list of strings; each element needs
      to be a valid variable name;  supersedes ``cluster_variable_prefix``

    - ``coefficient_prefix`` -- string (default: ``'y'``); it needs to be
      a valid variable name

    - ``coefficient_names`` -- list of strings; each element needs
      to be a valid variable name; supersedes ``cluster_variable_prefix``

    - ``principal_coefficients`` -- boolean (default: ``False``); supersedes any
      coefficient defined by ``data``

    ALGORITHM:

    The implementation is mainly based on [FZ2007]_ and [NZ2012]_.

    EXAMPLES::

        sage: B = matrix([(0, 1, 0, 0), (-1, 0, -1, 0), (0, 1, 0, 1), (0, 0, -2, 0), (-1, 0, 0, 0), (0, -1, 0, 0)])
        sage: A = ClusterAlgebra(B); A
        A Cluster Algebra with cluster variables x0, x1, x2, x3
         and coefficients y0, y1 over Integer Ring
        sage: A.gens()
        (x0, x1, x2, x3, y0, y1)
        sage: A = ClusterAlgebra(['A', 2]); A
        A Cluster Algebra with cluster variables x0, x1 and no coefficients
         over Integer Ring
        sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True); A.gens()
        (x0, x1, y0, y1)
        sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True, coefficient_prefix='x'); A.gens()
        (x0, x1, x2, x3)
        sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, cluster_variable_names=['a', 'b', 'c']); A.gens()
        (a, b, c, y0, y1, y2)
        sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, cluster_variable_names=['a', 'b'])
        Traceback (most recent call last):
        ...
        ValueError: cluster_variable_names should be an iterable of 3 valid variable names
        sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, coefficient_names=['a', 'b', 'c']); A.gens()
        (x0, x1, x2, a, b, c)
        sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, coefficient_names=['a', 'b'])
        Traceback (most recent call last):
        ...
        ValueError: coefficient_names should be an iterable of 3 valid variable names
    """

    @staticmethod
    def __classcall__(self, data, **kwargs):
        r"""
        Preparse input.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2]); A   # indirect doctest
            A Cluster Algebra with cluster variables x0, x1 and no coefficients
            over Integer Ring

        Check that :issue:`28176` is fixed::

            sage: A1 = ClusterAlgebra(['A',2])
            sage: A2 = ClusterAlgebra(['A',2], cluster_variable_prefix='x')
            sage: A1 is A2
            True
            sage: A3 = ClusterAlgebra(Matrix([[0,1],[-1,0]]))
            sage: A1 is A3
            True
            sage: A4 = ClusterAlgebra([[0,1]]) # built from a digraph
            sage: A1 is A4
            True
        """
        # Use ClusterQuiver to parse the input; eventually we may want to avoid this
        Q = ClusterQuiver(data)

        # Rank
        n = Q.n()

        # Exchange matrix
        B0 = Q.b_matrix()[:n, :]

        # Coefficient matrix
        if kwargs.pop('principal_coefficients', False):
            M0 = identity_matrix(n)
        else:
            M0 = Q.b_matrix()[n:, :]
        m = M0.nrows()

        B0 = block_matrix([[B0], [M0]])
        B0.set_immutable()

        # Determine the names of the initial cluster variables
        kwargs.setdefault('cluster_variable_prefix', 'x')
        kwargs['cluster_variable_names'] = tuple(kwargs.get('cluster_variable_names',
                [kwargs['cluster_variable_prefix'] + str(i) for i in range(n)]))
        if len(kwargs['cluster_variable_names']) != n:
            raise ValueError("cluster_variable_names should be an iterable of %d valid variable names" % n)

        # Determine the names of the coefficients
        coefficient_prefix = kwargs.pop('coefficient_prefix', 'y')
        offset = n if coefficient_prefix == kwargs['cluster_variable_prefix'] else 0
        kwargs['coefficient_names'] = tuple(kwargs.get('coefficient_names',
                [coefficient_prefix + str(i) for i in range(offset, m + offset)]))
        if len(kwargs['coefficient_names']) != m:
            raise ValueError("coefficient_names should be an iterable of %d valid variable names" % m)

        # Compute the next free index for new named variables
        # This is the first integer nfi such that for any j >= nfi
        # kwargs['cluster_variable_prefix']+str(j) is not the name of an
        # initial cluster variable nor a coefficient. This will be used in
        # mutate_initial to name new cluster variables.
        splitnames = (w.partition(kwargs['cluster_variable_prefix'])
                      for w in
                      kwargs['cluster_variable_names'] + kwargs['coefficient_names'])
        nfi = 1 + max((int(v) for u, _, v in splitnames
                       if u == '' and v.isdigit()), default=-1)
        kwargs.setdefault('next_free_index', nfi)

        # Determine scalars
        kwargs.setdefault('scalars', ZZ)

        return super().__classcall__(self, B0, **kwargs)

    def __init__(self, B, **kwargs):
        """
        Initialize ``self``.

        TESTS::

            sage: B = matrix([(0, 1, 0, 0), (-1, 0, -1, 0), (0, 1, 0, 1), (0, 0, -2, 0), (-1, 0, 0, 0), (0, -1, 0, 0)])
            sage: A = ClusterAlgebra(B)
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
            sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True)
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True, coefficient_prefix='x')
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
            sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, cluster_variable_names=['a','b','c'])
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
            sage: A = ClusterAlgebra(['A', 3], principal_coefficients=True, coefficient_names=['a','b','c'])
            sage: A.clear_computed_data()
            sage: TestSuite(A).run()
        """
        # Exchange matrix
        self._B0 = copy(B)

        # Rank
        self._n = B.ncols()

        M0 = B[self._n:, :]
        m = M0.nrows()

        # Ambient space for F-polynomials
        # NOTE: for speed purposes we need to have QQ here instead of the more
        # natural ZZ. The reason is that _mutated_F is faster if we do not cast
        # the result to polynomials but then we get "rational" coefficients
        self._U = PolynomialRing(QQ, ['u%s' % i for i in range(self._n)])

        # Setup infrastructure to store computed data
        self.clear_computed_data()

        # Data to build new named variables
        self._cluster_variable_prefix = kwargs['cluster_variable_prefix']
        self._next_free_index = kwargs['next_free_index']

        # Base ring
        base = LaurentPolynomialRing(kwargs['scalars'], kwargs['coefficient_names']) if m > 0 else kwargs['scalars']

        # Have we got principal coefficients?
        self.Element = PrincipalClusterAlgebraElement if M0 == identity_matrix(self._n) else ClusterAlgebraElement

        # Setup Parent and ambient
        names = kwargs['cluster_variable_names'] + kwargs['coefficient_names']
        self._ambient = LaurentPolynomialRing(kwargs['scalars'], names)
        Parent.__init__(self, base=base, category=Rings(kwargs['scalars']).Commutative().Subobjects(), names=names)

        # Data to compute cluster variables using separation of additions
        # NOTE: storing both _B0 as rectangular matrix and _yhat is redundant.
        # We keep both around for speed purposes.
        self._y = {self._U.gen(j): prod(self._base.gen(i) ** M0[i, j] for i in range(m))
                   for j in range(self._n)}
        self._yhat = {self._U.gen(j): prod(self._ambient.gen(i) ** self._B0[i, j]
                                           for i in range(self._n + m))
                      for j in range(self._n)}

        # Register embedding into self.ambient()
        embedding = SetMorphism(Hom(self, self.ambient()), lambda x: x.lift())
        self._populate_coercion_lists_(embedding=embedding)

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix(1), principal_coefficients=True); A
            A Cluster Algebra with cluster variable x0
             and coefficient y0 over Integer Ring
            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1
             and coefficients y0, y1 over Integer Ring
        """
        var_names = self.initial_cluster_variable_names()
        var_names_str = (" " if len(var_names) == 1 else "s ") + ", ".join(var_names)
        coeff_names = self.coefficient_names()
        coeff_prefix = " and" + (" " if len(coeff_names) > 0 else " no ") + "coefficient"
        coeff = coeff_prefix + (" " if len(coeff_names) == 1 else "s ") + ", ".join(coeff_names) + (" " if len(coeff_names) > 0 else "")
        return "A Cluster Algebra with cluster variable" + var_names_str + coeff + "over " + repr(self.scalars())

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.an_element()
            x0
        """
        return self.initial_cluster_variable(0)

    def _coerce_map_from_(self, other):
        r"""
        Test whether there is a coercion from ``other`` to ``self``.

        ALGORITHM:

        If ``other`` is an instance of :class:`ClusterAlgebra` then allow
        coercion if ``other.ambient()`` can be coerced into ``self.ambient()``
        and ``other`` can be obtained from ``self`` by permuting variables
        and coefficients and/or freezing some initial cluster variables.

        Otherwise allow anything that coerces into ``self.base()`` to coerce
        into ``self``.

        EXAMPLES::

            sage: B1 = matrix([(0, 1, 0, 0), (-1, 0, -1, 0), (0, 1, 0, 1), (0, 0, -2, 0), (-1, 0, 0, 0), (0, -1, 0, 0)])
            sage: B2 = B1.matrix_from_columns([0, 1, 2])
            sage: A1 = ClusterAlgebra(B1, coefficient_prefix='x')
            sage: A2 = ClusterAlgebra(B2, coefficient_prefix='x')
            sage: A1.has_coerce_map_from(A2)
            True
            sage: A2.has_coerce_map_from(A1)
            False
            sage: f = A1.coerce_map_from(A2)
            sage: seq = A2.find_g_vector((-1, 1, -1)); seq  # random
            [0, 2, 1]
            sage: S = A1.initial_seed(); S.mutate(seq)
            sage: S.cluster_variable(seq[-1]) == f(A2.cluster_variable((-1, 1, -1)))
            True
            sage: B3 = B1.matrix_from_columns([1, 2, 3]); B3
            [ 1  0  0]
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -2  0]
            [ 0  0  0]
            [-1  0  0]
            sage: G = PermutationGroup(['(1, 2, 3, 4)'])
            sage: B3.permute_rows(G.gen(0)); B3
            [ 0 -1  0]
            [ 1  0  1]
            [ 0 -2  0]
            [ 1  0  0]
            [ 0  0  0]
            [-1  0  0]
            sage: A3 = ClusterAlgebra(B3, cluster_variable_names=['x1', 'x2', 'x3'], coefficient_names=['x0', 'x4', 'x5'])
            sage: A1.has_coerce_map_from(A3)
            True
            sage: g = A1.coerce_map_from(A3)
            sage: seq1 = A3.find_g_vector((1, -2, 2))
            sage: seq2 = [G.gen(0)(x + 1) - 1 for x in seq1 ]
            sage: S = A1.initial_seed(); S.mutate(seq2)
            sage: S.cluster_variable(seq2[-1]) == g(A3.cluster_variable((1, -2, 2)))
            True

         Check that :issue:`23654` is fixed::

            sage: A = ClusterAlgebra(['A',2])
            sage: AA = ClusterAlgebra(['A',3])
            sage: A.has_coerce_map_from(AA)
            False
        """
        if isinstance(other, ClusterAlgebra):
            gen_s = self.gens()
            gen_o = other.gens()
            if len(gen_s) == len(gen_o):
                f = self.ambient().coerce_map_from(other.ambient())
                if f is not None:
                    perm = Permutation([gen_s.index(self(f(v))) + 1 for v in gen_o])
                    n = self.rank()
                    M0 = self._B0[n:, :]
                    m = M0.nrows()
                    B = block_matrix([[self.b_matrix(), -M0.transpose()], [M0, matrix(m)]])
                    B.permute_rows_and_columns(perm, perm)
                    return B[:, :other.rank()] == other._B0

        # everything that is in the base can be coerced to self
        return self.base().has_coerce_map_from(other)

    @cached_method
    def coxeter_element(self):
        r"""
        Return the Coxeter element associated to the initial exchange matrix, if acyclic.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,1,1],[-1,0,1],[-1,-1,0]]))
            sage: A.coxeter_element()
            [0, 1, 2]

        Raise an error if the initial exchange matrix is not acyclic::

            sage: A = ClusterAlgebra(matrix([[0,1,-1],[-1,0,1],[1,-1,0]]))
            sage: A.coxeter_element()
            Traceback (most recent call last):
            ...
            ValueError: the initial exchange matrix is not acyclic
        """
        dg = DiGraph(self.b_matrix().apply_map(lambda x: ZZ.zero() if x <= 0 else ZZ.one()))
        acyclic, coxeter = dg.is_directed_acyclic(certificate=True)
        if not acyclic:
            raise ValueError("the initial exchange matrix is not acyclic")
        return coxeter

    @cached_method
    def is_acyclic(self) -> bool:
        r"""
        Return ``True`` if the exchange matrix in the initial seed is acyclic, ``False`` otherwise.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,1,1],[-1,0,1],[-1,-1,0]]))
            sage: A.is_acyclic()
            True
            sage: A = ClusterAlgebra(matrix([[0,1,-1],[-1,0,1],[1,-1,0]]))
            sage: A.is_acyclic()
            False
        """
        dg = DiGraph(self.b_matrix().apply_map(lambda x: ZZ.zero() if x <= 0 else ZZ.one()))
        return dg.is_directed_acyclic()

    def rank(self):
        r"""
        Return the rank of ``self``, i.e. the number of cluster variables
        in any seed.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1
             and coefficients y0, y1 over Integer Ring
            sage: A.rank()
            2
        """
        return self._n

    def current_seed(self):
        r"""
        Return the current seed of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.current_seed()
            The initial seed of a Cluster Algebra with cluster variables x0, x1
             and no coefficients over Integer Ring
        """
        return self._seed

    def set_current_seed(self, seed):
        r"""
        Set the value reported by :meth:`current_seed` to ``seed``,
        if it makes sense.

        INPUT:

        - ``seed`` -- a :class:`ClusterAlgebraSeed`

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: S = copy(A.current_seed())
            sage: S.mutate([0, 1, 0])
            sage: A.current_seed() == S
            False
            sage: A.set_current_seed(S)
            sage: A.current_seed() == S
            True
            sage: A1 = ClusterAlgebra(['B', 2])
            sage: A.set_current_seed(A1.initial_seed())
            Traceback (most recent call last):
            ...
            ValueError: this is not a seed in this cluster algebra
        """
        if self.contains_seed(seed):
            self._seed = seed
        else:
            raise ValueError("this is not a seed in this cluster algebra")

    def reset_current_seed(self):
        r"""
        Reset the value reported by :meth:`current_seed`
        to :meth:`initial_seed`.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.current_seed().mutate([1, 0])
            sage: A.current_seed() == A.initial_seed()
            False
            sage: A.reset_current_seed()
            sage: A.current_seed() == A.initial_seed()
            True
        """
        self._seed = self.initial_seed()

    def clear_computed_data(self):
        r"""
        Clear the cache of computed g-vectors and F-polynomials
        and reset both the current seed and the exploring iterator.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: sorted(A.g_vectors_so_far())
            [(0, 1), (1, 0)]
            sage: A.current_seed().mutate([1, 0])
            sage: sorted(A.g_vectors_so_far())
            [(-1, 0), (0, -1), (0, 1), (1, 0)]
            sage: A.clear_computed_data()
            sage: sorted(A.g_vectors_so_far())
            [(0, 1), (1, 0)]
        """
        I = identity_matrix(self._n)
        self._path_dict = {v: [] for v in map(tuple, I.columns())}
        self._F_poly_dict = {v: self._U(1) for v in self._path_dict}
        self.reset_current_seed()
        self.reset_exploring_iterator()

    def contains_seed(self, seed) -> bool:
        r"""
        Test if ``seed`` is a seed of ``self``.

        INPUT:

        - ``seed`` -- a :class:`ClusterAlgebraSeed`

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
            sage: S = copy(A.current_seed())
            sage: A.contains_seed(S)
            True
        """
        computed_sd = self.initial_seed()
        computed_sd.mutate(seed._path, mutating_F=False)
        return computed_sd == seed

    def initial_seed(self):
        r"""
        Return the initial seed of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.initial_seed()
            The initial seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
        """
        n = self.rank()
        I = identity_matrix(n)
        return ClusterAlgebraSeed(self.b_matrix(), I, I, self)

    def b_matrix(self):
        r"""
        Return the initial exchange matrix of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.b_matrix()
            [ 0  1]
            [-1  0]
        """
        n = self.rank()
        return copy(self._B0[:n, :])

    def euler_matrix(self):
        r"""
        Return the Euler matrix associated to ``self``.

        ALGORITHM:

            This method returns the matrix of the bilinear form defined in Equation (2.1) of [ReSt2020]_ .

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,1,1],[-1,0,1],[-1,-1,0]]))
            sage: A.euler_matrix()
            [ 1  0  0]
            [-1  1  0]
            [-1 -1  1]

        Raise an error if the initial exchange matrix is not acyclic::

            sage: A = ClusterAlgebra(matrix([[0,1,-1],[-1,0,1],[1,-1,0]]))
            sage: A.euler_matrix()
            Traceback (most recent call last):
            ...
            ValueError: the initial exchange matrix is not acyclic
        """
        if not self.is_acyclic():
            raise ValueError("the initial exchange matrix is not acyclic")
        return 1 + self.b_matrix().apply_map(lambda x: min(ZZ.zero(), x))

    def d_vector_to_g_vector(self, d) -> tuple:
        r"""
        Return the g-vector of an element of ``self`` having d-vector ``d``

        INPUT:

        - ``d`` -- the d-vector

        ALGORITHM:

            This method implements the piecewise-linear map `\\nu_c` introduced in Section 9.1 of [ReSt2020]_.

        .. WARNING:

            This implementation works only when the initial exchange matrix is acyclic.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,1,1],[-1,0,1],[-1,-1,0]]))
            sage: A.d_vector_to_g_vector((1,0,-1))
            (-1, 1, 2)
        """
        dm = vector(x if x < 0 else 0 for x in d)
        dp = vector(d) - dm
        return tuple(- dm - self.euler_matrix() * dp)

    def g_vector_to_d_vector(self, g) -> tuple:
        r"""
        Return the d-vector of an element of ``self`` having g-vector ``g``

        INPUT:

        - ``g`` -- the g-vector

        ALGORITHM:

            This method implements the inverse of the piecewise-linear map `\\nu_c` introduced in Section 9.1 of [ReSt2020]_.

        .. WARNING:

            This implementation works only when the initial exchange matrix is acyclic.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,1,1],[-1,0,1],[-1,-1,0]]))
            sage: A.g_vector_to_d_vector((-1,1,2))
            (1, 0, -1)
        """
        E = -self.euler_matrix()
        c = self.coxeter_element()
        dp = vector(ZZ, self.rank())
        g = vector(g)
        for i in c:
            dp[i] = -min(g[i], 0)
            g += min(g[i], 0) * E.column(i)
        return tuple(-g + dp)

    def g_vectors(self, mutating_F=True):
        r"""
        Return an iterator producing all the g-vectors of ``self``.

        INPUT:

        - ``mutating_F`` -- boolean (default: ``True``); whether to compute
          F-polynomials; disable this for speed considerations

        ALGORITHM:

        This method does not use the caching framework provided by ``self``,
        but recomputes all the g-vectors from scratch. On the other hand it
        stores the results so that other methods like :meth:`g_vectors_so_far`
        can access them afterwards.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: len(list(A.g_vectors()))
            9
        """
        seeds = self.seeds(mutating_F=mutating_F)
        found_so_far = set()
        for g in next(seeds).g_vectors():
            found_so_far.add(g)
            yield g
        for S in seeds:
            j = S.path_from_initial_seed()[-1]
            g = S.g_vector(j)
            if g not in found_so_far:
                found_so_far.add(g)
                yield g

    def cluster_variables(self):
        r"""
        Return an iterator producing all the cluster variables of ``self``.

        ALGORITHM:

        This method does not use the caching framework provided by ``self``,
        but recomputes all the cluster variables from scratch. On the other
        hand it stores the results so that other methods like
        :meth:`cluster_variables_so_far` can access them afterwards.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: len(list(A.cluster_variables()))
            9
        """
        return map(self.cluster_variable, self.g_vectors())

    def F_polynomials(self):
        r"""
        Return an iterator producing all the F_polynomials of ``self``.

        ALGORITHM:

        This method does not use the caching framework provided by ``self``,
        but recomputes all the F-polynomials from scratch. On the other hand
        it stores the results so that other methods like
        :meth:`F_polynomials_so_far` can access them afterwards.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: len(list(A.F_polynomials()))
            9
        """
        return map(self.F_polynomial, self.g_vectors())

    def g_vectors_so_far(self) -> list:
        r"""
        Return a list of the g-vectors of cluster variables encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.current_seed().mutate(0)
            sage: sorted(A.g_vectors_so_far())
            [(-1, 1), (0, 1), (1, 0)]
        """
        return list(self._path_dict)

    def cluster_variables_so_far(self) -> list:
        r"""
        Return a list of the cluster variables encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.current_seed().mutate(0)
            sage: sorted(A.cluster_variables_so_far(), key=str)
            [(x1 + 1)/x0, x0, x1]
        """
        return [self.cluster_variable(v) for v in self.g_vectors_so_far()]

    def F_polynomials_so_far(self) -> list:
        r"""
        Return a list of the F-polynomials encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.current_seed().mutate(0)
            sage: sorted(A.F_polynomials_so_far(), key=str)
            [1, 1, u0 + 1]
        """
        return list(self._F_poly_dict.values())

    @cached_method(key=lambda a, b: tuple(b))
    def cluster_variable(self, g_vector):
        r"""
        Return the cluster variable with g-vector ``g_vector`` if it has
        been found.

        INPUT:

        - ``g_vector`` -- tuple; the g-vector of the cluster variable to return

        ALGORITHM:

        This function computes cluster variables from their g-vectors and
        F-polynomials using the "separation of additions" formula of
        Theorem 3.7 in [FZ2007]_.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.initial_seed().mutate(0)
            sage: A.cluster_variable((-1, 1))
            (x1 + 1)/x0
        """
        g_vector = tuple(g_vector)
        F = self.F_polynomial(g_vector)
        F_std = F.subs(self._yhat)
        g_mon = prod(self.ambient().gen(i) ** g_vector[i] for i in range(self.rank()))
        F_trop = self.ambient()(F.subs(self._y))._fraction_pair()[1]
        return self.retract(g_mon * F_std * F_trop)

    def F_polynomial(self, g_vector):
        r"""
        Return the F-polynomial with g-vector ``g_vector`` if it has
        been found.

        INPUT:

        - ``g_vector`` -- tuple; the g-vector of the F-polynomial to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.clear_computed_data()
            sage: A.F_polynomial((-1, 1))
            Traceback (most recent call last):
            ...
            KeyError: 'the g-vector (-1, 1) has not been found yet'
            sage: A.initial_seed().mutate(0, mutating_F=False)
            sage: A.F_polynomial((-1, 1))
            Traceback (most recent call last):
            ...
            KeyError: 'the F-polynomial with g-vector (-1, 1) has not been computed yet;
             you can compute it by mutating from the initial seed along the sequence [0]'
            sage: A.initial_seed().mutate(0)
            sage: A.F_polynomial((-1, 1))
            u0 + 1
        """
        g_vector = tuple(g_vector)
        try:
            return self._F_poly_dict[g_vector]
        except KeyError:
            if g_vector in self._path_dict:
                msg = "the F-polynomial with g-vector {} has not been computed yet; ".format(g_vector)
                msg += "you can compute it by mutating from the initial seed along the sequence "
                msg += str(self._path_dict[g_vector])
                raise KeyError(msg)
            else:
                raise KeyError("the g-vector %s has not been found yet" % str(g_vector))

    def find_g_vector(self, g_vector, depth=infinity):
        r"""
        Return a mutation sequence to obtain a seed containing the g-vector
        ``g_vector`` from the initial seed.

        INPUT:

        - ``g_vector`` -- tuple; the g-vector to find
        - ``depth`` -- positive integer or infinity (default: ``infinity``);
          the maximum distance from ``self.current_seed`` to reach

        OUTPUT:

        This function returns a list of integers if it can find ``g_vector``,
        otherwise it returns ``None``.  If the exploring iterator stops, it
        means that the algebra is of finite type and ``g_vector`` is not the
        g-vector of any cluster variable. In this case the function resets the
        iterator and raises an error.

        EXAMPLES::

            sage: A = ClusterAlgebra(['G', 2], principal_coefficients=True)
            sage: A.clear_computed_data()
            sage: A.find_g_vector((-2, 3), depth=2)
            sage: A.find_g_vector((-2, 3), depth=3)
            [0, 1, 0]
            sage: A.find_g_vector((1, 1), depth=3)
            sage: A.find_g_vector((1, 1), depth=4)
            Traceback (most recent call last):
            ...
            ValueError: (1, 1) is not the g-vector of any cluster variable of a
             Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1
             over Integer Ring
        """
        g_vector = tuple(g_vector)
        while g_vector not in self.g_vectors_so_far() and self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                if isinstance(seed, ClusterAlgebraSeed):
                    self._explored_depth = seed.depth()
                else:
                    # We got an exception because self._sd_iter caught a KeyboardInterrupt, let's raise it again
                    raise seed
            except StopIteration:
                # Unless self._sd_iter has been manually altered, we checked
                # all the seeds of self and did not find g_vector.
                # Do some house cleaning before failing
                self.reset_exploring_iterator()
                raise ValueError("%s is not the g-vector of any cluster variable of a %s" % (str(g_vector), str(self)[2:]))
        return copy(self._path_dict.get(g_vector, None))

    def ambient(self):
        r"""
        Return the Laurent polynomial ring containing ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.ambient()
            Multivariate Laurent Polynomial Ring in x0, x1, y0, y1 over Integer Ring
        """
        return self._ambient

    def scalars(self):
        r"""
        Return the ring of scalars over which ``self`` is defined.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.scalars()
            Integer Ring
        """
        return self.base().base()

    def lift(self, x):
        r"""
        Return ``x`` as an element of :meth:`ambient`.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: x = A.cluster_variable((1, 0))
            sage: A.lift(x).parent()
            Multivariate Laurent Polynomial Ring in x0, x1, y0, y1 over Integer Ring
        """
        return self.ambient()(x.value)

    def retract(self, x):
        r"""
        Return ``x`` as an element of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: L = A.ambient()
            sage: x = L.gen(0)
            sage: A.retract(x).parent()
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
        """
        return self(x)

    @cached_method
    def gens(self) -> tuple:
        r"""
        Return the list of initial cluster variables and coefficients of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.gens()
            (x0, x1, y0, y1)
            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True, coefficient_prefix='x')
            sage: A.gens()
            (x0, x1, x2, x3)
        """
        return tuple(map(self.retract, self.ambient().gens()))

    def coefficient(self, j):
        r"""
        Return the ``j``-th coefficient of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the coefficient to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.coefficient(0)
            y0
        """
        if not isinstance(self.base(), LaurentPolynomialRing_generic):
            raise ValueError("generator not defined")
        return self.retract(self.base().gen(j))

    @cached_method
    def coefficients(self) -> tuple:
        r"""
        Return the list of coefficients of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.coefficients()
            (y0, y1)
            sage: A1 = ClusterAlgebra(['B', 2])
            sage: A1.coefficients()
            ()
        """
        if isinstance(self.base(), LaurentPolynomialRing_generic):
            return tuple(map(self.retract, self.base().gens()))
        else:
            return ()

    def coefficient_names(self) -> tuple:
        r"""
        Return the list of coefficient names.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 3])
            sage: A.coefficient_names()
            ()
            sage: A1 = ClusterAlgebra(['B', 2], principal_coefficients=True)
            sage: A1.coefficient_names()
            ('y0', 'y1')
            sage: A2 = ClusterAlgebra(['C', 3], principal_coefficients=True, coefficient_prefix='x')
            sage: A2.coefficient_names()
            ('x3', 'x4', 'x5')
        """
        return self.variable_names()[self.rank():]

    def initial_cluster_variable(self, j):
        r"""
        Return the ``j``-th initial cluster variable of ``self``.

        INPUT:

        - ``j`` -- integer in ``range(self.parent().rank())``;
          the index of the cluster variable to return

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.initial_cluster_variable(0)
            x0
        """
        return self.retract(self.ambient().gen(j))

    @cached_method
    def initial_cluster_variables(self) -> tuple:
        r"""
        Return the list of initial cluster variables of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.initial_cluster_variables()
            (x0, x1)
        """
        return tuple(map(self.retract, self.ambient().gens()[:self.rank()]))

    def initial_cluster_variable_names(self) -> tuple:
        r"""
        Return the list of initial cluster variable names.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2], principal_coefficients=True)
            sage: A.initial_cluster_variable_names()
            ('x0', 'x1')
            sage: A1 = ClusterAlgebra(['B', 2], cluster_variable_prefix='a')
            sage: A1.initial_cluster_variable_names()
            ('a0', 'a1')
        """
        return self.variable_names()[:self.rank()]

    def seeds(self, **kwargs):
        r"""
        Return an iterator running over seeds of ``self``.

        INPUT:

        - ``from_current_seed`` -- boolean (default: ``False``); whether to start
          the iterator from :meth:`current_seed` or :meth:`initial_seed`

        - ``mutating_F`` -- boolean (default: ``True``); whether to compute
          F-polynomials also; disable this for speed considerations

        - ``allowed_directions`` -- iterable of integers
          (default: ``range(self.rank())``); the directions in which to mutate

        - ``depth`` -- positive integer or infinity (default: ``infinity``);
          the maximum depth at which to stop searching

        - ``catch_KeyboardInterrupt`` -- boolean (default: ``False``); whether to
          catch ``KeyboardInterrupt`` and return it rather then raising an
          exception -- this allows the iterator returned by this method to be
          resumed after being interrupted

        ALGORITHM:

        This function traverses the exchange graph in a breadth-first search.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 4])
            sage: A.clear_computed_data()
            sage: seeds = A.seeds(allowed_directions=[3, 0, 1])
            sage: _ = list(seeds)
            sage: sorted(A.g_vectors_so_far())
            [(-1, 0, 0, 0),
             (-1, 1, 0, 0),
             (0, -1, 0, 0),
             (0, 0, 0, -1),
             (0, 0, 0, 1),
             (0, 0, 1, 0),
             (0, 1, 0, 0),
             (1, 0, 0, 0)]
        """
        # should we begin from the current seed?
        if kwargs.get('from_current_seed', False):
            seed = copy(self.current_seed())
        else:
            seed = self.initial_seed()

        # yield first seed
        yield seed

        # keep track of depth
        depth_counter = 0

        # do we mutate F-polynomials?
        mutating_F = kwargs.get('mutating_F', True)

        # which directions are we allowed to mutate into
        allowed_dirs = sorted(kwargs.get('allowed_directions',
                                         range(self.rank())))

        # setup seeds storage
        cl = frozenset(seed.g_vectors())
        clusters = {}
        clusters[cl] = [seed, copy(allowed_dirs)]

        # ready, set, go!
        gets_bigger = True
        while gets_bigger and depth_counter < kwargs.get('depth', infinity):
            # remember if we got a new seed
            gets_bigger = False

            for cl in list(clusters):
                sd, directions = clusters[cl]
                while directions:
                    try:
                        # we can mutate in some direction
                        i = directions.pop()
                        new_sd = sd.mutate(i, inplace=False, mutating_F=mutating_F)
                        new_cl = frozenset(new_sd.g_vectors())
                        if new_cl in clusters:
                            # we already had new_sd, make sure it does not mutate to sd during next round
                            j = clusters[new_cl][0].g_vectors().index(new_sd.g_vector(i))
                            try:
                                clusters[new_cl][1].remove(j)
                            except ValueError:
                                pass
                        else:
                            # we got a new seed
                            gets_bigger = True
                            # next round do not mutate back to sd and make sure we only walk three sides of squares
                            new_directions = [j for j in allowed_dirs if j > i or new_sd.b_matrix()[j, i] != 0]
                            clusters[new_cl] = [new_sd, new_directions]
                            yield new_sd
                    except KeyboardInterrupt as e:
                        if kwargs.get('catch_KeyboardInterrupt', False):
                            print("caught a KeyboardInterrupt; cleaning up before returning")
                            # mutation in direction i was not completed; put it back in for next round
                            directions.append(i)
                            yield e
                            continue
                        else:
                            raise e
            # we went one step deeper
            depth_counter += 1

    def reset_exploring_iterator(self, mutating_F=True):
        r"""
        Reset the iterator used to explore ``self``.

        INPUT:

        - ``mutating_F`` -- boolean (default: ``True``); whether to also compute
          F-polynomials; disable this for speed considerations

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 4])
            sage: A.clear_computed_data()
            sage: A.reset_exploring_iterator(mutating_F=False)
            sage: A.explore_to_depth(infinity)
            sage: len(A.g_vectors_so_far())
            14
            sage: len(A.F_polynomials_so_far())
            4
        """
        self._sd_iter = self.seeds(mutating_F=mutating_F, catch_KeyboardInterrupt=True)
        self._explored_depth = 0

    def explore_to_depth(self, depth):
        r"""
        Explore the exchange graph of ``self`` up to distance ``depth``
        from the initial seed.

        INPUT:

        - ``depth`` -- positive integer or infinity; the maximum depth
          at which to stop searching

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 4])
            sage: A.explore_to_depth(infinity)
            sage: len(A.g_vectors_so_far())
            14
        """
        while self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                if isinstance(seed, ClusterAlgebraSeed):
                    self._explored_depth = seed.depth()
                else:
                    # We got an exception because self._sd_iter caught a KeyboardInterrupt, let's raise it again
                    raise seed
            except StopIteration:
                break

    def cluster_fan(self, depth=infinity):
        r"""
        Return the cluster fan (the fan of g-vectors) of ``self``.

        INPUT:

        - ``depth`` -- positive integer or infinity (default: ``infinity``);
          the maximum depth at which to compute

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', 2])
            sage: A.cluster_fan()                                                       # needs sage.geometry.polyhedron
            Rational polyhedral fan in 2-d lattice N
        """
        seeds = self.seeds(depth=depth, mutating_F=False)
        cones = [Cone(S.g_vectors()) for S in seeds]
        return Fan(cones)

    def mutate_initial(self, direction, **kwargs):
        r"""
        Return the cluster algebra obtained by mutating ``self`` at
        the initial seed.

        .. WARNING::

            This method is significantly slower than :meth:`ClusterAlgebraSeed.mutate`.
            It is therefore advisable to use the latter for exploration purposes.

        INPUT:

        - ``direction`` -- in which direction(s) to mutate, it can be:

          * an integer in ``range(self.rank())`` to mutate in one direction only
          * an iterable of such integers to mutate along a sequence
          * a string "sinks" or "sources" to mutate at all sinks or sources simultaneously

        - ``mutating_F`` -- boolean (default: ``True``); whether to compute
          F-polynomials while mutating

        .. NOTE::

            While knowing F-polynomials is essential to computing
            cluster variables, the process of mutating them is quite slow.
            If you care only about combinatorial data like g-vectors and
            c-vectors, setting ``mutating_F=False`` yields significant
            benefits in terms of speed.

        ALGORITHM:

        This function computes data for the new algebra from known data for
        the old algebra using Equation (4.2) from [NZ2012]_ for g-vectors, and
        Equation (6.21) from [FZ2007]_ for F-polynomials. The exponent `h`
        in the formula for F-polynomials is ``-min(0, old_g_vect[k])``
        due to [NZ2012]_ Proposition 4.2.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.explore_to_depth(infinity)
            sage: B = A.b_matrix()
            sage: B.mutate(0)
            sage: A1 = ClusterAlgebra(B)
            sage: A1.explore_to_depth(infinity)
            sage: A2 = A1.mutate_initial(0)
            sage: A2._F_poly_dict == A._F_poly_dict
            True

        Check that we did not mess up the original algebra because of :class:`UniqueRepresentation`::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.mutate_initial(0) is A
            False

        A faster example without recomputing F-polynomials::

            sage: A = ClusterAlgebra(matrix([[0,5],[-5,0]]))
            sage: A.mutate_initial([0,1]*10, mutating_F=False)
            A Cluster Algebra with cluster variables x20, x21 and no coefficients over Integer Ring

        Check that :issue:`28176` is fixed::

            sage: A = ClusterAlgebra( matrix(5,[0,1,-1,1,-1]), cluster_variable_names=['p13'], coefficient_names=['p12','p23','p34','p41']); A
            A Cluster Algebra with cluster variable p13 and coefficients p12, p23, p34, p41 over Integer Ring
            sage: A.mutate_initial(0)
            A Cluster Algebra with cluster variable x0 and coefficients p12, p23, p34, p41 over Integer Ring

            sage: A1 = ClusterAlgebra(['A',[2,1],1])
            sage: A2 = A1.mutate_initial([0,1,0])
            sage: len(A2.g_vectors_so_far()) == len(A2.F_polynomials_so_far())
            True
            sage: all(parent(f) == A2._U for f in A2.F_polynomials_so_far())
            True
            sage: A2.find_g_vector((0,0,1)) == []
            True
        """
        n = self.rank()

        # construct mutation sequence
        # if you change this be considerate and change also :class:`ClusterAlgebraSeed`.mutate
        if direction == "sinks":
            B = self.b_matrix()
            seq = [i for i in range(n) if all(x <= 0 for x in B.column(i))]
        elif direction == "sources":
            B = self.b_matrix()
            seq = [i for i in range(n) if all(x >= 0 for x in B.column(i))]
        else:
            try:
                seq = iter(direction)
            except TypeError:
                seq = iter((direction, ))

        # setup
        path_dict = copy(self._path_dict)
        path_to_current = copy(self.current_seed().path_from_initial_seed())
        B0 = copy(self._B0)
        cv_names = list(self.initial_cluster_variable_names())
        nfi = self._next_free_index
        I = identity_matrix(n)
        initial_g_vectors = frozenset(map(tuple, I.columns()))

        # go
        for k in seq:
            if k not in range(n):
                raise ValueError('cannot mutate in direction ' + str(k))

            # clear storage
            tmp_path_dict = {}

            # mutate B-matrix
            B0.mutate(k)

            for old_g_vect in path_dict:
                # compute new g-vector
                J = copy(I)
                old_g = vector(ZZ, old_g_vect)
                minus_eps = -old_g[k].sign()
                for j in range(n):
                    # here we have -eps*B0 rather than eps*B0
                    # because we want the k-th column of the old B0
                    J[j, k] += max(0, minus_eps * B0[j, k])
                J[k, k] = -1
                new_g_vect = tuple(J * old_g)

                # compute new path
                if new_g_vect in initial_g_vectors:
                    tmp_path_dict[new_g_vect] = []
                else:
                    new_path = path_dict[old_g_vect]
                    new_path = ([k] + new_path[:1] if new_path[:1] != [k] else []) + new_path[1:]
                    tmp_path_dict[new_g_vect] = new_path

            # update storage
            initial_g = (0,) * (k) + (1,) + (0,) * (n - k - 1)
            tmp_path_dict[initial_g] = []
            path_dict = tmp_path_dict
            path_to_current = ([k] + path_to_current[:1] if path_to_current[:1] != [k] else []) + path_to_current[1:]

            # name the new cluster variable
            cv_names[k] = self._cluster_variable_prefix + str(nfi)
            nfi += 1

        # create new algebra
        coeff_names = self.coefficient_names()
        scalars = self.scalars()
        A = ClusterAlgebra(B0, cluster_variable_names=cv_names,
                           next_free_index=nfi,
                           coefficient_names=coeff_names, scalars=scalars)

        # store computed data
        A._path_dict.update(path_dict)

        # reset self.current_seed() to the previous location
        S = A.initial_seed()
        S.mutate(path_to_current, mutating_F=False)
        A.set_current_seed(S)

        # recompute F-polynomials
        # We use forward mutation of F-polynomials because it is much faster
        # than backward mutation. Moreover the number of needed mutation is
        # linear in len(seq) rather than quadratic.
        if kwargs.get('mutating_F', True):
            for p in A._path_dict.values():
                A.initial_seed().mutate(p)

        return A

    def greedy_element(self, d_vector):
        r"""
        Return the greedy element with denominator vector ``d_vector``.

        INPUT:

        - ``d_vector`` -- tuple of 2 integers; the denominator vector of
          the element to compute

        ALGORITHM:

        This implements greedy elements of a rank 2 cluster algebra using
        Equation (1.5) from [LLZ2014]_.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', [1, 1], 1])
            sage: A.greedy_element((1, 1))
            (x0^2 + x1^2 + 1)/(x0*x1)
        """
        if self.rank() != 2:
            raise ValueError('greedy elements are only defined in rank 2')

        return self.theta_basis_element(self.d_vector_to_g_vector(d_vector))

    @cached_method(key=lambda a, b: tuple(b))
    def theta_basis_element(self, g_vector):
        r"""
        Return the element of the theta basis of ``self`` with g-vector ``g_vector``.

        INPUT:

        - ``g_vector`` -- tuple; the g-vector of the element to compute

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,-3],[2,0]]), principal_coefficients=True)
            sage: A.theta_basis_element((-1,-1))
            (x1^8*y0^4*y1 + 4*x1^6*y0^3*y1 + 6*x1^4*y0^2*y1 + x0^3*x1^2*y0 + 4*x1^2*y0*y1 + x0^3 + y1)/(x0^4*x1)

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.theta_basis_element((1, 0, 0, 0))
            Traceback (most recent call last):
            ...
            NotImplementedError: currently only implemented for cluster algebras of rank 2

        .. NOTE::

            Elements of the theta basis correspond with the associated cluster
            monomial only for appropriate coefficient choices. For example::

                sage: A = ClusterAlgebra(matrix([[0,-1],[1,0],[-1,0]]))
                sage: A.theta_basis_element((-1,0))
                (x1 + y0)/(x0*y0)

            while::

                sage: _ = A.find_g_vector((-1,0));
                sage: A.cluster_variable((-1,0))
                (x1 + y0)/x0

            In particular theta basis elements do not satisfy a separation of additions formula.

        .. WARNING::

            Currently only cluster algebras of rank 2 are supported

        .. SEEALSO::

            :meth:`sage.algebras.cluster_algebra.theta_basis_F_polynomial`
        """
        g_vector = tuple(g_vector)
        F = self.theta_basis_F_polynomial(g_vector).subs(self._yhat)
        g_mon = prod(self.ambient().gen(i) ** g_vector[i] for i in range(self.rank()))
        # we only return the monomal g_mon times the evaluated F-polynomial because this is how
        # theta basis elements behave.
        return self.retract(g_mon * F)

    @cached_method(key=lambda a, b: tuple(b))
    def theta_basis_F_polynomial(self, g_vector):
        r"""
        Return the F-polynomial of the element of the theta basis of ``self`` with g-vector ``g_vector``.

        INPUT:

        - ``g_vector`` -- tuple; the g-vector of the F-polynomial to compute

        .. WARNING::

            Elements of the theta basis do not satisfy a separation of additions formula.
            See the implementation of :meth:`sage.algebras.cluster_algebra.theta_basis_F_polynomial`
            for further details.

        ALGORITHM:

        This method uses the fact that the greedy basis and the theta basis
        coincide in rank 2 and uses the former defining recursion (Equation
        (1.5) from [LLZ2014]_) to compute.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix([[0,-3],[2,0]]), principal_coefficients=True)
            sage: A.theta_basis_F_polynomial((-1,-1))
            u0^4*u1 + 4*u0^3*u1 + 6*u0^2*u1 + 4*u0*u1 + u0 + u1 + 1

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.theta_basis_F_polynomial((1, 0, 0, 0))
            Traceback (most recent call last):
            ...
            NotImplementedError: currently only implemented for cluster algebras of rank 2
        """
        if self.rank() != 2:
            raise NotImplementedError("currently only implemented for cluster algebras of rank 2")

        # extract the part of g_vector not coming from the initial cluster
        d = tuple(max(x, 0) for x in self.g_vector_to_d_vector(g_vector))
        g = self.d_vector_to_g_vector(d)

        shifts = ((d[0] + g[0]) / self._B0[0, 1], (d[1] + g[1]) / self._B0[1, 0])
        signs = (self._B0[0, 1].sign(), self._B0[1, 0].sign())

        u = list(self._U.gens())
        output = self._U.zero()
        for p in range(d[1] + 1):
            for q in range(d[0] + 1):
                output += self._greedy_coefficient(d, p, q) * u[1] ** (signs[0] * p - shifts[0]) * u[0] ** (signs[1] * q - shifts[1])
        return output

    @cached_method
    def _greedy_coefficient(self, d_vector, p, q):
        r"""
        Return the coefficient of the monomial ``x1 ** (b * p) * x2 ** (c * q)``
        in the numerator of the greedy element with denominator vector ``d_vector``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A', [1, 1], 1])
            sage: A.greedy_element((1, 1))
            (x0^2 + x1^2 + 1)/(x0*x1)
            sage: A._greedy_coefficient((1, 1), 0, 0)
            1
            sage: A._greedy_coefficient((1, 1), 1, 0)
            1
        """
        b = abs(self._B0[0, 1])
        c = abs(self._B0[1, 0])
        a1, a2 = d_vector
        p = Integer(p)
        q = Integer(q)
        if p == 0 and q == 0:
            return ZZ.one()
        sum1 = 0
        for k in range(1, p + 1):
            bino = 0
            if a2 - c * q + k - 1 >= k:
                bino = binomial(a2 - c * q + k - 1, k)
            sum1 += (-1) ** (k - 1) * self._greedy_coefficient(d_vector, p - k, q) * bino
        sum2 = 0
        for l in range(1, q + 1):
            bino = 0
            if a1 - b * p + l - 1 >= l:
                bino = binomial(a1 - b * p + l - 1, l)
            sum2 += (-1) ** (l - 1) * self._greedy_coefficient(d_vector, p, q - l) * bino
        return Integer(max(sum1, sum2))

    # DESIDERATA
    # Some of these are probably unrealistic
    def upper_cluster_algebra(self):
        r"""
        Return the upper cluster algebra associated to ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.upper_cluster_algebra()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented yet
        """
        raise NotImplementedError("not implemented yet")

    def upper_bound(self):
        r"""
        Return the upper bound associated to ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.upper_bound()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented yet
        """
        raise NotImplementedError("not implemented yet")

    def lower_bound(self):
        r"""
        Return the lower bound associated to ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F', 4])
            sage: A.lower_bound()
            Traceback (most recent call last):
            ...
            NotImplementedError: not implemented yet
        """
        raise NotImplementedError("not implemented yet")
