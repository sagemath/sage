r"""
Catalog of discrete dynamical systems

This module contains constructors for several specific discrete
dynamical systems.
These are accessible through
:mod:`sage.dynamics.finite_dynamical_system_catalog. <sage.dynamics.finite_dynamical_system_catalog>`
or just through ``finite_dynamical_systems.``
(type either of these in Sage and hit ``tab`` for a list).

AUTHORS:

- Darij Grinberg, Tom Roby (2018): initial version

Functions
=========
"""
#*****************************************************************************
#       Copyright (C) 2018 Darij Grinberg <darijgrinberg@gmail.com>,
#                     2018 Tom Roby <tomrobyuconn@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.dynamics.finite_dynamical_system import DiscreteDynamicalSystem, \
        FiniteDynamicalSystem, InvertibleDiscreteDynamicalSystem, \
        InvertibleFiniteDynamicalSystem

def permutation(pi, invertible=True):
    r"""
    Return the invertible finite discrete dynamical system
    induced by the permutation ``pi`` of the set
    `\{1, 2, \ldots, n\}`.

    The permutation ``pi`` can be provided as an instance
    of :class:`Permutation`, but also as anything that can
    be cast into that class.

    See also :meth:`one_line` for a more general
    construction, which takes any map from
    `\{1, 2, \ldots, n\}` to `\{1, 2, \ldots, n\}` given in
    one-line notation, and builds a (not necessarily
    invertible) discrete dynamical system out of it.

    EXAMPLES::

        sage: F = finite_dynamical_systems.permutation([3, 5, 4, 1, 2])
        sage: F.verify_inverse_evolution()
        True
        sage: sorted(F.orbit_lengths())
        [2, 3]
        sage: F.orbit(3)
        [3, 4, 1]
        sage: F.is_homomesic(lambda x: 1)
        True
        sage: F.is_homomesic(lambda x: x)
        False
    """
    from sage.combinat.permutation import Permutation
    pi = Permutation(pi)
    n = len(pi)
    X = range(1, n+1)
    return InvertibleFiniteDynamicalSystem(X, pi, inverse=pi.inverse(), create_tuple=True)

def one_line(xs):
    r"""
    Return the finite discrete dynamical system
    with ground set `\{1, 2, \ldots, n\}` and evolution
    sending each `i` to `x_i`,
    where `(x_1, x_2, \ldots, x_n)` is the argument ``xs``
    provided.

    EXAMPLES::

        sage: F = finite_dynamical_systems.one_line([2, 2, 1, 2, 3])
        sage: F.orbit(3)
        [3, 1, 2]
        sage: F.orbit(5)
        [5, 3, 1, 2]
        sage: F.orbit(2)
        [2]
    """
    n = len(xs)
    X = range(1, n + 1)
    xs2 = tuple(xs)

    def pi(i):
        return xs2[i - 1]
    return FiniteDynamicalSystem(X, pi, create_tuple=True)


def bitstring_rotation(n, ones=None):
    r"""
    Return the invertible finite discrete dynamical system
    consisting of all bitstrings of size `n` (that is, of
    all `n`-tuples `(i_1, \ldots, i_n) \in \{0, 1\}^n`),
    evolving by cyclic rotation.

    If the optional parameter ``ones`` is provided, the
    system is restricted only to those bitstrings whose
    number of ones is the value of this parameter.

    EXAMPLES::

        sage: F = finite_dynamical_systems.bitstring_rotation(5)
        sage: sorted(F.orbit_lengths())
        [1, 1, 5, 5, 5, 5, 5, 5]
        sage: F.orbit((0, 1, 1, 0, 1))
        [(0, 1, 1, 0, 1),
         (1, 1, 0, 1, 0),
         (1, 0, 1, 0, 1),
         (0, 1, 0, 1, 1),
         (1, 0, 1, 1, 0)]
        sage: F.is_homomesic(lambda x: sum(1 for i in range(5) for j in range(i) if x[j] > x[i]))
        False
        sage: F.is_homomesic(lambda x: x[0])
        False
        sage: F = finite_dynamical_systems.bitstring_rotation(5, ones=3)
        sage: F.verify_inverse_evolution()
        True
        sage: sorted(F.orbit_lengths())
        [5, 5]
        sage: F.orbit((0, 1, 1, 0, 1))
        [(0, 1, 1, 0, 1),
         (1, 1, 0, 1, 0),
         (1, 0, 1, 0, 1),
         (0, 1, 0, 1, 1),
         (1, 0, 1, 1, 0)]
        sage: F.is_homomesic(lambda x: sum(1 for i in range(5) for j in range(i) if x[j] > x[i]))
        True
        sage: F.is_homomesic(lambda x: x[0])
        True
    """
    if ones is None:
        from sage.categories.cartesian_product import cartesian_product
        X = cartesian_product([[0,1]] * n)
    else:
        from itertools import combinations
        X = [tuple((1 if i in cs else 0) for i in range(n))
             for cs in combinations(range(n), ones)]
    if n == 0:
        phi = lambda x: x
        psi = phi
    else:
        phi = lambda x: x[1:] + (x[0],)
        psi = lambda x: (x[-1],) + x[:-1]
    return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)

def striker_sweep(E, predicate, elements, lazy=False):
    r"""
    Return the invertible finite discrete dynamical system
    on all subsets of a finite set ``E`` satisfying a
    boolean predicate ``predicate``, where evolution is
    the "Striker sweep" -- i.e., the composition of
    Striker toggles corresponding to the elements in the
    iterable ``elements`` (from first to last).

    Let `E` be a finite set.
    Let `\mathcal{L}` be a subset of the powerset of `E`.
    (In this implementation, `\mathcal{L}` should be
    specified via the boolean predicate ``predicate``,
    which takes a subset `F` of `E` and returns the truth
    value of `F \in \mathcal{L}`.)
    For any `e \in E`, the *Striker toggle* `t_e` is
    the involution of the set `\mathcal{L}` that sends
    each `F \in \mathcal{L}` to the symmetric difference
    `F \triangle \{ e \}` if this symmetric difference
    is in `\mathcal{L}`, and otherwise to `F` itself.
    If `(e_1, e_2, \ldots, e_k)` is a finite sequence of
    elements of `E` (to be provided as a tuple, via the
    argument ``elements``), then the *Striker sweep*
    corresponding to this sequence is the composition of
    maps
    `t_{e_k} \circ t_{e_{k-1}} \circ \cdots \circ t_{e_1}`.

    This generalizes classical constructions such as
    rowmotion on order ideals.

    The optional argument ``lazy`` can be set to
    ``True``; in that case, the ground set of the
    dynamical system will not be explicitly computed.

    .. WARNING::

        The sequence ``elements`` should be provided as
        a tuple, or as a list that is guaranteed not to
        be mutated (as otherwise, mutation will corrupt
        this DDS).

    .. TODO::

        Implement the ``lazy=True`` case.
        This should differ in that ``X`` is no longer a
        list, but an enumerated set.
        But how to build an enumerated set filtering a
        given enumerated set according to a predicate?

    EXAMPLES::

        sage: StS = finite_dynamical_systems.striker_sweep
        sage: E = range(1, 5)
        sage: lac = lambda S: all(s + 1 not in S for s in S) # lacunarity predicate
        sage: F = StS(E, lac, [1, 2, 3, 4])
        sage: F.ground_set()
        [{}, {1}, {2}, {3}, {4}, {1, 3}, {1, 4}, {2, 4}]
        sage: F.evolution()(Set([2, 4]))
        {}
        sage: F.evolution()(Set([]))
        {1, 3}
        sage: F.evolution()(Set([1, 3]))
        {4}
        sage: F.evolution()(Set([4]))
        {1}
        sage: F.inverse_evolution()(Set([1]))
        {4}
        sage: F.verify_inverse_evolution()
        True
        sage: sorted(F.orbit_lengths())
        [3, 5]
        sage: F.orbit(Set([2, 4]))
        [{2, 4}, {}, {1, 3}, {4}, {1}]
        sage: F.is_homomesic(lambda S: S.cardinality())
        False
        sage: F.is_homomesic(lambda S: bool(1 in S) - bool(4 in S), find_average=True)
        0
        sage: F.is_homomesic(lambda S: bool(2 in S) - bool(3 in S), find_average=True)
        0
        sage: F.is_homomesic(lambda S: bool(1 in S), find_average=True)
        False
    """
    from sage.combinat.subset import Subsets
    from sage.sets.set import Set
    X = [F for F in Subsets(E) if predicate(F)]

    def phi(F):
        for e in elements:
            G = F.symmetric_difference(Set([e]))
            if predicate(G):
                F = G
        return F

    def psi(F):
        for e in reversed(elements):
            G = F.symmetric_difference(Set([e]))
            if predicate(G):
                F = G
        return F
    return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)


def syt_promotion(lam):
    r"""
    Return the invertible finite discrete dynamical system
    consisting of all standard tableaux of shape ``lam`` (a
    given partition) and evolving according to promotion.

    EXAMPLES::

        sage: F = finite_dynamical_systems.syt_promotion([4, 4, 4])
        sage: sorted(F.orbit_lengths())
        [3, 3, 4, 4, 4, 6, 6, 6, 6, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12]
        sage: G = finite_dynamical_systems.syt_promotion([4, 3, 1])
        sage: sorted(G.orbit_lengths())
        [16, 22, 32]
        sage: G.verify_inverse_evolution()
        True
    """
    from sage.combinat.partition import Partition
    from sage.combinat.tableau import StandardTableaux
    lam = Partition(lam)
    X = StandardTableaux(lam)
    return InvertibleFiniteDynamicalSystem(X, lambda T: T.promotion(),
                                           inverse=lambda T: T.promotion_inverse())


def order_ideal_rowmotion(P):
    r"""
    Return the invertible finite discrete dynamical system
    consisting of all order ideals of the poset ``P``,
    evolving according to rowmotion.

    EXAMPLES::

        sage: P = RootSystem(["A", 6]).root_poset()
        sage: F = finite_dynamical_systems.order_ideal_rowmotion(P)
        sage: sorted(F.orbit_lengths())
        [2, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
         14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
        sage: F.is_homomesic(lambda I: len(I))
        False
        sage: F.is_homomesic(lambda I: sum((-1)**(P.rank(i)) for i in I))
        True

        sage: P = RootSystem(["A", 3]).root_poset()
        sage: F = finite_dynamical_systems.order_ideal_rowmotion(P)
        sage: F.verify_inverse_evolution()
        True
    """
    from sage.sets.set import Set
    X = [Set(P.order_ideal(A)) for A in P.antichains()]
    # Using P.order_ideals_lattice() instead causes intransparency issues:
    # sage can't always do P.rowmotion(I) when I is in P.order_ideals_lattice().
    # Bug in P.order_ideals_lattice() when P is facade?
    phi = P.rowmotion

    def psi(I):  # inverse of rowmotion
        result = I
        for i in P.linear_extension():
            result = P.order_ideal_toggle(result, i)
        return result
    return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)


def bulgarian_solitaire(n):
    r"""
    Return the finite discrete dynamical system defined
    by Bulgarian solitaire on partitions of size `n`.

    Let `n` be a nonnegative integer.
    Let `P` be the set of all integer partitions of
    size `n`.
    Let `B : P \to P` be the map that sends each
    partition
    `\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_k)`
    of `n` (with all the `\lambda_i` positive) to
    `(\lambda_1 - 1, \lambda_2 - 1, \ldots, \lambda_k - 1, k)`,
    where zero entries have been removed and the remaining
    entries sorted into decreasing order.
    (For example,
    `B(5, 4, 2, 2, 1, 1) = (6, 4, 3, 1, 1)`.)
    This method yields the finite DDS whose ground set
    is `P` and whose evolution is `B`.

    EXAMPLES::

        sage: BS = finite_dynamical_systems.bulgarian_solitaire
        sage: BS(3).evolution()(Partition([3]))
        [2, 1]
        sage: BS(3).evolution()(Partition([2, 1]))
        [2, 1]
        sage: BS(3).evolution()(Partition([1, 1, 1]))
        [3]
        sage: BS(4).evolution()(Partition([4]))
        [3, 1]
        sage: BS(4).orbit(Partition([4]))
        [[4], [3, 1], [2, 2], [2, 1, 1]]
        sage: BS(4).orbit(Partition([3, 1]))
        [[3, 1], [2, 2], [2, 1, 1]]
        sage: BS(7).orbit(Partition([6, 1]), preperiod=True)
        ([[6, 1], [5, 2], [4, 2, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1]], 2)
        sage: BS(6).is_homomesic(lambda lam: len(lam))
        True
        sage: BS(6).is_homomesic(lambda lam: lam[0])
        True
        sage: BS(6).is_homomesic(lambda lam: lam[-1])
        True
        sage: BS(8).is_homomesic(lambda lam: len(lam))
        True
        sage: BS(8).is_homomesic(lambda lam: lam[0])
        True
        sage: BS(8).is_homomesic(lambda lam: lam[-1])
        False
    """
    from sage.combinat.partition import Partition, Partitions
    X = Partitions(n)

    def phi(lam):
        mu = [p - 1 for p in lam if p > 0]
        nu = sorted(mu + [len(lam)], reverse=True)
        return Partition(nu)
    return FiniteDynamicalSystem(X, phi)
