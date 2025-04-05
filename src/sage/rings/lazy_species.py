r"""
Lazy Combinatorial Species

We regard a combinatorial species as a sequence of group actions of
the symmetric groups `\mathfrak S_n`, for `n\in\NN`.

Coefficients of lazy species are computed on demand.  They have
infinite precision, although equality can only be decided in special
cases.

AUTHORS:

- Mainak Roy, Martin Rubey, Travis Scrimshaw (2024-2025)

EXAMPLES::

    We can reproduce the molecular expansions from Appendix B in
    [GL2011]_ with little effort.  The molecular expansion of the
    species of point determining graphs can be computed as the
    species of graphs composed with the compositional inverse of the
    species of non-empty sets.  To make the result more readable, we
    provide a name for `E_2(X^2)`::

        sage: from sage.rings.lazy_species import LazySpecies
        sage: L.<X> = LazySpecies(QQ)
        sage: E = L.Sets()
        sage: E_2 = L(SymmetricGroup(2))
        sage: Ep = L.Sets().restrict(1)
        sage: G = L.Graphs()

    The molecular decomposition begins with::

        sage: P = G(Ep.revert())
        sage: P.truncate(6)
        1 + X + E_2 + (E_3+X*E_2) + (E_4+X*E_3+E_2(E_2)+X^2*E_2+E_2(X^2))
        + (E_5+E_2*E_3+X*E_4+X*E_2^2+X^2*E_3+2*X*E_2(E_2)+P_5+5*X*E_2(X^2)+3*X^3*E_2)

    Note that [GL2011]_ write `D_5` instead of `P_5`, and there is
    apparently a misprint: `X*E_2(E_2) + 4 X^3 E_2` should be `2 X
    E_2(E_2) + 3 X^3 E_2`.

    To compute the molecular decomposition of the species of
    connected graphs with no endpoints, we use Equation (3.3) in
    [GL2011]_.  Before that we need to define the species of
    connected graphs::

        sage: Gc = Ep.revert()(G-1)
        sage: Mc = Gc(X*E(-X)) + E_2(-X)
        sage: E(Mc).truncate(5)
        1 + X + E_2 + 2*E_3 + (2*E_4+E_2(E_2)+E_2^2+X*E_3)

    Note that [GL2011]_ apparently contains a misprint: `2 X E_3`
    should be `X E_3 + E_2^2`.  Indeed, the graphs on four vertices
    without endpoints are the complete graph and the empty graph, the
    square, the diamond graph and the triangle with an extra isolated
    vertex.


    To compute the molecular decomposition of the species of
    bi-point-determining graphs we use Corollary (4.6) in
    [GL2011]_::

        sage: B = G(2*Ep.revert() - X)
        sage: B.truncate(6)
        1 + X + E_2(X^2) + (P_5+5*X*E_2(X^2))

"""

from sage.misc.lazy_list import lazy_list
from sage.rings.integer_ring import ZZ
from sage.rings.lazy_series import LazyCompletionGradedAlgebraElement, LazyModuleElement
from sage.rings.lazy_series_ring import (LazyCompletionGradedAlgebra,
                                         LazyPowerSeriesRing,
                                         LazySymmetricFunctions)
from sage.rings.species import PolynomialSpecies, _label_sets
from sage.data_structures.stream import (Stream_zero,
                                         Stream_exact,
                                         Stream_truncated,
                                         Stream_function)
from sage.categories.tensor import tensor
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.subset import subsets
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition import Partitions
from sage.combinat.permutation import CyclicPermutations
from sage.combinat.set_partition import SetPartitions
from sage.graphs.graph_generators import graphs
from sage.groups.perm_gps.permgroup_named import SymmetricGroup, CyclicPermutationGroup
from sage.structure.element import parent
import itertools
from collections import defaultdict


def weighted_compositions(n, d, weight_multiplicities, _w0=0):
    r"""
    Return all compositions of `n` of weight `d`.

    The weight of a composition `n_1, n_2, \dots` is `\sum_i w_i n_i`.

    INPUT:

    - ``n`` -- a nonnegative integer, the sum of the parts
    - ``d`` -- a nonnegative integer, the total weight
    - ``weight_multiplicities`` -- an iterable,
      ``weight_multiplicities[i]`` is the number of positions with
      weight `i+1`.

    .. TODO::

        Possibly this could be merged with
        :class:`~sage.combinat.integer_vector_weighted.WeightedIntegerVectors`.
        However, that class does not support fixing the sum of the
        parts currently.

    EXAMPLES::

        sage: from sage.rings.lazy_species import weighted_compositions
        sage: list(weighted_compositions(1, 1, [2,1]))
        [[1, 0], [0, 1]]

        sage: list(weighted_compositions(2, 1, [2,1]))
        []

        sage: list(weighted_compositions(1, 2, [2,1,1]))
        [[0, 0, 1]]

        sage: list(weighted_compositions(3, 4, [2,2]))
        [[2, 0, 1, 0],
         [1, 1, 1, 0],
         [0, 2, 1, 0],
         [2, 0, 0, 1],
         [1, 1, 0, 1],
         [0, 2, 0, 1]]

    """
    # the empty composition exists if and only if n == d == 0
    if not n:
        if not d:
            yield []
        return
    if not d:
        return

    # otherwise we iterate over the possibilities for the first
    # weight_multiplicities[_w0] parts
    try:
        if _w0 >= len(weight_multiplicities):
            return
    except TypeError:
        pass
    if _w0 > d:
        return
    for s in range(n + 1):
        for c in weighted_compositions(n - s, d - s * (_w0 + 1), weight_multiplicities, _w0=_w0+1):
            m = weight_multiplicities[_w0]
            for v in map(list, IntegerVectors(s, length=m)):
                yield v + c


def weighted_vector_compositions(n_vec, d, weight_multiplicities_vec):
    r"""
    Return all compositions of the vector `n` of weight `d`.

    INPUT:

    - ``n_vec`` -- a `k`-tuple of non-negative integers.

    - ``d`` -- a non-negative integer, the total sum of the parts in
      all components

    - ``weight_multiplicities_vec`` -- `k`-tuple of iterables, an
      iterable, ``weight_multiplicities_vec[j][i]`` is the number of
      positions with weight `i+1` in the `j`-th component.

    EXAMPLES::

        sage: from sage.rings.lazy_species import weighted_vector_compositions
        sage: list(weighted_vector_compositions([1,1], 2, [[2,1,1], [1,1,1]]))
        [([1, 0], [1]), ([0, 1], [1])]

        sage: list(weighted_vector_compositions([3,1], 4, [[2,1,0,0,1], [2,1,0,0,1]]))
        [([3, 0], [1, 0]),
         ([3, 0], [0, 1]),
         ([2, 1], [1, 0]),
         ([2, 1], [0, 1]),
         ([1, 2], [1, 0]),
         ([1, 2], [0, 1]),
         ([0, 3], [1, 0]),
         ([0, 3], [0, 1])]
    """
    k = len(n_vec)
    for d_vec in IntegerVectors(d, length=k):
        yield from itertools.product(*map(weighted_compositions,
                                          n_vec, d_vec,
                                          weight_multiplicities_vec))

######################################################################


class LazySpeciesElement(LazyCompletionGradedAlgebraElement):
    r"""
    EXAMPLES:

    Compute the molecular expansion of `E(-X)`::

        sage: from sage.rings.lazy_species import LazySpecies
        sage: L = LazySpecies(ZZ, "X")
        sage: E = L(lambda n: SymmetricGroup(n))
        sage: E_inv = 1 / E
        sage: E_inv
        1 + (-X) + (-E_2+X^2) + (-E_3+2*X*E_2-X^3)
          + (-E_4+2*X*E_3+E_2^2-3*X^2*E_2+X^4)
          + (-E_5+2*X*E_4+2*E_2*E_3-3*X^2*E_3-3*X*E_2^2+4*X^3*E_2-X^5)
          + (-E_6+2*X*E_5+2*E_2*E_4-3*X^2*E_4+E_3^2-6*X*E_2*E_3+4*X^3*E_3-E_2^3+6*X^2*E_2^2-5*X^4*E_2+X^6)
          + O^7

    Compare with the explicit formula::

        sage: def coefficient(m):
        ....:     return sum((-1)^len(la) * multinomial((n := la.to_exp())) * prod(E[i]^ni for i, ni in enumerate(n, 1)) for la in Partitions(m))

        sage: all(coefficient(m) == E_inv[m] for m in range(10))
        True
    """
    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E.isotype_generating_series()
            1 + X + X^2 + X^3 + X^4 + X^5 + X^6 + O(X^7)

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: E(C).isotype_generating_series()
            1 + X + 2*X^2 + 3*X^3 + 5*X^4 + 7*X^5 + 11*X^6 + O(X^7)

            sage: L2.<X, Y> = LazySpecies(QQ)
            sage: E(X + Y).isotype_generating_series()
            1 + (X+Y) + (X^2+X*Y+Y^2) + (X^3+X^2*Y+X*Y^2+Y^3)
            + (X^4+X^3*Y+X^2*Y^2+X*Y^3+Y^4)
            + (X^5+X^4*Y+X^3*Y^2+X^2*Y^3+X*Y^4+Y^5)
            + (X^6+X^5*Y+X^4*Y^2+X^3*Y^3+X^2*Y^4+X*Y^5+Y^6)
            + O(X,Y)^7

            sage: C(X + Y).isotype_generating_series()
            (X+Y) + (X^2+X*Y+Y^2) + (X^3+X^2*Y+X*Y^2+Y^3)
            + (X^4+X^3*Y+2*X^2*Y^2+X*Y^3+Y^4)
            + (X^5+X^4*Y+2*X^3*Y^2+2*X^2*Y^3+X*Y^4+Y^5)
            + (X^6+X^5*Y+3*X^4*Y^2+4*X^3*Y^3+3*X^2*Y^4+X*Y^5+Y^6)
            + O(X,Y)^7
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        if P._arity == 1:
            def coefficient(n):
                return sum(c for M, c in self[n].monomial_coefficients().items())
        else:
            def coefficient(n):
                return sum(c * P.base_ring().prod(v ** d for v, d in zip(L.gens(), M.grade()))
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E.generating_series()
            1 + X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + O(X^7)

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: C.generating_series()
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4 + 1/5*X^5 + 1/6*X^6 + O(X^7)

            sage: L2.<X, Y> = LazySpecies(QQ)
            sage: E(X + Y).generating_series()
            1 + (X+Y) + (1/2*X^2+X*Y+1/2*Y^2)
            + (1/6*X^3+1/2*X^2*Y+1/2*X*Y^2+1/6*Y^3)
            + (1/24*X^4+1/6*X^3*Y+1/4*X^2*Y^2+1/6*X*Y^3+1/24*Y^4)
            + (1/120*X^5+1/24*X^4*Y+1/12*X^3*Y^2+1/12*X^2*Y^3+1/24*X*Y^4+1/120*Y^5)
            + (1/720*X^6+1/120*X^5*Y+1/48*X^4*Y^2+1/36*X^3*Y^3+1/48*X^2*Y^4+1/120*X*Y^5+1/720*Y^6)
            + O(X,Y)^7

            sage: C(X + Y).generating_series()
            (X+Y) + (1/2*X^2+X*Y+1/2*Y^2) + (1/3*X^3+X^2*Y+X*Y^2+1/3*Y^3)
            + (1/4*X^4+X^3*Y+3/2*X^2*Y^2+X*Y^3+1/4*Y^4)
            + (1/5*X^5+X^4*Y+2*X^3*Y^2+2*X^2*Y^3+X*Y^4+1/5*Y^5)
            + (1/6*X^6+X^5*Y+5/2*X^4*Y^2+10/3*X^3*Y^3+5/2*X^2*Y^4+X*Y^5+1/6*Y^6)
            + O(X,Y)^7
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        if P._arity == 1:
            def coefficient(n):
                return sum(c / M.permutation_group()[0].cardinality()
                           for M, c in self[n].monomial_coefficients().items())
        else:
            def coefficient(n):
                return sum(c / M.permutation_group()[0].cardinality()
                           * P.base_ring().prod(v ** d for v, d in zip(L.gens(), M.grade()))
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def cycle_index_series(self):
        r"""
        Return the cycle index series for this species.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: h = SymmetricFunctions(QQ).h()
            sage: LazySymmetricFunctions(h)(E.cycle_index_series())
            h[] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + O^7

            sage: s = SymmetricFunctions(QQ).s()
            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: s(C.cycle_index_series()[5])
            s[1, 1, 1, 1, 1] + s[2, 2, 1] + 2*s[3, 1, 1] + s[3, 2] + s[5]

            sage: L = LazySpecies(QQ, "X")
            sage: E = L.Sets()
            sage: L2.<X, Y> = LazySpecies(QQ)
            sage: E(X + Y).cycle_index_series()[3]
            1/6*p[] # p[1, 1, 1] + 1/2*p[] # p[2, 1] + 1/3*p[] # p[3]
            + 1/2*p[1] # p[1, 1] + 1/2*p[1] # p[2] + 1/2*p[1, 1] # p[1]
            + 1/6*p[1, 1, 1] # p[] + 1/2*p[2] # p[1] + 1/2*p[2, 1] # p[]
            + 1/3*p[3] # p[]
        """
        P = self.parent()
        p = SymmetricFunctions(P.base_ring().fraction_field()).p()
        if P._arity == 1:
            L = LazySymmetricFunctions(p)

            def coefficient(n):
                return sum(c * M.permutation_group()[0].cycle_index()
                           for M, c in self[n].monomial_coefficients().items())
        else:
            L = LazySymmetricFunctions(tensor([p for _ in range(P._arity)]))

            def coefficient(n):
                return sum(c * M.cycle_index()
                           for M, c in self[n].monomial_coefficients().items())

        return L(coefficient)

    def restrict(self, min_degree=None, max_degree=None):
        r"""
        Return the series obtained by keeping only terms of
        degree between ``min_degree`` and ``max_degree``.

        INPUT:

        - ``min_degree``, ``max_degree`` -- (optional) integers
          indicating which degrees to keep

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: G = L.Graphs()
            sage: list(G.isotypes(2))
            [Graph on 2 vertices, Graph on 2 vertices]

            sage: list(G.restrict(2, 2).isotypes(2))
            [Graph on 2 vertices, Graph on 2 vertices]
        """
        return RestrictedSpeciesElement(self, min_degree, max_degree)

    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        In particular, the method to obtain the structures is
        adapted.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: F = L(lambda n: SymmetricGroup(n))
            sage: list(E.structures([1,2,3]))
            [(E_3, ((1, 2, 3),))]
            sage: list((E+F).structures([1,2,3]))
            [(E_3, ((1, 2, 3),)), (E_3, ((1, 2, 3),))]

        """
        return SumSpeciesElement(self, other)

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        In particular, the method to obtain the structures is
        adapted.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: sorted((E^2).structures([1,2,3]))
            [((1, ()), (E_3, ((1, 2, 3),))),
             ((X, ((1,),)), (E_2, ((2, 3),))),
             ((X, ((2,),)), (E_2, ((1, 3),))),
             ((X, ((3,),)), (E_2, ((1, 2),))),
             ((E_2, ((1, 2),)), (X, ((3,),))),
             ((E_2, ((1, 3),)), (X, ((2,),))),
             ((E_2, ((2, 3),)), (X, ((1,),))),
             ((E_3, ((1, 2, 3),)), (1, ()))]
        """
        return ProductSpeciesElement(self, other)

    def structures(self, *labels):
        r"""
        Iterate over the structures on the given set of labels.

        Generically, this yields a list of relabelled representatives
        of the cosets of corresponding groups.

        The relabelling is such that the first few labels correspond
        to the first factor in the atomic decomposition, etc.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: list(E.structures([1,2,3]))
            [(E_3, ((1, 2, 3),))]

            sage: P = L(lambda n: CyclicPermutationGroup(n))
            sage: list(P.structures([1,2,3]))
            [(C_3, ((1, 2, 3),)), (C_3, ((1, 3, 2),))]

            sage: F = 1/(2-E)
            sage: sorted(F.structures([1,2,3]))
            [(E_3, ((1, 2, 3),)),
             (X*E_2, ((1,), (2, 3)), 0),
             (X*E_2, ((1,), (2, 3)), 1),
             (X*E_2, ((2,), (1, 3)), 0),
             (X*E_2, ((2,), (1, 3)), 1),
             (X*E_2, ((3,), (1, 2)), 0),
             (X*E_2, ((3,), (1, 2)), 1),
             (X^3, ((1,), (2,), (3,))),
             (X^3, ((1,), (3,), (2,))),
             (X^3, ((2,), (1,), (3,))),
             (X^3, ((2,), (3,), (1,))),
             (X^3, ((3,), (1,), (2,))),
             (X^3, ((3,), (2,), (1,)))]

            sage: from sage.rings.species import PolynomialSpecies
            sage: L = LazySpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: XY = L(P(PermutationGroup([], domain=[1, 2]), {0: [1], 1: [2]}))
            sage: list((XY).structures([1], ["a"]))
            [(X*Y, ((1,), ('a',)))]

            sage: sorted(E(XY).structures([1,2], [3, 4]))
            [((E_2, ((((1, 'X'), (3, 'Y')), ((2, 'X'), (4, 'Y'))),)),
              ((X*Y, ((1,), (3,))), (X*Y, ((2,), (4,))))),
             ((E_2, ((((1, 'X'), (4, 'Y')), ((2, 'X'), (3, 'Y'))),)),
              ((X*Y, ((1,), (4,))), (X*Y, ((2,), (3,)))))]

            sage: list(XY.structures([], [1, 2]))
            []
        """
        yield from self[sum(map(len, labels))].structures(*labels)

    def isotypes(self, labels):
        pass

    def polynomial(self, degree=None, names=None):
        r"""
        Return ``self`` as a polynomial if ``self`` is actually so.

        INPUT:

        - ``degree`` -- ``None`` or an integer
        - ``names`` -- names of the variables; if it is ``None``, the name of
          the variables of the series is used

        OUTPUT:

        If ``degree`` is not ``None``, the terms of the series of
        degree greater than ``degree`` are first truncated.  If
        ``degree`` is ``None`` and the series is not a polynomial
        polynomial, a ``ValueError`` is raised.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E.polynomial(3)
            1 + X + E_2 + E_3
        """
        S = self.parent()
        R = S._laurent_poly_ring

        if isinstance(self._coeff_stream, Stream_zero):
            return R.zero()

        if degree is None:
            if (isinstance(self._coeff_stream, Stream_exact)
                and not self._coeff_stream._constant):
                m = self._coeff_stream._degree
            else:
                raise ValueError("not a polynomial species")
        else:
            m = degree + 1

        return R.sum(self[:m])

    def __call__(self, *args):
        """

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E2 = L(SymmetricGroup(2))
            sage: E2(E2)
            E_2(E_2) + O^11

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(QQ, "X")
            sage: Gc = L(lambda n: sum(P(G.automorphism_group()) for G in graphs(n) if G.is_connected()) if n else 0)
            sage: E = L.Sets()
            sage: G = L.Graphs()
            sage: E(Gc) - G
            O^7

            sage: L.<X> = LazySpecies(QQ)
            sage: E = L.Sets()
            sage: A = L.undefined(1)
            sage: A.define(X*E(A))
            sage: A[5]
            X*E_4 + X^2*E_3 + 3*X^3*E_2 + X*E_2(X^2) + 3*X^5

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: F = E(C(A))
            sage: [sum(F[n].monomial_coefficients().values()) for n in range(1, 7)]
            [1, 3, 7, 19, 47, 130]
            sage: oeis(_)  # optional -- internet
            0: A001372: Number of unlabeled mappings (or mapping patterns)
             from n points to themselves; number of unlabeled endofunctions.

            sage: R.<q> = QQ[]
            sage: L = LazySpecies(R, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E1 = L(lambda n: SymmetricGroup(n) if n else 0)
            sage: E(q*E1)[4]
            (q^4+q)*E_4 + q^2*E_2(E_2) + q^2*X*E_3 + q^3*E_2^2

        TESTS::

            sage: L.<X> = LazySpecies(QQ)
            sage: E2 = L(SymmetricGroup(2))
            sage: X(X + E2)
            X + E_2 + O^8
            sage: E2(X + E2)
            E_2 + X*E_2 + E_2(E_2) + O^9

            sage: (1+E2)(X)
            1 + E_2 + O^7

            sage: L.<X,Y> = LazySpecies(QQ)
            sage: X(Y, 0)
            Y + O^8

            sage: L1 = LazySpecies(QQ, "X")
            sage: E = L1.Sets()
            sage: L.<X,Y> = LazySpecies(QQ)
            sage: E(X)
            1 + X + E_2(X) + E_3(X) + E_4(X) + E_5(X) + E_6(X) + O^7

        It would be extremely nice to allow the following, but this
        poses theoretical problems::

            sage: L.<X> = LazySpecies(QQ)
            sage: E1 = L.Sets().restrict(1)
            sage: Omega = L.undefined(1)
            sage: L.define_implicitly([Omega], [E1(Omega) - X])
            sage: Omega[1]  # not tested
        """
        fP = self.parent()
        if len(args) != fP._arity:
            raise ValueError("arity of must be equal to the number of arguments provided")
        # Find a good parent for the result
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        P = cm.common_parent(self.base_ring(), *[parent(g) for g in args])
        # f = 0
        if isinstance(self._coeff_stream, Stream_zero):
            return P.zero()

        # args = (0, ..., 0)
        if all((not isinstance(g, LazyModuleElement) and not g)
               or (isinstance(g, LazyModuleElement)
                   and isinstance(g._coeff_stream, Stream_zero))
               for g in args):
            return P(self[0])

        # f is a constant polynomial
        if (isinstance(self._coeff_stream, Stream_exact)
            and not self._coeff_stream._constant
            and self.polynomial().is_constant()):
            return P(self.polynomial())

        return CompositionSpeciesElement(self, *args)

    def revert(self):
        r"""
        Return the compositional inverse of ``self``.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L.<X> = LazySpecies(QQ)
            sage: E1 = L.Sets().restrict(1)
            sage: g = E1.revert()
            sage: g[:5]
            [X, -E_2, -E_3 + X*E_2, -E_4 + E_2(E_2) + X*E_3 - X^2*E_2]

            sage: E = L.Sets()
            sage: P = E(X*E1(-X))*(1+X) - 1
            sage: P.revert()[:5]
            [X, X^2, X*E_2 + 2*X^3, X*E_3 + 2*X^2*E_2 + E_2(X^2) + 5*X^4]

        TESTS::

            sage: (3 + 2*X).revert()
            (-3/2) + 1/2*X
        """
        P = self.parent()
        if P._arity != 1:
            raise ValueError("arity must be equal to 1")
        coeff_stream = self._coeff_stream
        if isinstance(coeff_stream, Stream_zero):
            raise ValueError("compositional inverse does not exist")
        R = P._laurent_poly_ring
        if (isinstance(coeff_stream, Stream_exact)
            and coeff_stream.order() >= 0
            and coeff_stream._degree == 2):
            # self = a + b * X; self.revert() = -a/b + 1/b * X
            a = coeff_stream[0]
            b = coeff_stream[1].coefficients()[0]
            X = R(SymmetricGroup(1))  # as a polynomial species
            coeff_stream = Stream_exact((-a/b, 1/b * X),
                                        order=0)
            return P.element_class(P, coeff_stream)

        # TODO: coefficients should not be checked here, it prevents
        # us from using self.define in some cases!
        if coeff_stream[0]:
            raise ValueError("cannot determine whether the compositional inverse exists")

        X_mol = P._laurent_poly_ring._indices.subset(1)[0]  # as a molecular species
        X = P(SymmetricGroup(1))  # as a lazy species

        def coefficient(n):
            if n:
                return 0
            c = coeff_stream[1].coefficient(X_mol)
            if c.is_unit():
                return ~c
            raise ValueError("compositional inverse does not exist")

        b = P(lambda n: 0 if n else coeff_stream[1].coefficient(X_mol))  # TODO: we want a lazy version of Stream_exact
        b_inv = P(coefficient)  # TODO: we want a lazy version of Stream_exact
        g = P.undefined(valuation=1)
        g.define(b_inv * (X - (self - b * X)(g)))
        return g

    compositional_inverse = revert


class SumSpeciesElement(LazySpeciesElement):
    def __init__(self, left, right):
        r"""
        Initialize the sum of two species.
        """
        F = super(LazySpeciesElement, type(left))._add_(left, right)
        super().__init__(F.parent(), F._coeff_stream)
        self._left = left
        self._right = right

    def structures(self, *labels):
        labels = _label_sets(self.parent()._arity, labels)
        yield from self._left.structures(*labels)
        yield from self._right.structures(*labels)


class ProductSpeciesElement(LazySpeciesElement):
    def __init__(self, left, right):
        r"""
        Initialize the product of two species.
        """
        F = super(LazySpeciesElement, type(left))._mul_(left, right)
        super().__init__(F.parent(), F._coeff_stream)
        self._left = left
        self._right = right

    def structures(self, *labels):
        """
        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: C = L.Cycles()
            sage: P = E * C
            sage: list(P.structures([1,2]))
            [((), [1, 2]), ((1,), [2]), ((2,), [1])]

            sage: P = E * E
            sage: list(P.structures([1,2]))
            [((), (1, 2)), ((1,), (2,)), ((2,), (1,)), ((1, 2), ())]

            sage: L.<X, Y> = LazySpecies(QQ)
            sage: list((X*Y).structures([1], [2]))
            [((X, ((1,),)), (Y, ((2,),)))]
        """
        def dissections(s):
            for subset in subsets(s):
                subset_set = set(subset)
                yield (subset, tuple([e for e in s if e not in subset_set]))

        labels = _label_sets(self.parent()._arity, labels)
        for d in itertools.product(*[dissections(u) for u in labels]):
            yield from itertools.product(self._left.structures(*[U for U, _ in d]),
                                         self._right.structures(*[V for _, V in d]))


class CompositionSpeciesElement(LazySpeciesElement):
    def __init__(self, left, *args):
        r"""
        Initialize the composition of species.

        TESTS::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: P.<X> = LazySpecies(QQ)
            sage: P.zero()(X)
            0
            sage: X(P.zero())
            0
            sage: (1+X)(P.zero())
            1
        """
        fP = left.parent()
        # Find a good parent for the result
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        P = cm.common_parent(left.base_ring(), *[parent(g) for g in args])

        args = [P(g) for g in args]

        for g in args:
            if g._coeff_stream._approximate_order == 0:
                if not g._coeff_stream.is_uninitialized() and g[0]:
                    raise ValueError("can only compose with a positive valuation series")
                g._coeff_stream._approximate_order = 1

        sorder = left._coeff_stream._approximate_order
        gv = min(g._coeff_stream._approximate_order for g in args)
        R = P._internal_poly_ring.base_ring()
        L = fP._internal_poly_ring.base_ring()

        def coeff(g, i):
            c = g._coeff_stream[i]
            if not isinstance(c, PolynomialSpecies.Element):
                return R(c)
            return c

        # args_flat and weights contain one list for each g
        weight_exp = [lazy_list(lambda j, g=g: len(coeff(g, j+1)))
                      for g in args]

        def flat(g):
            # function needed to work around python's scoping rules
            return itertools.chain.from_iterable((coeff(g, j) for j in itertools.count()))

        args_flat1 = [lazy_list(flat(g)) for g in args]

        def coefficient(n):
            if not n:
                if left[0]:
                    return R(list(left[0])[0][1])
                return R.zero()
            result = R.zero()
            for i in range(1, n // gv + 1):
                # skip i=0 because it produces a term only for n=0

                # compute homogeneous components
                lF = defaultdict(L)
                for M, c in left[i]:
                    lF[M.grade()] += L._from_dict({M: c})
                for mc, F in lF.items():
                    for degrees in weighted_vector_compositions(mc, n, weight_exp):
                        args_flat = [list(a[0:len(degrees[j])]) for j, a in enumerate(args_flat1)]
                        multiplicities = [c for alpha, g_flat in zip(degrees, args_flat)
                                          for d, (_, c) in zip(alpha, g_flat) if d]
                        molecules = [M for alpha, g_flat in zip(degrees, args_flat)
                                          for d, (M, _) in zip(alpha, g_flat) if d]
                        non_zero_degrees = [[d for d in alpha if d] for alpha in degrees]
                        names = ["X%s" % i for i in range(len(molecules))]
                        FX = F._compose_with_weighted_singletons(names,
                                                                 multiplicities,
                                                                 non_zero_degrees)
                        FG = [(M(*molecules), c) for M, c in FX]
                        result += R.sum_of_terms(FG)
            return result

        coeff_stream = Stream_function(coefficient, P._sparse, sorder * gv)
        super().__init__(P, coeff_stream)
        self._left = left
        self._args = args

    def structures(self, *labels):
        r"""

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E1 = L.Sets().restrict(1)
            sage: sorted(E(E1).structures([1,2,3]))
            [((((1, 'X'),), ((2, 'X'),), ((3, 'X'),)), ((1,), (2,), (3,))),
             ((((1, 'X'),), ((2, 'X'), (3, 'X'))), ((1,), (2, 3))),
             ((((1, 'X'), (2, 'X')), ((3, 'X'),)), ((1, 2), (3,))),
             ((((1, 'X'), (2, 'X'), (3, 'X')),), ((1, 2, 3),)),
             ((((1, 'X'), (3, 'X')), ((2, 'X'),)), ((1, 3), (2,)))]

            sage: C = L.Cycles()
            sage: L.<X, Y> = LazySpecies(QQ)
            sage: sum(1 for s in C(X*Y).structures([1,2,3], [1,2,3]))
            12

            sage: C(X*Y).generating_series()[6]
            1/3*X^3*Y^3

            sage: sum(1 for s in E(X*Y).structures([1,2,3], ["a", "b", "c"]))
            6
        """
        F = self._left
        G = self._args
        m = len(G)  # == F.parent()._arity
        k = self.parent()._arity  # == G[i].parent()._arity
        names = self.parent()._laurent_poly_ring._indices._indices._names
        labels = _label_sets(k, labels)
        # make label sets disjoint
        U = [(e, i) for l, i in zip(labels, names) for e in l]

        def split_set(C):
            C_split = defaultdict(list)
            for e, i in C:
                C_split[i].append(e)
            return [C_split[i] for i in names]

        Par_U = SetPartitions(U)
        for pi in Par_U:
            # Fix an arbitrary order of the blocks
            try:
                pi_list = sorted([sorted(b) for b in pi])
            except TypeError:
                pi_list = sorted([sorted(b, key=str) for b in pi], key=str)

            # Generate all functions chi from pi to {0, ..., m-1}
            for chi in itertools.product(range(m), repeat=len(pi_list)):
                chi_inv = defaultdict(list)
                for b, i in zip(pi_list, chi):
                    chi_inv[i].append(b)

                # The set of structures is the Cartesian product of
                # the structures in F[chi_inv[i] for i in range(m)]
                # and for each set C in chi_inv[i] the set of
                # structures in G_i[C]
                F_s = F.structures(*[[tuple(b) for b in chi_inv[i]] for i in range(m)])
                G_s = [G[i].structures(*split_set(C)) for i in range(m) for C in chi_inv[i]]
                yield from itertools.product(F_s, itertools.product(*G_s))


class LazySpecies(LazyCompletionGradedAlgebra):
    """
    The ring of combinatorial species.
    """
    Element = LazySpeciesElement

    @staticmethod
    def __classcall_private__(cls, base_ring, names, sparse=True):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: LazySpecies(QQ, "X")
            Lazy completion of Polynomial species in X over Rational Field
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, names, sparse)

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.combinat.CombinatorialFreeModule`, because we do
        not want to implement `gens`.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: P.<X, Y> = LazySpecies(QQ)  # indirect doctest
            sage: 1/(1-X-Y)
            1 + (X+Y) + (X^2+2*X*Y+Y^2) + (X^3+3*X^2*Y+3*X*Y^2+Y^3)
             + (X^4+4*X^3*Y+6*X^2*Y^2+4*X*Y^3+Y^4)
             + (X^5+5*X^4*Y+10*X^3*Y^2+10*X^2*Y^3+5*X*Y^4+Y^5)
             + (X^6+6*X^5*Y+15*X^4*Y^2+20*X^3*Y^3+15*X^2*Y^4+6*X*Y^5+Y^6) + O^7
        """
        return tuple([self(g) for g in self._laurent_poly_ring._first_ngens(n)])

    def __init__(self, base_ring, names, sparse):
        r"""
        Initialize the ring of lazy species.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: LazySpecies(QQ, "X, Y")
            Lazy completion of Polynomial species in X, Y over Rational Field

            sage: L = LazySpecies(QQ, "X")
            sage: G = L.Graphs()
            sage: P = L.SetPartitions()
            sage: S = L.Sets()
            sage: C = L.Cycles()

        TESTS::

            sage: LazySpecies(QQ, "X, Y, Z")._arity
            3
        """
        super().__init__(PolynomialSpecies(base_ring, names))
        self._arity = len(names)
        if self._arity == 1:
            self.Graphs = lambda: GraphSpecies(self)
            self.SetPartitions = lambda: SetPartitionSpecies(self)
            self.Sets = lambda: SetSpecies(self)
            self.Cycles = lambda: CycleSpecies(self)


class SetSpecies(LazySpeciesElement):
    def __init__(self, parent):
        r"""
        Initialize the species of sets.
        """
        S = parent(lambda n: SymmetricGroup(n))
        super().__init__(parent, S._coeff_stream)

    def structures(self, *labels):
        """

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: list(E.structures([1,2,3]))
            [(1, 2, 3)]
        """
        labels = _label_sets(self.parent()._arity, labels)
        yield labels[0]


class CycleSpecies(LazySpeciesElement):
    def __init__(self, parent):
        r"""
        Initialize the species of cycles.
        """
        S = parent(lambda n: CyclicPermutationGroup(n) if n else 0)
        super().__init__(parent, S._coeff_stream)

    def structures(self, *labels):
        """

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: C = L.Cycles()
            sage: list(C.structures([]))
            []
            sage: list(C.structures([1]))
            [[1]]
            sage: list(C.structures([1,2]))
            [[1, 2]]
            sage: list(C.structures([1,2,3]))
            [[1, 2, 3], [1, 3, 2]]
        """
        labels = _label_sets(self.parent()._arity, labels)
        yield from CyclicPermutations(labels[0])


class GraphSpecies(LazySpeciesElement):

    def __init__(self, parent):
        r"""
        Initialize the species of simple graphs.
        """
        P = parent._laurent_poly_ring
        S = parent(lambda n: sum(P(G.automorphism_group()) for G in graphs(n)))
        super().__init__(parent, S._coeff_stream)

    def isotypes(self, labels):
        if labels in ZZ:
            yield from graphs(labels)


class SetPartitionSpecies(LazySpeciesElement):
    def __init__(self, parent):
        r"""
        Initialize the species of set partitions.
        """
        E = parent.Sets()
        E1 = parent.Sets().restrict(1)
        super().__init__(parent, E(E1)._coeff_stream)

    def isotypes(self, labels):
        if labels in ZZ:
            yield from Partitions(labels)

    def structures(self, labels):
        labels = _label_sets(self.parent()._arity, labels)
        yield from SetPartitions(labels)


class RestrictedSpeciesElement(LazySpeciesElement):
    def __init__(self, F, min_degree, max_degree):
        r"""
        Initialize the restriction of a species to the given degrees.
        """
        self._F = F
        self._min = min_degree
        self._max = max_degree

        if max_degree is None and min_degree is None:
            coeff_stream = F._coeff_stream
        elif max_degree is None:
            v = max(F._coeff_stream._approximate_order, min_degree)
            coeff_stream = Stream_truncated(F._coeff_stream, 0, v)
        else:
            if min_degree is None:
                v = F._coeff_stream._approximate_order
            else:
                v = max(F._coeff_stream._approximate_order, min_degree)
            initial_coefficients = [F._coeff_stream[i] for i in range(v, max_degree + 1)]
            if not any(initial_coefficients):
                coeff_stream = Stream_zero()
            else:
                coeff_stream = Stream_exact(initial_coefficients, order=v)

        super().__init__(F.parent(), coeff_stream)

    def isotypes(self, labels):
        if (labels in ZZ
            and (self._min is None or self._min <= labels)
            and (self._max is None or labels <= self._max)):
            yield from self._F.isotypes(labels)

    def structures(self, *labels):
        n = sum(map(len, labels))
        if ((self._min is None or self._min <= n)
            and (self._max is None or n <= self._max)):
            yield from self._F.structures(*labels)
