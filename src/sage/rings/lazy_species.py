r"""
Lazy Combinatorial Species

We regard a combinatorial species as a sequence of group actions of
the symmetric groups `\mathfrak S_n`, for `n\in\NN`.

Coefficients of lazy species are computed on demand.  They have
infinite precision, although equality can only be decided in special
cases.

AUTHORS:

- Mainak Roy, Martin Rubey, Travis Scrimshaw (2024-2025)

EXAMPLES:

We can reproduce the molecular expansions from Appendix B in
[GL2011]_ with little effort.  The molecular expansion of the
species of point determining graphs can be computed as the
species of graphs composed with the compositional inverse of the
species of non-empty sets::

    sage: L.<X> = LazyCombinatorialSpecies(QQ)
    sage: E = L.Sets()
    sage: Ep = E.restrict(1)
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
    sage: E_2 = L(SymmetricGroup(2))
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
from sage.arith.misc import divisors, multinomial
from sage.functions.other import binomial, factorial
from sage.misc.lazy_list import lazy_list
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.lazy_series import (LazyCompletionGradedAlgebraElement,
                                    LazyModuleElement)
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
from sage.combinat.partition import _Partitions, Partitions
from sage.combinat.permutation import CyclicPermutations
from sage.combinat.set_partition import SetPartitions
from sage.graphs.graph_generators import graphs
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup_named import (AlternatingGroup,
                                                  CyclicPermutationGroup,
                                                  DihedralGroup,
                                                  SymmetricGroup)
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import parent
from sage.structure.unique_representation import UniqueRepresentation
import itertools
from collections import defaultdict


def weighted_compositions(n, d, weight_multiplicities, _w0=0):
    r"""
    Return all compositions of `n` of weight `d`.

    The weight of a composition `n_1, n_2, \dots` is `\sum_i w_i n_i`.

    INPUT:

    - ``n`` -- nonnegative integer; the sum of the parts
    - ``d`` -- nonnegative integer; the total weight
    - ``weight_multiplicities`` -- iterable;
      ``weight_multiplicities[i]`` is the number of positions with
      weight ``i+1``

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
    from sage.combinat.integer_lists.invlex import IntegerListsBackend_invlex
    for s in range(n + 1):
        for c in weighted_compositions(n - s, d - s * (_w0 + 1), weight_multiplicities, _w0=_w0+1):
            m = weight_multiplicities[_w0]
            for v in IntegerListsBackend_invlex(s, length=m)._iter():
                yield v + c


def weighted_vector_compositions(n_vec, d, weight_multiplicities_vec):
    r"""
    Return all compositions of the vector `n` of weight `d`.

    INPUT:

    - ``n_vec`` -- a `k`-tuple of non-negative integers

    - ``d`` -- a non-negative integer, the total sum of the parts in
      all components

    - ``weight_multiplicities_vec`` -- `k`-tuple of iterables, where
      ``weight_multiplicities_vec[j][i]`` is the number of
      positions with weight `i+1` in the `j`-th component

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
    from sage.combinat.integer_lists.invlex import IntegerListsBackend_invlex
    for d_vec in IntegerListsBackend_invlex(d, length=k)._iter():
        yield from itertools.product(*map(weighted_compositions,
                                          n_vec, d_vec,
                                          weight_multiplicities_vec))

######################################################################


class LazyCombinatorialSpeciesElement(LazyCompletionGradedAlgebraElement):
    r"""
    EXAMPLES:

    Compute the molecular expansion of `E(-X)`::

        sage: L = LazyCombinatorialSpecies(ZZ, "X")
        sage: E = L(SymmetricGroup)
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

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L(SymmetricGroup)
            sage: E.isotype_generating_series()
            1 + X + X^2 + X^3 + X^4 + X^5 + X^6 + O(X^7)

            sage: C = L(CyclicPermutationGroup, valuation=1)
            sage: E(C).isotype_generating_series()
            1 + X + 2*X^2 + 3*X^3 + 5*X^4 + 7*X^5 + 11*X^6 + O(X^7)

            sage: L2.<X, Y> = LazyCombinatorialSpecies(QQ)
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
                return sum(self[n].coefficients())
        else:
            def coefficient(n):
                return sum(c * P.base_ring().prod(v ** d for v, d in zip(L.gens(), M.grade()))
                           for M, c in self[n].monomial_coefficients().items())
        return L(coefficient)

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E.generating_series()
            1 + X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + O(X^7)

            sage: C = L.Cycles()
            sage: C.generating_series()
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4 + 1/5*X^5 + 1/6*X^6 + 1/7*X^7 + O(X^8)

            sage: L2.<X, Y> = LazyCombinatorialSpecies(QQ)
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
            + (1/7*X^7+X^6*Y+3*X^5*Y^2+5*X^4*Y^3+5*X^3*Y^4+3*X^2*Y^5+X*Y^6+1/7*Y^7)
            + O(X,Y)^8
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

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: h = SymmetricFunctions(QQ).h()
            sage: LazySymmetricFunctions(h)(E.cycle_index_series())
            h[] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + O^7

            sage: s = SymmetricFunctions(QQ).s()
            sage: C = L.Cycles()
            sage: s(C.cycle_index_series()[5])
            s[1, 1, 1, 1, 1] + s[2, 2, 1] + 2*s[3, 1, 1] + s[3, 2] + s[5]

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: L2.<X, Y> = LazyCombinatorialSpecies(QQ)
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

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
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

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L(SymmetricGroup)
            sage: list(E.structures([1,2,3]))
            [(E_3, ((1, 2, 3),))]
            sage: list((E+E).structures([1,2,3]))
            [((E_3, ((1, 2, 3),)), 'left'), ((E_3, ((1, 2, 3),)), 'right')]
        """
        return SumSpeciesElement(self, other)

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L(SymmetricGroup)
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

        Generically, this yields a list of pairs consisting of a
        molecular species and a relabelled representative of the
        cosets of corresponding groups.

        The relabelling is such that the first few labels correspond
        to the first factor in the atomic decomposition, etc.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L(SymmetricGroup)
            sage: list(E.structures([1,2,3]))
            [(E_3, ((1, 2, 3),))]

            sage: C = L(CyclicPermutationGroup, valuation=1)
            sage: list(C.structures([1,2,3]))
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
            sage: L = LazyCombinatorialSpecies(QQ, "X, Y")
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

    def _test_structures(self, tester=None, max_size=5, **options):
        r"""
        Check that structures and generating series are consistent.

        We check all structures with less than ``max_size`` labels.

        TESTS::

            sage: from sage.rings.species import PolynomialSpecies
            sage: L = LazyCombinatorialSpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: XY = L(P(PermutationGroup([], domain=[1, 2, 3]), {0: [1], 1: [2, 3]}))
            sage: XY._test_structures()
        """
        if tester is None:
            tester = self._tester(**options)
        P = self.parent()
        for n in range(max_size):
            if P._arity == 1:
                labels = list(range(n))
                s = list(self.structures(labels))
                tester.assertEqual(len(s), len(set(s)),
                                   f"structures for {labels} are {s}, which is not a set")
                coeff = self.generating_series()[n]
                tester.assertEqual(len(s) / factorial(n), coeff,
                                   f"the number of structures for {labels} is {len(s)}, but the generating series gives {coeff}")
            else:
                label_shapes = IntegerVectors(n, length=P._arity)
                for shape in label_shapes:
                    labels = [list(range(k)) for k in shape]
                    s = list(self.structures(*labels))
                    tester.assertEqual(len(s), len(set(s)), f"structures for {labels} are {s}, which is not a set")
                    coeff = self.generating_series()[n].coefficient(list(shape))
                    tester.assertEqual(len(s) / ZZ.prod(factorial(k) for k in shape),
                                       coeff,
                                       f"the number of structures for {labels} is {len(s)}, but the generating series gives {coeff}")

    def isotypes(self, *shape):
        r"""
        Iterate over the isotypes on the given list of sizes.

        Generically, this yields a list of tuples consisting of a
        molecular species and, if necessary, an index.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L(SymmetricGroup)
            sage: list(E.isotypes(3))
            [(E_3,)]

            sage: P = L(CyclicPermutationGroup, valuation=1)
            sage: list(P.isotypes(3))
            [(C_3,)]

            sage: F = 1/(2-E)
            sage: sorted(F.isotypes(3))
            [(E_3,), (X*E_2, 0), (X*E_2, 1), (X^3,)]

            sage: from sage.rings.species import PolynomialSpecies
            sage: L = LazyCombinatorialSpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: XY = L(P(PermutationGroup([], domain=[1, 2]), {0: [1], 1: [2]}))
            sage: list((XY).isotypes(1, 1))
            [(X*Y,)]

            sage: list(E(XY).isotypes(2, 2))
            [(E_2(X*Y),)]
        """
        multivariate = self.parent()._arity > 1
        shape = tuple(shape)
        if not all(e in ZZ for e in shape):
            raise NotImplementedError("isotypes with given labels are currently not supported")
        for M, c in self[sum(shape)]:
            if c not in ZZ or c < 0:
                raise NotImplementedError("only implemented for proper non-virtual species")

            if multivariate and tuple(M.grade()) != shape:
                continue

            if c == 1:
                yield tuple([M])
            else:
                for e in range(c):
                    yield (M, e)

    def _test_isotypes(self, tester=None, max_size=5, **options):
        r"""
        Check that isotypes and generating series are consistent.

        TESTS::

            sage: from sage.rings.species import PolynomialSpecies
            sage: L = LazyCombinatorialSpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: XY = L(P(PermutationGroup([], domain=[1, 2]), {0: [1], 1: [2]}))
            sage: XY._test_isotypes(max_size=5)
        """
        if tester is None:
            tester = self._tester(**options)
        P = self.parent()
        for n in range(max_size):
            if P._arity == 1:
                s = list(self.isotypes(n))
                tester.assertEqual(len(s), len(set(s)),
                                   f"isotypes for {n} are {s}, which is not a set")
                coeff = self.isotype_generating_series()[n]
                tester.assertEqual(len(s), coeff,
                                   f"the number of isotypes for {n} is {len(s)}, but the generating series gives {coeff}")
            else:
                shapes = IntegerVectors(n, length=P._arity)
                for shape in shapes:
                    s = list(self.isotypes(*shape))
                    tester.assertEqual(len(s), len(set(s)),
                                       f"isotypes for {shape} are {s}, which is not a set")
                    coeff = self.isotype_generating_series()[n].coefficient(list(shape))
                    tester.assertEqual(len(s), coeff,
                                       f"the number of isotypes for {shape} is {len(s)}, but the generating series gives {coeff}")

    def polynomial(self, degree=None, names=None):
        r"""
        Return ``self`` as a polynomial if ``self`` is actually so or up to
        specified degree.

        INPUT:

        - ``degree`` -- (optional) integer
        - ``names`` -- (default: name of the variables of the series) names of the variables

        OUTPUT:

        If ``degree`` is not ``None``, the terms of the series of
        degree greater than ``degree`` are first truncated.  If
        ``degree`` is ``None`` and the series is not a polynomial
        polynomial, a ``ValueError`` is raised.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L(SymmetricGroup)
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
        Evaluate ``self`` at ``*args``.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: E2 = L(SymmetricGroup(2))
            sage: E2(E2)
            E_2(E_2) + O^11

            sage: E = L.Sets()
            sage: A = L.undefined(1)
            sage: A.define(X*E(A))
            sage: A[5]  # random
            X*E_4 + X^2*E_3 + 3*X^3*E_2 + X*E_2(X^2) + 3*X^5
            sage: A[5] == X*E[4] + X^2*E[3] + 3*X^3*E[2] + X*E[2](X[1]^2) + 3*X^5
            True

            sage: C = L.Cycles()
            sage: F = E(C(A))
            sage: [sum(F[n].monomial_coefficients().values()) for n in range(1, 7)]
            [1, 3, 7, 19, 47, 130]
            sage: oeis(_)  # optional -- internet
            0: A001372: Number of unlabeled mappings (or mapping patterns)
             from n points to themselves; number of unlabeled endofunctions.

            sage: R.<q> = QQ[]
            sage: L = LazyCombinatorialSpecies(R, "X")
            sage: E = L.Sets()
            sage: E1 = E.restrict(1)
            sage: E(q*E1)[4]
            (q^4+q)*E_4 + q^2*E_2(E_2) + q^2*X*E_3 + q^3*E_2^2

        TESTS::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: E2 = L(SymmetricGroup(2))
            sage: X(X + E2)
            X + E_2 + O^8
            sage: E2(X + E2)
            E_2 + X*E_2 + E_2(E_2) + O^9

            sage: from sage.rings.species import PolynomialSpecies
            sage: P = PolynomialSpecies(QQ, "X")
            sage: Gc = L(lambda n: sum(P(G.automorphism_group()) for G in graphs(n) if G.is_connected()) if n else 0)
            sage: E = L.Sets()
            sage: G = L.Graphs()
            sage: E(Gc) - G
            O^7

            sage: (1+E2)(X)
            1 + E_2 + O^7

            sage: L.<X,Y> = LazyCombinatorialSpecies(QQ)
            sage: X(Y, 0)
            Y + O^8

            sage: L1 = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L1.Sets()
            sage: L.<X,Y> = LazyCombinatorialSpecies(QQ)
            sage: E(X)
            1 + X + E_2(X) + E_3(X) + E_4(X) + E_5(X) + E_6(X) + O^7

        It would be extremely nice to allow the following, but this
        poses theoretical problems::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
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
            and self._coeff_stream._degree == 1):
            c = self._coeff_stream[0]
            B = c.parent()
            if B is ZZ or B is QQ or B == self.base_ring():
                return P(c)
            return P(c.coefficients()[0])

        return CompositionSpeciesElement(self, *args)

    def revert(self):
        r"""
        Return the compositional inverse of ``self``.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
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

    def combinatorial_logarithm(self):
        r"""
        Return the combinatorial logarithm of ``self``.

        This is the series reversion of the species of non-empty sets
        applied to ``self - 1``.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: L.Sets().restrict(1).revert() - (1+X).combinatorial_logarithm()
            O^7

        This method is much faster, however::

            sage: (1+X).combinatorial_logarithm().generating_series()[10]
            -1/10
        """
        P = self.parent()
        if P._arity != 1:
            raise ValueError("arity must be equal to 1")
        log = self.log()
        P1 = P._laurent_poly_ring
        M1 = P1._indices
        A1 = M1._indices

        def E(mu):
            return M1({A1(SymmetricGroup(e)): a
                       for e, a in enumerate(mu.to_exp(), 1) if a})

        def pi(mu):
            return (-1)**(len(mu)-1) * multinomial(mu.to_exp()) / len(mu)

        F = P.undefined()

        def coefficient(n):
            if not n:
                return 0
            res = log[n].monomial_coefficients()
            for k in divisors(n):
                if k == 1:
                    continue
                for mu in Partitions(k):
                    for N, g_N in F[n / k].monomial_coefficients().items():
                        M = E(mu)(N)
                        res[M] = res.get(M, 0) - pi(mu) * g_N
            return P1._from_dict(res)

        F.define(P(coefficient))
        return F


class LazyCombinatorialSpeciesElementGeneratingSeriesMixin:
    r"""
    A lazy species element whose generating series are obtained
    by specializing the cycle index series rather than the molecular
    expansion.

    TESTS:

    We check that the series are correct even if the cycle index
    series are not in the powersum basis::

        sage: from sage.rings.lazy_species import LazyCombinatorialSpeciesElement, LazyCombinatorialSpeciesElementGeneratingSeriesMixin
        sage: class F(LazyCombinatorialSpeciesElementGeneratingSeriesMixin, LazyCombinatorialSpeciesElement):
        ....:     def __init__(self, parent):
        ....:         super().__init__(parent, parent(PermutationGroup([], domain=[1,2])))
        ....:     def cycle_index_series(self):
        ....:         s = SymmetricFunctions(QQ).s()
        ....:         L = LazySymmetricFunctions(s)
        ....:         return L(s[1, 1] + s[2])

        sage: L = LazyCombinatorialSpecies(QQ, "X")
        sage: F(L).generating_series()
        X^2 + O(X^7)

        sage: F(L).isotype_generating_series()
        X^2 + O(X^7)
        sage: TestSuite(F(L)).run(skip=['_test_category', '_test_pickling'])


        sage: class F(LazyCombinatorialSpeciesElementGeneratingSeriesMixin, LazyCombinatorialSpeciesElement):
        ....:     def __init__(self, parent):
        ....:         G = PermutationGroup([], domain=[1,2,3,4])
        ....:         pi = {0:[1,2],1:[3,4]}
        ....:         P = parent._laurent_poly_ring
        ....:         super().__init__(parent, parent(P(G, pi)))
        ....:     def cycle_index_series(self):
        ....:         s = SymmetricFunctions(QQ).s()
        ....:         L = LazySymmetricFunctions(tensor([s, s]))
        ....:         return L(self[4].support()[0].cycle_index())

        sage: L = LazyCombinatorialSpecies(QQ, "X, Y")
        sage: F(L).isotype_generating_series()
        X^2*Y^2 + O(X,Y)^7

        sage: F(L).generating_series()
        X^2*Y^2 + O(X,Y)^7

        sage: TestSuite(F(L)).run(skip=['_test_category', '_test_pickling'])
    """
    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        The series is obtained by applying the principal
        specialization of order `1` to the cycle index series, that
        is, setting `x_1 = x` and `x_k = 0` for `k > 1`.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Graphs().isotype_generating_series().truncate(8)
            1 + X + 2*X^2 + 4*X^3 + 11*X^4 + 34*X^5 + 156*X^6 + 1044*X^7

        TESTS::

            sage: L.Graphs().isotype_generating_series()[20]
            645490122795799841856164638490742749440
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        cis = self.cycle_index_series()
        one = ZZ.one()

        if P._arity == 1:
            return L(lambda n: cis[n].principal_specialization(one, one))

        vars = L._laurent_poly_ring.gens()
        parents = cis.parent()._laurent_poly_ring.tensor_factors()

        def coefficient(n):
            return sum(c * prod(S(la).principal_specialization(one, one)
                                * v**la.size()
                                for v, S, la in zip(vars, parents, M))
                       for M, c in cis[n].monomial_coefficients().items())

        return L(coefficient)

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        The series is obtained by applying the exponential
        specialization to the cycle index series.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Graphs().generating_series().truncate(7)
            1 + X + X^2 + 4/3*X^3 + 8/3*X^4 + 128/15*X^5 + 2048/45*X^6
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        cis = self.cycle_index_series()
        one = ZZ.one()

        if P._arity == 1:
            return L(lambda n: cis[n].exponential_specialization(one, one))

        vars = L._laurent_poly_ring.gens()
        parents = cis.parent()._laurent_poly_ring.tensor_factors()

        def coefficient(n):
            return sum(c * prod(S(la).exponential_specialization(one, one)
                                * v**la.size()
                                for v, S, la in zip(vars, parents, M))
                       for M, c in cis[n].monomial_coefficients().items())

        return L(coefficient)


class SumSpeciesElement(LazyCombinatorialSpeciesElement):
    def __init__(self, left, right):
        r"""
        Initialize the sum of two species.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: F = L.Sets() + L.SetPartitions()
            sage: TestSuite(F).run(skip=['_test_category', '_test_pickling'])
        """
        F = super(LazyCombinatorialSpeciesElement, type(left))._add_(left, right)
        super().__init__(F.parent(), F._coeff_stream)
        self._left = left
        self._right = right

    def structures(self, *labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: F = L.Sets() + L.SetPartitions()
            sage: list(F.structures([1,2,3]))
            [((1, 2, 3), 'left'),
             ({{1, 2, 3}}, 'right'),
             ({{1, 2}, {3}}, 'right'),
             ({{1, 3}, {2}}, 'right'),
             ({{1}, {2, 3}}, 'right'),
             ({{1}, {2}, {3}}, 'right')]
        """
        labels = _label_sets(self.parent()._arity, labels)
        yield from ((s, 'left') for s in self._left.structures(*labels))
        yield from ((s, 'right') for s in self._right.structures(*labels))

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: F = L.Sets() + L.SetPartitions()
            sage: F.generating_series()
            2 + 2*X + 3/2*X^2 + X^3 + 2/3*X^4 + 53/120*X^5 + 17/60*X^6 + O(X^7)

        TESTS::

            sage: F.generating_series()[20]
            3978781402721/187146308321280000
        """
        return self._left.generating_series() + self._right.generating_series()

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: F = L.Sets() + L.SetPartitions()
            sage: F.isotype_generating_series()
            2 + 2*X + 3*X^2 + 4*X^3 + 6*X^4 + 8*X^5 + 12*X^6 + O(X^7)

        TESTS::

            sage: F.isotype_generating_series()[20]
            628
        """
        return self._left.isotype_generating_series() + self._right.isotype_generating_series()


class ProductSpeciesElement(LazyCombinatorialSpeciesElement):
    def __init__(self, left, right):
        r"""
        Initialize the product of two species.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: F = L.Sets() * L.SetPartitions()
            sage: TestSuite(F).run(skip=['_test_category', '_test_pickling'])
        """
        F = super(LazyCombinatorialSpeciesElement, type(left))._mul_(left, right)
        super().__init__(F.parent(), F._coeff_stream)
        self._left = left
        self._right = right

    def structures(self, *labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: C = L.Cycles()
            sage: P = E * C
            sage: list(P.structures([1,2]))
            [((), (1, 2)), ((1,), (2,)), ((2,), (1,))]

            sage: P = E * E
            sage: list(P.structures([1,2]))
            [((), (1, 2)), ((1,), (2,)), ((2,), (1,)), ((1, 2), ())]

            sage: L.<X, Y> = LazyCombinatorialSpecies(QQ)
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

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: F = E*E
            sage: F.generating_series()
            1 + 2*X + 2*X^2 + 4/3*X^3 + 2/3*X^4 + 4/15*X^5 + 4/45*X^6 + O(X^7)
        """
        return self._left.generating_series() * self._right.generating_series()

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: F = E*E
            sage: F.isotype_generating_series()
            1 + 2*X + 3*X^2 + 4*X^3 + 5*X^4 + 6*X^5 + 7*X^6 + O(X^7)
        """
        return self._left.isotype_generating_series() * self._right.isotype_generating_series()


class CompositionSpeciesElement(LazyCombinatorialSpeciesElementGeneratingSeriesMixin,
                                LazyCombinatorialSpeciesElement):
    def __init__(self, left, *args):
        r"""
        Initialize the composition of species.

        TESTS::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: L.zero()(X)
            0
            sage: X(L.zero())
            0
            sage: (1+X)(L.zero())
            1

            sage: L2.<X,Y> = LazyCombinatorialSpecies(QQ)
            sage: F = L.Sets()(X + 2*Y)
            sage: TestSuite(F).run(skip=['_test_category', '_test_pickling'])
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
            return itertools.chain.from_iterable(coeff(g, j) for j in itertools.count())

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
                        args_flat = [list(a[0:len(degrees[j])])
                                     for j, a in enumerate(args_flat1)]
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
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E1 = L.Sets().restrict(1)
            sage: sorted(E(E1).structures([1,2,3]))
            [((((1, 'X'),), ((2, 'X'),), ((3, 'X'),)), ((1,), (2,), (3,))),
             ((((1, 'X'),), ((2, 'X'), (3, 'X'))), ((1,), (2, 3))),
             ((((1, 'X'), (2, 'X')), ((3, 'X'),)), ((1, 2), (3,))),
             ((((1, 'X'), (2, 'X'), (3, 'X')),), ((1, 2, 3),)),
             ((((1, 'X'), (3, 'X')), ((2, 'X'),)), ((1, 3), (2,)))]

            sage: C = L.Cycles()
            sage: L.<X, Y> = LazyCombinatorialSpecies(QQ)
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

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: F = E(E.restrict(1))
            sage: F.generating_series()
            1 + X + X^2 + 5/6*X^3 + 5/8*X^4 + 13/30*X^5 + 203/720*X^6 + O(X^7)
        """
        return self._left.generating_series()(*[G.generating_series() for G in self._args])

    def cycle_index_series(self):
        r"""
        Return the cycle index series for this species.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: F = E(E.restrict(1))
            sage: F.cycle_index_series()[5]
            h[2, 2, 1] - h[3, 1, 1] + 3*h[3, 2] + 2*h[4, 1] + 2*h[5]
        """
        return self._left.cycle_index_series()(*[G.cycle_index_series() for G in self._args])


class LazyCombinatorialSpecies(LazyCompletionGradedAlgebra):
    Element = LazyCombinatorialSpeciesElement

    @staticmethod
    def __classcall_private__(cls, base_ring, names, sparse=True):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: LazyCombinatorialSpecies(QQ, "X") is LazyCombinatorialSpecies(QQ, "X")
            True
        """
        from sage.structure.category_object import normalize_names
        names = normalize_names(-1, names)
        if len(names) == 1:
            return LazyCombinatorialSpeciesUnivariate(base_ring, names, sparse)
        return LazyCombinatorialSpeciesMultivariate(base_ring, names, sparse)

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.combinat.CombinatorialFreeModule`, because we do
        not want to implement `gens`.

        EXAMPLES::

            sage: L.<X, Y> = LazyCombinatorialSpecies(QQ)  # indirect doctest
            sage: 1/(1-X-Y)
            1 + (X+Y) + (X^2+2*X*Y+Y^2) + (X^3+3*X^2*Y+3*X*Y^2+Y^3)
             + (X^4+4*X^3*Y+6*X^2*Y^2+4*X*Y^3+Y^4)
             + (X^5+5*X^4*Y+10*X^3*Y^2+10*X^2*Y^3+5*X*Y^4+Y^5)
             + (X^6+6*X^5*Y+15*X^4*Y^2+20*X^3*Y^3+15*X^2*Y^4+6*X*Y^5+Y^6) + O^7
        """
        return tuple([self(g) for g in self._laurent_poly_ring._first_ngens(n)])

    def __init__(self, base_ring, names, sparse):
        r"""
        The ring of lazy species.

        EXAMPLES:

        We provide univariate and multivariate (mostly known as
        multisort) species::

            sage: LazyCombinatorialSpecies(QQ, "X")
            Lazy completion of Polynomial species in X over Rational Field

            sage: LazyCombinatorialSpecies(QQ, "X, Y")
            Lazy completion of Polynomial species in X, Y over Rational Field

        In the univariate case, several basic species are provided as
        methods::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Sets()
            Set species
            sage: L.Cycles()
            Cycle species
            sage: L.OrientedSets()
            Oriented Set species
            sage: L.Polygons()
            Polygon species
            sage: L.Graphs()
            Graph species
            sage: L.SetPartitions()
            Set Partition species

        TESTS::

            sage: LazyCombinatorialSpecies(QQ, "X, Y, Z")._arity
            3
        """
        super().__init__(PolynomialSpecies(base_ring, names))
        self._arity = len(names)


class LazyCombinatorialSpeciesUnivariate(LazyCombinatorialSpecies):
    def Sets(self):
        r"""
        Return the species of sets.

        This species corresponds to the sequence of trivial group
        actions.  Put differently, the stabilizers are the full
        symmetric groups.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Sets()
            sage: set(G.isotypes(4))
            {(E_4,)}
            sage: set(G.structures(["a", 1, x]))
            {(1, 'a', x)}
        """
        return SetSpecies(self)

    def Cycles(self):
        r"""
        Return the species of (oriented) cycles.

        This species corresponds to the sequence of group actions
        having the cyclic groups as stabilizers.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Cycles()
            sage: set(G.isotypes(4))
            {(C_4,)}
            sage: set(G.structures(["a", "b", "c"]))
            {('a', 'b', 'c'), ('a', 'c', 'b')}
        """
        return CycleSpecies(self)

    def Polygons(self):
        r"""
        Return the species of polygons.

        Polygons are cycles up to orientation.

        This species corresponds to the sequence of group actions
        having the dihedral groups as stabilizers.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Polygons()
            sage: set(G.isotypes(5))
            {(P_5,)}
            sage: set(G.structures(["a", 1, "b", 2]))
            {(E_2(E_2), ((1, 'a', 2, 'b'),)),
             (E_2(E_2), ((1, 'b', 2, 'a'),)),
             (E_2(E_2), ((1, 2, 'a', 'b'),))}
        """
        return PolygonSpecies(self)

    def OrientedSets(self):
        r"""
        Return the species of oriented sets.

        Oriented sets are total orders up to an even orientation.

        This species corresponds to the sequence of group actions
        having the alternating groups as stabilizers.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.OrientedSets()
            sage: set(G.isotypes(5))
            {(Eo_5,)}
            sage: set(G.structures(["a", 1, "b", 2]))
            {(Eo_4, ((1, 2, 'a', 'b'),)), (Eo_4, ((1, 2, 'b', 'a'),))}
        """
        return OrientedSetSpecies(self)

    def Chains(self):
        r"""
        Return the species of chains.

        Chains are linear orders up to reversal.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: Ch = L.Chains()
            sage: set(Ch.isotypes(4))
            {(E_2(X^2),)}
            sage: list(Ch.structures(["a", "b", "c"]))
            [('a', 'c', 'b'), ('a', 'b', 'c'), ('b', 'a', 'c')]
        """
        return ChainSpecies(self)

    def Graphs(self):
        r"""
        Return the species of vertex labelled simple graphs.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Graphs()
            sage: set(G.isotypes(2))
            {Graph on 2 vertices, Graph on 2 vertices}

            sage: G.isotype_generating_series()[20]
            645490122795799841856164638490742749440
        """
        return GraphSpecies(self)

    def SetPartitions(self):
        r"""
        Return the species of set partitions.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.SetPartitions()
            sage: set(G.isotypes(4))
            {[1, 1, 1, 1], [2, 1, 1], [2, 2], [3, 1], [4]}
            sage: list(G.structures(["a", "b", "c"]))
            [{{'a', 'b', 'c'}},
             {{'a', 'b'}, {'c'}},
             {{'a', 'c'}, {'b'}},
             {{'a'}, {'b', 'c'}},
             {{'a'}, {'b'}, {'c'}}]
        """
        return SetPartitionSpecies(self)


class LazyCombinatorialSpeciesMultivariate(LazyCombinatorialSpecies):
    pass


class SetSpecies(LazyCombinatorialSpeciesElement, UniqueRepresentation,
                 metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of sets.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: TestSuite(E).run(skip=['_test_category', '_test_pickling'])

            sage: E is L.Sets()
            True
        """
        S = parent(SymmetricGroup)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").Sets()  # indirect doctest
           Set species
        """
        return "Set species"

    def structures(self, labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: E = L.Sets()
            sage: list(E.structures([1,2,3]))
            [(1, 2, 3)]
        """
        labels = _label_sets(self.parent()._arity, [labels])
        yield labels[0]

    def generating_series(self):
        r"""
        Return the (exponential) generating series of the
        species of sets.

        This is the exponential.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: L.Sets().generating_series()
            1 + X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + O(X^7)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L.gen().exp()

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of the species of
        sets.

        This is the geometric series.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Sets().isotype_generating_series()
            1 + X + X^2 + O(X^3)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(constant=1)

    def cycle_index_series(self):
        r"""
        Return the cycle index series of the species of sets.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: L.Sets().cycle_index_series()
            h[] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + O^7
        """
        P = self.parent()
        h = SymmetricFunctions(P.base_ring().fraction_field()).h()
        L = LazySymmetricFunctions(h)
        return L(lambda n: h[n])


class CycleSpecies(LazyCombinatorialSpeciesElement, UniqueRepresentation,
                   metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of cycles.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: C = L.Cycles()
            sage: TestSuite(C).run(skip=['_test_category', '_test_pickling'])

            sage: C is L.Cycles()
            True
        """
        S = parent(CyclicPermutationGroup, valuation=1)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").Cycles()  # indirect doctest
           Cycle species
        """
        return "Cycle species"

    def structures(self, labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: C = L.Cycles()
            sage: list(C.structures([]))
            []
            sage: list(C.structures([1]))
            [(1,)]
            sage: list(C.structures([1,2]))
            [(1, 2)]
            sage: list(C.structures([1,2,3]))
            [(1, 2, 3), (1, 3, 2)]
        """
        labels = _label_sets(self.parent()._arity, [labels])
        # TODO: CyclicPermutations should yield hashable objects, not lists
        yield from map(tuple, CyclicPermutations(labels[0]))

    def generating_series(self):
        r"""
        Return the (exponential) generating series of the
        species of cycles.

        This is `-log(1-x)`.

        EXAMPLES::

            sage: L.<X> = LazyCombinatorialSpecies(QQ)
            sage: L.Cycles().generating_series()
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4 + 1/5*X^5 + 1/6*X^6 + 1/7*X^7 + O(X^8)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return -(L.one()-L.gen()).log()

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of the species of
        cycles.

        This is `x/(1-x)`.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Cycles().isotype_generating_series()
            X + X^2 + X^3 + O(X^4)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(constant=1, valuation=1)


class PolygonSpecies(LazyCombinatorialSpeciesElement, UniqueRepresentation,
                     metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of polygons.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: P = L.Polygons()
            sage: TestSuite(P).run(skip=['_test_category', '_test_pickling'])

            sage: P is L.Polygons()
            True
        """
        S = parent(DihedralGroup, valuation=3)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").Polygons()  # indirect doctest
           Polygon species
        """
        return "Polygon species"


class OrientedSetSpecies(LazyCombinatorialSpeciesElement, UniqueRepresentation,
                         metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of polygons.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: Eo = L.OrientedSets()
            sage: TestSuite(Eo).run(skip=['_test_category', '_test_pickling'])

            sage: Eo is L.OrientedSets()
            True
        """
        S = parent(AlternatingGroup, valuation=4)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").OrientedSets()  # indirect doctest
           Oriented Set species
        """
        return "Oriented Set species"


class ChainSpecies(LazyCombinatorialSpeciesElement, UniqueRepresentation,
                   metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of chains.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: Ch = L.Chains()
            sage: TestSuite(Ch).run(skip=['_test_category', '_test_pickling'])

            sage: Ch is L.Chains()
            True
        """
        P = parent._laurent_poly_ring

        def coefficient(n):
            if not n:
                return P.one()
            if n % 2:
                gen = [(i, i+1) for i in range(2, n+1, 2)]
            else:
                gen = [(i, i+1) for i in range(1, n+1, 2)]
            return P(PermutationGroup([gen]))

        S = parent(coefficient)
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").Chains()  # indirect doctest
           Chain species
        """
        return "Chain species"

    def structures(self, labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: Ch = L.Chains()
            sage: list(Ch.structures([1,2,3]))
            [(1, 3, 2), (1, 2, 3), (2, 1, 3)]
        """
        labels = _label_sets(self.parent()._arity, [labels])[0]
        n = len(labels)
        if not n:
            yield ()
        elif n == 1:
            yield labels
        else:
            for a, b in itertools.combinations(labels, 2):
                ia = labels.index(a)
                ib = labels.index(b)
                rest = labels[:ia] + labels[ia+1:ib] + labels[ib+1:]
                for pi in itertools.permutations(rest):
                    yield (a,) + pi + (b,)


class GraphSpecies(LazyCombinatorialSpeciesElementGeneratingSeriesMixin,
                   LazyCombinatorialSpeciesElement, UniqueRepresentation,
                   metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of simple graphs.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Graphs()
            sage: TestSuite(G).run(skip=['_test_category', '_test_pickling'])

            sage: G is L.Graphs()
            True
        """
        P = parent._laurent_poly_ring
        S = parent(lambda n: sum(P(G.automorphism_group()) for G in graphs(n)))
        super().__init__(parent, S._coeff_stream)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").Graphs()  # indirect doctest
           Graph species
        """
        return "Graph species"

    def isotypes(self, labels):
        r"""
        Iterate over the isotypes on the given list of sizes.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G = L.Graphs()
            sage: list(G.isotypes(2))
            [Graph on 2 vertices, Graph on 2 vertices]
        """
        if labels in ZZ:
            yield from (G.canonical_label().copy(immutable=True) for G in graphs(labels))
        else:
            raise NotImplementedError("isotypes with given labels are currently not supported")

    def generating_series(self):
        r"""
        Return the (exponential) generating series of the
        species of simple graphs.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Graphs().generating_series().truncate(7)
            1 + X + X^2 + 4/3*X^3 + 8/3*X^4 + 128/15*X^5 + 2048/45*X^6
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(lambda n: 2**binomial(n, 2) / factorial(n))

    def cycle_index_series(self):
        r"""
        Return the cycle index series of the species of simple graphs.

        The cycle index series is computed using Proposition 2.2.7 in
        [BLL1998]_.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: L.Graphs().cycle_index_series().truncate(4)
            p[] + p[1] + (p[1,1]+p[2]) + (4/3*p[1,1,1]+2*p[2,1]+2/3*p[3])

        Check that the number of isomorphism types is computed quickly::

            sage: L.Graphs().isotype_generating_series()[20]
            645490122795799841856164638490742749440
        """
        P = self.parent()
        p = SymmetricFunctions(P.base_ring().fraction_field()).p()
        L = LazySymmetricFunctions(p)

        def a(sigma):
            rho = sigma.to_exp()
            res1 = ZZ.sum(ZZ(i+1)._gcd(ZZ(j+1)) * rho[i] * rho[j]
                          for i in range(len(rho))
                          for j in range(i+1, len(rho)))
            res2 = ZZ.sum(ZZ(i+1) * rho[i]**2
                          for i in range(len(rho)))
            res3 = ZZ.sum(rho[::2])
            return ZZ(2) ** (res1 + (res2 - res3) / 2) / sigma.centralizer_size()

        def coefficient(n):
            return p._from_dict({sigma: a(sigma) for sigma in Partitions(n)})

        return L(coefficient)


class SetPartitionSpecies(CompositionSpeciesElement, UniqueRepresentation,
                          metaclass=InheritComparisonClasscallMetaclass):
    def __init__(self, parent):
        r"""
        Initialize the species of set partitions.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: p = L.SetPartitions()
            sage: TestSuite(p).run(skip=['_test_category', '_test_pickling'])

            sage: p is L.SetPartitions()
            True

            sage: p.generating_series()[20]
            263898766507/12412765347840000

            sage: SetPartitions(20).cardinality() / factorial(20)
            263898766507/12412765347840000

            sage: p.isotype_generating_series()[20]
            627

            sage: Partitions(20).cardinality()
            627
        """
        E = parent.Sets()
        E1 = parent.Sets().restrict(1)
        super().__init__(E, E1)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

           sage: LazyCombinatorialSpecies(QQ, "X").SetPartitions()  # indirect doctest
           Set Partition species
        """
        return "Set Partition species"

    def isotypes(self, labels):
        r"""
        Iterate over the isotypes on the given list of sizes.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: p = L.SetPartitions()
            sage: list(p.isotypes(3))
            [[3], [2, 1], [1, 1, 1]]
        """
        if labels in ZZ:
            yield from Partitions(labels)
        else:
            raise NotImplementedError("isotypes with given labels are currently not supported")

    def structures(self, labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: P = L.SetPartitions()
            sage: list(P.structures([]))
            [{}]
            sage: list(P.structures([1]))
            [{{1}}]
            sage: list(P.structures([1,2]))
            [{{1, 2}}, {{1}, {2}}]
            sage: list(P.structures([1,2,3]))
            [{{1, 2, 3}}, {{1, 2}, {3}}, {{1, 3}, {2}}, {{1}, {2, 3}}, {{1}, {2}, {3}}]
        """
        labels = _label_sets(self.parent()._arity, [labels])
        yield from SetPartitions(labels[0])

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: P = L.SetPartitions()
            sage: P.generating_series()
            1 + X + X^2 + 5/6*X^3 + 5/8*X^4 + 13/30*X^5 + 203/720*X^6 + O(X^7)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(lambda n: SetPartitions(n).cardinality() / factorial(n))

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: P = L.SetPartitions()
            sage: P.isotype_generating_series()
            1 + X + 2*X^2 + 3*X^3 + 5*X^4 + 7*X^5 + 11*X^6 + O(X^7)
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        return L(lambda n: Partitions(n).cardinality())


class RestrictedSpeciesElement(LazyCombinatorialSpeciesElement):
    def __init__(self, F, min_degree, max_degree):
        r"""
        Initialize the restriction of a species to the given degrees.

        TESTS::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: G3 = L.Graphs().restrict(3, 3)
            sage: TestSuite(G3).run(skip=['_test_category', '_test_pickling'])
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

    def isotypes(self, *shape):
        r"""
        Iterate over the isotypes on the given list of sizes.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: p = L.SetPartitions().restrict(2, 2)
            sage: list(p.isotypes(3))
            []
        """
        n = sum(shape)
        if ((self._min is None or self._min <= n)
            and (self._max is None or n <= self._max)):
            yield from self._F.isotypes(*shape)

    def structures(self, *labels):
        r"""
        Iterate over the structures on the given set of labels.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(ZZ, "X")
            sage: F = L.SetPartitions().restrict(3)
            sage: list(F.structures([1]))
            []
            sage: list(F.structures([1,2,3]))
            [{{1, 2, 3}}, {{1, 2}, {3}}, {{1, 3}, {2}}, {{1}, {2, 3}}, {{1}, {2}, {3}}]
        """
        n = sum(map(len, labels))
        if ((self._min is None or self._min <= n)
            and (self._max is None or n <= self._max)):
            yield from self._F.structures(*labels)

    def generating_series(self):
        r"""
        Return the (exponential) generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E.restrict(1, 5).generating_series()
            X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5
            sage: E.restrict(1).generating_series()
            X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + 1/5040*X^7 + O(X^8)
        """
        return self._F.generating_series().restrict(self._min, self._max)

    def isotype_generating_series(self):
        r"""
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E.restrict(1, 5).isotype_generating_series()
            X + X^2 + X^3 + X^4 + X^5

            sage: E.restrict(1).isotype_generating_series()
            X + X^2 + X^3 + O(X^4)
        """
        return self._F.isotype_generating_series().restrict(self._min, self._max)

    def cycle_index_series(self):
        r"""
        Return the cycle index series for this species.

        EXAMPLES::

            sage: L = LazyCombinatorialSpecies(QQ, "X")
            sage: E = L.Sets()
            sage: E.restrict(1, 5).cycle_index_series()
            h[1] + h[2] + h[3] + h[4] + h[5]

            sage: E.restrict(1).cycle_index_series()
            h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + h[7] + O^8
        """
        return self._F.cycle_index_series().restrict(self._min, self._max)
