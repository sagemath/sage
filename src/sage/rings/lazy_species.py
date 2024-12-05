from sage.rings.integer_ring import ZZ
from sage.rings.lazy_series import LazyCompletionGradedAlgebraElement, LazyModuleElement
from sage.rings.lazy_series_ring import (LazyCompletionGradedAlgebra,
                                         LazyPowerSeriesRing,
                                         LazySymmetricFunctions)
from sage.rings.species import PolynomialSpecies, _label_sets
from sage.data_structures.stream import (Stream_zero,
                                         Stream_exact,
                                         Stream_function)
from sage.categories.sets_cat import cartesian_product
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


def weighted_compositions(n, d, weights, offset=0):
    r"""
    Return all compositions of `n` of weight `d`.

    The weight of a composition `n_1, n_2, \dots` is `\sum_i w_i
    n_i`.

    ``weights`` is assumed to be a weakly increasing list of positive
    integers.

    EXAMPLES::

        sage: from sage.rings.lazy_species import weighted_compositions
        sage: list(weighted_compositions(1, 1, [1,1,2]))
        [[0, 1], [1]]

        sage: list(weighted_compositions(2, 1, [1,1,2]))
        []

        sage: list(weighted_compositions(1, 2, [1,1,2,3]))
        [[0, 0, 1]]

        sage: list(weighted_compositions(3, 4, [1,1,2,2,3,3,4,4,5]))
        [[0, 2, 0, 1], [0, 2, 1], [1, 1, 0, 1], [1, 1, 1], [2, 0, 0, 1], [2, 0, 1]]

    """
    # the empty composition exists if and only if n == d == 0
    if not n:
        if not d:
            yield []
        return
    if not d:
        return

    # otherwise we iterate over the possibilities for the first part
    if offset < len(weights):
        w0 = weights[offset]
    else:
        return
    if w0 > d:
        return
    for i in range(min(n, d // w0) + 1):
        for c in weighted_compositions(n - i, d - i * w0, weights, offset=offset+1):
            yield [i] + c


def weighted_vector_compositions(n_vec, d, weights_vec):
    r"""
    Return all compositions of the vector `n` of weight `d`.

    INPUT:

    - ``n_vec``, a `k`-tuple of non-negative integers.

    - ``d``, a non-negative integer.

    - ``weights_vec``, `k`-tuple of weakly increasing lists of
      positive integers.

    EXAMPLES::

        sage: from sage.rings.lazy_species import weighted_vector_compositions
        sage: list(weighted_vector_compositions([1,1], 2, [[1,1,2,3], [1,2,3]]))
        [([0, 1], [1]), ([1], [1])]

        sage: list(weighted_vector_compositions([3,1], 4, [[1,1,2,5], [1,1,2,5]]))
        [([0, 3], [0, 1]),
         ([0, 3], [1]),
         ([1, 2], [0, 1]),
         ([1, 2], [1]),
         ([2, 1], [0, 1]),
         ([2, 1], [1]),
         ([3], [0, 1]),
         ([3], [1])]

    """
    k = len(n_vec)
    for d_vec in IntegerVectors(d, length=k):
        yield from itertools.product(*map(weighted_compositions, n_vec, d_vec, weights_vec))

######################################################################


class LazySpeciesElement(LazyCompletionGradedAlgebraElement):
    r"""

    EXAMPLES:

    Compute the molecular expansion of `E(-X)`::

        sage: from sage.rings.lazy_species import LazySpecies
        sage: L = LazySpecies(ZZ, "X")
        sage: E = L(lambda n: SymmetricGroup(n))
        sage: 1 / E
        1 + (-X) + (-E_2+X^2) + (-E_3+2*X*E_2-X^3)
          + (-E_4+2*X*E_3+E_2^2-3*X^2*E_2+X^4)
          + (-E_5+2*X*E_4+2*E_2*E_3-3*X^2*E_3-3*X*E_2^2+4*X^3*E_2-X^5)
          + (-E_6+2*X*E_5+2*E_2*E_4-3*X^2*E_4+E_3^2-6*X*E_2*E_3+4*X^3*E_3-E_2^3+6*X^2*E_2^2-5*X^4*E_2+X^6)
          + O^7

    Compare with the explicit formula::

        sage: def coefficient(m):
        ....:     return sum((-1)^len(la) * multinomial((n := la.to_exp())) * prod(E[i]^ni for i, ni in enumerate(n, 1)) for la in Partitions(m))

        sage: all(coefficient(m) == (1/E)[m] for m in range(10))
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

            sage: from sage.rings.species import PolynomialSpecies
            sage: L2 = LazySpecies(QQ, "X, Y")
            sage: P2 = PolynomialSpecies(QQ, "X, Y")
            sage: X = L2(P2(SymmetricGroup(1), {0: [1]}))
            sage: Y = L2(P2(SymmetricGroup(1), {1: [1]}))
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
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E.generating_series()
            1 + X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + O(X^7)

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: C.generating_series()
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4 + 1/5*X^5 + 1/6*X^6 + O(X^7)

            sage: from sage.rings.species import PolynomialSpecies
            sage: L2 = LazySpecies(QQ, "X, Y")
            sage: P2 = PolynomialSpecies(QQ, "X, Y")
            sage: X = L2(P2(SymmetricGroup(1), {0: [1]}))
            sage: Y = L2(P2(SymmetricGroup(1), {1: [1]}))
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
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: h = SymmetricFunctions(QQ).h()
            sage: LazySymmetricFunctions(h)(E.cycle_index_series())
            h[] + h[1] + h[2] + h[3] + h[4] + h[5] + h[6] + O^7

            sage: s = SymmetricFunctions(QQ).s()
            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: s(C.cycle_index_series()[5])
            s[1, 1, 1, 1, 1] + s[2, 2, 1] + 2*s[3, 1, 1] + s[3, 2] + s[5]

            sage: from sage.rings.species import PolynomialSpecies
            sage: L = LazySpecies(QQ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: L2 = LazySpecies(QQ, "X, Y")
            sage: P2 = PolynomialSpecies(QQ, "X, Y")
            sage: X = L2(P2(SymmetricGroup(1), {0: [1]}))
            sage: Y = L2(P2(SymmetricGroup(1), {1: [1]}))
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
        labels = _label_sets(self.parent()._arity, labels)
        F = self[sum(map(len, labels))]
        for M, c in F.monomial_coefficients().items():
            if c not in ZZ or c < 0:
                raise NotImplementedError("only implemented for proper non-virtual species")
            if c == 1:
                for s in M.structures(*labels):
                    yield M, s
            else:
                for e, s in cartesian_product([range(c), M.structures(*labels)]):
                    yield M, s, e

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
            P_4 + O^11

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
            X*E_4 + X^2*E_3 + 3*X^3*E_2 + X*{((1,2)(3,4),)} + 3*X^5

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
            (q^4+q)*E_4 + q^2*P_4 + q^2*X*E_3 + q^3*E_2^2

        TESTS::

            sage: L.<X> = LazySpecies(QQ)
            sage: E2 = L(SymmetricGroup(2))
            sage: X(X + E2)
            X + E_2 + O^8
            sage: E2(X + E2)
            E_2 + X*E_2 + P_4 + O^9

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
        """
        return CompositionSpeciesElement(self, *args)


class SumSpeciesElement(LazySpeciesElement):
    def __init__(self, left, right):
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
        fP = parent(left)
        if len(args) != fP._arity:
            raise ValueError("arity of must be equal to the number of arguments provided")

        # Find a good parent for the result
        from sage.structure.element import get_coercion_model
        cm = get_coercion_model()
        P = cm.common_parent(left.base_ring(), *[parent(g) for g in args])

        # f = 0
        if isinstance(left._coeff_stream, Stream_zero):
            return P.zero()

        # args = (0, ..., 0)
        if all((not isinstance(g, LazyModuleElement) and not g)
               or (isinstance(g, LazyModuleElement)
                   and isinstance(g._coeff_stream, Stream_zero))
               for g in args):
            return P(left[0])

        # f is a constant polynomial
        if (isinstance(left._coeff_stream, Stream_exact)
            and not left._coeff_stream._constant
            and left.polynomial().is_constant()):
            return P(left.polynomial())

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

        def coefficient(n):
            if not n:
                if left[0]:
                    return R(list(left[0])[0][1])
                return R.zero()
            args_flat = [[(M, c) for i in range(n+1) for M, c in g[i]]
                         for g in args]
            weights = [[sum(M.grade()) for i in range(n+1) for M, _ in g[i]]
                       for g in args]
            result = R.zero()
            for i in range(n // gv + 1):
                # compute homogeneous components
                lF = defaultdict(L)
                for M, c in left[i]:
                    lF[M.grade()] += L._from_dict({M: c})
                for mc, F in lF.items():
                    for degrees in weighted_vector_compositions(mc, n, weights):
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
            sage: E1 = L.Sets(min=1)
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
        """
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
            self.Sets = lambda min=0: SetSpecies(self, min)
            self.Cycles = lambda: CycleSpecies(self)


class SetSpecies(LazySpeciesElement):
    def __init__(self, parent, min):
        S = parent(lambda n: SymmetricGroup(n) if n >= min else 0)
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
        P = parent._laurent_poly_ring
        S = parent(lambda n: sum(P(G.automorphism_group()) for G in graphs(n)))
        super().__init__(parent, S._coeff_stream)

    def isotypes(self, labels):
        if labels in ZZ:
            yield from graphs(labels)


class SetPartitionSpecies(LazySpeciesElement):
    def __init__(self, parent):
        E = parent.Sets()
        E1 = parent.Sets(min=1)
        super().__init__(parent, E(E1)._coeff_stream)

    def isotypes(self, labels):
        if labels in ZZ:
            yield from Partitions(labels)

    def structures(self, labels):
        labels = _label_sets(self.parent()._arity, labels)
        yield from SetPartitions(labels)
