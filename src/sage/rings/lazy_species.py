from sage.rings.lazy_series import LazyCompletionGradedAlgebraElement, LazyModuleElement
from sage.rings.lazy_series_ring import LazyCompletionGradedAlgebra
from sage.data_structures.stream import (Stream_zero,
                                         Stream_exact,
                                         Stream_function)
from sage.rings.species import PolynomialSpecies
from sage.libs.gap.libgap import libgap
from sage.categories.sets_cat import cartesian_product
from sage.combinat.integer_vector import IntegerVectors
from sage.structure.element import parent
import itertools
from collections import defaultdict
from sage.rings.lazy_series_ring import LazyPowerSeriesRing, LazySymmetricFunctions
from sage.combinat.sf.sf import SymmetricFunctions

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

            sage: L2 = LazySpecies(QQ, "X, Y")
            sage: P2 = PolynomialSpecies(QQ, "X, Y")
            sage: X = L2(P2(SymmetricGroup(1), {1: [1]}))
            sage: Y = L2(P2(SymmetricGroup(1), {2: [1]}))
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
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E.generating_series()
            1 + X + 1/2*X^2 + 1/6*X^3 + 1/24*X^4 + 1/120*X^5 + 1/720*X^6 + O(X^7)

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: C.generating_series()
            X + 1/2*X^2 + 1/3*X^3 + 1/4*X^4 + 1/5*X^5 + 1/6*X^6 + O(X^7)

            sage: L2 = LazySpecies(QQ, "X, Y")
            sage: P2 = PolynomialSpecies(QQ, "X, Y")
            sage: X = L2(P2(SymmetricGroup(1), {1: [1]}))
            sage: Y = L2(P2(SymmetricGroup(1), {2: [1]}))
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
            + (1/6*X^6+X^5*Y+3*X^4*Y^2+4*X^3*Y^3+3*X^2*Y^4+X*Y^5+1/6*Y^6)
            + O(X,Y)^7
        """
        P = self.parent()
        L = LazyPowerSeriesRing(P.base_ring().fraction_field(),
                                P._laurent_poly_ring._indices._indices.variable_names())
        if P._arity == 1:
            def coefficient(n):
                return sum(c / M._group.cardinality()
                           for M, c in self[n].monomial_coefficients().items())
        else:
            def coefficient(n):
                return sum(c / M._group.cardinality()
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
        """
        P = self.parent()
        p = SymmetricFunctions(P.base_ring().fraction_field()).p()
        if P._arity == 1:
            L = LazySymmetricFunctions(p)
            def coefficient(n):
                return sum(c * M._group.cycle_index()
                           for M, c in self[n].monomial_coefficients().items())
        else:
            raise NotImplementedError
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
            [((1, 2, 3), E_3)]
            sage: list((E+F).structures([1,2,3]))
            [((1, 2, 3), E_3), ((1, 2, 3), E_3)]

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
            sage: list((E^2).structures([1,2,3]))
            [(((), 1), ((1, 2, 3), E_3)),
             (((1,), X), ((2, 3), E_2)),
             (((2,), X), ((1, 3), E_2)),
             (((3,), X), ((1, 2), E_2)),
             (((1, 2), E_2), ((3,), X)),
             (((1, 3), E_2), ((2,), X)),
             (((2, 3), E_2), ((1,), X))]
        """
        return ProductSpeciesElement(self, other)

    def structures(self, labels):
        r"""
        Iterate over the structures on the given set of labels.

        Generically, this yields a list of relabelled representatives
        of the cosets of corresponding groups.

        The relabelling is such that the first few labels correspond
        to the first factor in the atomic decomposition, etc.

        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: L = LazySpecies(ZZ, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: list(E.structures([1,2,3]))
            [((1, 2, 3), E_3)]

            sage: P = L(lambda n: CyclicPermutationGroup(n))
            sage: list(P.structures([1,2,3]))
            [((1, 2, 3), C_3), ((2, 1, 3), C_3)]

            sage: F = 1/(2-E)
            sage: list(F.structures([1,2,3]))
            [((1, 2, 3), E_3),
             ((1, 2, 3), X*E_2, 0),
             ((3, 2, 1), X*E_2, 0),
             ((2, 3, 1), X*E_2, 0),
             ((1, 2, 3), X*E_2, 1),
             ((3, 2, 1), X*E_2, 1),
             ((2, 3, 1), X*E_2, 1),
             ((1, 2, 3), X^3),
             ((3, 2, 1), X^3),
             ((3, 1, 2), X^3),
             ((2, 1, 3), X^3),
             ((2, 3, 1), X^3),
             ((1, 3, 2), X^3)]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        n = len(labels)
        l = set(labels)
        assert len(l) == n, f"The argument labels must be a set, but {labels} has duplicates"
        S = SymmetricGroup(n)
        F = self[n]
        l = list(l)[::-1]
        for M, c in F.monomial_coefficients().items():
            if c < 0:
                raise NotImplementedError("only implemented for proper non-virtual species")
            types = [tuple(S(rep)._act_on_list_on_position(l))[::-1]
                     for rep in libgap.RightTransversal(S, M._group)]
            if c == 1:
                for s in types:
                    yield s, M
            else:
                for e, s in cartesian_product([range(c), types]):
                    yield s, M, e

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

            sage: P = PolynomialSpecies(QQ, "X")
            sage: Gc = L(lambda n: sum(P(G.automorphism_group()) for G in graphs(n) if G.is_connected()) if n else 0)
            sage: E = L(lambda n: SymmetricGroup(n))

            sage: G = L(lambda n: sum(P(G.automorphism_group()) for G in graphs(n)))
            sage: E(Gc) - G
            O^7

            sage: L = LazySpecies(QQ, "X")
            sage: X = L(SymmetricGroup(1))
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: A = L.undefined(1)
            sage: A.define(X*E(A))
            sage: A
            X + X^2 + (X^3+X*E_2) + (X^2*E_2+2*X^4+X*E_3)
             + (X^2*E_3+3*X^5+3*X^3*E_2+X*{((1,2)(3,4),)}+X*E_4)
             + (X^2*E_4+2*X^2*{((1,2)(3,4),))}+6*X^4*E_2+6*X^6
                +3*X^3*E_3+X^2*E_2^2+X*E_5)
             + (X^2*E_5+2*X^3*E_2^2+6*X^4*E_3+12*X^7+14*X^5*E_2
                +3*X^3*{((1,2)(3,4),))}+3*X^3*E_4
                +X*{((3,4),(1,2),(1,3)(2,4)(5,6)))}
                +X*{((1,2)(3,4)(5,6),)}
                +X*{((2,3)(4,5),(1,3)(5,6))}+2*X^2*E_2*E_3
                +E_2*{((1,2)(3,4),)}*X+X*E_6)
             + O^8

            sage: C = L(lambda n: CyclicPermutationGroup(n) if n else 0)
            sage: F = E(C(A))
            sage: [sum(F[n].monomial_coefficients().values()) for n in range(1, 7)]
            [1, 3, 7, 19, 47, 130]
            sage: oeis(_)
            0: A001372: Number of unlabeled mappings (or mapping patterns) from n points to themselves; number of unlabeled endofunctions.


            sage: R.<q> = QQ[]
            sage: L = LazySpecies(R, "X")
            sage: E = L(lambda n: SymmetricGroup(n))
            sage: E1 = L(lambda n: SymmetricGroup(n) if n else 0)
            sage: E(q*E1)
            1 + q*X + ((q^2+q)*E_2) + ((q^3+q)*E_3+q^2*X*E_2)
             + ((q^4+q)*E_4+q^2*P_4+q^2*X*E_3+q^3*E_2^2)
             + ((q^5+q)*E_5+(q^4+q^3+q^2)*E_3*E_2+q^2*X*E_4+q^3*P_4*X)
             + ((q^6+q)*E_6+q^2*{((1,2,3)(4,6),(1,4)(2,5)(3,6))}
                +(q^5+q^3+q^2)*E_4*E_2+q^2*X*E_5
                +q^3*{((1,4)(2,3),(1,6,3,2,5,4))}
                +q^3*E_3*E_2*X+q^4*E_2*P_4+q^4*E_3^2)
             + O^7

        TESTS::

            sage: L = LazySpecies(QQ, "X")
            sage: X = L(SymmetricGroup(1))
            sage: E2 = L(SymmetricGroup(2))
            sage: X(X + E2)
            X + E_2 + O^8
            sage: E2(X + E2)
            E_2 + X*E_2 + P_4 + O^9

            sage: (1+E2)(X)
            1 + E_2 + O^7

            sage: L = LazySpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: X = L(P(SymmetricGroup(1), {1: [1]}))
            sage: Y = L(P(SymmetricGroup(1), {2: [1]}))
            sage: X(Y, 0)
            Y + O^8

            sage: L1 = LazySpecies(QQ, "X")
            sage: L = LazySpecies(QQ, "X, Y")
            sage: P = PolynomialSpecies(QQ, "X, Y")
            sage: X = L(P(SymmetricGroup(1), {1:[1]}))
            sage: E = L1(lambda n: SymmetricGroup(n))
            sage: E(X)
            1 + X + E_2(X) + E_3(X) + E_4(X) + E_5(X) + E_6(X) + O^7
        """
        fP = parent(self)
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

        args = [P(g) for g in args]

        for g in args:
            if g._coeff_stream._approximate_order == 0:
                if not g._coeff_stream.is_uninitialized() and g[0]:
                    raise ValueError("can only compose with a positive valuation series")
                g._coeff_stream._approximate_order = 1

        sorder = self._coeff_stream._approximate_order
        gv = min(g._coeff_stream._approximate_order for g in args)
        R = P._internal_poly_ring.base_ring()

        def coefficient(n):
            if not n:
                if self[0]:
                    return R(list(self[0])[0][1])
                return R.zero()
            args_flat = [[(M, c) for i in range(n+1) for M, c in g[i]]
                         for g in args]
            weights = [[M._tc for i in range(n+1) for M, _ in g[i]]
                       for g in args]
            result = R.zero()
            for i in range(n // gv + 1):
                # compute homogeneous components
                lF = defaultdict(R)
                for M, c in self[i]:
                    lF[M._mc] += R._from_dict({M: c})
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
        return P.element_class(P, coeff_stream)

class SumSpeciesElement(LazySpeciesElement):
    def __init__(self, left, right):
        self._left = left
        self._right = right
        add = super(LazySpeciesElement, type(left))._add_(left, right)
        super().__init__(add.parent(), add._coeff_stream)

    def structures(self, labels):
        yield from self._left.structures(labels)
        yield from self._right.structures(labels)

class ProductSpeciesElement(LazySpeciesElement):
    def __init__(self, left, right):
        self._left = left
        self._right = right
        add = super(LazySpeciesElement, type(left))._mul_(left, right)
        super().__init__(add.parent(), add._coeff_stream)

    def structures(self, labels):
        n = len(labels)
        l = set(labels)
        assert len(l) == n, f"The argument labels must be a set, but {labels} has duplicates"
        for k in range(n):
            for U in itertools.combinations(l, k):
                V = l.difference(U)
                yield from itertools.product(self._left.structures(U),
                                             self._right.structures(V))

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

    def __init__(self, base_ring, names, sparse):
        """
        EXAMPLES::

            sage: from sage.rings.lazy_species import LazySpecies
            sage: LazySpecies(QQ, "X, Y")
            Lazy completion of Polynomial species in X, Y over Rational Field

        TESTS::

            sage: LazySpecies(QQ, "X, Y, Z")._arity
            3
        """
        super().__init__(PolynomialSpecies(base_ring, names))
        self._arity = len(names)
