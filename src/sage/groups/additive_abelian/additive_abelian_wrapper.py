r"""
Wrapper class for abelian groups

This class is intended as a template for anything in Sage that needs the
functionality of abelian groups. One can create an AdditiveAbelianGroupWrapper
object from any given set of elements in some given parent, as long as an
``_add_`` method has been defined.


EXAMPLES:

We create a toy example based on the Mordell-Weil group of an elliptic curve over `\QQ`::

    sage: E = EllipticCurve('30a2')
    sage: pts = [E(4,-7,1), E(7/4, -11/8, 1), E(3, -2, 1)]
    sage: M = AdditiveAbelianGroupWrapper(pts[0].parent(), pts, [3, 2, 2])
    sage: M
    Additive abelian group isomorphic to Z/3 + Z/2 + Z/2 embedded in Abelian
    group of points on Elliptic Curve defined by y^2 + x*y + y = x^3 - 19*x + 26
    over Rational Field
    sage: M.gens()
    ((4 : -7 : 1), (7/4 : -11/8 : 1), (3 : -2 : 1))
    sage: 3*M.0
    (0 : 1 : 0)
    sage: 3000000000000001 * M.0
    (4 : -7 : 1)
    sage: M == loads(dumps(M))  # known bug, see https://trac.sagemath.org/sage_trac/ticket/11599#comment:7
    True

We check that ridiculous operations are being avoided::

    sage: from sage.misc.verbose import set_verbose
    sage: set_verbose(2, 'additive_abelian_wrapper.py')
    sage: 300001 * M.0
    verbose 1 (...: additive_abelian_wrapper.py, _discrete_exp) Calling discrete exp on (1, 0, 0)
    (4 : -7 : 1)
    sage: set_verbose(0, 'additive_abelian_wrapper.py')


.. TODO::

    - Implement proper black-box discrete logarithm (using baby-step
      giant-step).  The discrete_exp function can also potentially be
      speeded up substantially via caching.

    - Think about subgroups and quotients, which probably won't work
      in the current implementation -- some fiddly adjustments will be
      needed in order to be able to pass extra arguments to the
      subquotient's init method.
"""

from . import additive_abelian_group as addgp
from sage.rings.integer_ring import ZZ
from sage.categories.morphism import Morphism
from sage.structure.element import parent
from sage.modules.free_module_element import vector


class UnwrappingMorphism(Morphism):
    r"""
    The embedding into the ambient group. Used by the coercion framework.
    """
    def __init__(self, domain):
        r"""
        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: F = QQbar.coerce_map_from(G); F
            Generic morphism:
              From: Additive abelian group isomorphic to Z + Z embedded in Algebraic Field
              To:   Algebraic Field
            sage: type(F)
            <class 'sage.groups.additive_abelian.additive_abelian_wrapper.UnwrappingMorphism'>
        """
        Morphism.__init__(self, domain.Hom(domain.universe()))

    def _call_(self, x):
        r"""
        TESTS::

            sage: E = EllipticCurve("65a1")
            sage: G = E.torsion_subgroup()
            sage: isinstance(G, sage.groups.additive_abelian.additive_abelian_wrapper.AdditiveAbelianGroupWrapper)
            True
            sage: P1 = E([1,-1,1])
            sage: P2 = E([0,1,0])
            sage: P1 in G # indirect doctest
            False
            sage: P2 in G
            True
            sage: (G(P2) + P1) in G
            False
            sage: (G(P2) + P1).parent()
            Abelian group of points on Elliptic Curve defined by y^2 + x*y = x^3 - x over Rational Field
        """
        return self.codomain()(x.element())


class AdditiveAbelianGroupWrapperElement(addgp.AdditiveAbelianGroupElement):
    """
    An element of an :class:`AdditiveAbelianGroupWrapper`.
    """

    def __init__(self, parent, vector, element=None, check=False):
        r"""
        EXAMPLES::

            sage: from sage.groups.additive_abelian.additive_abelian_wrapper import AdditiveAbelianGroupWrapper
            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: G.0 # indirect doctest
            1.414213562373095?
        """
        addgp.AdditiveAbelianGroupElement.__init__(self, parent, vector, check)
        if element is not None:
            element = self.parent().universe()(element)
        self._element = element

    def element(self):
        r"""
        Return the underlying object that this element wraps.

        EXAMPLES::

            sage: T = EllipticCurve('65a').torsion_subgroup().gen(0)
            sage: T; type(T)
            (0 : 0 : 1)
            <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup_with_category.element_class'>
            sage: T.element(); type(T.element())
            (0 : 0 : 1)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_number_field'>
        """
        if self._element is None:
            self._element = self.parent()._discrete_exp(self._hermite_lift())
        return self._element

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: T = EllipticCurve('65a').torsion_subgroup().gen(0)
            sage: repr(T) # indirect doctest
            '(0 : 0 : 1)'
        """
        return repr(self.element())


class AdditiveAbelianGroupWrapper(addgp.AdditiveAbelianGroup_fixed_gens):
    """
    The parent of :class:`AdditiveAbelianGroupWrapperElement`
    """

    Element = AdditiveAbelianGroupWrapperElement

    def __init__(self, universe, gens, invariants):
        r"""
        EXAMPLES::

            sage: AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0]) # indirect doctest
            Additive abelian group isomorphic to Z + Z embedded in Algebraic Field
        """
        self._universe = universe
        self._gen_elements = tuple(universe(x) for x in gens)
        self._gen_orders = invariants
        cover,rels = addgp.cover_and_relations_from_invariants(invariants)
        addgp.AdditiveAbelianGroup_fixed_gens.__init__(self, cover, rels, cover.gens())
        self._unset_coercions_used()
        self.register_embedding(UnwrappingMorphism(self))

    def universe(self):
        r"""
        The ambient group in which this abelian group lives.

        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: G.universe()
            Algebraic Field
        """
        return self._universe

    def generator_orders(self):
        r"""
        The orders of the generators with which this group was initialised.
        (Note that these are not necessarily a minimal set of generators.)
        Generators of infinite order are returned as 0. Compare
        ``self.invariants()``, which returns the orders of a minimal set of
        generators.

        EXAMPLES::

            sage: V = Zmod(6)**2
            sage: G = AdditiveAbelianGroupWrapper(V, [2*V.0, 3*V.1], [3, 2])
            sage: G.generator_orders()
            (3, 2)
            sage: G.invariants()
            (6,)
        """
        return tuple(self._gen_orders)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: repr(G) # indirect doctest
            'Additive abelian group isomorphic to Z + Z embedded in Algebraic Field'
        """
        return addgp.AdditiveAbelianGroup_fixed_gens._repr_(self) + " embedded in " + self.universe()._repr_()

    def _discrete_exp(self, v):
        r"""
        Given a list (or other iterable) of length equal to the number of
        generators of this group, compute the element of the ambient group with
        those exponents in terms of the generators of self.

        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(QQbar(2)), -1], [0, 0])
            sage: v = G._discrete_exp([3, 5]); v
            -0.7573593128807148?
            sage: v.parent() is QQbar
            True
        """
        from sage.misc.verbose import verbose
        v = self.V()(v)
        verbose("Calling discrete exp on %s" % v)
        # DUMB IMPLEMENTATION!
        return sum([self._gen_elements[i] * ZZ(v[i]) for i in range(len(v))], self.universe()(0))

    def _discrete_log_pgroup(self, p, aa, b):
        r"""
        Attempt to express an element of p-power order in terms of
        generators of a p-subgroup of this group.

        Used as a subroutine in the _discrete_log() method.

        ALGORITHM:

        This implements a basic version of the recursive algorithm
        from [Suth2008]_.
        The base cases are handled using a variant of Shanks'
        baby-step giant-step algorithm for products of cyclic groups.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([5, 5**2, 5**4, 5**4])
            sage: (a, b, c, d) = gs = G.gens()
            sage: A = AdditiveAbelianGroupWrapper(a.parent(), gs, [g.order() for g in gs])
            sage: A._discrete_log_pgroup(5, gs, a + 17 * b + 123 * c + 456 * d)
            (1, 17, 123, 456)
        """
        from sage.arith.misc import valuation
        from sage.functions.other import ceil
        from sage.misc.functional import sqrt
        from itertools import product as iproduct

        vals = [valuation(a.order(), p) for a in aa]
        qq = lambda j, k: vector(p ** (j + max(0, v - k)) for a, v in zip(aa, vals))
        subbasis = lambda j, k: [q * a for q, a in zip(qq(j, k), aa)]
        dotprod = lambda xs, ys: sum(x * y for x, y in zip(xs, ys))

        def _base(j, k, c):

            assert k - j == 1
            aajk = subbasis(j, k)
            assert all(a.order() in (1, p) for a in aajk)
            idxs = [i for i, a in enumerate(aajk) if a.order() == p]

            rs = [([0], [0]) for i in range(len(aajk))]
            for i in range(len(idxs)):
                rs[idxs[i]] = (range(p), [0]) if i % 2 else ([0], range(p))
            if len(idxs) % 2:
                m = ceil(sqrt(p))
                rs[idxs[-1]] = range(0, p, m), range(m)

            tab = {}
            for x in iproduct(*(r for r, _ in rs)):
                key = dotprod(x, aajk)
                if hasattr(key, 'set_immutable'):
                    key.set_immutable()
                tab[key] = vector(x)
            for y in iproduct(*(r for _, r in rs)):
                key = c - dotprod(y, aajk)
                if hasattr(key, 'set_immutable'):
                    key.set_immutable()
                if key in tab:
                    return tab[key] + vector(y)

            raise TypeError('Not in group')

        def _rec(j, k, c):

            assert 0 <= j < k

            if k - j <= 1: # base case
                return _base(j, k, c)

            w = 2
            js = list(range(j, k, (k-j+w-1) // w)) + [k]
            assert len(js) == w + 1

            x = vector([0] * len(aa))
            for i in reversed(range(w)):

                gamma = p ** (js[i] - j) * c - dotprod(x, subbasis(js[i], k))

                v = _rec(js[i], js[i+1], gamma)

                assert not any(q1 % q2 for q1, q2 in zip(qq(js[i], js[i+1]), qq(js[i], k)))
                x += vector(q1 // q2 * r for q1, q2, r in zip(qq(js[i], js[i+1]), qq(js[i], k), v))

            return x

        return _rec(0, max(vals), b)

    def _discrete_log(self, x, gens=None):
        r"""
        Given an element of the ambient group, attempt to express it in terms
        of the generators of this group or the given generators of a subgroup.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2, 2*3, 2*3*5, 2*3*5*7, 2*3*5*7*11])
            sage: A = AdditiveAbelianGroupWrapper(G.0.parent(), G.gens(), [g.order() for g in G.gens()])
            sage: A._discrete_log(G.0 + 5 * G.1 + 23 * G.2 + 127 * G.3 + 539 * G.4)
            (1, 5, 23, 127, 539)
            sage: V = Zmod(8)**2; G = AdditiveAbelianGroupWrapper(V, [[2,2],[4,0]], [4, 2])
            sage: G._discrete_log(V([6, 2]))
            (1, 1)
            sage: G._discrete_log(V([6, 4]))
            Traceback (most recent call last):
            ...
            TypeError: Not in group
            sage: F.<t> = GF(1009**2, modulus=x**2+11); E = EllipticCurve(j=F(940))
            sage: P, Q = E(900*t + 228, 974*t + 185), E(1007*t + 214, 865*t + 802)
            sage: E.abelian_group()._discrete_log(123 * P + 777 * Q, [P, Q])
            (123, 777)
            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(2)], [0])
            sage: G._discrete_log(QQbar(2*sqrt(2)))
            Traceback (most recent call last):
            ...
            NotImplementedError: No black-box discrete log for infinite abelian groups
        """
        from sage.arith.misc import CRT_list
        from sage.rings.infinity import Infinity

        if self.order() == Infinity:
            raise NotImplementedError("No black-box discrete log for infinite abelian groups")

        if gens is None:
            gens = self.gens()

        gens = [g if parent(g) is self.universe() else g.element() for g in gens]
        x = x if parent(x) is self.universe() else x.element()

        crt_data = [[] for _ in gens]
        for p, e in self.order().factor():
            cofactor = self.order() // p ** e
            pgens = [cofactor * g for g in gens]
            y = cofactor * x

            plog = self._discrete_log_pgroup(p, pgens, y)

            for i, (r, g) in enumerate(zip(plog, pgens)):
                crt_data[i].append((r, ZZ(g.order())))

        res = vector(CRT_list(*map(list, zip(*l))) for l in crt_data)
        assert x == sum(r * g for r, g in zip(res, gens))
        return res

    def _element_constructor_(self, x, check=False):
        r"""
        Create an element from x. This may be either an element of self, an element of the
        ambient group, or an iterable (in which case the result is the corresponding
        product of the generators of self).

        EXAMPLES::

            sage: V = Zmod(8)**2; G = AdditiveAbelianGroupWrapper(V, [[2,2],[4,0]], [4, 2])
            sage: G(V([6,2]))
            (6, 2)
            sage: G([1,1])
            (6, 2)
            sage: G(G([1,1]))
            (6, 2)
        """
        if parent(x) is self.universe():
            return self.element_class(self, self._discrete_log(x), element = x)
        return addgp.AdditiveAbelianGroup_fixed_gens._element_constructor_(self, x, check)

