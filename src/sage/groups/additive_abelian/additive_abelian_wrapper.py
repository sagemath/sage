# sage_setup: distribution = sagemath-modules
r"""
Wrapper class for abelian groups

This class is intended as a template for anything in Sage that needs the
functionality of abelian groups. One can create an ``AdditiveAbelianGroupWrapper``
object from any given set of elements in some given parent, as long as an
``_add_`` method has been defined.

EXAMPLES:

We create a toy example based on the Mordell-Weil group of an elliptic curve over `\QQ`::

    sage: # needs sage.schemes
    sage: E = EllipticCurve('30a2')
    sage: pts = [E(4,-7,1), E(7/4, -11/8, 1), E(3, -2, 1)]
    sage: M = AdditiveAbelianGroupWrapper(pts[0].parent(), pts, [3, 2, 2]); M
    Additive abelian group isomorphic to Z/3 + Z/2 + Z/2 embedded in Abelian
    group of points on Elliptic Curve defined by y^2 + x*y + y = x^3 - 19*x + 26
    over Rational Field
    sage: M.gens()
    ((4 : -7 : 1), (7/4 : -11/8 : 1), (3 : -2 : 1))
    sage: 3*M.0
    (0 : 1 : 0)
    sage: 3000000000000001 * M.0
    (4 : -7 : 1)
    sage: M == loads(dumps(M))  # known bug (https://github.com/sagemath/sage/issues/11599#comment:7)
    True

TESTS:

We check that ridiculous operations are being avoided::

    sage: from sage.misc.verbose import set_verbose
    sage: set_verbose(2, 'additive_abelian_wrapper.py')
    sage: 300001 * M.0                                                                  # needs sage.schemes
    verbose 1 (...: additive_abelian_wrapper.py, discrete_exp) Calling discrete exp on (1, 0, 0)
    (4 : -7 : 1)
    sage: set_verbose(0, 'additive_abelian_wrapper.py')


.. TODO::

    Think about subgroups and quotients, which probably won't work
    in the current implementation -- some fiddly adjustments will be
    needed in order to be able to pass extra arguments to the
    subquotient's init method.

AUTHORS:

- David Loeffler (2010)
- Lorenz Panny (2017): :meth:`AdditiveAbelianGroupWrapper.discrete_log`
- Lorenz Panny (2023): :meth:`AdditiveAbelianGroupWrapper.from_generators`
"""

# ****************************************************************************
#       Copyright (C) 2010 David Loeffler
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import additive_abelian_group as addgp
from sage.rings.integer_ring import ZZ
from sage.categories.morphism import Morphism
from sage.structure.element import parent
from sage.structure.sequence import Sequence
from sage.modules.free_module_element import vector

from sage.misc.superseded import deprecated_function_alias


class UnwrappingMorphism(Morphism):
    r"""
    The embedding into the ambient group. Used by the coercion framework.
    """
    def __init__(self, domain):
        r"""
        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar,                                # needs sage.rings.number_field
            ....:                                 [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: F = QQbar.coerce_map_from(G); F                                       # needs sage.rings.number_field
            Generic morphism:
              From: Additive abelian group isomorphic to Z + Z embedded in Algebraic Field
              To:   Algebraic Field
            sage: type(F)                                                               # needs sage.rings.number_field
            <class 'sage.groups.additive_abelian.additive_abelian_wrapper.UnwrappingMorphism'>
        """
        Morphism.__init__(self, domain.Hom(domain.universe()))

    def _call_(self, x):
        r"""
        TESTS::

            sage: # needs sage.schemes
            sage: E = EllipticCurve("65a1")
            sage: G = E.torsion_subgroup()
            sage: isinstance(G, sage.groups.additive_abelian.additive_abelian_wrapper.AdditiveAbelianGroupWrapper)
            True
            sage: P1 = E([1,-1,1])
            sage: P2 = E([0,1,0])
            sage: P1 in G  # indirect doctest
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
            sage: G = AdditiveAbelianGroupWrapper(QQbar,                                # needs sage.rings.number_field
            ....:                                 [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: G.0  # indirect doctest                                               # needs sage.rings.number_field
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

            sage: T = EllipticCurve('65a').torsion_subgroup().gen(0)                    # needs sage.schemes
            sage: T; type(T)                                                            # needs sage.schemes
            (0 : 0 : 1)
            <class 'sage.schemes.elliptic_curves.ell_torsion.EllipticCurveTorsionSubgroup_with_category.element_class'>
            sage: T.element(); type(T.element())                                        # needs sage.schemes
            (0 : 0 : 1)
            <class 'sage.schemes.elliptic_curves.ell_point.EllipticCurvePoint_number_field'>
        """
        if self._element is None:
            self._element = self.parent().discrete_exp(self._hermite_lift())
        return self._element

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: T = EllipticCurve('65a').torsion_subgroup().gen(0)                    # needs sage.schemes
            sage: repr(T)  # indirect doctest                                           # needs sage.schemes
            '(0 : 0 : 1)'
        """
        return repr(self.element())


class AdditiveAbelianGroupWrapper(addgp.AdditiveAbelianGroup_fixed_gens):
    """
    This class is used to wrap a subgroup of an existing
    additive abelian group as a new additive abelian group.

    EXAMPLES::

        sage: G2 = AdditiveAbelianGroupWrapper(Zmod(42), [2], [21]); G2
        Additive abelian group isomorphic to Z/21 embedded in Ring of integers modulo 42
        sage: G6 = AdditiveAbelianGroupWrapper(Zmod(42), [6], [7]); G6
        Additive abelian group isomorphic to Z/7 embedded in Ring of integers modulo 42
        sage: G = AdditiveAbelianGroupWrapper(Zmod(42), [21,14,6], [2,3,7]); G
        Additive abelian group isomorphic to Z/2 + Z/3 + Z/7 embedded in
         Ring of integers modulo 42
        sage: G.invariants()
        (42,)

    ::

        sage: AdditiveAbelianGroupWrapper(QQbar, [sqrt(2), sqrt(3)], [0, 0])            # needs sage.rings.number_field sage.symbolic
        Additive abelian group isomorphic to Z + Z embedded in Algebraic Field

    ::

        sage: EllipticCurve(GF(419**2), [1,0]).abelian_group()  # indirect doctest      # needs sage.rings.finite_rings sage.schemes
        Additive abelian group isomorphic to Z/420 + Z/420 embedded in
         Abelian group of points on Elliptic Curve
          defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2
    """

    Element = AdditiveAbelianGroupWrapperElement

    def __init__(self, universe, gens, invariants):
        r"""
        EXAMPLES::

            sage: AdditiveAbelianGroupWrapper(QQbar,  # indirect doctest                # needs sage.rings.number_field
            ....:                             [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            Additive abelian group isomorphic to Z + Z embedded in Algebraic Field
        """
        self._universe = universe
        self._gen_elements = tuple(universe(x) for x in gens)
        self._gen_orders = invariants
        cover, rels = addgp.cover_and_relations_from_invariants(invariants)
        addgp.AdditiveAbelianGroup_fixed_gens.__init__(self, cover, rels, cover.gens())
        self._unset_coercions_used()
        self.register_embedding(UnwrappingMorphism(self))

    def universe(self):
        r"""
        The ambient group in which this abelian group lives.

        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar,                                # needs sage.rings.number_field
            ....:                                 [sqrt(QQbar(2)), sqrt(QQbar(3))],
            ....:                                 [0, 0])
            sage: G.universe()                                                          # needs sage.rings.number_field
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

            sage: G = AdditiveAbelianGroupWrapper(QQbar,                                # needs sage.rings.number_field
            ....:                                 [sqrt(QQbar(2)), sqrt(QQbar(3))], [0, 0])
            sage: repr(G)  # indirect doctest                                           # needs sage.rings.number_field
            'Additive abelian group isomorphic to Z + Z embedded in Algebraic Field'
        """
        return addgp.AdditiveAbelianGroup_fixed_gens._repr_(self) + " embedded in " + self.universe()._repr_()

    def _element_constructor_(self, x, check=False):
        r"""
        Create an element from ``x``.

        This may be either an element of self, an element of the ambient
        group, or an iterable (in which case the result is the corresponding
        product of the generators of self).

        EXAMPLES::

            sage: V = Zmod(8)**2
            sage: G = AdditiveAbelianGroupWrapper(V, [[2,2],[4,0]], [4, 2])
            sage: G(V([6,2]))
            (6, 2)
            sage: G([1,1])
            (6, 2)
            sage: G(G([1,1]))
            (6, 2)
        """
        if parent(x) is self.universe():
            return self.element_class(self, self.discrete_log(x), element=x)
        return addgp.AdditiveAbelianGroup_fixed_gens._element_constructor_(self, x, check)

    def discrete_exp(self, v):
        r"""
        Given a list (or other iterable) of length equal to the number of
        generators of this group, compute the element of the ambient group
        with those exponents in terms of the generators of ``self``.

        EXAMPLES::

            sage: G = AdditiveAbelianGroupWrapper(QQbar,                                # needs sage.rings.number_field
            ....:                                 [sqrt(QQbar(2)), -1], [0, 0])
            sage: v = G.discrete_exp([3, 5]); v                                         # needs sage.rings.number_field
            -0.7573593128807148?
            sage: v.parent() is QQbar                                                   # needs sage.rings.number_field
            True

        This method is an inverse of :meth:`discrete_log`::

            sage: orders = [2, 2*3, 2*3*5, 2*3*5*7, 2*3*5*7*11]
            sage: G = AdditiveAbelianGroup(orders)
            sage: A = AdditiveAbelianGroupWrapper(G.0.parent(), G.gens(), orders)
            sage: el = A.random_element()
            sage: A.discrete_exp(A.discrete_log(el)) == el
            True

        TESTS:

        Check that :meth:`_discrete_exp` still works (for now)::

            sage: A._discrete_exp(list(range(1,6)))
            doctest:warning ...
            DeprecationWarning: _discrete_exp is deprecated. ...
            (1, 2, 3, 4, 5)
        """
        from sage.misc.verbose import verbose
        v = self.V()(v)
        verbose("Calling discrete exp on %s" % v)
        # DUMB IMPLEMENTATION!
        return sum([self._gen_elements[i] * ZZ(v[i]) for i in range(len(v))], self.universe()(0))

    _discrete_exp = deprecated_function_alias(32384, discrete_exp)

    def discrete_log(self, x, gens=None):
        r"""
        Given an element of the ambient group, attempt to express it in terms
        of the generators of this group or the given generators of a subgroup.

        ALGORITHM:

        This reduces to p-groups, then calls :func:`_discrete_log_pgroup` which
        implements a basic version of the recursive algorithm from [Suth2008]_.

        AUTHORS:

        - Lorenz Panny (2017)

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([2, 2*3, 2*3*5, 2*3*5*7, 2*3*5*7*11])
            sage: A = AdditiveAbelianGroupWrapper(G.0.parent(), G.gens(),
            ....:                                 [g.order() for g in G.gens()])
            sage: A.discrete_log(A.discrete_exp([1,5,23,127,539]))
            (1, 5, 23, 127, 539)

        ::

            sage: x = polygen(ZZ, 'x')
            sage: F.<t> = GF(1009**2, modulus=x**2+11); E = EllipticCurve(j=F(940))     # needs sage.rings.finite_rings sage.schemes
            sage: P, Q = E(900*t + 228, 974*t + 185), E(1007*t + 214, 865*t + 802)      # needs sage.rings.finite_rings sage.schemes
            sage: E.abelian_group().discrete_log(123 * P + 777 * Q, [P, Q])             # needs sage.rings.finite_rings sage.schemes
            (123, 777)

        ::

            sage: V = Zmod(8)**2
            sage: G = AdditiveAbelianGroupWrapper(V, [[2,2],[4,0]], [4, 2])
            sage: G.discrete_log(V([6, 2]))
            (1, 1)
            sage: G.discrete_log(V([6, 4]))
            Traceback (most recent call last):
            ...
            ValueError: not in group

        ::

            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(2)], [0])                # needs sage.rings.number_field sage.symbolic
            sage: G.discrete_log(QQbar(2*sqrt(2)))                                      # needs sage.rings.number_field sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: No black-box discrete log for infinite abelian groups

        TESTS:

        Check that :meth:`_discrete_log` still works (for now)::

            sage: orders = [2, 2*3, 2*3*5, 2*3*5*7, 2*3*5*7*11]
            sage: G = AdditiveAbelianGroup(orders)
            sage: A = AdditiveAbelianGroupWrapper(G.0.parent(), G.gens(), orders)
            sage: A._discrete_log(sum(i*g for i,g in enumerate(G.gens(),1)))
            doctest:warning ...
            DeprecationWarning: _discrete_log is deprecated. ...
            (1, 2, 3, 4, 5)
        """
        from sage.arith.misc import CRT_list
        from sage.rings.infinity import Infinity

        if self.order() == Infinity:
            raise NotImplementedError("No black-box discrete log for infinite abelian groups")

        if gens is None:
            gens = self.gens()
            ords = self.generator_orders()
        else:
            ords = [g.order() for g in gens]

        gens = [self._universe(g.element() if parent(g) is self else g) for g in gens]
        x = self._universe(x.element() if parent(x) is self else x)

        crt_data = [[] for _ in gens]
        for p in self.exponent().prime_factors():
            cofactor = self.exponent().prime_to_m_part(p)
            pgens = [cofactor * g for g in gens]
            y = cofactor * x

            pvals = [o.valuation(p) for o in ords]
            if not any(pvals):
                continue

            plog = _discrete_log_pgroup(p, pvals, pgens, y)

            for i, (r, v) in enumerate(zip(plog, pvals)):
                crt_data[i].append((r, p**v))

        res = vector(CRT_list(*map(list, zip(*l))) for l in crt_data)
        assert x == sum(r * g for r, g in zip(res, gens))
        return res

    _discrete_log = deprecated_function_alias(32384, discrete_log)

    def torsion_subgroup(self, n=None):
        r"""
        Return the `n`-torsion subgroup of this additive abelian group
        when `n` is given, and the torsion subgroup otherwise.

        The [`n`-]torsion subgroup consists of all elements whose order
        is finite [and divides `n`].

        EXAMPLES::

            sage: ords = [2, 2*3, 2*3*5, 0, 2*3*5*7, 2*3*5*7*11]
            sage: G = AdditiveAbelianGroup(ords)
            sage: A = AdditiveAbelianGroupWrapper(G.0.parent(), G.gens(), ords)
            sage: T = A.torsion_subgroup(5)
            sage: T
            Additive abelian group isomorphic to Z/5 + Z/5 + Z/5 embedded in
             Additive abelian group isomorphic to Z/2 + Z/6 + Z/30 + Z + Z/210 + Z/2310
            sage: T.gens()
            ((0, 0, 6, 0, 0, 0), (0, 0, 0, 0, 42, 0), (0, 0, 0, 0, 0, 462))

        ::

            sage: # needs sage.rings.finite_rings sage.schemes
            sage: E = EllipticCurve(GF(487^2), [311,205])
            sage: T = E.abelian_group().torsion_subgroup(42); T
            Additive abelian group isomorphic to Z/42 + Z/6 embedded in
             Abelian group of points on Elliptic Curve
              defined by y^2 = x^3 + 311*x + 205 over Finite Field in z2 of size 487^2
            sage: [P.order() for P in T.gens()]
            [42, 6]

        ::

            sage: # needs sage.schemes
            sage: E = EllipticCurve('574i1')
            sage: pts = [E(103,172), E(61,18)]
            sage: assert pts[0].order() == 7 and pts[1].order() == infinity
            sage: M = AdditiveAbelianGroupWrapper(pts[0].parent(), pts, [7,0]); M
            Additive abelian group isomorphic to Z/7 + Z embedded in
             Abelian group of points on Elliptic Curve defined by
              y^2 + x*y + y = x^3 - x^2 - 19353*x + 958713 over Rational Field
            sage: M.torsion_subgroup()
            Additive abelian group isomorphic to Z/7 embedded in
             Abelian group of points on Elliptic Curve defined by
              y^2 + x*y + y = x^3 - x^2 - 19353*x + 958713 over Rational Field
            sage: M.torsion_subgroup(7)
            Additive abelian group isomorphic to Z/7 embedded in
             Abelian group of points on Elliptic Curve defined by
              y^2 + x*y + y = x^3 - x^2 - 19353*x + 958713 over Rational Field
            sage: M.torsion_subgroup(5)
            Trivial group embedded in Abelian group of points on Elliptic Curve
             defined by y^2 + x*y + y = x^3 - x^2 - 19353*x + 958713 over Rational Field

        AUTHORS:

        - Lorenz Panny (2022)
        """
        genords = zip(self._gen_elements, self._gen_orders)
        if n is None:
            gens, ords = zip(*(t for t in genords if t[1]))
        else:
            n = ZZ(n)
            if n <= 0:
                raise ValueError('n must be a positive integer')
            gens, ords = [], []
            for g, o in genords:
                if not o:
                    continue
                d = n.gcd(o)
                if d == 1:
                    continue
                gens.append(o // d * g)
                ords.append(d)
        return AdditiveAbelianGroupWrapper(self.universe(), gens, ords)

    @staticmethod
    def from_generators(gens, universe=None):
        r"""
        This method constructs the subgroup generated by a sequence
        of *finite-order* elements in an additive abelian group.

        The elements need not be independent, hence this can be used
        to perform tasks such as finding relations between some given
        elements of an abelian group, computing the structure of the
        generated subgroup, enumerating all elements of the subgroup,
        and solving discrete-logarithm problems.

        EXAMPLES::

            sage: G = AdditiveAbelianGroup([15, 30, 45])
            sage: gs = [G((1,2,3)), G((4,5,6)), G((7,7,7)), G((3,2,1))]
            sage: H = AdditiveAbelianGroupWrapper.from_generators(gs); H
            Additive abelian group isomorphic to Z/90 + Z/15 embedded in
             Additive abelian group isomorphic to Z/15 + Z/30 + Z/45
            sage: H.gens()
            ((12, 13, 14), (1, 26, 21))

        TESTS:

        Random testing::

            sage: invs = []
            sage: while not 1 < prod(invs) < 10^4:
            ....:     invs = [randrange(1,100) for _ in range(randrange(1,20))]
            sage: G = AdditiveAbelianGroup(invs)
            sage: gs = [G.random_element() for _ in range(randrange(1,10))]
            sage: H = AdditiveAbelianGroupWrapper.from_generators(gs)
            sage: os = H.generator_orders()
            sage: vecs = cartesian_product_iterator(list(map(range, os)))
            sage: els = {sum(i*g for i,g in zip(vec, H.gens())) for vec in vecs}
            sage: len(els) == prod(os)
            True
        """
        if not gens:
            if universe is None:
                raise ValueError('need universe if no generators are given')
            return AdditiveAbelianGroupWrapper(universe, [], [])

        if universe is None:
            universe = Sequence(gens).universe()

        basis, ords = basis_from_generators(gens)
        return AdditiveAbelianGroupWrapper(universe, basis, ords)


def _discrete_log_pgroup(p, vals, aa, b):
    r"""
    Attempt to express an element of p-power order in terms of
    generators of a nontrivial p-subgroup of this group.

    Used as a subroutine in :meth:`discrete_log`.

    ALGORITHM:

    This implements a basic version of the recursive algorithm
    from [Suth2008]_.
    The base cases are handled using a variant of Shanks'
    baby-step giant-step algorithm for products of cyclic groups.

    EXAMPLES::

        sage: G = AdditiveAbelianGroup([5, 5**2, 5**4, 5**4])
        sage: (a, b, c, d) = gs = G.gens()
        sage: A = AdditiveAbelianGroupWrapper(a.parent(), gs, [g.order() for g in gs])
        sage: from sage.groups.additive_abelian.additive_abelian_wrapper import _discrete_log_pgroup
        sage: _discrete_log_pgroup(5, [1,2,4,4], gs, a + 17*b + 123*c + 456*d)
        (1, 17, 123, 456)

    TESTS:

    Check for :issue:`34716`::

        sage: # needs sage.rings.finite_rings sage.schemes
        sage: E = EllipticCurve(GF(487^2), [311,205])
        sage: G = E.abelian_group().torsion_subgroup(42)
        sage: G.invariants()
        (6, 42)
        sage: P, Q = G.torsion_subgroup(6).gens()
        sage: G.discrete_log(2*P + 3*Q, [P, Q])  # indirect doctest                     # needs sage.groups
        (2, 3)
    """
    from itertools import product as iproduct

    qq = lambda j, k: vector(p ** (j + max(0, v - k)) for a, v in zip(aa, vals))
    subbasis = lambda j, k: [q * a for q, a in zip(qq(j, k), aa)]
    dotprod = lambda xs, ys: sum(x * y for x, y in zip(xs, ys))

    def _base(j, k, c):

        assert k - j == 1
        aajk = subbasis(j, k)
        assert not any(p*a for a in aajk)  # orders are in {1,p}
        idxs = [i for i, a in enumerate(aajk) if a]

        rs = [([0], [0]) for i in range(len(aajk))]
        for i in range(len(idxs)):
            rs[idxs[i]] = (range(p), [0]) if i % 2 else ([0], range(p))
        if len(idxs) % 2:
            m = p.isqrt() + 1  # hence m^2 >= p
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

        raise ValueError('not in group')

    def _rec(j, k, c):

        assert 0 <= j < k

        if k - j <= 1:  # base case
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


def _expand_basis_pgroup(p, alphas, vals, beta, h, rel):
    r"""
    Given a basis of a `p`-subgroup of a finite abelian group
    and an element lying outside the subgroup, extend the basis
    to the subgroup spanned jointly by the original subgroup and
    the new element.

    Used as a subroutine in :func:`basis_from_generators`.

    This function modifies ``alphas`` and ``vals`` in place.

    ALGORITHM: [Suth2007]_, Algorithm 9.2

    INPUT:

    - ``p`` -- prime integer `p`
    - ``alphas`` -- list; basis for a `p`-subgroup of an abelian group
    - ``vals`` -- list; valuation at `p` of the orders of the ``alphas``
    - ``beta`` -- element of the same abelian group as the ``alphas``
    - ``h`` -- integer; valuation at `p` of the order of ``beta``
    - ``rel`` -- list of integers; relation on ``alphas + [beta]``

    OUTPUT: basis of the subgroup generated by ``alphas + [beta]``

    EXAMPLES::

        sage: from sage.groups.additive_abelian.additive_abelian_wrapper import _expand_basis_pgroup
        sage: A = AdditiveAbelianGroup([9,3])
        sage: alphas = [A((5,2))]
        sage: beta = A((1,0))
        sage: vals = [2]
        sage: rel = next([ZZ(r),ZZ(s)] for s in range(9) for r in range(9) if s > 1 and not r*alphas[0] + s*beta)
        sage: _expand_basis_pgroup(3, alphas, vals, beta, 2, rel)
        sage: alphas
        [(5, 2), (6, 2)]
        sage: vals
        [2, 1]
        sage: len({i*alphas[0] + j*alphas[1] for i in range(3^2) for j in range(3^1)})
        27
    """
    # The given assertions should hold, but were commented out for speed.

    k = len(rel)
    if not (isinstance(alphas, list) and isinstance(vals, list)):
        raise TypeError('alphas and vals must be lists for mutability')
    if not len(alphas) == len(vals) == k - 1:
        raise ValueError('alphas and/or vals have incorrect length')
    #    assert not sum(r*a for r,a in zip(rel, alphas+[beta]))
    #    assert all(a.order() == p**v for a,v in zip(alphas,vals))

    if rel[-1] < 0:
        raise ValueError('rel must have nonnegative entries')

    # step 1
    min_r = rel[-1] or float('inf')
    for i in range(k-1):
        if not rel[i]:
            continue
        if rel[i] < 0:
            raise ValueError('rel must have nonnegative entries')
        q = rel[i].p_primary_part(p)
        alphas[i] *= rel[i] // q
        rel[i] = q
        min_r = min(q, min_r)
    if min_r == float('inf'):
        raise ValueError('rel must have at least one nonzero entry')
    val_rlast = rel[-1].valuation(p)
#    assert rel[-1] == p ** val_rlast
#    assert not sum(r*a for r,a in zip(rel, alphas+[beta]))

    # step 2
    if rel[-1] == min_r:
        for i in range(k-1):
            beta += alphas[i] * (rel[i]//rel[-1])
        alphas.append(beta)
        vals.append(val_rlast)
#        assert alphas[-1].order() == p**vals[-1]
        return

    # step 3
    j = next(j for j, r in enumerate(rel) if r == min_r)
    alphas[j] = sum(a * (r // rel[j]) for a, r in zip(alphas + [beta], rel))

    # step 4
    if not alphas[j]:
        del alphas[j], vals[j]
        if not alphas:
            alphas.append(beta)
            vals.append(val_rlast)
#            assert alphas[-1].order() == p**vals[-1]
            return

    # step 5
    beta_q = beta
    for v in range(1, h):
        beta_q *= p
        try:
            e = _discrete_log_pgroup(p, vals, alphas, -beta_q)
        except ValueError:
            continue
        # step 6
        _expand_basis_pgroup(p, alphas, vals, beta, h, list(e) + [p**v])
        break
    else:
        alphas.append(beta)
        vals.append(h)
    #    assert alphas[-1].order() == p**vals[-1]


def basis_from_generators(gens, ords=None):
    r"""
    Given a generating set of some finite abelian group
    (additively written), compute and return a basis of
    the group.

    .. NOTE::

        A *basis* of a finite abelian group is a generating
        set `\{g_1, \ldots, g_n\}` such that each element of the
        group can be written as a unique linear combination
        `\alpha_1 g_1 + \cdots + \alpha_n g_n` with each
        `\alpha_i \in \{0, \ldots, \mathrm{ord}(g_i)-1\}`.

    ALGORITHM: [Suth2007]_, Algorithm 9.1 & Remark 9.1

    EXAMPLES::

        sage: # needs sage.groups sage.rings.finite_rings
        sage: from sage.groups.additive_abelian.additive_abelian_wrapper import basis_from_generators
        sage: E = EllipticCurve(GF(31337^6,'a'), j=37)
        sage: E.order()
        946988065073788930380545280
        sage: (R,S), (ordR,ordS) = basis_from_generators(E.gens())
        sage: ordR, ordS
        (313157428926517503432720, 3024)
        sage: R.order() == ordR
        True
        sage: S.order() == ordS
        True
        sage: ordR * ordS == E.order()
        True
        sage: R.weil_pairing(S, ordR).multiplicative_order() == ordS
        True
        sage: E.abelian_group().invariants()
        (3024, 313157428926517503432720)
    """
    if not gens:
        return [], []
    if ords is None:
        ords = [g.order() for g in gens]

    from sage.arith.functions import lcm
    lam = lcm(ords)
    ps = sorted(lam.prime_factors(), key=lam.valuation)

    gammas = []
    ms = []
    for p in ps:
        pgens = [(o.prime_to_m_part(p) * g, o.p_primary_part(p))
                 for g, o in zip(gens, ords) if not o % p]
        assert pgens
        pgens.sort(key=lambda tup: tup[1])

        alpha, ord_alpha = pgens.pop()
        vals = [ord_alpha.valuation(p)]
        alphas = [alpha]

        while pgens:
            beta, ord_beta = pgens.pop()
            try:
                _ = _discrete_log_pgroup(p, vals, alphas, beta)
            except ValueError:
                pass
            else:
                continue

            # step 4
            val_beta = ord_beta.valuation(p)
            beta_q = beta
            for v in range(1, val_beta):
                beta_q *= p
#                assert beta_q == beta * p**v
                try:
                    e = _discrete_log_pgroup(p, vals, alphas, -beta_q)
                except ValueError:
                    continue
                _expand_basis_pgroup(p, alphas, vals, beta, val_beta, list(e) + [p**v])
#                assert all(a.order() == p**v for a,v in zip(alphas, vals))
                break
            else:
                alphas.append(beta)
                vals.append(val_beta)

        for i, (v, a) in enumerate(sorted(zip(vals, alphas), reverse=True)):
            if i < len(gammas):
                gammas[i] += a
                ms[i] *= p ** v
            else:
                gammas.append(a)
                ms.append(p ** v)

#    assert len({sum(i*g for i,g in zip(vec,gammas))
#                for vec in __import__('itertools').product(*map(range,ms))}) \
#               == __import__('sage').misc.misc_c.prod(ms)

    return gammas, ms
