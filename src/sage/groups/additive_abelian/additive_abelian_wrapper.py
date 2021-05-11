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
from sage.rings.all import ZZ
from sage.categories.morphism import Morphism
from sage.structure.element import parent


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

    def _discrete_log(self,x):
        r"""
        Given an element of the ambient group, attempt to express it in terms of the
        generators of self.

        EXAMPLES::

            sage: V = Zmod(8)**2; G = AdditiveAbelianGroupWrapper(V, [[2,2],[4,0]], [4, 2])
            sage: G._discrete_log(V([6, 2]))
            (1, 1)
            sage: G._discrete_log(V([6, 4]))
            Traceback (most recent call last):
            ...
            TypeError: Not in group
            sage: G = AdditiveAbelianGroupWrapper(QQbar, [sqrt(2)], [0])
            sage: G._discrete_log(QQbar(2*sqrt(2)))
            Traceback (most recent call last):
            ...
            NotImplementedError: No black-box discrete log for infinite abelian groups
        """
        # EVEN DUMBER IMPLEMENTATION!
        from sage.rings.infinity import Infinity
        if self.order() == Infinity:
            raise NotImplementedError("No black-box discrete log for infinite abelian groups")
        u = [y for y in self.list() if y.element() == x]
        if len(u) == 0: raise TypeError("Not in group")
        if len(u) > 1: raise NotImplementedError
        return u[0].vector()

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

