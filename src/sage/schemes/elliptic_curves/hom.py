"""
Elliptic-curve morphisms

This class serves as a common parent for various specializations of
morphisms between elliptic curves, with the aim of providing a common
interface regardless of implementation details.

Current implementations of elliptic-curve morphisms (child classes):

- :class:`~sage.schemes.elliptic_curves.ell_curve_isogeny.EllipticCurveIsogeny`
- :class:`~sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism`
- :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`
- :class:`~sage.schemes.elliptic_curves.hom_velusqrt.EllipticCurveHom_velusqrt`

AUTHORS:

- See authors of :class:`EllipticCurveIsogeny`. Some of the code
  in this class was lifted from there.

- Lorenz Panny (2021): Refactor isogenies and isomorphisms into
  the common :class:`EllipticCurveHom` interface.
"""

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp_not_equal, richcmp, op_EQ, op_NE

from sage.categories.morphism import Morphism

from sage.arith.misc import integer_floor

from sage.rings.finite_rings import finite_field_base
from sage.rings.number_field import number_field_base


class EllipticCurveHom(Morphism):
    """
    Base class for elliptic-curve morphisms.
    """
    def _repr_type(self):
        """
        Return a textual representation of what kind of morphism
        this is. Used by :meth:`Morphism._repr_`.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom._repr_type(None)
            'Elliptic-curve'
        """
        return 'Elliptic-curve'

    @staticmethod
    def _composition_impl(left, right):
        """
        Called by :meth:`_composition_`.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom._composition_impl(None, None)
            NotImplemented
        """
        return NotImplemented

    def _composition_(self, other, homset):
        """
        Return the composition of this elliptic-curve morphism
        with another elliptic-curve morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(19), [1,0])
            sage: phi = E.isogeny(E(0,0))
            sage: iso = E.change_weierstrass_model(5,0,0,0).isomorphism_to(E)
            sage: phi * iso
            Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 9*x over Finite Field of size 19 to Elliptic Curve defined by y^2 = x^3 + 15*x over Finite Field of size 19
            sage: phi.dual() * phi
            Composite map:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
              Defn:   Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19 to Elliptic Curve defined by y^2 = x^3 + 15*x over Finite Field of size 19
                    then
                      Isogeny of degree 2 from Elliptic Curve defined by y^2 = x^3 + 15*x over Finite Field of size 19 to Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 19
        """
        if not isinstance(self, EllipticCurveHom) or not isinstance(other, EllipticCurveHom):
            raise TypeError(f'cannot compose {type(self)} with {type(other)}')

        ret = self._composition_impl(self, other)
        if ret is not NotImplemented:
            return ret

        ret = other._composition_impl(self, other)
        if ret is not NotImplemented:
            return ret

        # fall back to generic formal composite map
        return Morphism._composition_(self, other, homset)


    @staticmethod
    def _comparison_impl(left, right, op):
        """
        Called by :meth:`_richcmp_`.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom._comparison_impl(None, None, None)
            NotImplemented
        """
        return NotImplemented

    def _richcmp_(self, other, op):
        r"""
        Compare :class:`EllipticCurveHom` objects.

        ALGORITHM:

        The method first makes sure that domain, codomain and degree match.
        Then, it determines if there is a specialized comparison method by
        trying :meth:`_comparison_impl` on either input. If not, it falls
        back to comparing :meth:`rational_maps`.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi_v = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi_k = EllipticCurveIsogeny(E, [0,1])
            sage: phi_k == phi_v
            True
            sage: E_F17 = EllipticCurve(GF(17), [0,0,0,1,0])
            sage: phi_p = EllipticCurveIsogeny(E_F17, [0,1])
            sage: phi_p == phi_v
            False
            sage: E = EllipticCurve('11a1')
            sage: phi = E.isogeny(E(5,5))
            sage: phi == phi
            True
            sage: phi == -phi
            False
            sage: psi = E.isogeny(phi.kernel_polynomial())
            sage: phi == psi
            True
            sage: phi.dual() == psi.dual()
            True

        ::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve([9,9])
            sage: F = E.change_ring(GF(71))
            sage: wE = WeierstrassIsomorphism(E, (1,0,0,0))
            sage: wF = WeierstrassIsomorphism(F, (1,0,0,0))
            sage: mE = E.multiplication_by_m_isogeny(1)
            sage: mF = F.multiplication_by_m_isogeny(1)
            sage: [mE == wE, mF == wF]
            [True, True]
            sage: [a == b for a in (wE,mE) for b in (wF,mF)]
            [False, False, False, False]

        .. SEEALSO::

            - :meth:`_comparison_impl`
            - :func:`compare_via_evaluation`
        """
        if not isinstance(self, EllipticCurveHom) or not isinstance(other, EllipticCurveHom):
            raise TypeError(f'cannot compare {type(self)} to {type(other)}')

        if op == op_NE:
            return not self._richcmp_(other, op_EQ)

        # We first compare domain, codomain, and degree; cf. Trac #11327

        lx, rx = self.domain(), other.domain()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx, rx = self.codomain(), other.codomain()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        lx, rx = self.degree(), other.degree()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        # Do self or other have specialized comparison methods?

        ret = self._comparison_impl(self, other, op)
        if ret is not NotImplemented:
            return ret

        ret = other._comparison_impl(self, other, op)
        if ret is not NotImplemented:
            return ret

        # If not, fall back to comparing rational maps; cf. Trac #11327

        return richcmp(self.rational_maps(), other.rational_maps(), op)


    def degree(self):
        r"""
        Return the degree of this elliptic-curve morphism.

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi.degree()
            2
            sage: phi = EllipticCurveIsogeny(E, [0,1,0,1])
            sage: phi.degree()
            4

            sage: E = EllipticCurve(GF(31), [1,0,0,1,2])
            sage: phi = EllipticCurveIsogeny(E, [17, 1])
            sage: phi.degree()
            3

        Degrees are multiplicative, so the degree of a composite isogeny
        is the product of the degrees of the individual factors::

            sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
            doctest:warning ...
            sage: E = EllipticCurve(GF(419), [1,0])
            sage: P, = E.gens()
            sage: phi = EllipticCurveHom_composite(E, P+P)
            sage: phi.degree()
            210
            sage: phi.degree() == prod(f.degree() for f in phi.factors())
            True

        Isomorphisms always have degree `1` by definition::

            sage: E1 = EllipticCurve([1,2,3,4,5])
            sage: E2 = EllipticCurve_from_j(E1.j_invariant())
            sage: E1.isomorphism_to(E2).degree()
            1

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.degree(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        try:
            return self._degree
        except AttributeError:
            raise NotImplementedError('children must implement')

    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this elliptic-curve morphism.

        Implemented by child classes. For examples, see:

        - :meth:`EllipticCurveIsogeny.kernel_polynomial`
        - :meth:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism.kernel_polynomial`
        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite.kernel_polynomial`

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.kernel_polynomial(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        raise NotImplementedError('children must implement')

    def dual(self):
        r"""
        Return the dual of this elliptic-curve morphism.

        Implemented by child classes. For examples, see:

        - :meth:`EllipticCurveIsogeny.dual`
        - :meth:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism.dual`
        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite.dual`

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.dual(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        raise NotImplementedError('children must implement')

    def rational_maps(self):
        r"""
        Return the pair of explicit rational maps defining this
        elliptic-curve morphism as fractions of bivariate
        polynomials in `x` and `y`.

        Implemented by child classes. For examples, see:

        - :meth:`EllipticCurveIsogeny.rational_maps`
        - :meth:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism.rational_maps`
        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite.rational_maps`

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.rational_maps(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        raise NotImplementedError('children must implement')

    def x_rational_map(self):
        r"""
        Return the `x`-coordinate rational map of this elliptic-curve
        morphism as a univariate rational expression in `x`.

        Implemented by child classes. For examples, see:

        - :meth:`EllipticCurveIsogeny.x_rational_map`
        - :meth:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism.x_rational_map`
        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite.x_rational_map`

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.x_rational_map(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        #TODO: could have a default implementation that simply
        #      returns the first component of rational_maps()
        raise NotImplementedError('children must implement')


    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        elliptic-curve morphism.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        Implemented by child classes. For examples, see:

        - :meth:`EllipticCurveIsogeny.scaling_factor`
        - :meth:`sage.schemes.elliptic_curves.weierstrass_morphism.WeierstrassIsomorphism.scaling_factor`
        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite.scaling_factor`

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: EllipticCurveHom.scaling_factor(None)
            Traceback (most recent call last):
            ...
            NotImplementedError: ...
        """
        #TODO: could have a default implementation that simply
        #      returns .formal()[1], but it seems safer to fail
        #      visibly to make sure we would notice regressions
        raise NotImplementedError('children must implement')

    def formal(self, prec=20):
        r"""
        Return the formal isogeny associated to this elliptic-curve
        morphism as a power series in the variable `t=-x/y` on the
        domain curve.

        INPUT:

        - ``prec`` -- (default: 20), the precision with which the
          computations in the formal group are carried out.

        EXAMPLES::

            sage: E = EllipticCurve(GF(13),[1,7])
            sage: phi = E.isogeny(E(10,4))
            sage: phi.formal()
            t + 12*t^13 + 2*t^17 + 8*t^19 + 2*t^21 + O(t^23)

        ::

            sage: E = EllipticCurve([0,1])
            sage: phi = E.isogeny(E(2,3))
            sage: phi.formal(prec=10)
            t + 54*t^5 + 255*t^7 + 2430*t^9 + 19278*t^11 + O(t^13)

        ::

            sage: E = EllipticCurve('11a2')
            sage: R.<x> = QQ[]
            sage: phi = E.isogeny(x^2 + 101*x + 12751/5)
            sage: phi.formal(prec=7)
            t - 2724/5*t^5 + 209046/5*t^7 - 4767/5*t^8 + 29200946/5*t^9 + O(t^10)
        """
        Eh = self._domain.formal()
        f, g = self.rational_maps()
        xh = Eh.x(prec=prec)
        assert xh.valuation() == -2, f"xh has valuation {xh.valuation()} (should be -2)"
        yh = Eh.y(prec=prec)
        assert yh.valuation() == -3, f"yh has valuation {yh.valuation()} (should be -3)"
        fh = f(xh,yh)
        assert fh.valuation() == -2, f"fh has valuation {fh.valuation()} (should be -2)"
        gh = g(xh,yh)
        assert gh.valuation() == -3, f"gh has valuation {gh.valuation()} (should be -3)"
        th = -fh/gh
        assert th.valuation() == +1, f"th has valuation {th.valuation()} (should be +1)"
        return th


    def is_normalized(self):
        r"""
        Determine whether this morphism is a normalized isogeny.

        .. NOTE::

            An isogeny `\varphi\colon E_1\to E_2` between two given
            Weierstrass equations is said to be *normalized* if the
            `\varphi^*(\omega_2) = \omega_1`, where `\omega_1` and
            `\omega_2` are the invariant differentials on `E_1` and
            `E_2` corresponding to the given equation.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
            sage: R.<x> = GF(7)[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (3, 0, 0, 0))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (5, 0, 0, 0))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True

        ::

            sage: F = GF(2^5, 'alpha'); alpha = F.gen()
            sage: E = EllipticCurve(F, [1,0,1,1,1])
            sage: R.<x> = F[]
            sage: phi = EllipticCurveIsogeny(E, x+1)
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (alpha, 0, 0, 0))
            sage: phi.is_normalized()
            True
            sage: phi = isom * phi
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1/alpha, 0, 0, 0))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True

        ::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^3 - x^2 - 10*x - 79/4
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (2, 0, 0, 0))
            sage: phi.is_normalized()
            True
            sage: phi = isom * phi
            sage: phi.is_normalized()
            False
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1/2, 0, 0, 0))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True
            sage: isom = WeierstrassIsomorphism(phi.codomain(), (1, 1, 1, 1))
            sage: phi = isom * phi
            sage: phi.is_normalized()
            True

        ALGORITHM: We check if :meth:`scaling_factor` returns `1`.
        """
        return self.scaling_factor() == 1


    def is_separable(self):
        r"""
        Determine whether or not this morphism is separable.

        .. NOTE::

            This method currently always returns ``True`` as Sage does
            not yet implement inseparable isogenies. This will probably
            change in the future.

        EXAMPLES::

            sage: E = EllipticCurve(GF(17), [0,0,0,3,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.is_separable()
            True

        ::

            sage: E = EllipticCurve('11a1')
            sage: phi = EllipticCurveIsogeny(E, E.torsion_points())
            sage: phi.is_separable()
            True
        """
        return True

    def is_surjective(self):
        r"""
        Determine whether or not this morphism is surjective.

        .. NOTE::

            This method currently always returns ``True``, since a
            non-constant map of algebraic curves must be surjective,
            and Sage does not yet implement the constant zero map.
            This will probably change in the future.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 + x - 29/5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.is_surjective()
            True

        ::

            sage: E = EllipticCurve(GF(7), [0,0,0,1,0])
            sage: phi = EllipticCurveIsogeny(E,  E((0,0)))
            sage: phi.is_surjective()
            True

        ::

            sage: F = GF(2^5, 'omega')
            sage: E = EllipticCurve(j=F(0))
            sage: R.<x> = F[]
            sage: phi = EllipticCurveIsogeny(E, x)
            sage: phi.is_surjective()
            True
        """
        return bool(self.degree())

    def is_injective(self):
        r"""
        Determine whether or not this morphism has trivial kernel.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: R.<x> = QQ[]
            sage: f = x^2 + x - 29/5
            sage: phi = EllipticCurveIsogeny(E, f)
            sage: phi.is_injective()
            False
            sage: phi = EllipticCurveIsogeny(E, R(1))
            sage: phi.is_injective()
            True

        ::

            sage: F = GF(7)
            sage: E = EllipticCurve(j=F(0))
            sage: phi = EllipticCurveIsogeny(E, [ E((0,-1)), E((0,1))])
            sage: phi.is_injective()
            False
            sage: phi = EllipticCurveIsogeny(E, E(0))
            sage: phi.is_injective()
            True
        """
        if not self.is_separable():
            #TODO: should implement .separable_degree() or similar
            raise NotImplementedError
        return self.degree() == 1

    def is_zero(self):
        """
        Check whether this elliptic-curve morphism is the zero map.

        .. NOTE::

            This function currently always returns ``True`` as Sage
            does not yet implement the constant zero morphism. This
            will probably change in the future.

        EXAMPLES::

            sage: E = EllipticCurve(j=GF(7)(0))
            sage: phi = EllipticCurveIsogeny(E, [E(0,1), E(0,-1)])
            sage: phi.is_zero()
            False
        """
        return not self.degree()

    def __neg__(self):
        r"""
        Return the negative of this elliptic-curve morphism. In other
        words, return `[-1]\circ\varphi` where `\varphi` is ``self``
        and `[-1]` is the negation automorphism on the codomain curve.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom import EllipticCurveHom
            sage: E = EllipticCurve(GF(1019), [5,5])
            sage: phi = E.isogeny(E.lift_x(73))
            sage: f,g = phi.rational_maps()
            sage: psi = EllipticCurveHom.__neg__(phi)
            sage: psi.rational_maps() == (f, -g)
            True
        """
        import sage.schemes.elliptic_curves.weierstrass_morphism as wm
        a1,_,a3,_,_ = self.codomain().a_invariants()
        return wm.WeierstrassIsomorphism(self.codomain(), (-1,0,-a1,-a3)) * self

    @cached_method
    def __hash__(self):
        """
        Return a hash value for this elliptic-curve morphism.

        ALGORITHM:

        Hash a tuple containing the domain, codomain, and kernel
        polynomial of this morphism. (The base field is factored
        into the computation as part of the (co)domain hashes.)

        EXAMPLES::

            sage: E = EllipticCurve(QQ, [0,0,0,1,0])
            sage: phi_v = EllipticCurveIsogeny(E, E((0,0)))
            sage: phi_k = EllipticCurveIsogeny(E, [0,1])
            sage: phi_k.__hash__() == phi_v.__hash__()
            True
            sage: E_F17 = EllipticCurve(GF(17), [0,0,0,1,1])
            sage: phi_p = EllipticCurveIsogeny(E_F17, E_F17([0,1]))
            sage: phi_p.__hash__() == phi_v.__hash__()
            False

        ::

            sage: E = EllipticCurve('49a3')
            sage: R.<X> = QQ[]
            sage: EllipticCurveIsogeny(E,X^3-13*X^2-58*X+503,check=False)
            Isogeny of degree 7 from Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 107*x + 552 over Rational Field to Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 5252*x - 178837 over Rational Field
        """
        return hash((self.domain(), self.codomain(), self.kernel_polynomial()))


def compare_via_evaluation(left, right):
    r"""
    Test if two elliptic-curve morphisms are equal by evaluating
    them at enough points.

    INPUT:

    - ``left``, ``right`` -- :class:`EllipticCurveHom` objects

    ALGORITHM:

    We use the fact that two isogenies of equal degree `d` must be
    the same if and only if they behave identically on more than
    `4d` points. (It suffices to check this on a few points that
    generate a large enough subgroup.)

    If the domain curve does not have sufficiently many rational
    points, the base field is extended first: Taking an extension
    of degree `O(\log(d))` suffices.

    EXAMPLES::

        sage: E = EllipticCurve(GF(83), [1,0])
        sage: phi = E.isogeny(12*E.0, model='montgomery'); phi
        Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 83 to Elliptic Curve defined by y^2 = x^3 + 70*x^2 + x over Finite Field of size 83
        sage: psi = phi.dual(); psi
        Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 70*x^2 + x over Finite Field of size 83 to Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 83
        sage: from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
        sage: mu = EllipticCurveHom_composite.from_factors([phi, psi])
        sage: from sage.schemes.elliptic_curves.hom import compare_via_evaluation
        sage: compare_via_evaluation(mu, E.multiplication_by_m_isogeny(7))
        True

    .. SEEALSO::

        - :meth:`sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite._richcmp_`
    """
    if left.domain() != right.domain():
        return False
    if left.codomain() != right.codomain():
        return False
    if left.degree() != right.degree():
        return False

    E = left.domain()
    F = E.base_ring()

    if isinstance(F, finite_field_base.FiniteField):
        q = F.cardinality()
        d = left.degree()
        e = integer_floor(1 + 2 * (2*d.sqrt() + 1).log(q))  # from Hasse bound
        e = next(i for i,n in enumerate(E.count_points(e+1), 1) if n > 4*d)
        EE = E.base_extend(F.extension(e))
        Ps = EE.gens()
        return all(left._eval(P) == right._eval(P) for P in Ps)

    elif isinstance(F, number_field_base.NumberField):
        for _ in range(100):
            P = E.lift_x(F.random_element(), extend=True)
            if not P.has_finite_order():
                return left._eval(P) == right._eval(P)
        else:
            assert False, "couldn't find a point of infinite order"

    else:
        raise NotImplementedError('not implemented for this base field')

