"""
Gauss valuations on polynomial rings

This file implements Gauss valuations for polynomial rings, i.e. discrete
valuations which assign to a polynomial the minimal valuation of its
coefficients.

AUTHORS:

- Julian Rüth (2013-04-15): initial version

EXAMPLES:

A Gauss valuation maps a polynomial to the minimal valuation of any of its
coefficients::

    sage: R.<x> = QQ[]
    sage: v0 = QQ.valuation(2)
    sage: v = GaussValuation(R, v0); v
    Gauss valuation induced by 2-adic valuation
    sage: v(2*x + 2)
    1

Gauss valuations can also be defined iteratively based on valuations over
polynomial rings::

    sage: v = v.augmentation(x, 1/4); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4 ]
    sage: v = v.augmentation(x^4+2*x^3+2*x^2+2*x+2, 4/3); v
    [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: S.<T> = R[]
    sage: w = GaussValuation(S, v); w
    Gauss valuation induced by [ Gauss valuation induced by 2-adic valuation, v(x) = 1/4, v(x^4 + 2*x^3 + 2*x^2 + 2*x + 2) = 4/3 ]
    sage: w(2*T + 1)
    0
"""
# ****************************************************************************
#       Copyright (C) 2013-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .inductive_valuation import NonFinalInductiveValuation

from sage.misc.cachefunc import cached_method
from sage.structure.factory import UniqueFactory


class GaussValuationFactory(UniqueFactory):
    r"""
    Create a Gauss valuation on ``domain``.

    INPUT:

    - ``domain`` -- a univariate polynomial ring

    - ``v`` -- a valuation on the base ring of ``domain``, the underlying
      valuation on the constants of the polynomial ring (if unspecified take
      the natural valuation on the valued ring ``domain``.)

    EXAMPLES:

    The Gauss valuation is the minimum of the valuation of the coefficients::

        sage: v = QQ.valuation(2)
        sage: R.<x> = QQ[]
        sage: w = GaussValuation(R, v)
        sage: w(2)
        1
        sage: w(x)
        0
        sage: w(x + 2)
        0
    """
    def create_key(self, domain, v=None):
        r"""
        Normalize and check the parameters to create a Gauss valuation.

        TESTS::

            sage: v = QQ.valuation(2)
            sage: R.<x> = ZZ[]
            sage: GaussValuation.create_key(R, v)
            Traceback (most recent call last):
            ...
            ValueError: the domain of v must be the base ring of domain but 2-adic valuation is not defined over Integer Ring but over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        if not isinstance(domain, PolynomialRing_general):
            raise TypeError("GaussValuations can only be created over polynomial rings but %r is not a polynomial ring" % (domain,))
        if not domain.ngens() == 1:
            raise NotImplementedError("domain must be univariate but %r is not univariate" % (domain,))

        if v is None:
            v = domain.base_ring().valuation()

        if v.domain() is not domain.base_ring():
            raise ValueError("the domain of v must be the base ring of domain but %r is not defined over %r but over %r" % (v, domain.base_ring(), v.domain()))
        if not v.is_discrete_valuation():
            raise ValueError("v must be a discrete valuation but %r is not" % (v,))

        return (domain, v)

    def create_object(self, version, key, **extra_args):
        r"""
        Create a Gauss valuation from normalized parameters.

        TESTS::

            sage: v = QQ.valuation(2)
            sage: R.<x> = QQ[]
            sage: GaussValuation.create_object(0, (R, v))
            Gauss valuation induced by 2-adic valuation
        """
        domain, v = key
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(domain)
        return parent.__make_element_class__(GaussValuation_generic)(parent, v)


GaussValuation = GaussValuationFactory("sage.rings.valuation.gauss_valuation.GaussValuation")


class GaussValuation_generic(NonFinalInductiveValuation):
    """
    A Gauss valuation on a polynomial ring ``domain``.

    INPUT:

    - ``domain`` -- a univariate polynomial ring over a valued ring `R`

    - ``v`` -- a discrete valuation on `R`

    EXAMPLES::

        sage: R = Zp(3,5)
        sage: S.<x> = R[]                                                               # needs sage.libs.ntl
        sage: v0 = R.valuation()
        sage: v = GaussValuation(S, v0); v                                              # needs sage.libs.ntl
        Gauss valuation induced by 3-adic valuation

        sage: S.<x> = QQ[]
        sage: v = GaussValuation(S, QQ.valuation(5)); v
        Gauss valuation induced by 5-adic valuation

    TESTS::

        sage: TestSuite(v).run()                # long time                             # needs sage.geometry.polyhedron
    """
    def __init__(self, parent, v):
        """
        TESTS::

            sage: from sage.rings.valuation.gauss_valuation import GaussValuation_generic
            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(5))
            sage: isinstance(v, GaussValuation_generic)
            True
        """
        NonFinalInductiveValuation.__init__(self, parent, parent.domain().gen())

        self._base_valuation = v

    def value_group(self):
        """
        Return the value group of this valuation.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(5))
            sage: v.value_group()
            Additive Abelian Group generated by 1
        """
        return self._base_valuation.value_group()

    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(5))
            sage: v.value_semigroup()
            Additive Abelian Semigroup generated by -1, 1
        """
        return self._base_valuation.value_semigroup()

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(5))
            sage: v # indirect doctest
            Gauss valuation induced by 5-adic valuation
        """
        return "Gauss valuation induced by %r" % self._base_valuation

    @cached_method
    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., a uniformizer of the
        valuation of the base ring.

        EXAMPLES::

            sage: S.<x> = QQ[]
            sage: v = GaussValuation(S, QQ.valuation(5))
            sage: v.uniformizer()
            5
            sage: v.uniformizer().parent() is S
            True
        """
        return self.domain()(self._base_valuation.uniformizer())

    def valuations(self, f, coefficients=None, call_error=False):
        r"""
        Return the valuations of the `f_i\phi^i` in the expansion `f=\sum f_i\phi^i`.

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        - ``coefficients`` -- the coefficients of ``f`` as produced by
          :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.coefficients` or ``None`` (default: ``None``); this can be
          used to speed up the computation when the expansion of ``f`` is
          already known from a previous computation.

        - ``call_error`` -- whether or not to speed up the computation by
          assuming that the result is only used to compute the valuation of
          ``f`` (default: ``False``)

        OUTPUT: list, each entry a rational numbers or infinity, the valuations
        of `f_0, f_1\phi, \dots`

        EXAMPLES::

            sage: R = ZZ
            sage: S.<x> = R[]
            sage: v = GaussValuation(S, R.valuation(2))
            sage: f = x^2 + 2*x + 16
            sage: list(v.valuations(f))
            [4, 1, 0]
        """
        f = self.domain().coerce(f)

        if f.is_constant():
            yield self._base_valuation(f[0])
            return

        from sage.rings.infinity import infinity
        from sage.rings.rational_field import QQ
        if f == self.domain().gen():
            yield infinity
            yield QQ(0)
            return

        if call_error:
            lowest_valuation = infinity
        for c in coefficients or f.coefficients(sparse=False):
            if call_error:
                if lowest_valuation is not infinity:
                    v = self._base_valuation.lower_bound(c)
                    if v is infinity or v >= lowest_valuation:
                        yield infinity
                        continue
            ret = self._base_valuation(c)
            if call_error:
                if ret is not infinity and (lowest_valuation is infinity or ret < lowest_valuation):
                    lowest_valuation = ret
            yield ret

    @cached_method
    def residue_ring(self):
        """
        Return the residue ring of this valuation, i.e., the elements of
        valuation zero module the elements of positive valuation.

        EXAMPLES::

            sage: S.<x> = Qp(2,5)[]                                                     # needs sage.libs.ntl
            sage: v = GaussValuation(S)                                                 # needs sage.libs.ntl
            sage: v.residue_ring()                                                      # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using ...)
        """
        return self.domain().change_ring(self._base_valuation.residue_ring())

    def reduce(self, f, check=True, degree_bound=None, coefficients=None, valuations=None):
        """
        Return the reduction of ``f`` modulo this valuation.

        INPUT:

        - ``f`` -- an integral element of the domain of this valuation

        - ``check`` -- whether or not to check whether ``f`` has nonnegative
          valuation (default: ``True``)

        - ``degree_bound`` -- an a-priori known bound on the degree of the
          result which can speed up the computation (default: not set)

        - ``coefficients`` -- the coefficients of ``f`` as produced by
          :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.coefficients` or ``None`` (default: ``None``); ignored

        - ``valuations`` -- the valuations of ``coefficients`` or ``None``
          (default: ``None``); ignored

        OUTPUT: a polynomial in the :meth:`residue_ring` of this valuation

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: S.<x> = Qp(2,5)[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 16
            sage: v.reduce(f)
            x^2
            sage: v.reduce(f).parent() is v.residue_ring()
            True

        The reduction is only defined for integral elements::

            sage: f = x^2/2                                                             # needs sage.libs.ntl
            sage: v.reduce(f)                                                           # needs sage.libs.ntl
            Traceback (most recent call last):
            ...
            ValueError: reduction not defined for non-integral elements and (2^-1 + O(2^4))*x^2 is not integral over Gauss valuation induced by 2-adic valuation

        .. SEEALSO::

            :meth:`lift`
        """
        f = self.domain().coerce(f)

        if degree_bound is not None:
            f = f.truncate(degree_bound + 1)

        try:
            return f.map_coefficients(self._base_valuation.reduce, self._base_valuation.residue_field())
        except Exception:
            if check and not all(v >= 0 for v in self.valuations(f)):
                raise ValueError("reduction not defined for non-integral elements and %r is not integral over %r" % (f, self))
            raise

    def lift(self, F):
        """
        Return a lift of ``F``.

        INPUT:

        - ``F`` -- a polynomial over the :meth:`residue_ring` of this valuation

        OUTPUT:

        a (possibly non-monic) polynomial in the domain of this valuation which
        reduces to ``F``

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: f = x^2 + 2*x + 16
            sage: F = v.reduce(f); F
            x^2 + 2*x + 1
            sage: g = v.lift(F); g
            (1 + O(3^5))*x^2 + (2 + O(3^5))*x + 1 + O(3^5)
            sage: v.is_equivalent(f,g)
            True
            sage: g.parent() is v.domain()
            True

        .. SEEALSO::

            :meth:`reduce`
        """
        F = self.residue_ring().coerce(F)
        return F.map_coefficients(self._base_valuation.lift,
                                  self._base_valuation.domain())

    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` from the :meth:`residue_ring` to
        a key polynomial over this valuation.

        INPUT:

        - ``F`` -- an irreducible non-constant monic polynomial in
          :meth:`residue_ring` of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which, for a suitable equivalence
        unit `R`, satisfies that the reduction of `Rf` is ``F``

        EXAMPLES::

            sage: R.<u> = QQ
            sage: S.<x> = R[]
            sage: v = GaussValuation(S, QQ.valuation(2))
            sage: y = v.residue_ring().gen()
            sage: f = v.lift_to_key(y^2 + y + 1); f
            x^2 + x + 1
        """
        F = self.residue_ring().coerce(F)

        if F.is_constant():
            raise ValueError("F must not be constant but %r is constant" % (F,))
        if not F.is_monic():
            raise ValueError("F must be monic but %r is not monic" % (F,))
        if not F.is_irreducible():
            raise ValueError("F must be irreducible but %r factors" % (F,))

        return self.lift(F)

    @cached_method
    def equivalence_unit(self, s, reciprocal=False):
        """
        Return an equivalence unit of valuation ``s``.

        INPUT:

        - ``s`` -- an element of the :meth:`value_group`

        - ``reciprocal`` -- boolean (default: ``False``); whether or not to
          return the equivalence unit as the :meth:`~sage.rings.valuation.inductive_valuation.InductiveValuation.equivalence_reciprocal` of
          the equivalence unit of valuation ``-s``

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_unit(2)
            3^2 + O(3^7)
            sage: v.equivalence_unit(-2)
            3^-2 + O(3^3)
        """
        if reciprocal:
            return self.equivalence_reciprocal(self.equivalence_unit(-s))

        ret = self._base_valuation.element_with_valuation(s)
        return self.domain()(ret)

    def element_with_valuation(self, s):
        r"""
        Return a polynomial of minimal degree with valuation ``s``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v.element_with_valuation(-2)
            1/4
        """
        return self.equivalence_unit(s)

    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1
        """
        from sage.rings.integer_ring import ZZ
        return ZZ.one()

    def F(self):
        """
        Return the degree of the residue field extension of this valuation
        over the Gauss valuation, i.e., 1.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1
        """
        from sage.rings.integer_ring import ZZ
        return ZZ.one()

    def change_domain(self, ring):
        r"""
        Return this valuation as a valuation over ``ring``.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.change_domain(QQ['x'])
            Gauss valuation induced by 2-adic valuation
        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        if isinstance(ring, PolynomialRing_general) and ring.ngens() == 1:
            base_valuation = self._base_valuation.change_domain(ring.base_ring())
            return GaussValuation(self.domain().change_ring(ring.base_ring()), base_valuation)
        return super().change_domain(ring)

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.extensions(GaussianIntegers()['x'])                                 # needs sage.rings.number_field
            [Gauss valuation induced by 2-adic valuation]
        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        if isinstance(ring, PolynomialRing_general) and ring.ngens() == 1:
            if self.domain().is_subring(ring):
                return [GaussValuation(ring, w) for w in self._base_valuation.extensions(ring.base_ring())]
        return super().extensions(ring)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.restriction(ZZ)
            2-adic valuation
        """
        if ring.is_subring(self.domain().base_ring()):
            return self._base_valuation.restriction(ring)
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        if isinstance(ring, PolynomialRing_general) and ring.ngens() == 1:
            if ring.base().is_subring(self.domain().base()):
                return GaussValuation(ring, self._base_valuation.restriction(ring.base()))
        return super().restriction(ring)

    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True
        """
        return True

    def augmentation_chain(self):
        r"""
        Return a list with the chain of augmentations down to the underlying
        Gauss valuation.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.augmentation_chain()
            [Gauss valuation induced by 2-adic valuation]
        """
        return [self]

    def is_trivial(self):
        r"""
        Return whether this is a trivial valuation (sending everything but zero
        to zero.)

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v.is_trivial()
            True
        """
        return self._base_valuation.is_trivial()

    def monic_integral_model(self, G):
        r"""
        Return a monic integral irreducible polynomial which defines the same
        extension of the base ring of the domain as the irreducible polynomial
        ``G`` together with maps between the old and the new polynomial.

        EXAMPLES::

            sage: R.<x> = Qp(2, 5)[]                                                    # needs sage.libs.ntl
            sage: v = GaussValuation(R)                                                 # needs sage.libs.ntl
            sage: v.monic_integral_model(5*x^2 + 1/2*x + 1/4)                           # needs sage.libs.ntl
            (Ring endomorphism of Univariate Polynomial Ring in x over 2-adic Field with capped relative precision 5
               Defn: (1 + O(2^5))*x |--> (2^-1 + O(2^4))*x,
             Ring endomorphism of Univariate Polynomial Ring in x over 2-adic Field with capped relative precision 5
               Defn: (1 + O(2^5))*x |--> (2 + O(2^6))*x,
             (1 + O(2^5))*x^2 + (1 + 2^2 + 2^3 + O(2^5))*x + 1 + 2^2 + 2^3 + O(2^5))
        """
        if not G.is_monic():
            # this might fail if the base ring is not a field
            G = G / G.leading_coefficient()

        x = G.parent().gen()
        u = self._base_valuation.uniformizer()

        factor = 1
        substitution = x
        H = G
        while self(H) < 0:
            # this might fail if the base ring is not a field
            factor *= u
            substitution = x / factor
            H = G(substitution) * (factor ** G.degree())

        assert H.is_monic()
        return H.parent().hom(substitution, G.parent()), G.parent().hom(x / substitution[1], H.parent()), H

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal to ``other``
        everywhere.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: w = GaussValuation(R, QQ.valuation(3))
            sage: v >= w
            False
            sage: w >= v
            False
        """
        if isinstance(other, GaussValuation_generic):
            return self._base_valuation >= other._base_valuation
        from .augmented_valuation import AugmentedValuation_base
        if isinstance(other, AugmentedValuation_base):
            return False
        if other.is_trivial():
            return other.is_discrete_valuation()
        return super()._ge_(other)

    def scale(self, scalar):
        r"""
        Return this valuation scaled by ``scalar``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: 3*v # indirect doctest
            Gauss valuation induced by 3 * 2-adic valuation
        """
        from sage.rings.rational_field import QQ
        if scalar in QQ and scalar > 0 and scalar != 1:
            return GaussValuation(self.domain(), self._base_valuation.scale(scalar))
        return super().scale(scalar)

    def _relative_size(self, f):
        r"""
        Return an estimate on the coefficient size of ``f``.

        The number returned is an estimate on the factor between the number of
        Bits used by ``f`` and the minimal number of bits used by an element
        Congruent to ``f``.

        This is used by :meth:`simplify` to decide whether simplification of
        Coefficients is going to lead to a significant shrinking of the
        Coefficients of ``f``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v._relative_size(x + 1024)
            6

        For performance reasons, only the constant coefficient is considered.
        (In common applications, the constant coefficient shows the most
        critical coefficient growth)::

            sage: v._relative_size(1024*x + 1)
            1
        """
        return self._base_valuation._relative_size(f[0])

    def simplify(self, f, error=None, force=False, size_heuristic_bound=32, effective_degree=None, phiadic=True):
        r"""
        Return a simplified version of ``f``.

        Produce an element which differs from ``f`` by an element of valuation
        strictly greater than the valuation of ``f`` (or strictly greater than
        ``error`` if set.)

        INPUT:

        - ``f`` -- an element in the domain of this valuation

        - ``error`` -- a rational, infinity, or ``None`` (default: ``None``),
          the error allowed to introduce through the simplification

        - ``force`` -- whether or not to simplify ``f`` even if there is
          heuristically no change in the coefficient size of ``f`` expected
          (default: ``False``)

        - ``effective_degree`` -- when set, assume that coefficients beyond
          ``effective_degree`` can be safely dropped (default: ``None``)

        - ``size_heuristic_bound`` -- when ``force`` is not set, the expected
          factor by which the coefficients need to shrink to perform an actual
          simplification (default: 32)

        - ``phiadic`` -- whether to simplify in the `x`-adic expansion; the
          parameter is ignored as no other simplification is implemented

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = x^10/2 + 1
            sage: v.simplify(f)
            (2^-1 + O(2^4))*x^10 + 1 + O(2^5)
        """
        f = self.domain().coerce(f)

        if effective_degree is not None:
            if effective_degree < f.degree():
                f = f.truncate(effective_degree + 1)

        if not force and self._relative_size(f) < size_heuristic_bound:
            return f

        if error is None:
            # if the caller was sure that we should simplify, then we should try to do the best simplification possible
            error = self(f) if force else self.uppper_bound(f)

        return f.map_coefficients(lambda c: self._base_valuation.simplify(c, error=error, force=force))

    def lower_bound(self, f):
        r"""
        Return a lower bound of this valuation at ``f``.

        Use this method to get an approximation of the valuation of ``f``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.lower_bound(1024*x + 2)
            1
            sage: v(1024*x + 2)
            1
        """
        from sage.rings.infinity import infinity
        coefficients = f.coefficients(sparse=True)
        coefficients.reverse()
        ret = infinity
        for c in coefficients:
            v = self._base_valuation.lower_bound(c)
            if c is not infinity and (ret is infinity or v < ret):
                ret = v
        return ret

    def upper_bound(self, f):
        r"""
        Return an upper bound of this valuation at ``f``.

        Use this method to get an approximation of the valuation of ``f``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.upper_bound(1024*x + 1)
            10
            sage: v(1024*x + 1)
            0
        """
        f = self.domain().coerce(f)
        coefficients = f.coefficients(sparse=True)
        if not coefficients:
            from sage.rings.infinity import infinity
            return infinity
        else:
            return self._base_valuation.upper_bound(coefficients[-1])
