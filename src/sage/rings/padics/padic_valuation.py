r"""
`p`-adic Valuations on Number Fields and Their Subrings and Completions

EXAMPLES::

    sage: ZZ.valuation(2)
    2-adic valuation
    sage: QQ.valuation(3)
    3-adic valuation
    sage: CyclotomicField(5).valuation(5)                                               # needs sage.rings.number_field
    5-adic valuation
    sage: GaussianIntegers().valuation(7)                                               # needs sage.rings.number_field
    7-adic valuation
    sage: Zp(11).valuation()
    11-adic valuation

These valuations can then, e.g., be used to compute approximate factorizations
in the completion of a ring::

    sage: v = ZZ.valuation(2)
    sage: R.<x> = ZZ[]
    sage: f = x^5 + x^4 + x^3 + x^2 + x - 1
    sage: v.montes_factorization(f, required_precision=20)                              # needs sage.geometry.polyhedron
    (x + 676027) * (x^4 + 372550*x^3 + 464863*x^2 + 385052*x + 297869)

AUTHORS:

- Julian Rüth (2013-03-16): initial version

REFERENCES:

The theory used here was originally developed in [Mac1936I]_ and [Mac1936II]_. An
overview can also be found in Chapter 4 of [Rüt2014]_.
"""
#*****************************************************************************
#       Copyright (C) 2013-2020 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.valuation.valuation import DiscreteValuation
from sage.rings.valuation.value_group import DiscreteValueSemigroup
from sage.rings.valuation.mapped_valuation import FiniteExtensionFromLimitValuation
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method

from sage.rings.infinity import infinity


class PadicValuationFactory(UniqueFactory):
    r"""
    Create a ``prime``-adic valuation on ``R``.

    INPUT:

    - ``R`` -- a subring of a number field or a subring of a local field in
      characteristic zero

    - ``prime`` -- a prime that does not split, a discrete (pseudo-)valuation,
      a fractional ideal, or ``None`` (default: ``None``)

    EXAMPLES:

    For integers and rational numbers, ``prime`` is just a prime of the
    integers::

        sage: valuations.pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: valuations.pAdicValuation(QQ, 3)
        3-adic valuation

    ``prime`` may be ``None`` for local rings::

        sage: valuations.pAdicValuation(Qp(2))
        2-adic valuation

        sage: valuations.pAdicValuation(Zp(2))
        2-adic valuation

    But it must be specified in all other cases::

        sage: valuations.pAdicValuation(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: prime must be specified for this ring

    It can sometimes be beneficial to define a number field extension as a
    quotient of a polynomial ring (since number field extensions always compute
    an absolute polynomial defining the extension which can be very costly)::

        sage: # needs sage.rings.number_field
        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^2 + 1)
        sage: R.<x> = K[]
        sage: L.<b> = R.quo(x^2 + a)
        sage: valuations.pAdicValuation(L, 2)
        2-adic valuation

    .. SEEALSO::

        :meth:`NumberField_generic.valuation() <sage.rings.number_field.number_field.NumberField_generic.valuation>`,
        :meth:`Order.valuation() <sage.rings.number_field.order.Order.valuation>`,
        :meth:`pAdicGeneric.valuation() <sage.rings.padics.padic_generic.pAdicGeneric.valuation>`,
        :meth:`RationalField.valuation() <sage.rings.rational_field.RationalField.valuation>`,
        :meth:`IntegerRing_class.valuation() <sage.rings.integer_ring.IntegerRing_class.valuation>`.
    """
    def create_key_and_extra_args(self, R, prime=None, approximants=None):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: QQ.valuation(2)  # indirect doctest
            2-adic valuation
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.number_field.number_field_base import NumberField
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic

        if R.characteristic() != 0:
            # We do not support equal characteristic yet
            raise ValueError("R must be a ring of characteristic zero.")

        if R is ZZ or R is QQ:
            return self.create_key_for_integers(R, prime), {}
        elif isinstance(R, pAdicGeneric):
            return self.create_key_for_local_ring(R, prime), {}
        elif isinstance(R.fraction_field(), NumberField) or isinstance(R, PolynomialQuotientRing_generic):
            return self.create_key_and_extra_args_for_number_field(R, prime, approximants=approximants)
        else:
            raise NotImplementedError("p-adic valuations not implemented for %r" % (R,))

    def create_key_for_integers(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: QQ.valuation(2)  # indirect doctest
            2-adic valuation
        """
        from sage.rings.integer_ring import ZZ
        if prime is None:
            raise ValueError("prime must be specified for this ring")
        from sage.rings.valuation.valuation import DiscretePseudoValuation
        if isinstance(prime, DiscretePseudoValuation):
            prime = prime.uniformizer()
        if prime not in ZZ or not ZZ(prime).is_prime():
            raise ValueError("prime must be a prime in the integers but %s is not" % (prime,))
        return R, prime

    def create_key_for_local_ring(self, R, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: Qp(2).valuation()  # indirect doctest
            2-adic valuation
        """
        # We do not care much about the value of prime since there is only one
        # reasonable p-adic valuation here
        if prime is not None:
            if prime in R:
                if R(prime).valuation() <= 0:
                    raise ValueError("prime must be an element of positive valuation")
            elif prime(R.prime()) <= 0:
                raise ValueError("prime must be an element of positive valuation")

        return (R,)

    def create_key_and_extra_args_for_number_field(self, R, prime, approximants):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``prime``.

        EXAMPLES::

            sage: GaussianIntegers().valuation(2)  # indirect doctest                   # needs sage.rings.number_field
            2-adic valuation
        """
        K, L, G = self._normalize_number_field_data(R)

        from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
        from sage.rings.valuation.valuation import DiscretePseudoValuation
        if isinstance(prime, DiscretePseudoValuation):
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, prime, prime, approximants=approximants)
        elif prime in K:
            return self.create_key_and_extra_args_for_number_field_from_valuation(R, K.valuation(prime), prime, approximants=approximants)
        elif prime in L or isinstance(prime, NumberFieldFractionalIdeal):
            return self.create_key_and_extra_args_for_number_field_from_ideal(R, L.fractional_ideal(prime), prime)
        else:
            raise ValueError("prime must be a discrete pseudo-valuation, a prime in the base ring, or a fractional ideal")

    def create_key_and_extra_args_for_number_field_from_valuation(self, R, v, prime, approximants):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``v``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args`, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: GaussianIntegers().valuation(ZZ.valuation(2))  # indirect doctest     # needs sage.rings.number_field
            2-adic valuation

        TESTS:

        We can extend to the field of fractions of a quotient ring::

            sage: R.<x> = ZZ[]
            sage: S = R.quo(x^2 + 1)
            sage: v = valuations.pAdicValuation(S, 2)                                   # needs sage.geometry.polyhedron
            sage: R.<x> = QQ[]
            sage: S = R.quo(x^2 + 1)
            sage: v = valuations.pAdicValuation(S, v)                                   # needs sage.geometry.polyhedron
        """
        K, L, G = self._normalize_number_field_data(R)

        if v.domain().is_subring(G.parent()):
            # v is defined on a subring of K[x].
            # We try to lift v to a pseudo-valuation on K[x].
            if _fraction_field(v.domain()) is not _fraction_field(G.parent()):
                # First, we lift valuations defined on subrings of K to
                # valuations on K[x].
                if v.domain().is_subring(K):
                    if v.domain() is not K:
                        v = K.valuation(v)
                    from sage.rings.valuation.gauss_valuation import GaussValuation
                    v = GaussValuation(G.parent(), v)
            if v.domain() != G.parent():
                # Then, we lift valuations defined on polynomial rings which are
                # subrings of K[x] to K[x]
                v = v.extension(G.parent())
        elif _fraction_field(v.domain()) == L:
            # v is defined on a ring whose field of fractions is L
            v = v._base_valuation._initial_approximation.change_domain(G.parent())
        else:
            raise NotImplementedError("cannot rewrite %r which is defined on %r as a pseudo-valuation on %r" % (v, v.domain(), G.parent()))

        assert (v.domain() is G.parent())

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        vK = v.restriction(v.domain().base_ring()).extension(K)
        if approximants is None:
            approximants = vK.mac_lane_approximants(G, require_incomparability=True)
        approximants = [approximant.extension(v.domain()) for approximant in approximants]
        approximant = vK.mac_lane_approximant(G, v, approximants=tuple(approximants))

        return (R, approximant), {'approximants': approximants}

    def create_key_and_extra_args_for_number_field_from_ideal(self, R, I, prime):
        r"""
        Create a unique key identifying the valuation of ``R`` with respect to
        ``I``.

        .. NOTE::

            ``prime``, the original parameter that was passed to
            :meth:`create_key_and_extra_args`, is only used to provide more
            meaningful error messages

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: GaussianIntegers().valuation(GaussianIntegers().number_field().fractional_ideal(2))  # indirect doctest
            2-adic valuation

        TESTS:

        Verify that :issue:`28976` has been resolved::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^6 - 18*x^4 - 24*x^3 + 27*x^2 + 36*x - 6)
            sage: I = K.fractional_ideal((2, -7/44*a^5 + 19/44*a^4 + 87/44*a^3 - 87/44*a^2 - 5/2*a + 39/22))
            sage: I.norm()
            2
            sage: I in K.primes_above(2)
            True
            sage: K.valuation(I)
            [ 2-adic valuation, v(x + 1) = 1/2 ]-adic valuation

        ::

            sage: # needs sage.rings.number_field
            sage: K.<a, b> = NumberField([x^2 - 2, x^2 + x + 1])
            sage: K.valuation(2)
            2-adic valuation
        """
        K, L, G = self._normalize_number_field_data(R)

        # To obtain uniqueness of p-adic valuations, we need a canonical
        # description of v. We consider all extensions of vK to L and select
        # the one approximated by v.
        # Of course, this only works if I comes from a single prime downstairs.
        p = I.relative_norm()
        F = p.factor()
        if len(F) != 1:
            raise ValueError("%r does not lie over a single prime of %r" % (I, K))
        vK = K.valuation(F[0][0])
        approximants = vK.mac_lane_approximants(G, require_incomparability=True)

        candidates = approximants[:]

        # The correct approximant has v(g) > 0 for all g in the ideal.
        # Unfortunately, the generators of I, even though defined over K have
        # their polynomial() defined over the rationals so we need to turn them
        # into polynomials over K[x] explicitly.
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        gens = I.gens()
        gens = [PolynomialRing(K, 'x')(list(g.vector())) for g in gens]

        # Refine candidates until we can detect which valuation corresponds to the ideal I
        while True:
            assert any(candidates), "the defining polynomial of the extension factored but we still could not figure out which valuation corresponds to the given ideal"

            match = [i for (i, v) in enumerate(candidates) if v and all(v(g) > 0 for g in gens)]

            if len(match) > 1:
                raise ValueError("%s does not single out a unique extension of %s to %s" % (prime, vK, L))
            if len(match) == 1:
                return (R, approximants[match[0]]), {'approximants': approximants}

            # We refine candidates which increases v(g) for all g in I;
            # however, we cannot augment the valuations which are already at
            # v(G) = +∞ which we ignore by setting them to None.
            candidates = [v.mac_lane_step(G)[0] if v and v.is_discrete_valuation() else None for v in candidates]

    def _normalize_number_field_data(self, R):
        r"""
        Helper method which returns the defining data of the number field
        ``R``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K = R.quo(x^2 + 1)
            sage: valuations.pAdicValuation._normalize_number_field_data(K)
            (Rational Field,
             Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 1,
             x^2 + 1)
        """
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        from sage.rings.number_field.number_field_base import NumberField
        if isinstance(R.fraction_field(), NumberField):
            L = R.fraction_field()
            G = L.relative_polynomial()
            K = L.base_ring()
        elif isinstance(R, PolynomialQuotientRing_generic):
            from sage.categories.number_fields import NumberFields
            if R.base_ring().fraction_field() not in NumberFields():
                raise NotImplementedError("cannot normalize quotients over %r" % (R.base_ring(),))
            L = R.fraction_field()
            K = R.base_ring().fraction_field()
            G = R.modulus().change_ring(K)
        else:
            raise NotImplementedError("cannot normalize %r" % (R,))

        return K, L, G

    def create_object(self, version, key, **extra_args):
        r"""
        Create a `p`-adic valuation from ``key``.

        EXAMPLES::

            sage: ZZ.valuation(5)  # indirect doctest
            5-adic valuation
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        from sage.rings.padics.padic_generic import pAdicGeneric
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        from sage.rings.number_field.number_field_base import NumberField
        R = key[0]
        parent = DiscretePseudoValuationSpace(R)
        if isinstance(R, pAdicGeneric):
            assert (len(key) == 1)
            return parent.__make_element_class__(pAdicValuation_padic)(parent)
        elif R is ZZ or R is QQ:
            prime = key[1]
            assert (len(key) == 2)
            return parent.__make_element_class__(pAdicValuation_int)(parent, prime)
        else:
            v = key[1]
            approximants = extra_args['approximants']
            parent = DiscretePseudoValuationSpace(R)
            K = R.fraction_field()
            if isinstance(K, NumberField):
                G = K.relative_polynomial()
            elif isinstance(R, PolynomialQuotientRing_generic):
                G = R.modulus()
            else:
                raise NotImplementedError
            return parent.__make_element_class__(pAdicFromLimitValuation)(parent, v, G.change_ring(R.base_ring()), approximants)


pAdicValuation = PadicValuationFactory("sage.rings.padics.padic_valuation.pAdicValuation")


class pAdicValuation_base(DiscreteValuation):
    r"""
    Abstract base class for `p`-adic valuations.

    INPUT:

    - ``ring`` -- an integral domain

    - ``p`` -- a rational prime over which this valuation lies, not
      necessarily a uniformizer for the valuation

    EXAMPLES::

        sage: ZZ.valuation(3)
        3-adic valuation

        sage: QQ.valuation(5)
        5-adic valuation

     For `p`-adic rings, ``p`` has to match the `p` of the ring. ::

        sage: v = valuations.pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        ValueError: prime must be an element of positive valuation

    TESTS::

        sage: TestSuite(ZZ.valuation(3)).run()  # long time                             # needs sage.geometry.polyhedron
        sage: TestSuite(QQ.valuation(5)).run()  # long time                             # needs sage.geometry.polyhedron
        sage: TestSuite(Zp(5).valuation()).run()        # long time                     # needs sage.geometry.polyhedron
    """
    def __init__(self, parent, p):
        r"""
        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_base
            sage: isinstance(ZZ.valuation(3), pAdicValuation_base)
            True
        """
        DiscreteValuation.__init__(self, parent)

        from sage.rings.integer_ring import ZZ
        self._p = ZZ(p)

    def p(self):
        r"""
        Return the `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: GaussianIntegers().valuation(2).p()                                   # needs sage.rings.number_field
            2
        """
        return self._p

    def reduce(self, x):
        r"""
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT: an element of the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`

        EXAMPLES::

            sage: v = ZZ.valuation(3)
            sage: v.reduce(4)
            1
        """
        x = self.domain().coerce(x)

        if self(x) < 0:
            raise ValueError("reduction is only defined for elements of nonnegative valuation")

        return self.residue_field()(x)

    def lift(self, x):
        r"""
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`

        EXAMPLES::

            sage: v = ZZ.valuation(3)
            sage: xbar = v.reduce(4)
            sage: v.lift(xbar)
            1
        """
        x = self.residue_field().coerce(x)

        return self.domain()(x)

    def is_unramified(self, G, include_steps=False, assume_squarefree=False):
        r"""
        Return whether ``G`` defines a single unramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- boolean (default: ``False``); whether to
          include the approximate valuations that were used to determine the
          result in the return value

        - ``assume_squarefree`` -- boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        EXAMPLES:

        We consider an extension as unramified if its ramification index is 1.
        Hence, a trivial extension is unramified::

            sage: R.<x> = QQ[]
            sage: v = QQ.valuation(2)
            sage: v.is_unramified(x)
            True

        If ``G`` remains irreducible in reduction, then it defines an
        unramified extension::

            sage: v.is_unramified(x^2 + x + 1)
            True

        However, even if ``G`` factors, it might define an unramified
        extension::

            sage: v.is_unramified(x^2 + 2*x + 4)                                        # needs sage.geometry.polyhedron
            True
        """
        R = G.parent()

        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        if not isinstance(R, PolynomialRing_generic) or R.base_ring() is not self.domain() or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.valuation.gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.E() > 1:
                ret = False
                break
            if v.F() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = True
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next) > 1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def is_totally_ramified(self, G, include_steps=False, assume_squarefree=False):
        r"""
        Return whether ``G`` defines a single totally ramified extension of the
        completion of the domain of this valuation.

        INPUT:

        - ``G`` -- a monic squarefree polynomial over the domain of this valuation

        - ``include_steps`` -- boolean (default: ``False``); where to include
          the valuations produced during the process of checking whether ``G``
          is totally ramified in the return value

        - ``assume_squarefree`` -- boolean (default: ``False``); whether to
          assume that ``G`` is square-free over the completion of the domain of
          this valuation. Setting this to ``True`` can significantly improve
          the performance.

        ALGORITHM:

        This is a simplified version of :meth:`sage.rings.valuation.valuation.DiscreteValuation.mac_lane_approximants`.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: k = Qp(5,4)
            sage: v = k.valuation()
            sage: R.<x> = k[]
            sage: G = x^2 + 1
            sage: v.is_totally_ramified(G)                                              # needs sage.geometry.polyhedron
            False
            sage: G = x + 1
            sage: v.is_totally_ramified(G)
            True
            sage: G = x^2 + 2
            sage: v.is_totally_ramified(G)
            False
            sage: G = x^2 + 5
            sage: v.is_totally_ramified(G)                                              # needs sage.geometry.polyhedron
            True
            sage: v.is_totally_ramified(G, include_steps=True)                          # needs sage.geometry.polyhedron
            (True, [Gauss valuation induced by 5-adic valuation, [ Gauss valuation induced by 5-adic valuation, v((1 + O(5^4))*x) = 1/2 ]])

        We consider an extension as totally ramified if its ramification index
        matches the degree. Hence, a trivial extension is totally ramified::

            sage: R.<x> = QQ[]
            sage: v = QQ.valuation(2)
            sage: v.is_totally_ramified(x)
            True

        TESTS:

        An example that Sebastian Pauli used at Sage Days 87::

            sage: R = ZpFM(3, 20)
            sage: S.<x> = R[]
            sage: f = x^9 + 9*x^2 + 3
            sage: R.valuation().is_totally_ramified(f)                                  # needs sage.geometry.polyhedron
            True
        """
        R = G.parent()

        from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
        if not isinstance(R, PolynomialRing_generic) or R.base_ring() is not self.domain() or not G.is_monic():
            raise ValueError("G must be a monic univariate polynomial over the domain of this valuation")
        if not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.valuation.gauss_valuation import GaussValuation

        steps = [ GaussValuation(R, self) ]
        while True:
            v = steps[-1]
            if v.F() > 1:
                ret = False
                break
            if v.E() == G.degree():
                ret = True
                break

            assert v(G) is not infinity
            if v.is_key(G):
                ret = False
                break

            next = v.mac_lane_step(G, assume_squarefree=True)
            if len(next) > 1:
                ret = False
                break
            steps.append(next[0])

        if include_steps:
            return ret, steps
        else:
            return ret

    def change_domain(self, ring):
        r"""
        Change the domain of this valuation to ``ring`` if possible.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.change_domain(QQ).domain()
            Rational Field
        """
        return pAdicValuation(ring, self.p())

    def _extensions_to_quotient(self, ring, approximants=None):
        r"""
        Return the extensions of this valuation to an integral quotient over
        the domain of this valuation.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: QQ.valuation(2)._extensions_to_quotient(R.quo(x^2 + x + 1))
            [2-adic valuation]
        """
        approximants = approximants or self.mac_lane_approximants(ring.modulus().change_ring(self.domain()), assume_squarefree=True, require_incomparability=True)
        return [pAdicValuation(ring, approximant, approximants) for approximant in approximants]

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.extensions(GaussianIntegers())                                      # needs sage.rings.number_field
            [2-adic valuation]

        TESTS::

            sage: # needs sage.rings.number_field
            sage: R.<a> = QQ[]
            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = QQ.extension(x^3 - 2)
            sage: R.<b> = L[]
            sage: M.<b> = L.extension(b^2 + 2*b + a)
            sage: M.valuation(2)
            2-adic valuation

        Check that we can extend to a field written as a quotient::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K.<a> = QQ.extension(x^2 + 1)
            sage: R.<y> = K[]
            sage: L.<b> = R.quo(x^2 + a)
            sage: QQ.valuation(2).extensions(L)
            [2-adic valuation]

        A case where there was at some point an internal error in the
        approximants code::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: L.<a> = NumberField(x^4 + 2*x^3 + 2*x^2 + 8)
            sage: QQ.valuation(2).extensions(L)
            [[ 2-adic valuation, v(x + 2) = 3/2 ]-adic valuation,
             [ 2-adic valuation, v(x) = 1/2 ]-adic valuation]

        A case where the extension was incorrect at some point::

            sage: # needs sage.rings.number_field
            sage: v = QQ.valuation(2)
            sage: L.<a> = NumberField(x^2 + 2)
            sage: M.<b> = L.extension(x^2 + 1)
            sage: w = v.extension(L).extension(M)
            sage: w(w.uniformizer())
            1/4

        A case where the extensions could not be separated at some point::

            sage: # needs sage.rings.number_field
            sage: v = QQ.valuation(2)
            sage: R.<x> = QQ[]
            sage: F = x^48 + 120*x^45 + 56*x^42 + 108*x^36 + 32*x^33 + 40*x^30 + 48*x^27 + 80*x^24 + 112*x^21 + 96*x^18 + 96*x^15 + 24*x^12 + 96*x^9 + 16*x^6 + 96*x^3 + 68
            sage: L.<a> = QQ.extension(F)
            sage: v.extensions(L)
            [[ 2-adic valuation, v(x) = 1/24, v(x^24 + 4*x^18 + 10*x^12 + 12*x^6 + 8*x^3 + 6) = 29/8 ]-adic valuation,
             [ 2-adic valuation, v(x) = 1/24, v(x^24 + 4*x^18 + 2*x^12 + 12*x^6 + 8*x^3 + 6) = 29/8 ]-adic valuation]
        """
        if self.domain() is ring:
            return [self]
        domain_fraction_field = _fraction_field(self.domain())
        if domain_fraction_field is not self.domain():
            if domain_fraction_field.is_subring(ring):
                return pAdicValuation(domain_fraction_field, self).extensions(ring)
        if self.domain().is_subring(ring):
            from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
            if isinstance(ring, PolynomialQuotientRing_generic):
                if isinstance(self.domain(), PolynomialQuotientRing_generic):
                    if self.domain().modulus() == ring.modulus():
                        base_extensions = self._base_valuation.extensions(self._base_valuation.domain().change_ring(self._base_valuation.domain().base_ring().fraction_field()))
                        return [pAdicValuation(ring, base._initial_approximation) for base in base_extensions]
                if ring.base_ring() is self.domain():
                    from sage.categories.integral_domains import IntegralDomains
                    if ring in IntegralDomains():
                        return self._extensions_to_quotient(ring)
                elif self.domain().is_subring(ring.base_ring()):
                    return sum([w.extensions(ring) for w in self.extensions(ring.base_ring())], [])
            from sage.rings.number_field.number_field_base import NumberField
            if isinstance(ring.fraction_field(), NumberField):
                if ring.base_ring().fraction_field() is self.domain().fraction_field():
                    approximants = self.mac_lane_approximants(ring.fraction_field().relative_polynomial().change_ring(self.domain()), assume_squarefree=True, require_incomparability=True)
                    return [pAdicValuation(ring, approximant, approximants) for approximant in approximants]
                if ring.base_ring() is not ring and self.domain().is_subring(ring.base_ring()):
                    return sum([w.extensions(ring) for w in self.extensions(ring.base_ring())], [])
        return super().extensions(ring)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: v = GaussianIntegers().valuation(2)                                   # needs sage.rings.number_field
            sage: v.restriction(ZZ)                                                     # needs sage.rings.number_field
            2-adic valuation
        """
        if ring is self.domain():
            return self

        if not ring.is_subring(self.domain()):
            raise ValueError("ring must be a subring of the domain of this valuation but %r is not a subring of %r" % (ring, self.domain()))

        return pAdicValuation(ring, self.p())

    @cached_method
    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        EXAMPLES::

            sage: v = GaussianIntegers().valuation(2)                                   # needs sage.rings.number_field
            sage: v.value_semigroup()                                                   # needs sage.rings.number_field
            Additive Abelian Semigroup generated by 1/2
        """
        from sage.categories.fields import Fields
        v = self(self.uniformizer())
        if self.domain() in Fields():
            return DiscreteValueSemigroup([-v,v])
        else:
            return DiscreteValueSemigroup([v])


class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a complete `p`-adic ring.

    INPUT:

    - ``R`` -- a `p`-adic ring

    EXAMPLES::

        sage: v = Qp(2).valuation(); v  # indirect doctest
        2-adic valuation

    TESTS::

        sage: TestSuite(v).run()                # long time                             # needs sage.geometry.polyhedron
    """
    def __init__(self, parent):
        """
        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_padic
            sage: isinstance(Qp(2).valuation(), pAdicValuation_padic)
            True
        """
        pAdicValuation_base.__init__(self, parent, parent.domain().prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT: an element of the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`

        EXAMPLES::

            sage: R = Zp(3)
            sage: Zp(3).valuation().reduce(R(4))
            1
        """
        x = self.domain().coerce(x)
        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field` to the domain of this
        valuation.

        INPUT:

        - ``x`` -- an element of the residue field of this valuation

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = R.valuation()
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)
        """
        x = self.residue_field().coerce(x)
        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation.

        EXAMPLES::

            sage: v = Zp(3).valuation()
            sage: v.uniformizer()
            3 + O(3^21)
        """
        return self.domain().uniformizer()

    def element_with_valuation(self, v):
        """
        Return an element of valuation ``v``.

        INPUT:

        - ``v`` -- an element of the :meth:`pAdicValuation_base.value_semigroup` of this valuation

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = R.valuation()
            sage: v.element_with_valuation(3)
            3^3 + O(3^23)

            sage: # needs sage.libs.ntl
            sage: K = Qp(3)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 + 3*y + 3)
            sage: L.valuation().element_with_valuation(3/2)
            y^3 + O(y^43)
        """
        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        v = QQ(v)
        if v not in self.value_semigroup():
            raise ValueError("%r is not in the value semigroup of %r" % (v, self))
        v = ZZ(v * self.domain().absolute_e())
        return self.domain().one() << v

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: ZZ.valuation(3)._repr_()
            '3-adic valuation'
        """
        return "%s-adic valuation" % (self.p())

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<y> = K[]                                                           # needs sage.libs.ntl
            sage: L.<y> = K.extension(y^2 - 3)                                          # needs sage.libs.ntl
            sage: L.valuation()(3)                                                      # needs sage.libs.ntl
            1
        """
        return x.ordp()

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: Qq(9, names='a').valuation().residue_ring()                           # needs sage.libs.ntl
            Finite Field in a0 of size 3^2
        """
        return self.domain().residue_field()

    def shift(self, x, s):
        r"""
        Shift ``x`` in its expansion with respect to :meth:`uniformizer` by
        ``s`` "digits".

        For nonnegative ``s``, this just returns ``x`` multiplied by a
        power of the uniformizer `\pi`.

        For negative ``s``, it does the same but when not over a field, it
        drops coefficients in the `\pi`-adic expansion which have negative
        valuation.

        EXAMPLES::

            sage: R = ZpCA(2)
            sage: v = R.valuation()
            sage: v.shift(R.one(), 1)
            2 + O(2^20)
            sage: v.shift(R.one(), -1)
            O(2^19)

            sage: # needs sage.libs.ntl sage.rings.padics
            sage: S.<y> = R[]
            sage: S.<y> = R.extension(y^3 - 2)
            sage: v = S.valuation()
            sage: v.shift(1, 5)
            y^5 + O(y^60)
        """
        x = self.domain().coerce(x)
        s = self.value_group()(s)
        return x << s

    def simplify(self, x, error=None, force=False):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        - ``error`` -- a rational, infinity, or ``None`` (default: ``None``),
          the error allowed to introduce through the simplification

        - ``force`` -- ignored

        EXAMPLES::

            sage: R = Zp(2)
            sage: v = R.valuation()
            sage: v.simplify(6)
            2 + O(2^21)
            sage: v.simplify(6, error=0)
            0
        """
        x = self.domain().coerce(x)

        if error is None:
            error = self(x)
        from sage.rings.infinity import infinity
        if error is infinity:
            return x
        # we need to scale by the ramification index because p-adics use a
        # different normalization
        normalized_error = (error / self.value_group().gen()).ceil()
        return x.add_bigoh(normalized_error + 1).lift_to_precision()


class pAdicValuation_int(pAdicValuation_base):
    r"""
    A `p`-adic valuation on the integers or the rationals.

    EXAMPLES::

        sage: v = ZZ.valuation(3); v
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run()                # long time                             # needs sage.geometry.polyhedron
    """
    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: ZZ.valuation(3)._repr_()
            '3-adic valuation'
        """
        return "%s-adic valuation" % (self.p())

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        EXAMPLES::

            sage: ZZ.valuation(3)(9)
            2
        """
        if x.is_zero():
            # x.valuation() is a factor 10 slower when computing the valuation
            # of a rational zero than when computing the valuation of another
            # small rational. Special casing this is a factor 100 faster.
            return infinity
        return x.valuation(self._p)

    def uniformizer(self):
        """
        Return a uniformizer of this `p`-adic valuation, i.e., `p` as an
        element of the domain.

        EXAMPLES::

            sage: v = ZZ.valuation(3)
            sage: v.uniformizer()
            3
        """
        return self.domain()(self.p())

    def residue_ring(self):
        """
        Return the residue field of this valuation.

        EXAMPLES::

            sage: v = ZZ.valuation(3)
            sage: v.residue_ring()
            Finite Field of size 3
        """
        from sage.rings.finite_rings.finite_field_constructor import GF
        return GF(self.p())

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater than or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: w = valuations.TrivialValuation(ZZ)
            sage: v >= w
            True
        """
        if other.is_trivial():
            return other.is_discrete_valuation()
        if isinstance(other, pAdicValuation_int):
            return self.p() == other.p()
        return super(pAdicValuation_base, self)._ge_(other)

    def _relative_size(self, x):
        r"""
        Return an estimate on the coefficient size of ``x``.

        The number returned is an estimate on the factor between the number of
        bits used by ``x`` and the minimal number of bits used by an element
        congruent to ``x``.

        This is used by :meth:`simplify` to decide whether simplification of
        coefficients is going to lead to a significant shrinking of the
        coefficients of ``x``.

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v._relative_size(2)
            1
            sage: v._relative_size(2**20)
            11
        """
        x = self.domain().coerce(x)
        return (x.numerator().nbits() + x.denominator().nbits())//self.p().nbits()

    def simplify(self, x, error=None, force=False, size_heuristic_bound=32):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        - ``error`` -- a rational, infinity, or ``None`` (default: ``None``),
          the error allowed to introduce through the simplification

        - ``force`` -- ignored

        - ``size_heuristic_bound`` -- when ``force`` is not set, the expected
          factor by which the ``x`` need to shrink to perform an actual
          simplification (default: 32)

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: v.simplify(6, force=True)
            2
            sage: v.simplify(6, error=0, force=True)
            0

        In this example, the usual rational reconstruction misses a good answer
        for some moduli (because the absolute value of the numerator is not
        bounded by the square root of the modulus)::

            sage: v = QQ.valuation(2)
            sage: v.simplify(110406, error=16, force=True)
            562/19
            sage: Qp(2, 16)(110406).rational_reconstruction()
            Traceback (most recent call last):
            ...
            ArithmeticError: rational reconstruction of 55203 (mod 65536) does not exist
        """
        if not force and self._relative_size(x) <= size_heuristic_bound:
            return x

        x = self.domain().coerce(x)

        v = self(x)
        if error is None:
            error = v
        from sage.rings.infinity import infinity
        if error is infinity:
            return x
        if error < v:
            return self.domain().zero()

        from sage.rings.rational_field import QQ
        from sage.rings.padics.factory import Qp
        precision_ring = Qp(self.p(), QQ(error).floor() + 1 - v)
        reduced = precision_ring(x)
        lift = (reduced >> v).lift()
        best = self.domain()(lift) * self.p()**v

        if self._relative_size(x) < self._relative_size(best):
            best = x

        # We implement a modified version of the usual rational reconstruction
        # algorithm (based on the extended Euclidean algorithm) here. We do not
        # get the uniqueness properties but we do not need them actually.
        # This is certainly slower than the implementation in Cython.
        from sage.categories.fields import Fields
        m = self.p()**(QQ(error).floor() + 1 - v)
        if self.domain() in Fields():
            r = (m, lift)
            s = (0, 1)
            while r[1]:
                qq, rr = r[0].quo_rem(r[1])
                r = r[1], rr
                s = s[1], s[0] - qq*s[1]
                from sage.arith.misc import GCD as gcd
                if s[1] != 0 and gcd(s[1], r[1]) == 1:
                    rational = self.domain()(r[1]) / self.domain()(s[1]) * self.p()**v
                    if self._relative_size(rational) < self._relative_size(best):
                        best = rational

        assert (self(x-best) > error)

        return best

    def inverse(self, x, precision):
        r"""
        Return an approximate inverse of ``x``.

        The element returned is such that the product differs from 1 by an
        element of valuation at least ``precision``.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        - ``precision`` -- a rational or infinity

        EXAMPLES::

            sage: v = ZZ.valuation(2)
            sage: x = 3
            sage: y = v.inverse(3, 2); y
            3
            sage: x*y - 1
            8

        This might not be possible for elements of positive valuation::

            sage: v.inverse(2, 2)
            Traceback (most recent call last):
            ...
            ValueError: element has no approximate inverse in this ring

        Unless the precision is very small::

            sage: v.inverse(2, 0)
            1
        """
        if not x.is_zero():
            y = ~x
            if y in self.domain():
                return self.domain()(y)
        if precision <= 0:
            return self.domain().one()

        from sage.rings.infinity import infinity
        if self(x) > 0 or precision is infinity:
            raise ValueError("element has no approximate inverse in this ring")

        from sage.rings.integer_ring import ZZ
        from sage.rings.rational_field import QQ
        return self.domain()(ZZ(x).inverse_mod(self.p() ** QQ(precision).ceil()))


class pAdicFromLimitValuation(FiniteExtensionFromLimitValuation, pAdicValuation_base):
    r"""
    A `p`-adic valuation on a number field or a subring thereof, i.e., a
    valuation that extends the `p`-adic valuation on the integers.

    EXAMPLES::

        sage: v = GaussianIntegers().valuation(3); v                                    # needs sage.rings.number_field
        3-adic valuation

    TESTS::

        sage: TestSuite(v).run(skip='_test_shift')      # long time                     # needs sage.rings.number_field

    The ``_test_shift`` test fails because the parent of the shift is
    incorrect, see :issue:`23971`::

        sage: v.shift(1, -1).parent()                                                   # needs sage.rings.number_field
        Number Field in I with defining polynomial x^2 + 1 with I = 1*I
    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        TESTS::

            sage: v = GaussianIntegers().valuation(3)                                   # needs sage.rings.number_field
            sage: from sage.rings.padics.padic_valuation import pAdicFromLimitValuation
            sage: isinstance(v, pAdicFromLimitValuation)                                # needs sage.rings.number_field
            True
        """
        FiniteExtensionFromLimitValuation.__init__(self, parent, approximant, G, approximants)
        pAdicValuation_base.__init__(self, parent, approximant.restriction(approximant.domain().base_ring()).p())

    def _to_base_domain(self, f):
        r"""
        Return ``f``, an element of the underlying limit valuation, as an
        element of the domain of this valuation.

        EXAMPLES::

            sage: v = GaussianIntegers().valuation(3)                                   # needs sage.rings.number_field
            sage: I = GaussianIntegers().fraction_field().gen()                         # needs sage.rings.number_field
            sage: v._to_base_domain(I)                                                  # needs sage.rings.number_field
            x

        TESTS:

        Check that this also works for relative extensions::

            sage: # needs sage.rings.number_field
            sage: v = QQ.valuation(2)
            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = NumberField(x^2 + 2)
            sage: M.<b> = L.extension(x^2 + 1)
            sage: w = v.extension(L).extension(M)
            sage: w._to_base_domain(b)
            x
        """
        polynomial = f.lift()
        return polynomial(self._base_valuation.domain().gen())

    def _from_base_domain(self, f):
        r"""
        Return ``f``, an element of the domain of this valuation, as an element
        of the domain of the underlying limit valuation.

        EXAMPLES::

            sage: v = GaussianIntegers().valuation(3)                                   # needs sage.rings.number_field
            sage: v._from_base_domain(v._base_valuation.domain().gen())                 # needs sage.rings.number_field
            I
        """
        return self.domain()(f)

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: v = GaussianIntegers().valuation(3)                                   # needs sage.rings.number_field
            sage: v.extensions(v.domain().fraction_field())                             # needs sage.rings.number_field
            [3-adic valuation]
        """
        if ring is self.domain().fraction_field():
            if self.domain() is not self.domain().fraction_field():
                G = ring.relative_polynomial()
                approximant = self._base_valuation.change_domain(G.parent())._initial_approximation
                return [pAdicValuation(ring, approximant)]
        return super().extensions(ring)


def _fraction_field(ring):
    r"""
    Return a fraction field of ``ring``.

    EXAMPLES:

    This works around some annoyances with ``ring.fraction_field()``::

        sage: R.<x> = ZZ[]
        sage: S = R.quo(x^2 + 1)
        sage: S.fraction_field()
        Fraction Field of Univariate Quotient Polynomial Ring in xbar over Integer Ring with modulus x^2 + 1

        sage: from sage.rings.padics.padic_valuation import _fraction_field
        sage: _fraction_field(S)
        Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 1
    """
    from sage.categories.fields import Fields
    if ring in Fields():
        return ring

    from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
    if isinstance(ring, PolynomialQuotientRing_generic):
        from sage.categories.integral_domains import IntegralDomains
        if ring in IntegralDomains():
            return ring.base().change_ring(ring.base_ring().fraction_field()).quo(ring.modulus())
    return ring.fraction_field()
