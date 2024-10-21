# sage_setup: distribution = sagemath-pari
r"""
Inductive valuations on polynomial rings

This module provides functionality for inductive valuations, i.e., finite
chains of :mod:`augmented valuations <sage.rings.valuation.augmented_valuation>` on top of a :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>`.

AUTHORS:

- Julian Rüth (2016-11-01): initial version

EXAMPLES:

A :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>` is an example of an inductive valuation::

    sage: R.<x> = QQ[]
    sage: v = GaussValuation(R, QQ.valuation(2))

Generally, an inductive valuation is an augmentation of an inductive valuation,
i.e., a valuation that was created from a Gauss valuation in a finite number of
augmentation steps::

    sage: w = v.augmentation(x, 1)
    sage: w = w.augmentation(x^2 + 2*x + 4, 3)

REFERENCES:

Inductive valuations are originally discussed in [Mac1936I]_ and [Mac1936II]_.
An introduction is also given in Chapter 4 of [Rüt2014]_.
"""
# ****************************************************************************
#       Copyright (C) 2016-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .valuation import DiscreteValuation, InfiniteDiscretePseudoValuation
from .developing_valuation import DevelopingValuation

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method


class InductiveValuation(DevelopingValuation):
    r"""
    Abstract base class for iterated :mod:`augmented valuations <sage.rings.valuation.augmented_valuation>` on top of a :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>`.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, QQ.valuation(5))

    TESTS::

        sage: TestSuite(v).run()                # long time                             # needs sage.geometry.polyhedron
    """
    def is_equivalence_unit(self, f, valuations=None):
        r"""
        Return whether the polynomial ``f`` is an equivalence unit, i.e., an
        element of :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.effective_degree`
        zero (see [Mac1936II]_ p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R = Zp(2,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_unit(x)
            False
            sage: v.is_equivalence_unit(S.zero())
            False
            sage: v.is_equivalence_unit(2*x + 1)
            True
        """
        f = self.domain().coerce(f)

        if f.is_zero():
            return False
        return self.effective_degree(f, valuations=valuations) == 0

    def equivalence_reciprocal(self, f, coefficients=None, valuations=None, check=True):
        r"""
        Return an equivalence reciprocal of ``f``.

        An equivalence reciprocal of `f` is a polynomial `h` such that `f\cdot
        h` is equivalent to 1 modulo this valuation (see [Mac1936II]_ p.497.)

        INPUT:

        - ``f`` -- a polynomial in the domain of this valuation which is an
          :meth:`equivalence_unit`

        - ``coefficients`` -- the coefficients of ``f`` in the :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.phi`-adic
          expansion if known (default: ``None``)

        - ``valuations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``check`` -- whether or not to check the validity of ``f`` (default:
          ``True``)

        .. WARNING::

            This method may not work over `p`-adic rings due to problems with
            the xgcd implementation there.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R = Zp(3,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: f = 3*x + 2
            sage: h = v.equivalence_reciprocal(f); h
            2 + O(3^5)
            sage: v.is_equivalent(f*h, 1)
            True

        In an extended valuation over an extension field::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: f = 2*x + u
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        Extending the valuation once more::

            sage: # needs sage.libs.ntl
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: h = v.equivalence_reciprocal(f); h
            (u + 1) + O(2^5)
            sage: v.is_equivalent(f*h, 1)
            True

        TESTS:

        A case that caused problems at some point::

            sage: # needs sage.libs.ntl
            sage: K = Qp(2, 4)
            sage: R.<x> = K[]
            sage: L.<a> = K.extension(x^4 + 4*x^3 + 6*x^2 + 4*x + 2)
            sage: R.<t> = L[]
            sage: v = GaussValuation(R)
            sage: w = v.augmentation(t + 1, 5/16)
            sage: w = w.augmentation(t^4 + (a^8 + a^12 + a^14 + a^16 + a^17 + a^19 + a^20 + a^23)*t^3 + (a^6 + a^9 + a^13 + a^15 + a^18 + a^19 + a^21)*t^2 + a^10*t + 1 + a^4 + a^5 + a^8 + a^13 + a^14 + a^15, 17/8)
            sage: f = a^-15*t^2 + (a^-11 + a^-9 + a^-6 + a^-5 + a^-3 + a^-2)*t + a^-15
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1
            sage: f = f.parent()([f[0], f[1].add_bigoh(1), f[2]])
            sage: f_ = w.equivalence_reciprocal(f)
            sage: w.reduce(f*f_)
            1
        """
        f = self.domain().coerce(f)

        if check:
            if coefficients is None:
                coefficients = list(self.coefficients(f))
            if valuations is None:
                valuations = list(self.valuations(f, coefficients=coefficients))
            if not self.is_equivalence_unit(f, valuations=valuations):
                raise ValueError("f must be an equivalence unit but %r is not" % (f,))

        if coefficients is None:
            e0 = next(self.coefficients(f))
        else:
            e0 = coefficients[0]

        # f is an equivalence unit, its valuation is given by the constant coefficient
        if valuations is None:
            vf = self(e0)
        else:
            vf = valuations[0]

        e0 = self.simplify(e0, error=vf)
        s_ = self.equivalence_unit(-vf)
        residue = self.reduce(e0 * s_)
        if not isinstance(self, FinalInductiveValuation):
            assert residue.is_constant()
            residue = residue[0]
        h = self.lift(~residue) * s_

        h = self.simplify(h, -vf)

        # it might be the case that f*h has nonzero valuation because h has
        # insufficient precision, so we must not assert that here but only
        # until we lifted to higher precision

        # We do not actually need g*phi + h*e0 = 1, it is only important that
        # the RHS is 1 in reduction.
        # This allows us to do two things:
        # - we may lift h to arbitrary precision
        # - we can add anything which times e0 has positive valuation, e.g., we
        # may drop coefficients of positive valuation
        return h.map_coefficients(_lift_to_maximal_precision)

    @cached_method
    def mu(self):
        r"""
        Return the valuation of :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.phi`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v.mu()
            0
        """
        return self(self.phi())

    @abstract_method
    def equivalence_unit(self, s, reciprocal=False):
        """
        Return an equivalence unit of valuation ``s``.

        INPUT:

        - ``s`` -- an element of the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.value_group`

        - ``reciprocal`` -- boolean (default: ``False``); whether or not to
          return the equivalence unit as the :meth:`equivalence_reciprocal` of
          the equivalence unit of valuation ``-s``.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: S.<x> = Qp(3,5)[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_unit(2)
            3^2 + O(3^7)
            sage: v.equivalence_unit(-2)
            3^-2 + O(3^3)

        Note that this might fail for negative ``s`` if the domain is not
        defined over a field::

            sage: v = ZZ.valuation(2)
            sage: R.<x> = ZZ[]
            sage: w = GaussValuation(R, v)
            sage: w.equivalence_unit(1)
            2
            sage: w.equivalence_unit(-1)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value semigroup of this valuation
            but -1 is not in Additive Abelian Semigroup generated by 1
        """

    @abstract_method
    def augmentation_chain(self):
        r"""
        Return a list with the chain of augmentations down to the underlying
        :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>`.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.augmentation_chain()
            [Gauss valuation induced by 2-adic valuation]
        """

    @abstract_method
    def is_gauss_valuation(self):
        r"""
        Return whether this valuation is a Gauss valuation over the domain.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_gauss_valuation()
            True
        """

    @abstract_method
    def E(self):
        """
        Return the ramification index of this valuation over its underlying
        Gauss valuation.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.E()
            1
        """

    @abstract_method
    def F(self):
        """
        Return the residual degree of this valuation over its Gauss extension.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.F()
            1
        """

    @abstract_method
    def monic_integral_model(self, G):
        r"""
        Return a monic integral irreducible polynomial which defines the same
        extension of the base ring of the domain as the irreducible polynomial
        ``G`` together with maps between the old and the new polynomial.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v.monic_integral_model(5*x^2 + 1/2*x + 1/4)
            (Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
               Defn: x |--> 1/2*x,
             Ring endomorphism of Univariate Polynomial Ring in x over Rational Field
               Defn: x |--> 2*x,
             x^2 + 1/5*x + 1/5)
        """

    @abstract_method
    def element_with_valuation(self, s):
        r"""
        Return a polynomial of minimal degree with valuation ``s``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v.element_with_valuation(-2)
            1/4

        Depending on the base ring, an element of valuation ``s`` might not
        exist::

            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, ZZ.valuation(2))
            sage: v.element_with_valuation(-2)
            Traceback (most recent call last):
            ...
            ValueError: s must be in the value semigroup of this valuation
            but -2 is not in Additive Abelian Semigroup generated by 1
        """

    def _test_element_with_valuation_inductive_valuation(self, **options):
        r"""
        Test the correctness of :meth:`element_with_valuation`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v._test_element_with_valuation_inductive_valuation()
        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for s in tester.some_elements(self.value_group().some_elements()):
            try:
                R = self.element_with_valuation(s)
            except (ValueError, NotImplementedError):
                # this is often not possible unless the underlying ring of
                # constants is a field
                from sage.categories.fields import Fields
                if self.domain().base() not in Fields():
                    continue
                raise
            tester.assertEqual(self(R), s)
            if chain != [self]:
                base = chain[1]
                if s in base.value_group():
                    S = base.element_with_valuation(s)
                    tester.assertEqual(self(S), s)
                    tester.assertGreaterEqual(S.degree(), R.degree())

    def _test_EF(self, **options):
        r"""
        Test the correctness of :meth:`E` and :meth:`F`.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v._test_EF()
        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        for w, v in zip(chain, chain[1:]):
            from sage.rings.infinity import infinity
            from sage.rings.integer_ring import ZZ
            if w(w.phi()) is infinity:
                tester.assertEqual(w.E(), v.E())
            tester.assertIn(w.E(), ZZ)
            tester.assertIn(w.F(), ZZ)
            tester.assertGreaterEqual(w.E(), v.E())
            tester.assertGreaterEqual(w.F(), v.F())

    def _test_augmentation_chain(self, **options):
        r"""
        Test the correctness of :meth:`augmentation_chain`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_augmentation_chain()
        """
        tester = self._tester(**options)
        chain = self.augmentation_chain()
        tester.assertIs(chain[0], self)
        tester.assertTrue(chain[-1].is_gauss_valuation())
        for w, v in zip(chain, chain[1:]):
            tester.assertGreaterEqual(w, v)

    def _test_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_equivalence_unit()
        """
        tester = self._tester(**options)

        if self.is_gauss_valuation():
            value_group = self.value_group()
        else:
            value_group = self.augmentation_chain()[1].value_group()

        for s in tester.some_elements(value_group.some_elements()):
            try:
                R = self.equivalence_unit(s)
            except (ValueError, NotImplementedError):
                # this is often not possible unless the underlying ring of
                # constants is a field
                from sage.categories.fields import Fields
                if self.domain().base() not in Fields():
                    continue
                raise
            tester.assertIs(R.parent(), self.domain())
            tester.assertEqual(self(R), s)
            tester.assertTrue(self.is_equivalence_unit(R))

    def _test_is_equivalence_unit(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_unit`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_is_equivalence_unit()
        """
        tester = self._tester(**options)
        tester.assertFalse(self.is_equivalence_unit(self.phi()))

    def _test_equivalence_reciprocal(self, **options):
        r"""
        Test the correctness of :meth:`equivalence_reciprocal`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_equivalence_reciprocal()
        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if self.is_equivalence_unit(f):
                try:
                    g = self.equivalence_reciprocal(f)
                except (ValueError, NotImplementedError):
                    # this is often not possible unless the underlying ring of
                    # constants is a field
                    from sage.categories.fields import Fields
                    if self.domain().base() not in Fields():
                        continue
                    raise
                tester.assertEqual(self.reduce(f * g), 1)

    def _test_inductive_valuation_inheritance(self, **options):
        r"""
        Test that every instance that is a :class:`InductiveValuation` is
        either a :class:`FiniteInductiveValuation` or a
        :class:`InfiniteInductiveValuation`. Same for
        :class:`FinalInductiveValuation` and
        :class:`NonFinalInductiveValuation`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_inductive_valuation_inheritance()
        """
        tester = self._tester(**options)
        tester.assertNotEqual(isinstance(self, InfiniteInductiveValuation),
                              isinstance(self, FiniteInductiveValuation))
        tester.assertNotEqual(isinstance(self, FinalInductiveValuation),
                              isinstance(self, NonFinalInductiveValuation))


class FiniteInductiveValuation(InductiveValuation, DiscreteValuation):
    r"""
    Abstract base class for iterated :mod:`augmented valuations <sage.rings.valuation.augmented_valuation>`
    on top of a :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>` which is a discrete valuation,
    i.e., the last key polynomial has finite valuation.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: from sage.rings.valuation.inductive_valuation import FiniteInductiveValuation
            sage: isinstance(v, FiniteInductiveValuation)
            True
        """
        InductiveValuation.__init__(self, parent, phi)
        DiscreteValuation.__init__(self, parent)

    def extensions(self, other):
        r"""
        Return the extensions of this valuation to ``other``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(ZZ))
            sage: K.<x> = FunctionField(QQ)
            sage: v.extensions(K)
            [Trivial valuation on Rational Field]
        """
        from sage.categories.function_fields import FunctionFields
        if other in FunctionFields() and other.ngens() == 1:
            # extend to K[x] and from there to K(x)
            v = self.extension(self.domain().change_ring(self.domain().base().fraction_field()))
            return [other.valuation(v)]
        return super().extensions(other)


class NonFinalInductiveValuation(FiniteInductiveValuation, DiscreteValuation):
    r"""
    Abstract base class for iterated :mod:`augmented valuations <sage.rings.valuation.augmented_valuation>`
    on top of a :mod:`Gauss valuation <sage.rings.valuation.gauss_valuation>` which can be extended further
    through :meth:`augmentation`.

    EXAMPLES::

        sage: # needs sage.libs.ntl
        sage: R.<u> = Qq(4,5)
        sage: S.<x> = R[]
        sage: v = GaussValuation(S)
        sage: v = v.augmentation(x^2 + x + u, 1)
    """
    def __init__(self, parent, phi):
        r"""
        TESTS::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: from sage.rings.valuation.inductive_valuation import NonFinalInductiveValuation
            sage: isinstance(v, NonFinalInductiveValuation)
            True
        """
        FiniteInductiveValuation.__init__(self, parent, phi)
        DiscreteValuation.__init__(self, parent)

    def augmentation(self, phi, mu, check=True):
        r"""
        Return the inductive valuation which extends this valuation by mapping
        ``phi`` to ``mu``.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation; this must be
          a key polynomial, see :meth:`is_key` for properties of key
          polynomials.

        - ``mu`` -- a rational number or infinity, the valuation of ``phi`` in
          the extended valuation

        - ``check`` -- boolean (default: ``True``); whether or not to check
          the correctness of the parameters

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v
            [ Gauss valuation induced by 2-adic valuation,
              v((1 + O(2^5))*x^2 + (1 + O(2^5))*x + u + O(2^5)) = 1,
              v((1 + O(2^5))*x^4
                + (2^2 + O(2^6))*x^3
                + (1 + (u + 1)*2 + O(2^5))*x^2
                + ((u + 1)*2^2 + O(2^6))*x
                + (u + 1) + (u + 1)*2 + (u + 1)*2^2 + (u + 1)*2^3 + (u + 1)*2^4 + O(2^5)) = 3 ]

        TESTS:

        Make sure that we do not make the assumption that the degrees of the
        key polynomials are strictly increasing::

            sage: v_K = QQ.valuation(3)
            sage: A.<t> = QQ[]
            sage: v0 = GaussValuation(A,v_K)

            sage: v1 = v0.augmentation(t, 1/12)
            sage: v2 = v1.augmentation(t^12 + 3, 7/6)
            sage: v3 = v2.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v4 = v1.augmentation(t^12 + 3*t^2 + 3, 9/4)
            sage: v3 <= v4 and v3 >= v4
            True

        .. SEEALSO::

            :mod:`~sage.rings.valuation.augmented_valuation`
        """
        from .augmented_valuation import AugmentedValuation
        return AugmentedValuation(self, phi, mu, check)

    def mac_lane_step(self, G, principal_part_bound=None, assume_squarefree=False, assume_equivalence_irreducible=False, report_degree_bounds_and_caches=False, coefficients=None, valuations=None, check=True, allow_equivalent_key=True):
        r"""
        Perform an approximation step towards the squarefree monic non-constant
        integral polynomial ``G`` which is not an :meth:`equivalence unit <InductiveValuation.is_equivalence_unit>`.

        This performs the individual steps that are used in
        :meth:`~sage.rings.valuation.valuation.DiscreteValuation.mac_lane_approximants`.

        INPUT:

        - ``G`` -- a squarefree monic non-constant integral polynomial ``G``
          which is not an :meth:`equivalence unit <InductiveValuation.is_equivalence_unit>`

        - ``principal_part_bound`` -- integer or ``None`` (default:
          ``None``), a bound on the length of the principal part, i.e., the
          section of negative slope, of the Newton polygon of ``G``

        - ``assume_squarefree`` -- whether or not to assume that ``G`` is
          squarefree (default: ``False``)

        - ``assume_equivalence_irreducible`` -- whether or not to assume that
          ``G`` is equivalence irreducible (default: ``False``)

        - ``report_degree_bounds_and_caches`` -- whether or not to include internal state with the returned value (used by :meth:`~sage.rings.valuation.valuation.DiscreteValuation.mac_lane_approximants` to speed up sequential calls)

        - ``coefficients`` -- the coefficients of ``G`` in the :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.phi`-adic expansion if known (default: ``None``)

        - ``valuations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``check`` -- whether to check that ``G`` is a squarefree monic
          non-constant  integral polynomial and not an :meth:`equivalence unit <InductiveValuation.is_equivalence_unit>`
          (default: ``True``)

        - ``allow_equivalent_key`` -- whether to return valuations which end in
          essentially the same key polynomial as this valuation but have a
          higher valuation assigned to that key polynomial (default: ``True``)

        EXAMPLES:

        We can use this method to perform the individual steps of
        :meth:`~sage.rings.valuation.valuation.DiscreteValuation.mac_lane_approximants`::

            sage: R.<x> = QQ[]
            sage: v = QQ.valuation(2)
            sage: f = x^36 + 1160/81*x^31 + 9920/27*x^30 + 1040/81*x^26 + 52480/81*x^25 + 220160/81*x^24 - 5120/81*x^21 - 143360/81*x^20 - 573440/81*x^19 + 12451840/81*x^18 - 266240/567*x^16 - 20316160/567*x^15 - 198737920/189*x^14 - 1129840640/81*x^13 - 1907359744/27*x^12 + 8192/81*x^11 + 655360/81*x^10 + 5242880/21*x^9 + 2118123520/567*x^8 + 15460204544/567*x^7 + 6509559808/81*x^6 - 16777216/567*x^2 - 268435456/567*x - 1073741824/567
            sage: v.mac_lane_approximants(f)                                            # needs sage.geometry.polyhedron
            [[ Gauss valuation induced by 2-adic valuation, v(x + 2056) = 23/2 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 11/9 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 2/5, v(x^5 + 4) = 7/2 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^10 + 8*x^5 + 64) = 7 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^5 + 8) = 5 ]]

        Starting from the Gauss valuation, a MacLane step branches off with
        some linear key polynomials in the above example::

            sage: v0 = GaussValuation(R, v)
            sage: V1 = sorted(v0.mac_lane_step(f)); V1                                  # needs sage.geometry.polyhedron
            [[ Gauss valuation induced by 2-adic valuation, v(x) = 2/5 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3/5 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 11/9 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3 ]]

        The computation of MacLane approximants would now perform a MacLane
        step on each of these branches, note however, that a direct call to
        this method might produce some unexpected results::

            sage: V1[1].mac_lane_step(f)                                                # needs sage.geometry.polyhedron
            [[ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^5 + 8) = 5 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^10 + 8*x^5 + 64) = 7 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 11/9 ]]

        Note how this detected the two augmentations of ``V1[1]`` but also two
        other valuations that we had seen in the previous step and that are
        greater than ``V1[1]``. To ignore such trivial augmentations, we can
        set ``allow_equivalent_key``::

            sage: V1[1].mac_lane_step(f, allow_equivalent_key=False)                    # needs sage.geometry.polyhedron
            [[ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^5 + 8) = 5 ],
             [ Gauss valuation induced by 2-adic valuation, v(x) = 3/5, v(x^10 + 8*x^5 + 64) = 7 ]]

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: S.<y> = K[]
            sage: F = y^2 - x^2 - x^3 - 3
            sage: v0 = GaussValuation(K._ring, QQ.valuation(3))
            sage: v1 = v0.augmentation(K._ring.gen(), 1/3)
            sage: mu0 = K.valuation(v1)
            sage: eta0 = GaussValuation(S, mu0)
            sage: eta1 = eta0.mac_lane_step(F)[0]                                       # needs sage.geometry.polyhedron
            sage: eta2 = eta1.mac_lane_step(F)[0]                                       # needs sage.geometry.polyhedron
            sage: eta2                                                                  # needs sage.geometry.polyhedron
            [ Gauss valuation induced by Valuation on rational function field induced by [ Gauss valuation induced by 3-adic valuation, v(x) = 1/3 ], v(y + x) = 2/3 ]

        Check that :issue:`26066` has been resolved::

            sage: R.<x> = QQ[]
            sage: v = QQ.valuation(2)
            sage: v = GaussValuation(R, v).augmentation(x+1, 1/2)
            sage: f = x^4 - 30*x^2 - 75
            sage: v.mac_lane_step(f)                                                    # needs sage.geometry.polyhedron
            [[ Gauss valuation induced by 2-adic valuation, v(x + 1) = 3/4 ]]
        """
        G = self.domain().coerce(G)

        if G.is_constant():
            raise ValueError("G must not be constant")

        from itertools import islice
        from sage.misc.verbose import verbose
        verbose("Augmenting %s towards %s" % (self, G), level=10)

        if not G.is_monic():
            raise ValueError("G must be monic")

        if coefficients is None:
            coefficients = self.coefficients(G)
            if principal_part_bound:
                coefficients = islice(coefficients, 0,
                                      int(principal_part_bound) + 1, 1)
            coefficients = list(coefficients)
        if valuations is None:
            valuations = self.valuations(G, coefficients=coefficients)
            if principal_part_bound:
                valuations = islice(valuations, 0,
                                    int(principal_part_bound) + 1, 1)
            valuations = list(valuations)

        if check and min(valuations) < 0:
            raise ValueError("G must be integral")

        if check and self.is_equivalence_unit(G, valuations=valuations):
            raise ValueError("G must not be an equivalence-unit")

        if check and not assume_squarefree and not G.is_squarefree():
            raise ValueError("G must be squarefree")

        from sage.rings.infinity import infinity
        assert self(G) is not infinity  # this is a valuation and G is nonzero

        ret = []

        F = self.equivalence_decomposition(G, assume_not_equivalence_unit=True, coefficients=coefficients, valuations=valuations, compute_unit=False, degree_bound=principal_part_bound)
        assert len(F), "%s equivalence-decomposes as an equivalence-unit %s" % (G, F)
        if len(F) == 1 and F[0][1] == 1 and F[0][0].degree() == G.degree():
            assert self.is_key(G, assume_equivalence_irreducible=assume_equivalence_irreducible)
            ret.append((self.augmentation(G, infinity, check=False), G.degree(), principal_part_bound, None, None))
        else:
            for phi, e in F:
                if G == phi:
                    # Something strange happened here:
                    # G is not a key (we checked that before) but phi==G is; so phi must have less precision than G
                    # this can happen if not all coefficients of G have the same precision
                    # if we drop some precision of G then it will be a key (but is
                    # that really what we should do?)
                    assert not G.base_ring().is_exact()
                    prec = min([c.precision_absolute() for c in phi.list()])
                    g = G.map_coefficients(lambda c: c.add_bigoh(prec))
                    assert self.is_key(g)
                    ret.append((self.augmentation(g, infinity, check=False), g.degree(), principal_part_bound, None, None))
                    assert len(F) == 1
                    break

                if not allow_equivalent_key and self.phi().degree() == phi.degree():
                    # We ignore augmentations that could have been detected in
                    # the previous MacLane step, see [Rüt2014, Theorem 4.33],
                    # i.e., we ignore key polynomials that are equivalent to
                    # the current key in the sense of that theorem.
                    if self.is_equivalent(self.phi(), phi):
                        continue

                verbose("Determining the augmentation of %s for %s" % (self, phi), level=11)

                base = self
                if phi.degree() == base.phi().degree():
                    # very frequently, the degree of the key polynomials
                    # stagnate for a bit while the valuation of the key
                    # polynomial is slowly increased.
                    # In this case, we can drop previous key polynomials
                    # of the same degree. (They have no influence on the
                    # phi-adic expansion.)
                    if not base.is_gauss_valuation():
                        base = base._base_valuation
                old_mu = self(phi)
                w = base.augmentation(phi, old_mu, check=False)

                # we made some experiments here: instead of computing the
                # coefficients again from scratch, update the coefficients when
                # phi - self.phi() is a constant.
                # It turned out to be slightly slower than just recomputing the
                # coefficients. The main issue with the approach was that we
                # needed to keep track of all the coefficients and not just of
                # the coefficients up to principal_part_bound.

                w_coefficients = w.coefficients(G)
                if principal_part_bound:
                    w_coefficients = islice(w_coefficients, 0,
                                            int(principal_part_bound) + 1, 1)
                w_coefficients = list(w_coefficients)

                w_valuations = w.valuations(G, coefficients=w_coefficients)
                if principal_part_bound:
                    w_valuations = islice(w_valuations, 0,
                                          int(principal_part_bound) + 1, 1)
                w_valuations = list(w_valuations)

                from sage.geometry.newton_polygon import NewtonPolygon
                NP = NewtonPolygon(w.newton_polygon(G, valuations=w_valuations).vertices(), last_slope=0)

                verbose("Newton-Polygon for v(phi)=%s : %s" % (self(phi), NP), level=11)
                slopes = NP.slopes(repetition=True)
                multiplicities = {slope: len([s for s in slopes if s == slope]) for slope in slopes}
                slopes = list(multiplicities)
                if NP.vertices()[0][0] != 0:
                    slopes = [-infinity] + slopes
                    multiplicities[-infinity] = 1

                for i, slope in enumerate(slopes):
                    verbose("Slope = %s" % slope, level=12)
                    new_mu = old_mu - slope
                    new_valuations = [val - (j * slope if slope is not -infinity else (0 if j == 0 else -infinity))
                                      for j, val in enumerate(w_valuations)]
                    if phi.degree() == self.phi().degree():
                        assert new_mu > self(phi), "the valuation of the key polynomial must increase when the degree stagnates"
                    # phi has already been simplified internally by the
                    # equivalence_decomposition method but we can now possibly
                    # simplify it further as we know exactly up to which
                    # precision it needs to be defined.
                    phi = base.simplify(phi, new_mu, force=True)
                    w = base.augmentation(phi, new_mu, check=False)
                    verbose("Augmented %s to %s" % (self, w), level=13)
                    assert slope is -infinity or 0 in w.newton_polygon(G).slopes(repetition=False)

                    from sage.rings.integer_ring import ZZ
                    assert (phi.degree() / self.phi().degree()) in ZZ
                    degree_bound = multiplicities[slope] * phi.degree()
                    assert degree_bound <= G.degree()
                    assert degree_bound >= phi.degree()
                    ret.append((w, degree_bound, multiplicities[slope], w_coefficients, new_valuations))

        if len(ret) == 0:
            assert not allow_equivalent_key, "a MacLane step produced no augmentation"
            assert 0 not in self.newton_polygon(G).slopes(), "a MacLane step produced no augmentation but the valuation given to the key polynomial was correct, i.e., it appears to come out of a call to mac_lane_approximants"

        assert ret, "a MacLane step produced no augmentations"
        if not report_degree_bounds_and_caches:
            ret = [v for v, _, _, _, _ in ret]
        return ret

    def is_key(self, phi, explain=False, assume_equivalence_irreducible=False):
        r"""
        Return whether ``phi`` is a key polynomial for this valuation, i.e.,
        whether it is monic, whether it :meth:`is_equivalence_irreducible`, and
        whether it is :meth:`is_minimal`.

        INPUT:

        - ``phi`` -- a polynomial in the domain of this valuation

        - ``explain`` -- boolean (default: ``False``); if ``True``, return a
          string explaining why ``phi`` is not a key polynomial

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_key(x)
            True
            sage: v.is_key(2*x, explain=True)
            (False, 'phi must be monic')
            sage: v.is_key(x^2, explain=True)
            (False, 'phi must be equivalence irreducible')
            sage: w = v.augmentation(x, 1)
            sage: w.is_key(x + 1, explain = True)
            (False, 'phi must be minimal')
        """
        phi = self.domain().coerce(phi)

        reason = None

        if not phi.is_monic():
            reason = "phi must be monic"
        elif not assume_equivalence_irreducible and not self.is_equivalence_irreducible(phi):
            reason = "phi must be equivalence irreducible"
        elif not self.is_minimal(phi, assume_equivalence_irreducible=True):
            reason = "phi must be minimal"

        if explain:
            return reason is None, reason
        else:
            return reason is None

    def is_minimal(self, f, assume_equivalence_irreducible=False):
        r"""
        Return whether the polynomial ``f`` is minimal with respect to this
        valuation.

        A polynomial `f` is minimal with respect to `v` if it is not a constant
        and any nonzero polynomial `h` which is `v`-divisible by `f` has at
        least the degree of `f`.

        A polynomial `h` is `v`-divisible by `f` if there is a polynomial `c`
        such that `fc` :meth:`~sage.rings.valuation.valuation.DiscretePseudoValuation.is_equivalent` to `h`.

        ALGORITHM:

        Based on Theorem 9.4 of [Mac1936II]_.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4, 5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_minimal(x + 1)
            True
            sage: w = v.augmentation(x, 1)
            sage: w.is_minimal(x + 1)
            False

        TESTS::

            sage: # needs sage.libs.ntl
            sage: K = Qp(2, 10)
            sage: R.<x> = K[]
            sage: vp = K.valuation()
            sage: v0 = GaussValuation(R, vp)
            sage: v1 = v0.augmentation(x, 1/4)
            sage: v2 = v1.augmentation(x^4 + 2, 5/4)
            sage: v2.is_minimal(x^5 + x^4 + 2)
            False

        Polynomials which are equivalent to the key polynomial are minimal if
        and only if they have the same degree as the key polynomial::

            sage: v2.is_minimal(x^4 + 2)                                                # needs sage.libs.ntl
            True
            sage: v2.is_minimal(x^4 + 4)                                                # needs sage.libs.ntl
            False
        """
        f = self.domain().coerce(f)

        if f.is_constant():
            return False

        if not assume_equivalence_irreducible and not self.is_equivalence_irreducible(f):
            # any factor divides f with respect to this valuation
            return False

        if not f.is_monic():
            # divide out the leading factor, it does not change minimality
            v = self
            if not self.domain().base_ring().is_field():
                domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
                v = self.extension(domain)
                f = domain(f)
            return v.is_minimal(f / f.leading_coefficient())

        if self.is_gauss_valuation():
            if self(f) == 0:
                F = self.reduce(f, check=False)
                assert not F.is_constant()
                return F.is_irreducible()
            else:
                assert (self(f) <= 0)  # f is monic
                # f is not minimal:
                # Let g be f stripped of its leading term, i.e., g = f - x^n.
                # Then g and f are equivalent with respect to this valuation
                # and in particular g divides f with respect to this valuation
                return False

        if self.is_equivalent(self.phi(), f):
            assert f.degree() >= self.phi().degree()
            # If an h divides f with respect to this valuation, then it also divides phi:
            # v(f - c*h) > v(f) = v(c*h) => v(phi - c*h) = v((phi - f) + (f - c*h)) > v(phi) = v(c*h)
            # So if f were not minimal then phi would not be minimal but it is.
            return f.degree() == self.phi().degree()

        else:
            tau = self.value_group().index(self._base_valuation.value_group())
            # see Theorem 9.4 of [Mac1936II]
            return list(self.valuations(f))[-1] == self(f) and \
                   list(self.coefficients(f))[-1].is_constant() and \
                   list(self.valuations(f))[0] == self(f) and \
                   tau.divides(len(list(self.coefficients(f))) - 1)

    def _equivalence_reduction(self, f, coefficients=None, valuations=None, degree_bound=None):
        r"""
        Helper method for :meth:`is_equivalence_irreducible` and
        :meth:`equivalence_decomposition` which essentially returns the
        reduction of ``f`` after multiplication with an ``R`` which
        :meth:`is_equivalence_unit`.

        This only works when ``f`` is not divisible by :meth:`phi` with respect
        to this valuation. Therefore, we also return the number of times that
        we took out :meth:`phi` of ``f`` before we computed the reduction.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: v._equivalence_reduction(2*x^6 + 4*x^5 + 2*x^4 + 8)
            (1, 4, x^2 + 1)
        """
        f = self.domain().coerce(f)

        # base change from R[x] to K[x], so divisions work and sufficient
        # elements of negative valuation exist
        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            assert self.residue_ring() is v.residue_ring()
            return v._equivalence_reduction(f)

        if coefficients is None:
            coefficients = list(self.coefficients(f))
        if valuations is None:
            valuations = list(self.valuations(f, coefficients=coefficients))
        valuation = min(valuations)
        for phi_divides in range(len(valuations)):
            # count how many times phi divides f
            if valuations[phi_divides] <= valuation:
                break

        if phi_divides:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(f.parent(), 'phi')
            f = R(coefficients[phi_divides:])(self.phi())
        valuations = [vv - self.mu() * phi_divides
                      for vv in valuations[phi_divides:]]
        coefficients = coefficients[phi_divides:]
        valuation = min(valuations)

        R = self.equivalence_unit(-valuation)
        R = next(self.coefficients(R))
        fR_valuations = [vv - valuation for vv in valuations]
        from sage.rings.infinity import infinity
        fR_coefficients = [next(self.coefficients(c * R))
                           if vv is not infinity and vv == 0 else 0
                           for c, vv in zip(coefficients, fR_valuations)]

        return valuation, phi_divides, self.reduce(f * R, check=False, degree_bound=degree_bound, coefficients=fR_coefficients, valuations=fR_valuations)

    def is_equivalence_irreducible(self, f, coefficients=None, valuations=None):
        r"""
        Return whether the polynomial ``f`` is equivalence-irreducible, i.e.,
        whether its :meth:`equivalence_decomposition` is trivial.

        ALGORITHM:

        We use the same algorithm as in :meth:`equivalence_decomposition` we
        just do not lift the result to key polynomials.

        INPUT:

        - ``f`` -- a non-constant polynomial in the domain of this valuation

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,5)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.is_equivalence_irreducible(x)
            True
            sage: v.is_equivalence_irreducible(x^2)
            False
            sage: v.is_equivalence_irreducible(x^2 + 2)
            False
        """
        f = self.domain().coerce(f)

        if not self.domain().base_ring().is_field():
            domain = self.domain().change_ring(self.domain().base_ring().fraction_field())
            v = self.extension(domain)
            return v.is_equivalence_irreducible(v.domain()(f))

        if f.is_constant():
            raise ValueError("f must not be constant")

        _, phi_divides, F = self._equivalence_reduction(f, coefficients=coefficients, valuations=valuations)
        if phi_divides == 0:
            return F.is_constant() or F.is_irreducible()
        if phi_divides == 1:
            return F.is_constant()
        if phi_divides > 1:
            return False

    def equivalence_decomposition(self, f, assume_not_equivalence_unit=False, coefficients=None, valuations=None, compute_unit=True, degree_bound=None):
        r"""
        Return an equivalence decomposition of ``f``, i.e., a polynomial
        `g(x)=e(x)\prod_i \phi_i(x)` with `e(x)` an :meth:`equivalence unit
        <InductiveValuation.is_equivalence_unit>` and the `\phi_i` :meth:`key
        polynomials <is_key>` such that ``f`` :meth:`~sage.rings.valuation.valuation.DiscretePseudoValuation.is_equivalent` to `g`.

        INPUT:

        - ``f`` -- a nonzero polynomial in the domain of this valuation

        - ``assume_not_equivalence_unit`` -- whether or not to assume that
          ``f`` is not an :meth:`equivalence unit <InductiveValuation.is_equivalence_unit>`
          (default: ``False``)

        - ``coefficients`` -- the coefficients of ``f`` in the
          :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.phi`-adic
          expansion if known (default: ``None``)

        - ``valuations`` -- the valuations of ``coefficients`` if known
          (default: ``None``)

        - ``compute_unit`` -- whether or not to compute the unit part of the
          decomposition (default: ``True``)

        - ``degree_bound`` -- a bound on the degree of the
          :meth:`_equivalence_reduction` of ``f`` (default: ``None``)

        ALGORITHM:

        We use the algorithm described in Theorem 4.4 of [Mac1936II]_. After
        removing all factors `\phi` from a polynomial `f`, there is an
        equivalence unit `R` such that `Rf` has valuation zero. Now `Rf` can be
        factored as `\prod_i \alpha_i` over the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`. Lifting
        all `\alpha_i` to key polynomials `\phi_i` gives `Rf=\prod_i R_i f_i`
        for suitable equivalence units `R_i` (see :meth:`lift_to_key`). Taking
        `R'` an :meth:`~InductiveValuation.equivalence_reciprocal` of `R`, we have `f` equivalent
        to `(R'\prod_i R_i)\prod_i\phi_i`.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.equivalence_decomposition(S.zero())
            Traceback (most recent call last):
            ...
            ValueError: equivalence decomposition of zero is not defined
            sage: v.equivalence_decomposition(S.one())
            1 + O(2^10)
            sage: v.equivalence_decomposition(x^2+2)
            ((1 + O(2^10))*x)^2
            sage: v.equivalence_decomposition(x^2+1)
            ((1 + O(2^10))*x + 1 + O(2^10))^2

        A polynomial that is an equivalence unit, is returned as the unit part
        of a :class:`~sage.structure.factorization.Factorization`, leading to a unit
        non-minimal degree::

            sage: w = v.augmentation(x, 1)                                              # needs sage.libs.ntl
            sage: F = w.equivalence_decomposition(x^2+1); F                             # needs sage.libs.ntl
            (1 + O(2^10))*x^2 + 1 + O(2^10)
            sage: F.unit()                                                              # needs sage.libs.ntl
            (1 + O(2^10))*x^2 + 1 + O(2^10)

        However, if the polynomial has a non-unit factor, then the unit might
        be replaced by a factor of lower degree::

            sage: f = x * (x^2 + 1)                                                     # needs sage.libs.ntl
            sage: F = w.equivalence_decomposition(f); F                                 # needs sage.libs.ntl
            (1 + O(2^10))*x
            sage: F.unit()                                                              # needs sage.libs.ntl
            1 + O(2^10)

        Examples over an iterated unramified extension::

            sage: # needs sage.libs.ntl
            sage: v = v.augmentation(x^2 + x + u, 1)
            sage: v = v.augmentation((x^2 + x + u)^2 + 2*x*(x^2 + x + u) + 4*x, 3)
            sage: v.equivalence_decomposition(x)
            (1 + O(2^10))*x
            sage: F = v.equivalence_decomposition( v.phi() )
            sage: len(F)
            1
            sage: F = v.equivalence_decomposition( v.phi() * (x^4 + 4*x^3 + (7 + 2*u)*x^2 + (8 + 4*u)*x + 1023 + 3*u) )
            sage: len(F)
            2

        TESTS::

            sage: # needs sage.geometry.polyhedron sage.groups sage.rings.number_field
            sage: R.<x> = QQ[]
            sage: K1.<pi> = NumberField(x^3 - 2)
            sage: K.<alpha> = K1.galois_closure()
            sage: R.<x> = K[]
            sage: vp = QQ.valuation(2)
            sage: vp = vp.extension(K)
            sage: v0 = GaussValuation(R, vp)
            sage: G = x^36 + 36*x^35 + 630*x^34 + 7144*x^33 + 59055*x^32 + 379688*x^31 +1978792*x^30 + 8604440*x^29 + 31895428*x^28 + 102487784*x^27 + 289310720*x^26 + 725361352*x^25 + 1629938380*x^24 + 3307417800*x^23 + 6098786184*x^22+10273444280*x^21 + 15878121214*x^20 + 22596599536*x^19 + 29695703772*x^18 +36117601976*x^17 + 40722105266*x^16 + 42608585080*x^15 + 41395961848*x^14 +37344435656*x^13 + 31267160756*x^12 + 24271543640*x^11 + 17439809008*x^10 + 11571651608*x^9 + 7066815164*x^8 + 3953912472*x^7 + 2013737432*x^6 + 925014888*x^5 + 378067657*x^4 + 134716588*x^3 + 40441790*x^2 + 9532544*x + 1584151
            sage: v1 = v0.mac_lane_step(G)[0]
            sage: V = v1.mac_lane_step(G)
            sage: v2 = V[0]
            sage: F = v2.equivalence_decomposition(G); F
            (x^4 + 2*alpha + 1)^3 * (x^4 + 1/2*alpha^4 + alpha + 1)^3 * (x^4 + 1/2*alpha^4 + 3*alpha + 1)^3
            sage: v2.is_equivalent(F.prod(), G)
            True

        Check that :issue:`33422` is fixed::

            sage: R.<x> = QQ[]
            sage: v_7 = QQ.valuation(7)
            sage: v0 = GaussValuation(R, v_7)
            sage: v1 = v0.augmentation(x, 3/2)
            sage: v2 = v1.augmentation(x^2-686, 7/2)
            sage: f = x^4 - 8001504*x^2 - 592815428352
            sage: F = v2.equivalence_decomposition(f); F
            x^4 - 343/2*x^2 + 1294139
            sage: v2.is_equivalent(F.prod(), f)
            True
        """
        f = self.domain().coerce(f)

        if f.is_zero():
            raise ValueError("equivalence decomposition of zero is not defined")

        from sage.structure.factorization import Factorization
        if not assume_not_equivalence_unit and self.is_equivalence_unit(f):
            return Factorization([], unit=f, sort=False)

        if not self.domain().base_ring().is_field():
            nonfractions = self.domain().base_ring()
            domain = self.domain().change_ring(nonfractions.fraction_field())
            v = self.extension(domain)
            ret = v.equivalence_decomposition(v.domain()(f))
            return Factorization([(self._eliminate_denominators(g), e)
                                  for (g, e) in ret], unit=self._eliminate_denominators(ret.unit()), sort=False)

        valuation, phi_divides, F = self._equivalence_reduction(f, coefficients=coefficients, valuations=valuations, degree_bound=degree_bound)
        F = F.factor()
        from sage.misc.verbose import verbose
        verbose("%s factors as %s = %s in reduction" % (f, F.prod(), F), level=20)

        unit = self.domain().one()
        if compute_unit:
            R_ = self.equivalence_unit(valuation, reciprocal=True)
            unit = self.lift(self.residue_ring()(F.unit())) * R_
        F = list(F)

        if compute_unit:
            from sage.misc.misc_c import prod
            unit *= self.lift(self.residue_ring()(prod([psi.leading_coefficient()**e for psi, e in F])))
            if not self.is_gauss_valuation():
                unit *= prod([self._Q_reciprocal(e*psi.degree()) for psi, e in F])
            unit = self.simplify(unit, effective_degree=0, force=True)

        # A potential speedup that we tried to implement here:
        # When F factors as T^n - a, then instead of using any lift of T^n - a
        # we tried to take a lift that approximates well an n-th root of the
        # constant coefficient of f[0]. Doing so saved a few invocations of
        # mac_lane_step but in the end made hardly any difference.

        F = [(self.lift_to_key(psi / psi.leading_coefficient()), e)
             for psi, e in F]

        if phi_divides:
            for i, (g, e) in enumerate(F):
                if g == self.phi():
                    F[i] = (self.phi(), e + phi_divides)
                    break
            else:
                F.append((self.phi(), phi_divides))

        ret = Factorization(F, unit=unit, sort=False)

        if compute_unit:
            assert self.is_equivalent(ret.prod(), f)  # this might fail because of leading zeros in inexact rings
            assert self.is_equivalence_unit(ret.unit())

        return ret

    def minimal_representative(self, f):
        r"""
        Return a minimal representative for ``f``, i.e., a pair `e, a` such
        that ``f`` :meth:`~sage.rings.valuation.valuation.DiscretePseudoValuation.is_equivalent` to `e a`, `e` is an
        :meth:`equivalence unit <InductiveValuation.is_equivalence_unit>`, and `a` :meth:`is_minimal` and monic.

        INPUT:

        - ``f`` -- a nonzero polynomial which is not an equivalence unit

        OUTPUT: a factorization which has `e` as its unit and `a` as its unique factor

        ALGORITHM:

        We use the algorithm described in the proof of Lemma 4.1 of [Mac1936II]_.
        In the expansion `f=\sum_i f_i\phi^i` take `e=f_i` for the largest `i`
        with `f_i\phi^i` minimal (see :meth:`~sage.rings.valuation.developing_valuation.DevelopingValuation.effective_degree`).
        Let `h` be the :meth:`~InductiveValuation.equivalence_reciprocal` of `e` and take `a` given
        by the terms of minimal valuation in the expansion of `e f`.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x

            sage: # needs sage.libs.ntl
            sage: v = v.augmentation(x, 1)
            sage: v.minimal_representative(x + 2)
            (1 + O(2^10))*x + 2 + O(2^11)
            sage: f = x^3 + 6*x + 4
            sage: F = v.minimal_representative(f); F
            (2 + 2^2 + O(2^11)) * ((1 + O(2^10))*x + 2 + O(2^11))
            sage: v.is_minimal(F[0][0])
            True
            sage: v.is_equivalent(F.prod(), f)
            True
        """
        f = self.domain().coerce(f)

        from sage.categories.fields import Fields
        if self.domain().base_ring() not in Fields():
            raise NotImplementedError("only implemented for polynomial rings over fields")

        if f.is_zero():
            raise ValueError("zero has no minimal representative")

        degree = self.effective_degree(f)
        if degree == 0:
            raise ValueError("equivalence units cannot have a minimal representative")

        e = list(self.coefficients(f))[degree]
        h = self.equivalence_reciprocal(e).map_coefficients(_lift_to_maximal_precision)
        g = h * f
        vg = self(g)

        coeffs = [c if v == vg else c.parent().zero()
                  for v, c in zip(self.valuations(g), self.coefficients(g))]
        coeffs[degree] = self.domain().base_ring().one()
        ret = sum([c * self._phi**i for i, c in enumerate(coeffs)])

        assert self.effective_degree(ret) == degree
        assert ret.is_monic()
        assert self.is_minimal(ret)

        from sage.structure.factorization import Factorization
        ret = Factorization([(ret, 1)], unit=e, sort=False)

        assert self.is_equivalent(ret.prod(), f)  # this might fail because of leading zeros
        return ret

    @abstract_method
    def lift_to_key(self, F):
        """
        Lift the irreducible polynomial ``F`` from the
        :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`
        to a key polynomial over this valuation.

        INPUT:

        - ``F`` -- an irreducible non-constant monic polynomial in
          :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`
          of this valuation

        OUTPUT:

        A polynomial `f` in the domain of this valuation which is a key
        polynomial for this valuation and which is such that an
        :meth:`augmentation` with this polynomial adjoins a root of ``F`` to
        the resulting :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`.

        More specifically, if ``F`` is not the generator of the residue ring,
        then multiplying ``f`` with the :meth:`~InductiveValuation.equivalence_reciprocal` of the
        :meth:`~InductiveValuation.equivalence_unit` of the valuation of ``f``, produces a unit
        which reduces to ``F``.

        EXAMPLES::

            sage: # needs sage.libs.ntl
            sage: R.<u> = Qq(4,10)
            sage: S.<x> = R[]
            sage: v = GaussValuation(S)
            sage: y = v.residue_ring().gen()
            sage: u0 = v.residue_ring().base_ring().gen()
            sage: f = v.lift_to_key(y^2 + y + u0); f
            (1 + O(2^10))*x^2 + (1 + O(2^10))*x + u + O(2^10)
        """

    def _eliminate_denominators(self, f):
        r"""
        Return a polynomial in the domain of this valuation that
        :meth:`is_equivalent` to ``f``.

        INPUT:

        - ``f`` -- a polynomial with coefficients in the fraction field of the
          base ring of the domain of this valuation

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, ZZ.valuation(2))
            sage: v._eliminate_denominators(x/3)
            x

        In general such a polynomial may not exist::

            sage: w = v.augmentation(x, 1)
            sage: w._eliminate_denominators(x/2)
            Traceback (most recent call last):
            ...
            ValueError: element has no approximate inverse in this ring

        In general it exists iff the coefficients of minimal valuation in the
        `\phi`-adic expansion of ``f`` do not have denominators of positive
        valuation and if the same is true for these coefficients in their
        expansion; at least if the coefficient ring's residue ring is already a
        field::

            sage: w._eliminate_denominators(x^3/2 + x)
            x
        """
        if f in self.domain():
            return self.domain()(f)

        nonfractions = self.domain().base_ring()
        fractions = nonfractions.fraction_field()

        extended_domain = self.domain().change_ring(fractions)

        g = extended_domain.coerce(f)

        w = self.extension(extended_domain)
        # drop coefficients whose valuation is not minimal (recursively)
        valuation = w(g)
        g = w.simplify(g, error=valuation, force=True, phiadic=True)

        if g in self.domain():
            return self.domain()(g)

        nonfraction_valuation = self.restriction(nonfractions)
        # if this fails then there is no equivalent polynomial in the domain of this valuation
        ret = g.map_coefficients(
                lambda c: c.numerator() * nonfraction_valuation.inverse(c.denominator(),
                        valuation
                        + nonfraction_valuation(c.denominator())
                        - nonfraction_valuation(c.numerator())
                        + nonfraction_valuation.value_group().gen()),
                nonfractions)
        assert w.is_equivalent(f, ret)
        return ret

    def _test_eliminate_denominators(self, **options):
        r"""
        Test the correctness of :meth:`_eliminate_denominators`.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: v = GaussValuation(R, ZZ.valuation(2))
            sage: v._test_eliminate_denominators()
        """
        tester = self._tester(**options)

        nonfractions = self.domain().base_ring()
        fractions = nonfractions.fraction_field()
        extended_domain = self.domain().change_ring(fractions)
        w = self.extension(extended_domain)

        S = tester.some_elements(w.domain().some_elements())
        for f in S:
            try:
                g = self._eliminate_denominators(f)
            except ValueError:
                continue
            tester.assertIs(g.parent(), self.domain())
            tester.assertTrue(w.is_equivalent(f, g))

    def _test_lift_to_key(self, **options):
        r"""
        Test the correctness of :meth:`lift_to_key`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_lift_to_key()                                                 # needs sage.rings.number_field
        """
        tester = self._tester(**options)

        try:
            self.residue_ring()
        except NotImplementedError:
            from sage.categories.fields import Fields
            if self.domain().base() in Fields():
                raise
            return

        S = tester.some_elements(self.residue_ring().some_elements())
        for F in S:
            if F.is_monic() and not F.is_constant() and F.is_irreducible():
                try:
                    f = self.lift_to_key(F)
                except NotImplementedError:
                    from sage.categories.fields import Fields
                    if self.domain().base() in Fields():
                        raise
                    continue
                tester.assertIs(f.parent(), self.domain())
                tester.assertTrue(self.is_key(f))

                # check that augmentation produces a valuation with roots of F
                # in the residue ring
                from sage.rings.infinity import infinity
                w = self.augmentation(f, infinity)
                F = F.change_ring(w.residue_ring())
                roots = F.roots(multiplicities=False)
                tester.assertGreaterEqual(len(roots), 1)

                # check that f has the right reduction
                if F == F.parent().gen():
                    tester.assertTrue(self.is_equivalent(f, self.phi()))
                else:
                    tester.assertEqual(self.reduce(f * self.equivalence_reciprocal(self.equivalence_unit(self(f)))), F)

    def _test_is_equivalence_irreducible(self, **options):
        r"""
        Test the correctness of :meth:`is_equivalence_irreducible`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
            sage: v._test_is_equivalence_irreducible()
        """
        tester = self._tester(**options)
        S = tester.some_elements(self.domain().some_elements())
        for f in S:
            if f.is_constant():
                continue
            is_equivalence_irreducible = self.is_equivalence_irreducible(f)
            F = self.equivalence_decomposition(f)
            tester.assertEqual(is_equivalence_irreducible, len(F) == 0 or (len(F) == 1 and F[0][1] == 1))
            if self.is_equivalence_unit(f):
                tester.assertTrue(f.is_constant() or self.is_equivalence_irreducible(f))

        tester.assertTrue(self.is_equivalence_irreducible(self.phi()))
        tester.assertTrue(self.is_equivalence_irreducible(-self.phi()))
        tester.assertFalse(self.is_equivalence_irreducible(self.phi() ** 2))


class FinalInductiveValuation(InductiveValuation):
    r"""
    Abstract base class for an inductive valuation which cannot be augmented further.

    TESTS::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, valuations.TrivialValuation(QQ))
        sage: w = v.augmentation(x^2 + x + 1, infinity)
        sage: from sage.rings.valuation.inductive_valuation import FinalInductiveValuation
        sage: isinstance(w, FinalInductiveValuation)
        True
    """


class InfiniteInductiveValuation(FinalInductiveValuation, InfiniteDiscretePseudoValuation):
    r"""
    Abstract base class for an inductive valuation which is not discrete, i.e.,
    which assigns infinite valuation to its last key polynomial.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: v = GaussValuation(R, QQ.valuation(2))
        sage: w = v.augmentation(x^2 + x + 1, infinity)
    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: from sage.rings.valuation.inductive_valuation import InfiniteInductiveValuation
            sage: isinstance(w, InfiniteInductiveValuation)
            True
        """
        FinalInductiveValuation.__init__(self, parent, base_valuation)
        InfiniteDiscretePseudoValuation.__init__(self, parent)

    def change_domain(self, ring):
        r"""
        Return this valuation over ``ring``.

        EXAMPLES:

        We can turn an infinite valuation into a valuation on the quotient::

            sage: R.<x> = QQ[]
            sage: v = GaussValuation(R, QQ.valuation(2))
            sage: w = v.augmentation(x^2 + x + 1, infinity)
            sage: w.change_domain(R.quo(x^2 + x + 1))
            2-adic valuation
        """
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        if isinstance(ring, PolynomialQuotientRing_generic) and ring.base() is self.domain() and ring.modulus() == self.phi():
            return self.restriction(self.domain().base())._extensions_to_quotient(ring, approximants=[self])[0]
        return super().change_domain(ring)


def _lift_to_maximal_precision(c):
    r"""
    Lift ``c`` to maximal precision if the parent is not exact.

    EXAMPLES::

        sage: R = Zp(2,5)
        sage: x = R(1,2); x
        1 + O(2^2)
        sage: from sage.rings.valuation.inductive_valuation import _lift_to_maximal_precision
        sage: _lift_to_maximal_precision(x)
        1 + O(2^5)

        sage: x = 1
        sage: _lift_to_maximal_precision(x)
        1
    """
    return c if c.parent().is_exact() else c.lift_to_precision()
