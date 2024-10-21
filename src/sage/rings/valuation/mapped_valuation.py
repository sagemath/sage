# sage_setup: distribution = sagemath-pari
r"""
Valuations which are implemented through a map to another valuation

EXAMPLES:

Extensions of valuations over finite field extensions `L=K[x]/(G)` are realized
through an infinite valuation on `K[x]` which maps `G` to infinity::

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(QQ)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x)

    sage: v = K.valuation(0)                                                            # needs sage.rings.function_field
    sage: w = v.extension(L); w                                                         # needs sage.rings.function_field
    (x)-adic valuation

    sage: w._base_valuation                                                             # needs sage.rings.function_field
    [ Gauss valuation induced by (x)-adic valuation, v(y) = 1/2 , … ]

AUTHORS:

- Julian Rüth (2016-11-10): initial version
"""
# ****************************************************************************
#       Copyright (C) 2016-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .valuation import DiscreteValuation, DiscretePseudoValuation
from sage.misc.abstract_method import abstract_method


class MappedValuation_base(DiscretePseudoValuation):
    r"""
    A valuation which is implemented through another proxy "base" valuation.

    EXAMPLES::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = K.valuation(0)
        sage: w = v.extension(L); w
        (x)-adic valuation

    TESTS::

        sage: TestSuite(w).run()                # long time                             # needs sage.rings.function_field
    """
    def __init__(self, parent, base_valuation):
        r"""
        .. TODO::

            It is annoying that we have to wrap any possible method on
            ``base_valuation`` in this class. It would be nice if this would
            somehow be done automagically, e.g., by adding annotations to the
            methods in ``base_valuation`` that explain which parameters and
            return values need to be mapped and how.

        TESTS::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^2 + 1)
            sage: v = K.valuation(0)
            sage: w = v.extension(L); w
            (x)-adic valuation
            sage: from sage.rings.valuation.mapped_valuation import MappedValuation_base
            sage: isinstance(w, MappedValuation_base)
            True
        """
        DiscretePseudoValuation.__init__(self, parent)

        self._base_valuation = base_valuation

    @abstract_method
    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        Subclasses must override this method.

        EXAMPLES::

            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)                                          # needs sage.rings.number_field
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: v.extension(L)  # indirect doctest                                    # needs sage.rings.number_field
            2-adic valuation
        """

    def residue_ring(self):
        r"""
        Return the residue ring of this valuation.

        EXAMPLES::

            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)                                          # needs sage.rings.number_field
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: v.extension(L).residue_ring()                                         # needs sage.rings.number_field
            Finite Field of size 2
        """
        return self._base_valuation.residue_ring()

    def uniformizer(self):
        r"""
        Return a uniformizing element of this valuation.

        EXAMPLES::

            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)                                          # needs sage.rings.number_field
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: v.extension(L).uniformizer()                                          # needs sage.rings.number_field
            t + 1
        """
        return self._from_base_domain(self._base_valuation.uniformizer())

    def _to_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of ``_base_valuation``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]
            sage: w._to_base_domain(y).parent()
            Univariate Polynomial Ring in y over Rational function field in x over Rational Field
        """
        return self._base_valuation.domain().coerce(f)

    def _from_base_domain(self, f):
        r"""
        Return ``f`` as an element in the domain of this valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w._from_base_domain(w._base_valuation.domain().gen()).parent()
            Function field in y defined by y^2 - x
        """
        return self.domain().coerce(f)

    def _call_(self, f):
        r"""
        Evaluate this valuation at ``f``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w(y)  # indirect doctest
            1/2
        """
        return self._base_valuation(self._to_base_domain(f))

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` in the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field` of this valuation.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - (x - 2))
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: w.reduce(y)
            u1
        """
        return self._from_base_residue_ring(self._base_valuation.reduce(self._to_base_domain(f)))

    def lift(self, F):
        r"""
        Lift ``F`` from the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`
        of this valuation into its domain.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(2)
            sage: w = v.extension(L)
            sage: w.lift(w.residue_field().gen())
            y
        """
        F = self.residue_ring().coerce(F)
        F = self._to_base_residue_ring(F)
        f = self._base_valuation.lift(F)
        return self._from_base_domain(f)

    def simplify(self, x, error=None, force=False):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        If ``force`` is not set, then expensive simplifications may be avoided.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]

        As :meth:`_relative_size` misses the bloated term ``x^32``, the
        following term does not get simplified::

            sage: w.simplify(y + x^32)                                                  # needs sage.rings.function_field
            y + x^32

        In this case the simplification can be forced but this should not
        happen as a default as the recursive simplification can be quite
        costly::

            sage: w.simplify(y + x^32, force=True)                                      # needs sage.rings.function_field
            y
        """
        return self._from_base_domain(self._base_valuation.simplify(self._to_base_domain(x), error=error, force=force))

    def _relative_size(self, x):
        r"""
        Return an estimate on the coefficient size of ``x``.

        The number returned is an estimate on the factor between the number of
        bits used by ``x`` and the minimal number of bits used by an element
        congruent to ``x``.

        This can be used by :meth:`simplify` to decide whether simplification
        of coefficients is going to lead to a significant shrinking of the
        coefficients of ``x``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]

        In this example, the method misses the size of the bloated term
        ``x^32``::

            sage: w._relative_size(y + x^32)                                            # needs sage.rings.function_field
            1
        """
        return self._base_valuation._relative_size(self._to_base_domain(x))

    def _to_base_residue_ring(self, F):
        r"""
        Return ``F``, an element of :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`,
        as an element of the residue ring of the ``_base_valuation``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]
            sage: w._to_base_residue_ring(1)
            1
        """
        return self._base_valuation.residue_ring().coerce(F)

    def _from_base_residue_ring(self, F):
        r"""
        Return ``F``, an element of the residue ring of ``_base_valuation``, as
        an element of this valuation's :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_ring`.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]
            sage: w._from_base_residue_ring(1)
            1
        """
        return self.residue_ring().coerce(F)

    def element_with_valuation(self, s):
        r"""
        Return an element with valuation ``s``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 5)
            sage: u,uu = v.extensions(L)
            sage: u.element_with_valuation(1)
            5
        """
        return self._from_base_domain(self._base_valuation.element_with_valuation(s))

    def _test_to_from_base_domain(self, **options):
        r"""
        Check the correctness of :meth:`to_base_domain` and
        :meth:`from_base_domain`.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]
            sage: w._test_to_from_base_domain()
        """
        tester = self._tester(**options)

        for x in tester.some_elements(self.domain().some_elements()):
            tester.assertEqual(x, self._from_base_domain(self._to_base_domain(x)))
            # note that the converse might not be true

    def _test_to_from_base_residue_ring(self, **options):
        r"""
        Check the correctness of :meth:`to_base_residue_ring` and
        :meth:`from_base_residue_ring`.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extensions(L)[0]
            sage: w._test_to_from_base_residue_ring()
        """
        tester = self._tester(**options)

        for x in tester.some_elements(self.residue_ring().some_elements()):
            tester.assertEqual(x, self._from_base_residue_ring(self._to_base_residue_ring(x)))
        for x in tester.some_elements(self._base_valuation.residue_ring().some_elements()):
            tester.assertEqual(x, self._to_base_residue_ring(self._from_base_residue_ring(x)))


class FiniteExtensionFromInfiniteValuation(MappedValuation_base, DiscreteValuation):
    r"""
    A valuation on a quotient of the form `L=K[x]/(G)` with an irreducible `G`
    which is internally backed by a pseudo-valuations on `K[x]` which sends `G`
    to infinity.

    INPUT:

    - ``parent`` -- the containing valuation space (usually the space of
      discrete valuations on `L`)

    - ``base_valuation`` -- an infinite valuation on `K[x]` which takes `G` to
      infinity

    EXAMPLES::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = K.valuation(0)
        sage: w = v.extension(L); w
        (x)-adic valuation
    """
    def __init__(self, parent, base_valuation):
        r"""
        TESTS::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(0)
            sage: w = v.extension(L)
            sage: from sage.rings.valuation.mapped_valuation import FiniteExtensionFromInfiniteValuation
            sage: isinstance(w, FiniteExtensionFromInfiniteValuation)
            True
            sage: TestSuite(w).run()            # long time
        """
        MappedValuation_base.__init__(self, parent, base_valuation)
        DiscreteValuation.__init__(self, parent)

    def _eq_(self, other):
        r"""
        Return whether this valuation is indistinguishable from ``other``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: ww = v.extension(L)
            sage: w == ww # indirect doctest
            True
        """
        return (isinstance(other, FiniteExtensionFromInfiniteValuation)
                and self._base_valuation == other._base_valuation)

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: w.restriction(K) is v
            True
        """
        if ring.is_subring(self._base_valuation.domain().base()):
            return self._base_valuation.restriction(ring)
        return super().restriction(ring)

    def _weakly_separating_element(self, other):
        r"""
        Return an element in the domain of this valuation which has
        positive valuation with respect to this valuation and higher
        valuation with respect to this valuation than with respect to
        ``other``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: v = valuations.pAdicValuation(QQ, 5)
            sage: u,uu = v.extensions(L)
            sage: u.separating_element([w,uu]) # indirect doctest
            1/20*t + 7/20
        """
        if isinstance(other, FiniteExtensionFromInfiniteValuation):
            return self.domain()(self._base_valuation._weakly_separating_element(other._base_valuation))
        super()._weakly_separating_element(other)

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

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 2)
            sage: w = v.extension(L)
            sage: w._relative_size(1024*t + 1024)
            6
        """
        return self._base_valuation._relative_size(self._to_base_domain(x))

    def simplify(self, x, error=None, force=False):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 5)
            sage: u,uu = v.extensions(L)
            sage: f = 125*t + 1
            sage: u.simplify(f, error=u(f), force=True)
            1
        """
        x = self.domain().coerce(x)

        if error is None:
            error = self.upper_bound(x)

        return self._from_base_domain(self._base_valuation.simplify(self._to_base_domain(x), error, force=force))

    def lower_bound(self, x):
        r"""
        Return a lower bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 5)
            sage: u,uu = v.extensions(L)
            sage: u.lower_bound(t + 2)
            0
            sage: u(t + 2)
            1
        """
        x = self.domain().coerce(x)
        return self._base_valuation.lower_bound(self._to_base_domain(x))

    def upper_bound(self, x):
        r"""
        Return an upper bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K = QQ
            sage: R.<t> = K[]
            sage: L.<t> = K.extension(t^2 + 1)
            sage: v = valuations.pAdicValuation(QQ, 5)
            sage: u,uu = v.extensions(L)
            sage: u.upper_bound(t + 2) >= 1
            True
            sage: u(t + 2)
            1
        """
        x = self.domain().coerce(x)
        return self._base_valuation.upper_bound(self._to_base_domain(x))


class FiniteExtensionFromLimitValuation(FiniteExtensionFromInfiniteValuation):
    r"""
    An extension of a valuation on a finite field extensions `L=K[x]/(G)` which
    is induced by an infinite limit valuation on `K[x]`.

    EXAMPLES::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: v = K.valuation(1)
        sage: w = v.extensions(L); w
        [[ (x - 1)-adic valuation, v(y + 1) = 1 ]-adic valuation,
         [ (x - 1)-adic valuation, v(y - 1) = 1 ]-adic valuation]

    TESTS::

        sage: TestSuite(w[0]).run()             # long time                             # needs sage.rings.function_field
        sage: TestSuite(w[1]).run()             # long time                             # needs sage.rings.function_field
    """
    def __init__(self, parent, approximant, G, approximants):
        r"""
        EXAMPLES:

        Note that this implementation is also used when the underlying limit is
        only taken over a finite sequence of valuations::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: v = K.valuation(2)
            sage: w = v.extension(L); w
            (x - 2)-adic valuation
            sage: from sage.rings.valuation.mapped_valuation import FiniteExtensionFromLimitValuation
            sage: isinstance(w, FiniteExtensionFromLimitValuation)
            True
        """
        # keep track of all extensions to this field extension so we can print
        # this valuation nicely, dropping any unnecessary information
        self._approximants = approximants

        from .limit_valuation import LimitValuation
        limit = LimitValuation(approximant, G)
        FiniteExtensionFromInfiniteValuation.__init__(self, parent, limit)

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: valuations.pAdicValuation(GaussianIntegers().fraction_field(), 2)  # indirect doctest                 # needs sage.rings.number_field
            2-adic valuation
        """
        from .limit_valuation import MacLaneLimitValuation
        if isinstance(self._base_valuation, MacLaneLimitValuation):
            # print the minimal information that singles out this valuation from all approximants
            assert (self._base_valuation._initial_approximation in self._approximants)
            approximants = [v.augmentation_chain()[::-1] for v in self._approximants]
            augmentations = self._base_valuation._initial_approximation.augmentation_chain()[::-1]
            unique_approximant = None
            for l in range(len(augmentations)):
                if len([a for a in approximants if a[:l + 1] == augmentations[:l + 1]]) == 1:
                    unique_approximant = augmentations[:l + 1]
                    break
            assert (unique_approximant is not None)
            if unique_approximant[0].is_gauss_valuation():
                unique_approximant[0] = unique_approximant[0].restriction(unique_approximant[0].domain().base_ring())
            if len(unique_approximant) == 1:
                return repr(unique_approximant[0])
            from .augmented_valuation import AugmentedValuation_base
            return "[ %s ]-adic valuation" % (", ".join("v(%r) = %r" % (v._phi, v._mu) if (isinstance(v, AugmentedValuation_base) and v.domain() == self._base_valuation.domain()) else repr(v) for v in unique_approximant))
        return "%s-adic valuation" % (self._base_valuation)
