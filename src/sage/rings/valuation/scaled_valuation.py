# sage_setup: distribution = sagemath-pari
r"""
Valuations which are scaled versions of another valuation

EXAMPLES::

    sage: 3*ZZ.valuation(3)
    3 * 3-adic valuation

AUTHORS:

- Julian Rüth (2016-11-10): initial version
"""
# ****************************************************************************
#       Copyright (C) 2016-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.factory import UniqueFactory

from .valuation import DiscreteValuation


class ScaledValuationFactory(UniqueFactory):
    r"""
    Return a valuation which scales the valuation ``base`` by the factor ``s``.

    EXAMPLES::

        sage: 3*ZZ.valuation(2) # indirect doctest
        3 * 2-adic valuation
    """
    def create_key(self, base, s):
        r"""
        Create a key which uniquely identifies a valuation.

        TESTS::

            sage: 3*ZZ.valuation(2) is 2*(3/2*ZZ.valuation(2)) # indirect doctest
            True
        """
        from sage.rings.infinity import infinity
        from sage.rings.rational_field import QQ
        if s is infinity or s not in QQ or s <= 0:
            # for these values we can not return a TrivialValuation() in
            # create_object() because that would override that instance's
            # _factory_data and lead to pickling errors
            raise ValueError("s must be a positive rational")
        if base.is_trivial():
            # for the same reason we can not accept trivial valuations here
            raise ValueError("base must not be trivial")
        s = QQ.coerce(s)
        if s == 1:
            # we would override the _factory_data of base if we just returned
            # it in create_object() so we just refuse to do so
            raise ValueError("s must not be 1")

        if isinstance(base, ScaledValuation_generic):
            return self.create_key(base._base_valuation, s * base._scale)

        return base, s

    def create_object(self, version, key):
        r"""
        Create a valuation from ``key``.

        TESTS::

            sage: 3*ZZ.valuation(2) # indirect doctest
            3 * 2-adic valuation
        """
        base, s = key

        assert not isinstance(base, ScaledValuation_generic)

        from .valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(base.domain())
        return parent.__make_element_class__(ScaledValuation_generic)(parent, base, s)


ScaledValuation = ScaledValuationFactory("sage.rings.valuation.scaled_valuation.ScaledValuation")


class ScaledValuation_generic(DiscreteValuation):
    r"""
    A valuation which scales another ``base_valuation`` by a finite positive factor ``s``.

    EXAMPLES::

        sage: v = 3*ZZ.valuation(3); v
        3 * 3-adic valuation

    TESTS::

        sage: TestSuite(v).run()                # long time                             # needs sage.geometry.polyhedron
    """
    def __init__(self, parent, base_valuation, s):
        r"""
        .. TODO::

            It is annoying that we have to wrap any possible method on
            ``base_valuation`` in this class. It would be nice if this would
            somehow be done automagically, e.g., by adding annotations to the
            methods in ``base_valuation`` that explain which parameters and
            return values need to be scaled.

        TESTS::

            sage: v = 3*ZZ.valuation(2)
            sage: from sage.rings.valuation.scaled_valuation import ScaledValuation_generic
            sage: isinstance(v, ScaledValuation_generic)
            True
        """
        DiscreteValuation.__init__(self, parent)
        self._base_valuation = base_valuation
        self._scale = s

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: 3*ZZ.valuation(2) # indirect doctest
            3 * 2-adic valuation
        """
        return "%r * %r" % (self._scale, self._base_valuation)

    def residue_ring(self):
        r"""
        Return the residue field of this valuation.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(2)
            sage: v.residue_ring()
            Finite Field of size 2
        """
        return self._base_valuation.residue_ring()

    def uniformizer(self):
        r"""
        Return a uniformizing element of this valuation.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(2)
            sage: v.uniformizer()
            2
        """
        return self._base_valuation.uniformizer()

    def _call_(self, f):
        r"""
        Evaluate this valuation at ``f``.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(2)
            sage: v(2)
            3
        """
        return self._scale * self._base_valuation(f)

    def reduce(self, f):
        r"""
        Return the reduction of ``f`` in the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field` of this valuation.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(2)
            sage: v.reduce(1)
            1
        """
        return self._base_valuation.reduce(f)

    def lift(self, F):
        r"""
        Lift ``F`` from the :meth:`~sage.rings.valuation.valuation_space.DiscretePseudoValuationSpace.ElementMethods.residue_field`
        of this valuation into its
        domain.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(2)
            sage: v.lift(1)
            1
        """
        return self._base_valuation.lift(F)

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: v = 3*ZZ.valuation(5)
            sage: v.extensions(GaussianIntegers().fraction_field())                     # needs sage.rings.number_field
            [3 * [ 5-adic valuation, v(x + 2) = 1 ]-adic valuation,
             3 * [ 5-adic valuation, v(x + 3) = 1 ]-adic valuation]
        """
        return [ScaledValuation(w, self._scale) for w in self._base_valuation.extensions(ring)]

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: v = 3*QQ.valuation(5)
            sage: v.restriction(ZZ)
            3 * 5-adic valuation
        """
        restricted = self._base_valuation.restriction(ring)
        if restricted.is_trivial():
            # trivial valuations can not be scaled explicitly for caching reasons
            return restricted
        return ScaledValuation(restricted, self._scale)

    def _strictly_separating_element(self, other):
        r"""
        Return an element in the domain of this valuation which has positive
        valuation with respect to this valuation but negative valuation with
        respect to ``other``.

        EXAMPLES::

            sage: v2 = QQ.valuation(2)
            sage: v3 = 12 * QQ.valuation(3)
            sage: v2._strictly_separating_element(v3)
            2/3
        """
        return self._base_valuation._strictly_separating_element(other)

    def _weakly_separating_element(self, other):
        r"""
        Return an element in the domain of this valuation which has
        positive valuation with respect to this valuation and higher
        valuation with respect to this valuation than with respect to
        ``other``.

        EXAMPLES::

            sage: v2 = QQ.valuation(2)
            sage: v3 = 12 * QQ.valuation(3)
            sage: v2._weakly_separating_element(v3)
            2
        """
        return self._base_valuation._weakly_separating_element(other)

    def _ge_(self, other):
        r"""
        Return whether this valuation is greater or equal to ``other``, a
        valuation on the same domain.

        EXAMPLES::

            sage: v2 = QQ.valuation(2)
            sage: 2*v2 >= v2
            True
            sage: v2/2 >= 2*v2
            False
            sage: 3*v2 > 2*v2
            True

        Test that non-scaled valuations call through to this method to resolve
        the scaling::

            sage: v2 > v2/2
            True
        """
        if self == other:
            return True
        if isinstance(other, ScaledValuation_generic):
            return (self._scale / other._scale) * self._base_valuation >= other._base_valuation
        if self._scale >= 1:
            if self._base_valuation >= other:
                return True
        else:
            assert not self.is_trivial()
            if self._base_valuation <= other:
                return False
        return super()._ge_(other)

    def _le_(self, other):
        r"""
        Return whether this valuation is smaller or equal to ``other``, a
        valuation on the same domain.

        EXAMPLES::

            sage: v2 = QQ.valuation(2)
            sage: 2*v2 <= v2
            False
            sage: v2/2 <= 2*v2
            True
            sage: 3*v2 < 2*v2
            False

        Test that non-scaled valuations call through to this method to resolve
        the scaling::

            sage: v2 < v2/2
            False
        """
        return other / self._scale >= self._base_valuation

    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        EXAMPLES::

            sage: v2 = QQ.valuation(2)
            sage: (2*v2).value_semigroup()
            Additive Abelian Semigroup generated by -2, 2
        """
        return self._scale * self._base_valuation.value_semigroup()
