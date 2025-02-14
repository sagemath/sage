"""
Valuation for power series and Laurent series.
"""

# ****************************************************************************
#       Copyright (C) 2025 Xavier Caruso <xavier@caruso.ovh>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.fields import Fields
from sage.structure.factory import UniqueFactory
from sage.rings.valuation.valuation import DiscreteValuation

class SeriesValuationFactory(UniqueFactory):
    r"""
    Create a valuation over a power series ring or a Laurent series ring.

    INPUT:

    - ``R`` -- the ring

    EXAMPLES:

        sage: R.<t> = PowerSeriesRing(QQ)
        sage: v = R.valuation()
        sage: v
        t-adic valuation
    """
    def create_key_and_extra_args(self, R):
        r"""
        TESTS::

            sage: R.<t> = LaurentSeriesRing(ZZ)
            sage: R.valuation()
            Traceback (most recent call last):
            ...
            NotImplementedError: valuation only implemented when the base ring is a field
        """
        if R.base_ring() not in Fields():
            raise NotImplementedError("valuation only implemented when the base ring is a field")
        return R, {}

    def create_object(self, version, key, **extra_args):
        r"""
        Create the valuation identified by the given key.

        EXAMPLES::

            sage: R.<t> = PowerSeriesRing(GF(5))
            sage: R.valuation()  # indirect doctest
            t-adic valuation
        """
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(key)
        return parent.__make_element_class__(SeriesValuation_generic)(parent)

SeriesValuation = SeriesValuationFactory("sage.rings.series_valuation.SeriesValuation")


class SeriesValuation_generic(DiscreteValuation):
    r"""
    A class for valuation on power series and Laurent Series

    TESTS::

        sage: K.<t> = LaurentSeriesRing(QQ)
        sage: v = K.valuation()
        sage: TestSuite(v).run()
    """
    def _repr_(self):
        r"""
        Return a string representation of this valuation.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: v  # indirect doctest
            t-adic valuation
        """
        return "%s-adic valuation" % self.uniformizer()

    def _call_(self, x):
        r"""
        Return the valuation of ``x``.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: v(t)  # indirect doctest
            1
        """
        return x.valuation()

    def uniformizer(self):
        r"""
        Return the uniformizer of the underlying ring.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: v.uniformizer()
            t
        """
        return self.domain().uniformizer()

    def residue_ring(self):
        r"""
        Return the residue field of the underlying ring.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: v.residue_field()
            Rational Field
        """
        return self.domain().residue_field()

    residue_field = residue_ring

    def reduce(self, x):
        r"""
        Return the image of ``x`` in the residue field.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: v.reduce(t)
            0
            sage: v.reduce(t + 1)
            1

        When ``x`` has negative valuation, an error is raised::

            sage: v.reduce(1/t)
            Traceback (most recent call last):
            ...
            ValueError: reduction is only defined for elements of nonnegative valuation
        """
        x = self.domain().coerce(x)
        if x.valuation() < 0:
            raise ValueError("reduction is only defined for elements of nonnegative valuation")
        return x[0]

    def lift(self, x):
        r"""
        Return a lift of ``x`` in the underlying ring.

        EXAMPLES::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: x = v.lift(1)
            sage: x
            1
            sage: x.parent()
            Laurent Series Ring in t over Rational Field
        """
        x = self.residue_field().coerce(x)
        return self.domain()(x)
