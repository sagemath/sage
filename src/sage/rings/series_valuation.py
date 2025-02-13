"""
Valuation for power series and Laurent series.
"""

from sage.rings.valuation.valuation import DiscreteValuation

class SeriesValuation(DiscreteValuation):
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
        return "%s-adic valuation" % self.domain().uniformizer()

    def __reduce__(self):
        r"""
        TESTS::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: loads(dumps(v)) == v
            True
        """
        return (self.__class__, (self.parent(),))

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

    def _eq_(self, other):
        r"""
        Return whether ``self`` and ``other`` define the same
        valuation.

        TESTS::

            sage: K.<t> = LaurentSeriesRing(QQ)
            sage: v = K.valuation()
            sage: w = K.valuation()
            sage: v == w
            True

        ::

            sage: R = K.power_series_ring()
            sage: w = R.valuation()
            sage: v == w
            False
        """
        t = self.domain().uniformizer()
        return other(t) == 1
