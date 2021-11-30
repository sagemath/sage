r"""
Precision management for non-exact objects

Manage the default precision for non-exact objects such as power series rings or
laurent series rings.

EXAMPLES::

    sage: R.<x> = PowerSeriesRing(QQ)
    sage: R.default_prec()
    20
    sage: cos(x)
    1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6 + 1/40320*x^8 - 1/3628800*x^10 +
     1/479001600*x^12 - 1/87178291200*x^14 + 1/20922789888000*x^16 -
     1/6402373705728000*x^18 + O(x^20)

::

    sage: R.<x> = PowerSeriesRing(QQ, default_prec=10)
    sage: R.default_prec()
    10
    sage: cos(x)
    1 - 1/2*x^2 + 1/24*x^4 - 1/720*x^6 + 1/40320*x^8 + O(x^10)

.. NOTE::

    Subclasses of :class:`Nonexact` which require to change the default
    precision should implement a method `set_default_prec`.

"""

from sage.rings.integer import Integer

class Nonexact:
    r"""
    A non-exact object with default precision.

    INPUT:

    - ``prec`` -- a non-negative integer representing the default precision of
      ``self`` (default: ``20``)

    """
    def __init__(self, prec=20):
        if prec < 0:
            raise ValueError(f"prec (= {prec}) must be non-negative")
        self._default_prec = Integer(prec)

    def default_prec(self):
        r"""
        Return the default precision for ``self``.

        EXAMPLES::

            sage: R = QQ[[x]]
            sage: R.default_prec()
            20

        ::

            sage: R.<x> = PowerSeriesRing(QQ, default_prec=10)
            sage: R.default_prec()
            10

        """
        try:
            return self._default_prec
        except AttributeError:
            self._default_prec = 20
            return self._default_prec

    def set_default_prec(self, prec):
        r"""
        Set the default precision for ``self``

        .. WARNING::

            This method is outdated. If a subclass of class:`Nonexact` requires
            this method, please overload it instead.

        """
        # TODO: remove in Sage 9.4
        from sage.misc.superseded import deprecation
        msg = "The method set_default_prec() is deprecated and will be removed "
        msg += "in a future version of Sage. The default precision is set "
        msg += "during construction."
        deprecation(18416, msg)
        self._default_prec = Integer(prec)


