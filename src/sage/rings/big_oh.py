"""
Big O for various types (power series, `p`-adics, etc.)

.. SEEALSO::

    - `asymptotic expansions <../../../asymptotic/index.html>`_
    - `p-adic numbers <../../../padics/index.html>`_
    - `power series <../../../power_series/index.html>`_
    - `polynomials <../../../polynomial_rings/index.html>`_
"""

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.padics.factory', ['Qp', 'Zp'])
from sage.rings.polynomial.polynomial_element import Polynomial

try:
    from .puiseux_series_ring_element import PuiseuxSeries
except ImportError:
    PuiseuxSeries = ()

from sage.rings import (
    multi_power_series_ring_element,
    power_series_ring_element,
)
from sage.rings.integer import Integer
from sage.rings.rational import Rational


def O(*x, **kwds):
    """
    Big O constructor for various types.

    EXAMPLES:

    This is useful for writing power series elements::

        sage: R.<t> = ZZ[['t']]
        sage: (1+t)^10 + O(t^5)
        1 + 10*t + 45*t^2 + 120*t^3 + 210*t^4 + O(t^5)

    A power series ring is created implicitly if a polynomial
    element is passed::

        sage: R.<x> = QQ['x']
        sage: O(x^100)
        O(x^100)
        sage: 1/(1+x+O(x^5))
        1 - x + x^2 - x^3 + x^4 + O(x^5)
        sage: R.<u,v> = QQ[[]]
        sage: 1 + u + v^2 + O(u, v)^5
        1 + u + v^2 + O(u, v)^5

    This is also useful to create `p`-adic numbers::

        sage: O(7^6)                                                                    # needs sage.rings.padics
        O(7^6)
        sage: 1/3 + O(7^6)                                                              # needs sage.rings.padics
        5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + O(7^6)

    It behaves well with respect to adding negative powers of `p`::

        sage: a = O(11^-32); a                                                          # needs sage.rings.padics
        O(11^-32)
        sage: a.parent()                                                                # needs sage.rings.padics
        11-adic Field with capped relative precision 20

    There are problems if you add a rational with very negative
    valuation to an `O`-Term::

        sage: 11^-12 + O(11^15)                                                         # needs sage.rings.padics
        11^-12 + O(11^8)

    The reason that this fails is that the constructor doesn't know
    the right precision cap to use. If you cast explicitly or use
    other means of element creation, you can get around this issue::

        sage: # needs sage.rings.padics
        sage: K = Qp(11, 30)
        sage: K(11^-12) + O(11^15)
        11^-12 + O(11^15)
        sage: 11^-12 + K(O(11^15))
        11^-12 + O(11^15)
        sage: K(11^-12, absprec=15)
        11^-12 + O(11^15)
        sage: K(11^-12, 15)
        11^-12 + O(11^15)

    We can also work with `asymptotic expansions`_::

        sage: A.<n> = AsymptoticRing(growth_group='QQ^n * n^QQ * log(n)^QQ',            # needs sage.symbolic
        ....:                        coefficient_ring=QQ); A
        Asymptotic Ring <QQ^n * n^QQ * log(n)^QQ * Signs^n> over Rational Field
        sage: O(n)                                                                      # needs sage.symbolic
        O(n)

    Application with Puiseux series::

        sage: P.<y> = PuiseuxSeriesRing(ZZ)
        sage: y^(1/5) + O(y^(1/3))
        y^(1/5) + O(y^(1/3))
        sage: y^(1/3) + O(y^(1/5))
        O(y^(1/5))


    TESTS::

        sage: var('x, y')                                                               # needs sage.symbolic
        (x, y)
        sage: O(x)                                                                      # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ArithmeticError: O(x) not defined
        sage: O(y)                                                                      # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ArithmeticError: O(y) not defined
        sage: O(x, y)
        Traceback (most recent call last):
        ...
        ArithmeticError: O(x, y) not defined
        sage: O(4, 2)
        Traceback (most recent call last):
        ...
        ArithmeticError: O(4, 2) not defined

    ::

        sage: R.<x> = QQ[]
        sage: O(2*x)
        Traceback (most recent call last):
        ...
        NotImplementedError: completion only currently defined for the maximal ideal (x)
        sage: R.<x> = LazyPowerSeriesRing(QQ)
        sage: O(x^5)
        O(x^5)
        sage: t = O(Zp(5)(2*5^2)); t
        O(5^2)
        sage: t.parent()
        5-adic Ring with capped relative precision 20
        sage: t = O(Qp(5)(2*5^-2)); t
        O(5^-2)
        sage: t.parent()
        5-adic Field with capped relative precision 20
        sage: O(-6)
        Traceback (most recent call last):
        ...
        ArithmeticError: x must be a prime power >= 2
        sage: O(6)
        Traceback (most recent call last):
        ...
        ArithmeticError: x must be prime power
        sage: O(11/2)
        Traceback (most recent call last):
        ...
        ArithmeticError: x must be prime power
        sage: O(Rational(8))
        O(2^3)
    """
    if len(x) > 1:
        if isinstance(x[0], multi_power_series_ring_element.MPowerSeries):
            return multi_power_series_ring_element.MO(x, **kwds)
        else:
            raise ArithmeticError("O(%s) not defined" %
                                  (', '.join(str(e) for e in x),))

    x = x[0]

    if isinstance(x, power_series_ring_element.PowerSeries):
        return x.parent()(0, x.degree(), **kwds)

    if isinstance(x, Polynomial):
        if x.parent().ngens() != 1:
            raise NotImplementedError("completion only currently defined "
                                      "for univariate polynomials")
        if not x.is_monomial():
            raise NotImplementedError("completion only currently defined "
                                      "for the maximal ideal (x)")

    if isinstance(x, (int, Integer, Rational)):
        # p-adic number
        if x <= 0:
            raise ArithmeticError("x must be a prime power >= 2")
        if isinstance(x, (int, Integer)):
            x = Integer(x)
            p, r = x.perfect_power()
        else:
            if x.denominator() == 1:
                p, r = x.numerator().perfect_power()
            elif x.numerator() == 1:
                p, r = x.denominator().perfect_power()
                r = -r
            else:
                raise ArithmeticError("x must be prime power")
        if not p.is_prime():
            raise ArithmeticError("x must be prime power")
        if r >= 0:
            return Zp(p, prec=max(r, 20),
                      type='capped-rel')(0, absprec=r, **kwds)
        else:
            return Qp(p, prec=max(r, 20),
                      type='capped-rel')(0, absprec=r, **kwds)

    if isinstance(x, PuiseuxSeries):
        # note that add_bigoh() of PuiseuxSeries adapts the precision
        # to the ramification index of the input, thus we cannot do
        # zero.add_bigoh() because zero has ramification index 1
        return x.add_bigoh(x.valuation(), **kwds)

    try:
        return x.parent().zero().add_bigoh(x.valuation(), **kwds)
    except AttributeError:
        if hasattr(x, 'O'):  # this case is used for AsymptoticRing
            return x.O(**kwds)
        raise ArithmeticError("O(%s) not defined" % (x,))
