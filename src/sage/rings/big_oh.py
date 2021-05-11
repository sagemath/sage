"""
Big O for various types (power series, p-adics, etc.)

.. SEEALSO::

    - `asymptotic expansions <../../../asymptotic/index.html>`_
    - `p-adic numbers <../../../padics/index.html>`_
    - `power series <../../../power_series/index.html>`_
    - `polynomials <../../../polynomial_rings/index.html>`_
"""

import sage.arith.all as arith
from . import laurent_series_ring_element
from sage.rings.puiseux_series_ring_element import PuiseuxSeries
import sage.rings.padics.factory as padics_factory
import sage.rings.padics.padic_generic_element as padic_generic_element
from . import power_series_ring_element
from . import integer
from . import rational
from sage.rings.polynomial.polynomial_element import Polynomial
from . import multi_power_series_ring_element


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

        sage: O(7^6)
        O(7^6)
        sage: 1/3 + O(7^6)
        5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + O(7^6)

    It behaves well with respect to adding negative powers of `p`::

        sage: a = O(11^-32); a
        O(11^-32)
        sage: a.parent()
        11-adic Field with capped relative precision 20

    There are problems if you add a rational with very negative
    valuation to an `O`-Term::

        sage: 11^-12 + O(11^15)
        11^-12 + O(11^8)

    The reason that this fails is that the constructor doesn't know
    the right precision cap to use. If you cast explicitly or use
    other means of element creation, you can get around this issue::

        sage: K = Qp(11, 30)
        sage: K(11^-12) + O(11^15)
        11^-12 + O(11^15)
        sage: 11^-12 + K(O(11^15))
        11^-12 + O(11^15)
        sage: K(11^-12, absprec = 15)
        11^-12 + O(11^15)
        sage: K(11^-12, 15)
        11^-12 + O(11^15)

    We can also work with `asymptotic expansions`_::

        sage: A.<n> = AsymptoticRing(growth_group='QQ^n * n^QQ * log(n)^QQ', coefficient_ring=QQ); A
        Asymptotic Ring <QQ^n * n^QQ * log(n)^QQ * Signs^n> over Rational Field
        sage: O(n)
        O(n)

    Application with Puiseux series::

        sage: P.<y> = PuiseuxSeriesRing(ZZ)
        sage: y^(1/5) + O(y^(1/3))
        y^(1/5) + O(y^(1/3))
        sage: y^(1/3) + O(y^(1/5))
        O(y^(1/5))


    TESTS::

        sage: var('x, y')
        (x, y)
        sage: O(x)
        Traceback (most recent call last):
        ...
        ArithmeticError: O(x) not defined
        sage: O(y)
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

    elif isinstance(x, Polynomial):
        if x.parent().ngens() != 1:
            raise NotImplementedError("completion only currently defined "
                                      "for univariate polynomials")
        if not x.is_monomial():
            raise NotImplementedError("completion only currently defined "
                                      "for the maximal ideal (x)")
        return x.parent().completion(x.parent().gen())(0, x.degree(), **kwds)

    elif isinstance(x, laurent_series_ring_element.LaurentSeries):
        return laurent_series_ring_element.LaurentSeries(x.parent(), 0).\
            add_bigoh(x.valuation(), **kwds)

    elif isinstance(x, PuiseuxSeries):
        return x.add_bigoh(x.valuation(), **kwds)

    elif isinstance(x, (int, integer.Integer, rational.Rational)):
        # p-adic number
        if x <= 0:
            raise ArithmeticError("x must be a prime power >= 2")
        F = arith.factor(x)
        if len(F) != 1:
            raise ArithmeticError("x must be prime power")
        p, r = F[0]
        if r >= 0:
            return padics_factory.Zp(p, prec=max(r, 20),
                                     type='capped-rel')(0, absprec=r, **kwds)
        else:
            return padics_factory.Qp(p, prec=max(r, 20),
                                     type='capped-rel')(0, absprec=r, **kwds)

    elif isinstance(x, padic_generic_element.pAdicGenericElement):
         return x.parent()(0, absprec=x.valuation(), **kwds)
    elif hasattr(x, 'O'):
        return x.O(**kwds)
    raise ArithmeticError("O(%s) not defined" % (x,))
