r"""
MacMahon's Partition Analysis Omega Operator

This module implements :func:`MacMahon's Omega Operator <MacMahonOmega>`
[Mac1915]_, which takes a quotient of Laurent polynomials and
removes all negative exponents in the corresponding power series.


Examples
========

In the following example, all negative exponents of `\mu` are removed.
The formula

.. MATH::

    \Omega_{\ge} \frac{1}{(1 - x\mu) (1 - y/\mu)}
    = \frac{1}{(1 - x) (1 - xy)}

can be calculated and verified by
::

    sage: L.<mu, x, y> = LaurentPolynomialRing(ZZ)
    sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu])
    1 * (-x + 1)^-1 * (-x*y + 1)^-1


Various
=======

AUTHORS:

- Daniel Krenn (2016)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.

Functions
=========
"""
# ****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import operator
from sage.misc.cachefunc import cached_function


def MacMahonOmega(var, expression, denominator=None, op=operator.ge,
          Factorization_sort=False, Factorization_simplify=True):
    r"""
    Return `\Omega_{\mathrm{op}}` of ``expression`` with respect to ``var``.

    To be more precise, calculate

    .. MATH::

        \Omega_{\mathrm{op}} \frac{n}{d_1 \dots d_n}

    for the numerator `n` and the factors `d_1`, ..., `d_n` of
    the denominator, all of which are Laurent polynomials in ``var``
    and return a (partial) factorization of the result.

    INPUT:

    - ``var`` -- a variable or a representation string of a variable

    - ``expression`` -- a
      :class:`~sage.structure.factorization.Factorization`
      of Laurent polynomials or, if ``denominator`` is specified,
      a Laurent polynomial interpreted as the numerator of the
      expression

    - ``denominator`` -- a Laurent polynomial or a
      :class:`~sage.structure.factorization.Factorization` (consisting
      of Laurent polynomial factors) or a tuple/list of factors (Laurent
      polynomials)

    - ``op`` -- (default: ``operator.ge``) an operator

      At the moment only ``operator.ge`` is implemented.

    - ``Factorization_sort`` (default: ``False``) and
      ``Factorization_simplify`` (default: ``True``) -- are passed on to
      :class:`sage.structure.factorization.Factorization` when creating
      the result

    OUTPUT:

    A (partial) :class:`~sage.structure.factorization.Factorization`
    of the result whose factors are Laurent polynomials

    .. NOTE::

        The numerator of the result may not be factored.

    REFERENCES:

    - [Mac1915]_

    - [APR2001]_

    EXAMPLES::

        sage: L.<mu, x, y, z, w> = LaurentPolynomialRing(ZZ)

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu])
        1 * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu, 1 - z/mu])
        1 * (-x + 1)^-1 * (-x*y + 1)^-1 * (-x*z + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu])
        (-x*y*z + 1) * (-x + 1)^-1 * (-y + 1)^-1 * (-x*z + 1)^-1 * (-y*z + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu^2])
        1 * (-x + 1)^-1 * (-x^2*y + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu^2, 1 - y/mu])
        (x*y + 1) * (-x + 1)^-1 * (-x*y^2 + 1)^-1

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu^2])
        (-x^2*y*z - x*y^2*z + x*y*z + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-x^2*z + 1)^-1 * (-y^2*z + 1)^-1

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu^3])
        1 * (-x + 1)^-1 * (-x^3*y + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y/mu^4])
        1 * (-x + 1)^-1 * (-x^4*y + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu^3, 1 - y/mu])
        (x*y^2 + x*y + 1) * (-x + 1)^-1 * (-x*y^3 + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu^4, 1 - y/mu])
        (x*y^3 + x*y^2 + x*y + 1) * (-x + 1)^-1 * (-x*y^4 + 1)^-1

        sage: MacMahonOmega(mu, 1, [1 - x*mu^2, 1 - y/mu, 1 - z/mu])
        (x*y*z + x*y + x*z + 1) *
        (-x + 1)^-1 * (-x*y^2 + 1)^-1 * (-x*z^2 + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu^2, 1 - y*mu, 1 - z/mu])
        (-x*y*z^2 - x*y*z + x*z + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-x*z^2 + 1)^-1 * (-y*z + 1)^-1

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z*mu, 1 - w/mu])
        (x*y*z*w^2 + x*y*z*w - x*y*w - x*z*w - y*z*w + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-z + 1)^-1 *
        (-x*w + 1)^-1 * (-y*w + 1)^-1 * (-z*w + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu, 1 - w/mu])
        (x^2*y*z*w + x*y^2*z*w - x*y*z*w - x*y*z - x*y*w + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 *
        (-x*z + 1)^-1 * (-x*w + 1)^-1 * (-y*z + 1)^-1 * (-y*w + 1)^-1

        sage: MacMahonOmega(mu, mu^-2, [1 - x*mu, 1 - y/mu])
        x^2 * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: MacMahonOmega(mu, mu^-1, [1 - x*mu, 1 - y/mu])
        x * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: MacMahonOmega(mu, mu, [1 - x*mu, 1 - y/mu])
        (-x*y + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: MacMahonOmega(mu, mu^2, [1 - x*mu, 1 - y/mu])
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

    We demonstrate the different allowed input variants::

        sage: MacMahonOmega(mu,
        ....:     Factorization([(mu, 2), (1 - x*mu, -1), (1 - y/mu, -1)]))
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: MacMahonOmega(mu, mu^2,
        ....:     Factorization([(1 - x*mu, 1), (1 - y/mu, 1)]))
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: MacMahonOmega(mu, mu^2, [1 - x*mu, 1 - y/mu])
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: MacMahonOmega(mu, mu^2, (1 - x*mu)*(1 - y/mu))  # not tested because not fully implemented
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: MacMahonOmega(mu, mu^2 / ((1 - x*mu)*(1 - y/mu)))  # not tested because not fully implemented
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

    TESTS::

        sage: MacMahonOmega(mu, 1, [1 - x*mu])
        1 * (-x + 1)^-1
        sage: MacMahonOmega(mu, 1, [1 - x/mu])
        1
        sage: MacMahonOmega(mu, 0, [1 - x*mu])
        0
        sage: MacMahonOmega(mu, L(1), [])
        1
        sage: MacMahonOmega(mu, L(0), [])
        0
        sage: MacMahonOmega(mu, 2, [])
        2
        sage: MacMahonOmega(mu, 2*mu, [])
        2
        sage: MacMahonOmega(mu, 2/mu, [])
        0

    ::

        sage: MacMahonOmega(mu, Factorization([(1/mu, 1), (1 - x*mu, -1),
        ....:                                  (1 - y/mu, -2)], unit=2))
        2*x * (-x + 1)^-1 * (-x*y + 1)^-2
        sage: MacMahonOmega(mu, Factorization([(mu, -1), (1 - x*mu, -1),
        ....:                                  (1 - y/mu, -2)], unit=2))
        2*x * (-x + 1)^-1 * (-x*y + 1)^-2
        sage: MacMahonOmega(mu, Factorization([(mu, -1), (1 - x, -1)]))
        0
        sage: MacMahonOmega(mu, Factorization([(2, -1)]))
        1 * 2^-1

    ::

        sage: MacMahonOmega(mu, 1, [1 - x*mu, 1 - z, 1 - y/mu])
        1 * (-z + 1)^-1 * (-x + 1)^-1 * (-x*y + 1)^-1

    ::

        sage: MacMahonOmega(mu, 1, [1 - x*mu], op=operator.lt)
        Traceback (most recent call last):
        ...
        NotImplementedError: At the moment, only Omega_ge is implemented.

        sage: MacMahonOmega(mu, 1, Factorization([(1 - x*mu, -1)]))
        Traceback (most recent call last):
        ...
        ValueError: Factorization (-mu*x + 1)^-1 of the denominator
        contains negative exponents.

        sage: MacMahonOmega(2*mu, 1, [1 - x*mu])
        Traceback (most recent call last):
        ...
        ValueError: 2*mu is not a variable.

        sage: MacMahonOmega(mu, 1, Factorization([(0, 2)]))
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Denominator contains a factor 0.

        sage: MacMahonOmega(mu, 1, [2 - x*mu])
        Traceback (most recent call last):
        ...
        NotImplementedError: Factor 2 - x*mu is not normalized.

        sage: MacMahonOmega(mu, 1, [1 - x*mu - mu^2])
        Traceback (most recent call last):
        ...
        NotImplementedError: Cannot handle factor 1 - x*mu - mu^2.

    ::

        sage: L.<mu, x, y, z, w> = LaurentPolynomialRing(QQ)
        sage: MacMahonOmega(mu, 1/mu,
        ....:     Factorization([(1 - x*mu, 1), (1 - y/mu, 2)], unit=2))
        1/2*x * (-x + 1)^-1 * (-x*y + 1)^-2
    """
    from sage.arith.misc import factor
    from sage.misc.misc_c import prod
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring \
        import LaurentPolynomialRing, LaurentPolynomialRing_univariate
    from sage.structure.factorization import Factorization

    if op != operator.ge:
        raise NotImplementedError('At the moment, only Omega_ge is implemented.')

    if denominator is None:
        if isinstance(expression, Factorization):
            numerator = expression.unit() * \
                        prod(f**e for f, e in expression if e > 0)
            denominator = tuple(f for f, e in expression if e < 0
                                for _ in range(-e))
        else:
            numerator = expression.numerator()
            denominator = expression.denominator()
    else:
        numerator = expression
    # at this point we have numerator/denominator

    if isinstance(denominator, (list, tuple)):
        factors_denominator = denominator
    else:
        if not isinstance(denominator, Factorization):
            denominator = factor(denominator)
        if not denominator.is_integral():
            raise ValueError('Factorization {} of the denominator '
                             'contains negative exponents.'.format(denominator))
        numerator *= ZZ(1) / denominator.unit()
        factors_denominator = tuple(factor
                                    for factor, exponent in denominator
                                    for _ in range(exponent))
    # at this point we have numerator/factors_denominator

    P = var.parent()
    if isinstance(P, LaurentPolynomialRing_univariate) and P.gen() == var:
        L = P
        L0 = L.base_ring()
    elif var in P.gens():
        var = repr(var)
        L0 = LaurentPolynomialRing(
            P.base_ring(), tuple(v for v in P.variable_names() if v != var))
        L = LaurentPolynomialRing(L0, var)
        var = L.gen()
    else:
        raise ValueError('{} is not a variable.'.format(var))

    other_factors = []
    to_numerator = []
    decoded_factors = []
    for factor in factors_denominator:
        factor = L(factor)
        D = factor.dict()
        if not D:
            raise ZeroDivisionError('Denominator contains a factor 0.')
        elif len(D) == 1:
            exponent, coefficient = next(iter(D.items()))
            if exponent == 0:
                other_factors.append(L0(factor))
            else:
                to_numerator.append(factor)
        elif len(D) == 2:
            if D.get(0, 0) != 1:
                raise NotImplementedError('Factor {} is not normalized.'.format(factor))
            D.pop(0)
            exponent, coefficient = next(iter(D.items()))
            decoded_factors.append((-coefficient, exponent))
        else:
            raise NotImplementedError('Cannot handle factor {}.'.format(factor))
    numerator = L(numerator) / prod(to_numerator)

    result_numerator, result_factors_denominator = \
        _Omega_(numerator.dict(), decoded_factors)
    if result_numerator == 0:
        return Factorization([], unit=result_numerator)

    return Factorization([(result_numerator, 1)] +
                         list((f, -1) for f in other_factors) +
                         list((1-f, -1) for f in result_factors_denominator),
                         sort=Factorization_sort,
                         simplify=Factorization_simplify)


def _simplify_(numerator, terms):
    r"""
    Cancels common factors of numerator and denominator.

    INPUT:

    - ``numerator`` -- a Laurent polynomial

    - ``terms`` -- tuple or other iterable of Laurent polynomials

      The denominator is the product of factors `1 - t` for each
      `t` in ``terms``.

    OUTPUT:

    A pair of a Laurent polynomial and a tuple of Laurent polynomials
    representing numerator and denominator as described in the INPUT-section.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import _simplify_
        sage: L.<x, y> = LaurentPolynomialRing(ZZ)
        sage: _simplify_(1-x^2, (x, y))
        (x + 1, (y,))

    TESTS::

        sage: _simplify_(1-x^2, (x, -x))
        (1, ())
        sage: _simplify_(1-x^2, (y^2, y))
        (-x^2 + 1, (y^2, y))
        sage: _simplify_(1-x^2, (x, L(2)))
        (x + 1, (2,))
    """
    new_terms = []
    for t in terms:
        if not t.is_constant():
            quo, rem = numerator.quo_rem(1 - t)
            if rem == 0:
                numerator = quo
                continue
        new_terms.append(t)
    return numerator, tuple(new_terms)


def _Omega_(A, decoded_factors):
    r"""
    Helper function for :func:`MacMahonOmega` which accesses the low level functions
    and does the substituting.

    INPUT:

    - ``A`` -- dictionary mapping `a` to `c` representing a summand
      `c\mu^a` of the numerator

    - ``decoded_factors`` -- tuple or list of pairs `(z, e)` representing
      a factor `1 - z \mu^e`

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a Laurent polynomial, its second component a factorization
    of the denominator as a tuple of Laurent polynomials, where each
    Laurent polynomial `z` represents a factor `1 - z`.

    TESTS:

    Extensive testing of this function is done in :func:`MacMahonOmega`.

    ::

        sage: L.<mu, x, y> = LaurentPolynomialRing(ZZ)
        sage: MacMahonOmega(mu, mu^-2, [1 - x*mu, 1 - y/mu])
        x^2 * (-x + 1)^-1 * (-x*y + 1)^-1

    internally calls
    ::

        sage: from sage.rings.polynomial.omega import _Omega_
        sage: _Omega_({-2: 1}, [(x, 1), (y, -1)])
        (x^2, (x, x*y))

    ::

        sage: _Omega_({0: 2, 1: 40, -1: -3}, [])
        (42, ())
        sage: _Omega_({-1: 42}, [])
        (0, ())


    ::

        sage: MacMahonOmega(mu, 1 - x^2, [1 - x*mu, 1 - y/mu])
        (x + 1) * (-x*y + 1)^-1
    """
    if not decoded_factors:
        return sum(c for a, c in A.items() if a >= 0), tuple()

    # Below we sort to make the caching more efficient. Doing this here
    # (in contrast to directly in Omega_ge) results in much cleaner
    # code and prevents an additional substitution or passing of a permutation.
    values, exponents = zip(*sorted(decoded_factors, key=lambda k: -k[1]))

    numerator = 0
    factors_denominator = None
    rules = None
    for a, c in A.items():
        n, fd = Omega_ge(a, exponents)
        if factors_denominator is None:
            factors_denominator = fd
        else:
            assert factors_denominator == fd
        if rules is None:
            rules = dict(zip(n.parent().gens(), values))
        numerator += c * n.subs(rules)

    if numerator == 0:
        factors_denominator = tuple()
    return _simplify_(numerator,
                      tuple(f.subs(rules) for f in factors_denominator))


@cached_function
def Omega_ge(a, exponents):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    To be more precise, calculate

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - z_0 \mu^{e_0}) \dots (1 - z_{n-1} \mu^{e_{n-1}})}

    and return its numerator and a factorization of its denominator.
    Note that `z_0`, ..., `z_{n-1}` only appear in the output, but not in the
    input.

    INPUT:

    - ``a`` -- integer

    - ``exponents`` -- tuple of integers

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a Laurent polynomial, its second component a factorization
    of the denominator as a tuple of Laurent polynomials, where each
    Laurent polynomial `z` represents a factor `1 - z`.

    The parents of these Laurent polynomials is always a
    Laurent polynomial ring in `z_0`, ..., `z_{n-1}` over `\ZZ`, where
    `n` is the length of ``exponents``.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import Omega_ge
        sage: Omega_ge(0, (1, -2))
        (1, (z0, z0^2*z1))
        sage: Omega_ge(0, (1, -3))
        (1, (z0, z0^3*z1))
        sage: Omega_ge(0, (1, -4))
        (1, (z0, z0^4*z1))

        sage: Omega_ge(0, (2, -1))
        (z0*z1 + 1, (z0, z0*z1^2))
        sage: Omega_ge(0, (3, -1))
        (z0*z1^2 + z0*z1 + 1, (z0, z0*z1^3))
        sage: Omega_ge(0, (4, -1))
        (z0*z1^3 + z0*z1^2 + z0*z1 + 1, (z0, z0*z1^4))

        sage: Omega_ge(0, (1, 1, -2))
        (-z0^2*z1*z2 - z0*z1^2*z2 + z0*z1*z2 + 1, (z0, z1, z0^2*z2, z1^2*z2))
        sage: Omega_ge(0, (2, -1, -1))
        (z0*z1*z2 + z0*z1 + z0*z2 + 1, (z0, z0*z1^2, z0*z2^2))
        sage: Omega_ge(0, (2, 1, -1))
        (-z0*z1*z2^2 - z0*z1*z2 + z0*z2 + 1, (z0, z1, z0*z2^2, z1*z2))

    ::

        sage: Omega_ge(0, (2, -2))
        (-z0*z1 + 1, (z0, z0*z1, z0*z1))
        sage: Omega_ge(0, (2, -3))
        (z0^2*z1 + 1, (z0, z0^3*z1^2))
        sage: Omega_ge(0, (3, 1, -3))
        (-z0^3*z1^3*z2^3 + 2*z0^2*z1^3*z2^2 - z0*z1^3*z2
         + z0^2*z2^2 - 2*z0*z2 + 1,
         (z0, z1, z0*z2, z0*z2, z0*z2, z1^3*z2))

    ::

        sage: Omega_ge(0, (3, 6, -1))
        (-z0*z1*z2^8 - z0*z1*z2^7 - z0*z1*z2^6 - z0*z1*z2^5 - z0*z1*z2^4 +
         z1*z2^5 - z0*z1*z2^3 + z1*z2^4 - z0*z1*z2^2 + z1*z2^3 -
         z0*z1*z2 + z0*z2^2 + z1*z2^2 + z0*z2 + z1*z2 + 1,
         (z0, z1, z0*z2^3, z1*z2^6))

    TESTS::

        sage: Omega_ge(0, (2, 2, 1, 1, 1, -1, -1))[0].number_of_terms()  # long time
        1695
        sage: Omega_ge(0, (2, 2, 1, 1, 1, 1, 1, -1, -1))[0].number_of_terms()  # not tested (too long, 1 min)
        27837

    ::

        sage: Omega_ge(1, (2,))
        (1, (z0,))
    """
    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_ge: a=%s, exponents=%s', a, exponents)

    from sage.arith.functions import lcm
    from sage.arith.srange import srange
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.rings.number_field.number_field import CyclotomicField

    if not exponents or any(e == 0 for e in exponents):
        raise NotImplementedError

    rou = sorted(set(abs(e) for e in exponents) - set([1]))
    ellcm = lcm(rou)
    B = CyclotomicField(ellcm, 'zeta')
    zeta = B.gen()
    z_names = tuple('z{}'.format(i) for i in range(len(exponents)))
    L = LaurentPolynomialRing(B, ('t',) + z_names, len(z_names) + 1)
    t = L.gens()[0]
    Z = LaurentPolynomialRing(ZZ, z_names, len(z_names))
    powers = {i: L(zeta**(ellcm//i)) for i in rou}
    powers[2] = L(-1)
    powers[1] = L(1)
    exponents_and_values = tuple(
        (e, tuple(powers[abs(e)]**j * z for j in srange(abs(e))))
        for z, e in zip(L.gens()[1:], exponents))
    x = tuple(v for e, v in exponents_and_values if e > 0)
    y = tuple(v for e, v in exponents_and_values if e < 0)

    def subs_power(expression, var, exponent):
        r"""
        Substitute ``var^exponent`` by ``var`` in ``expression``.

        It is assumed that ``var`` only occurs with exponents
        divisible by ``exponent``.
        """
        p = tuple(var.dict().popitem()[0]).index(1)  # var is the p-th generator

        def subs_e(e):
            e = list(e)
            assert e[p] % exponent == 0
            e[p] = e[p] // exponent
            return tuple(e)
        parent = expression.parent()
        result = parent({subs_e(e): c for e, c in expression.dict().items()})
        return result

    def de_power(expression):
        expression = Z(expression)
        for e, var in zip(exponents, Z.gens()):
            if abs(e) == 1:
                continue
            expression = subs_power(expression, var, abs(e))
        return expression

    logger.debug('Omega_ge: preparing denominator')
    factors_denominator = tuple(de_power(1 - factor)
                                for factor in _Omega_factors_denominator_(x, y))

    logger.debug('Omega_ge: preparing numerator')
    numerator = de_power(_Omega_numerator_(a, x, y, t))

    logger.info('Omega_ge: completed')
    return numerator, factors_denominator


def _Omega_numerator_(a, x, y, t):
    r"""
    Return the numerator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, calculate

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)}

    and return its numerator.

    This function is meant to be a helper function of :func:`MacMahonOmega`.

    INPUT:

    - ``a`` -- integer

    - ``x``, ``y`` -- tuple of tuples of Laurent polynomials

      The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.
      The non-flatness of these parameters is to be interface-consistent
      with :func:`_Omega_factors_denominator_`.

    - ``t`` -- a temporary Laurent polynomial variable used for substituting

    OUTPUT: a Laurent polynomial

    The output is normalized such that the corresponding denominator
    (:func:`_Omega_factors_denominator_`) has constant term `1`.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import _Omega_numerator_, _Omega_factors_denominator_

        sage: L.<x0, x1, x2, x3, y0, y1, t> = LaurentPolynomialRing(ZZ)
        sage: _Omega_numerator_(0, ((x0,),), ((y0,),), t)
        1
        sage: _Omega_numerator_(0, ((x0,), (x1,)), ((y0,),), t)
        -x0*x1*y0 + 1
        sage: _Omega_numerator_(0, ((x0,),), ((y0,), (y1,)), t)
        1
        sage: _Omega_numerator_(0, ((x0,), (x1,), (x2,)), ((y0,),), t)
        x0*x1*x2*y0^2 + x0*x1*x2*y0 - x0*x1*y0 - x0*x2*y0 - x1*x2*y0 + 1
        sage: _Omega_numerator_(0, ((x0,), (x1,)), ((y0,), (y1,)), t)
        x0^2*x1*y0*y1 + x0*x1^2*y0*y1 - x0*x1*y0*y1 - x0*x1*y0 - x0*x1*y1 + 1

        sage: _Omega_numerator_(-2, ((x0,),), ((y0,),), t)
        x0^2
        sage: _Omega_numerator_(-1, ((x0,),), ((y0,),), t)
        x0
        sage: _Omega_numerator_(1, ((x0,),), ((y0,),), t)
        -x0*y0 + y0 + 1
        sage: _Omega_numerator_(2, ((x0,),), ((y0,),), t)
        -x0*y0^2 - x0*y0 + y0^2 + y0 + 1

    TESTS::

        sage: _Omega_factors_denominator_((), ())
        ()
        sage: _Omega_numerator_(0, (), (), t)
        1
        sage: _Omega_numerator_(+2, (), (), t)
        1
        sage: _Omega_numerator_(-2, (), (), t)
        0

        sage: _Omega_factors_denominator_(((x0,),), ())
        (-x0 + 1,)
        sage: _Omega_numerator_(0, ((x0,),), (), t)
        1
        sage: _Omega_numerator_(+2, ((x0,),), (), t)
        1
        sage: _Omega_numerator_(-2, ((x0,),), (), t)
        x0^2

        sage: _Omega_factors_denominator_((), ((y0,),))
        ()
        sage: _Omega_numerator_(0, (), ((y0,),), t)
        1
        sage: _Omega_numerator_(+2, (), ((y0,),), t)
        y0^2 + y0 + 1
        sage: _Omega_numerator_(-2, (), ((y0,),), t)
        0

    ::

        sage: L.<X, Y, t> = LaurentPolynomialRing(ZZ)
        sage: _Omega_numerator_(2, ((X,),), ((Y,),), t)
        -X*Y^2 - X*Y + Y^2 + Y + 1
    """
    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    x_flat = sum(x, tuple())
    y_flat = sum(y, tuple())
    n = len(x_flat)
    m = len(y_flat)
    xy = x_flat + y_flat

    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_numerator: a=%s, n=%s, m=%s', a, n, m)

    if m == 0:
        result = 1 - (prod(_Omega_factors_denominator_(x, y)) *
                      sum(homogeneous_symmetric_function(j, xy)
                          for j in srange(-a))
                      if a < 0 else 0)
    elif n == 0:
        result = sum(homogeneous_symmetric_function(j, xy)
                     for j in srange(a+1))
    else:
        result = _Omega_numerator_P_(a, x_flat[:-1], y_flat, t).subs({t: x_flat[-1]})
    L = t.parent()
    result = L(result)

    logger.info('_Omega_numerator_: %s terms', result.number_of_terms())
    return result


def _Omega_numerator_P_(a, x, y, t):
    r"""
    Helper function for :func:`_Omega_numerator_`.

    This is an implementation of the function `P` of [APR2001]_.

    INPUT:

    - ``a`` -- integer

    - ``x``, ``y`` -- tuple of Laurent polynomials

      The tuple ``x`` here is the flattened ``x`` of :func:`_Omega_numerator_`
      but without its last entry.

    - ``t`` -- a temporary Laurent polynomial variable

      In the (final) result, ``t`` has to be substituted by the last
      entry of the flattened ``x`` of :func:`_Omega_numerator_`.

    OUTPUT: a Laurent polynomial

    TESTS::

        sage: from sage.rings.polynomial.omega import _Omega_numerator_P_
        sage: L.<x0, x1, y0, y1, t> = LaurentPolynomialRing(ZZ)
        sage: _Omega_numerator_P_(0, (x0,), (y0,), t).subs({t: x1})
        -x0*x1*y0 + 1
    """
    # This function takes Laurent polynomials as inputs. It would
    # be possible to input only the sizes of ``x`` and ``y`` and
    # perform a substitution afterwards; in this way caching of this
    # function would make sense. However, the way it is now allows
    # automatic collection and simplification of the summands, which
    # makes it more efficient for higher powers at the input of
    # :func:`Omega_ge`.
    # Caching occurs in :func:`Omega_ge`.

    import logging
    logger = logging.getLogger(__name__)

    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    n = len(x)
    if n == 0:
        x0 = t
        result = x0**(-a) + \
            (prod(1 - x0*yy for yy in y) *
             sum(homogeneous_symmetric_function(j, y) * (1-x0**(j-a))
                 for j in srange(a))
             if a > 0 else 0)
    else:
        Pprev = _Omega_numerator_P_(a, x[:n-1], y, t)
        x2 = x[n-1]
        logger.debug('Omega_numerator: P(%s): substituting...', n)
        x1 = t
        p1 = Pprev
        p2 = Pprev.subs({t: x2})
        logger.debug('Omega_numerator: P(%s): preparing...', n)
        dividend = x1 * (1-x2) * prod(1 - x2*yy for yy in y) * p1 - \
                x2 * (1-x1) * prod(1 - x1*yy for yy in y) * p2
        logger.debug('Omega_numerator: P(%s): dividing...', n)
        q, r = dividend.quo_rem(x1 - x2)
        assert r == 0
        result = q
    logger.debug('Omega_numerator: P(%s) has %s terms', n, result.number_of_terms())
    return result


@cached_function
def _Omega_factors_denominator_(x, y):
    r"""
    Return the denominator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, calculate

    .. MATH::

        \Omega_{\ge} \frac{1}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)}

    and return a factorization of its denominator.

    This function is meant to be a helper function of :func:`MacMahonOmega`.

    INPUT:

    - ``x``, ``y`` -- tuple of tuples of Laurent polynomials

      The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.

    OUTPUT:

    A factorization of the denominator as
    a tuple of Laurent polynomials

    The output is normalized such that it has constant term `1`.

    .. NOTE::

        The assumption is that the ``x`` and ``y`` are collected in
        such a way that one entry of ``x`` corresponds to the orbit of
        some ``x_j`` under multiplication by `d`-th roots of unity and that
        the output is collected in a corresponding way.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import _Omega_factors_denominator_

        sage: L.<x0, x1, x2, x3, y0, y1> = LaurentPolynomialRing(ZZ)
        sage: _Omega_factors_denominator_(((x0,),), ((y0,),))
        (-x0 + 1, -x0*y0 + 1)
        sage: _Omega_factors_denominator_(((x0,),), ((y0,), (y1,)))
        (-x0 + 1, -x0*y0 + 1, -x0*y1 + 1)
        sage: _Omega_factors_denominator_(((x0,), (x1,)), ((y0,),))
        (-x0 + 1, -x1 + 1, -x0*y0 + 1, -x1*y0 + 1)
        sage: _Omega_factors_denominator_(((x0,), (x1,), (x2,)), ((y0,),))
        (-x0 + 1, -x1 + 1, -x2 + 1, -x0*y0 + 1, -x1*y0 + 1, -x2*y0 + 1)
        sage: _Omega_factors_denominator_(((x0,), (x1,)), ((y0,), (y1,)))
        (-x0 + 1, -x1 + 1, -x0*y0 + 1, -x0*y1 + 1, -x1*y0 + 1, -x1*y1 + 1)

    ::

        sage: # needs sage.rings.number_field
        sage: B.<zeta> = ZZ.extension(cyclotomic_polynomial(3))
        sage: L.<x, y> = LaurentPolynomialRing(B)
        sage: _Omega_factors_denominator_(((x, -x),), ((y,),))
        (-x^2 + 1, -x^2*y^2 + 1)
        sage: _Omega_factors_denominator_(((x, -x),), ((y, zeta*y, zeta^2*y),))
        (-x^2 + 1, -x^6*y^6 + 1)
        sage: _Omega_factors_denominator_(((x, -x),), ((y, -y),))
        (-x^2 + 1, -x^2*y^2 + 1, -x^2*y^2 + 1)

    TESTS::

        sage: L.<x0, y0> = LaurentPolynomialRing(ZZ)
        sage: _Omega_factors_denominator_((), ())
        ()
        sage: _Omega_factors_denominator_(((x0,),), ())
        (-x0 + 1,)
        sage: _Omega_factors_denominator_((), ((y0,),))
        ()
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.misc.misc_c import prod

    result = tuple(prod(1 - xx for xx in gx) for gx in x) + \
             sum(((prod(1 - xx*yy for xx in gx for yy in gy),)
                  if len(gx) != len(gy)
                  else tuple(prod(1 - xx*yy for xx in gx) for yy in gy)
                  for gx in x for gy in y),
                 tuple())

    logger.info('Omega_denominator: %s factors', len(result))
    return result


def partition(items, predicate=bool):
    r"""
    Split ``items`` into two parts by the given ``predicate``.

    INPUT:

    - ``item`` -- an iterator

    - ``predicate`` -- a function

    OUTPUT:

    A pair of iterators; the first contains the elements not satisfying
    the ``predicate``, the second the elements satisfying the ``predicate``.

    ALGORITHM:

    Source of the code:
    `http://nedbatchelder.com/blog/201306/filter_a_list_into_two_parts.html
    <http://nedbatchelder.com/blog/201306/filter_a_list_into_two_parts.html>`_

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import partition
        sage: E, O = partition(srange(10), is_odd)
        sage: tuple(E), tuple(O)
        ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))
    """
    from itertools import tee
    a, b = tee((predicate(item), item) for item in items)
    return ((item for pred, item in a if not pred),
            (item for pred, item in b if pred))


def homogeneous_symmetric_function(j, x):
    r"""
    Return a complete homogeneous symmetric polynomial
    (:wikipedia:`Complete_homogeneous_symmetric_polynomial`).

    INPUT:

    - ``j`` -- the degree as a nonnegative integer

    - ``x`` -- an iterable of variables

    OUTPUT: a polynomial of the common parent of all entries of ``x``

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import homogeneous_symmetric_function
        sage: P = PolynomialRing(ZZ, 'X', 3)
        sage: homogeneous_symmetric_function(0, P.gens())
        1
        sage: homogeneous_symmetric_function(1, P.gens())
        X0 + X1 + X2
        sage: homogeneous_symmetric_function(2, P.gens())
        X0^2 + X0*X1 + X1^2 + X0*X2 + X1*X2 + X2^2
        sage: homogeneous_symmetric_function(3, P.gens())
        X0^3 + X0^2*X1 + X0*X1^2 + X1^3 + X0^2*X2 +
        X0*X1*X2 + X1^2*X2 + X0*X2^2 + X1*X2^2 + X2^3
    """
    from sage.combinat.integer_vector import IntegerVectors
    from sage.misc.misc_c import prod

    return sum(prod(xx**pp for xx, pp in zip(x, p))
               for p in IntegerVectors(j, length=len(x)))
