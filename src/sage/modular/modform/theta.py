# sage_setup: distribution = sagemath-schemes
"""
`q`-expansions of theta series

AUTHOR:

- William Stein
"""
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing

from math import sqrt


def theta2_qexp(prec=10, var='q', K=ZZ, sparse=False):
    r"""
    Return the `q`-expansion of the series `\theta_2 = \sum_{n \text{ odd}} q^{n^2}`.

    INPUT:

    - ``prec`` -- integer; the absolute precision of the output
    - ``var`` -- (default: ``'q'``) variable name
    - ``K`` -- (default: ZZ) base ring of answer

    OUTPUT: a power series over K

    EXAMPLES::

        sage: theta2_qexp(18)
        q + q^9 + O(q^18)
        sage: theta2_qexp(49)
        q + q^9 + q^25 + O(q^49)
        sage: theta2_qexp(100, 'q', QQ)
        q + q^9 + q^25 + q^49 + q^81 + O(q^100)
        sage: f = theta2_qexp(100, 't', GF(3)); f
        t + t^9 + t^25 + t^49 + t^81 + O(t^100)
        sage: parent(f)
        Power Series Ring in t over Finite Field of size 3
        sage: theta2_qexp(200)
        q + q^9 + q^25 + q^49 + q^81 + q^121 + q^169 + O(q^200)
        sage: f = theta2_qexp(20,sparse=True); f
        q + q^9 + O(q^20)
        sage: parent(f)
        Sparse Power Series Ring in q over Integer Ring
    """
    prec = Integer(prec)
    if prec <= 0:
        raise ValueError("prec must be positive")
    if sparse:
        v = {}
    else:
        v = [Integer(0)] * prec
    one = Integer(1)
    n = int(sqrt(prec))
    if n*n < prec:
        n += 1
    for m in range(1, n, 2):
        v[m*m] = one
    R = PowerSeriesRing(K, sparse=sparse, names=var)
    return R(v, prec=prec)


def theta_qexp(prec=10, var='q', K=ZZ, sparse=False):
    r"""
    Return the `q`-expansion of the standard `\theta` series
    `\theta = 1 + 2\sum_{n=1}^{\infty} q^{n^2}`.

    INPUT:

    - ``prec`` -- integer; the absolute precision of the output
    - ``var`` -- (default: ``'q'``) variable name
    - ``K`` -- (default: ZZ) base ring of answer

    OUTPUT: a power series over K

    EXAMPLES::

        sage: theta_qexp(25)
        1 + 2*q + 2*q^4 + 2*q^9 + 2*q^16 + O(q^25)
        sage: theta_qexp(10)
        1 + 2*q + 2*q^4 + 2*q^9 + O(q^10)
        sage: theta_qexp(100)
        1 + 2*q + 2*q^4 + 2*q^9 + 2*q^16 + 2*q^25 + 2*q^36 + 2*q^49 + 2*q^64 + 2*q^81 + O(q^100)
        sage: theta_qexp(100, 't')
        1 + 2*t + 2*t^4 + 2*t^9 + 2*t^16 + 2*t^25 + 2*t^36 + 2*t^49 + 2*t^64 + 2*t^81 + O(t^100)
        sage: theta_qexp(100, 't', GF(2))
        1 + O(t^100)
        sage: f = theta_qexp(20,sparse=True); f
        1 + 2*q + 2*q^4 + 2*q^9 + 2*q^16 + O(q^20)
        sage: parent(f)
        Sparse Power Series Ring in q over Integer Ring
    """
    prec = Integer(prec)
    if prec <= 0:
        raise ValueError("prec must be positive")
    if sparse:
        v = {}
    else:
        v = [Integer(0)] * prec
    v[0] = Integer(1)
    two = Integer(2)
    n = int(sqrt(prec))
    if n*n != prec:
        n += 1
    for m in range(1, n):
        v[m*m] = two

    R = PowerSeriesRing(K, sparse=sparse, names=var)
    return R(v, prec=prec)
