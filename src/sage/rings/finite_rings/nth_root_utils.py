# sage.doctest: needs sage.rings.finite_rings
r"""
Utilities for :meth:`~sage.rings.finite_rings.element_base.FiniteRingElement.nth_root`.

AUTHORS:

- SageMath developers (2026)
"""
# ****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.structure.factorization import Factorization


def gcd_factorization_from_order_factorization(order_factorization, n_exp, q):
    r"""
    Return the factorization of `\gcd(n_{\text{exp}}, q-1)` from a factorization
    of the multiplicative group order `q-1`.

    INPUT:

    - ``order_factorization`` -- iterable of pairs `(p, e)` with `p` prime and
      `e \geq 1`, whose product is `q-1`. A :class:`~sage.structure.factorization.Factorization`
      object is allowed.

    - ``n_exp`` -- integer (the original exponent passed to :meth:`nth_root`)

    - ``q`` -- integer (the cardinality of the base field)

    The factorization of `\gcd(n_{\text{exp}}, q-1)` is obtained from
    the factorization of `q-1` by taking `\min(v_p(n_{\text{exp}}), e)` for
    each prime power `p^e` in the factorization of `q-1`.

    EXAMPLES::

        sage: from sage.rings.finite_rings.nth_root_utils import gcd_factorization_from_order_factorization
        sage: gcd_factorization_from_order_factorization([(2, 4)], 4, 17)
        2^4
        sage: gcd_factorization_from_order_factorization([(2, 4)], 8, 17)
        2^4
        sage: gcd_factorization_from_order_factorization([(2, 4)], 3, 17)
        1

    TESTS::

        sage: from sage.rings.finite_rings.nth_root_utils import gcd_factorization_from_order_factorization
        sage: gcd_factorization_from_order_factorization(factor(16), 4, 17)
        2^4
        sage: gcd_factorization_from_order_factorization([(2, 4)], 4, 17) == 16.factor()
        True
        sage: gcd_factorization_from_order_factorization([(2, 3)], 4, 17)
        Traceback (most recent call last):
        ...
        ValueError: order_factorization does not multiply to q-1
        sage: gcd_factorization_from_order_factorization([(4, 1)], 4, 17)
        Traceback (most recent call last):
        ...
        ValueError: 4 is not prime
        sage: gcd_factorization_from_order_factorization([(2, 3), (2, 1)], 4, 17)
        Traceback (most recent call last):
        ...
        ValueError: duplicate prime in order_factorization
    """
    n_exp = ZZ(n_exp)
    q = ZZ(q)
    q1 = q - 1
    if q1 <= 0:
        raise ValueError("q must be at least 2")

    pairs = []
    prod = ZZ(1)
    for item in order_factorization:
        if len(item) != 2:
            raise ValueError("order_factorization must be a list of (prime, exponent) pairs")
        p, e = ZZ(item[0]), ZZ(item[1])
        if p < 2:
            raise ValueError("primes in order_factorization must be at least 2")
        if not p.is_prime():
            raise ValueError(f"{p} is not prime")
        if e < 1:
            raise ValueError("exponents in order_factorization must be positive")
        if prod % p == 0:
            raise ValueError("duplicate prime in order_factorization")
        prod *= p**e
        pairs.append((p, e))

    if prod != q1:
        raise ValueError("order_factorization does not multiply to q-1")

    gcd_fac = []
    for p, e in pairs:
        v = min(n_exp.valuation(p), e)
        if v > 0:
            gcd_fac.append((p, v))

    return Factorization(gcd_fac) if gcd_fac else ZZ(1).factor()
