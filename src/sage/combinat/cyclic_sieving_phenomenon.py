r"""
Cyclic sieving phenomenon

Implementation of the Cyclic Sieving Phenomenon as described by
Reiner, Stanton, and White in [RSW2004]_.

We define the :func:`CyclicSievingPolynomial` of a finite set `S`
together with cyclic action ``cyc_act`` (of order `n`) to be the
unique polynomial ``P(q)`` of order < `n` such that the triple (`S`,
``cyc_act``, ``P(q)``) exhibits the cyclic sieving phenomenon.

AUTHORS:

- Christian Stump

REFERENCES:

.. [RSW2004] Reiner, Stanton, White - *The cyclic sieving phenomenon*,
             Journal of Combinatorial Theory A 108 (2004)
"""
# ****************************************************************************
#       Copyright (C) 2010 Christian Stump christian.stump@univie.ac.at
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from sage.rings.integer_ring import ZZ
from sage.arith.all import lcm
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


def CyclicSievingPolynomial(L, cyc_act=None, order=None, get_order=False):
    """
    Return the unique polynomial ``p`` of degree smaller than ``order`` such
    that the triple ``(L, cyc_act, p)`` exhibits the Cyclic Sieving Phenomenon.

    If ``cyc_act`` is None, ``L`` is expected to contain the orbit lengths.

    INPUT:

    - ``L`` -- if ``cyc_act`` is ``None``: list of orbit sizes,
      otherwise list of objects

    - ``cyc_act`` -- (default:``None``) bijective function from ``L`` to ``L``

    - ``order`` -- (default:``None``) if set to an integer, this
        cyclic order of ``cyc_act`` is used (must be an integer multiple
        of the order of ``cyc_act``) otherwise, the order of ``cyc_action`` is
        used

    - ``get_order`` -- (default:``False``) if ``True``, a tuple ``[p,n]``
      is returned where ``p`` is as above, and ``n`` is the order

    EXAMPLES::

        sage: from sage.combinat.cyclic_sieving_phenomenon import CyclicSievingPolynomial
        sage: S42 = Subsets([1,2,3,4], 2)
        sage: def cyc_act(S): return Set(i.mod(4) + 1 for i in S)
        sage: cyc_act([1,3])
        {2, 4}
        sage: cyc_act([1,4])
        {1, 2}
        sage: CyclicSievingPolynomial(S42, cyc_act)
        q^3 + 2*q^2 + q + 2
        sage: CyclicSievingPolynomial(S42, cyc_act, get_order=True)
        [q^3 + 2*q^2 + q + 2, 4]
        sage: CyclicSievingPolynomial(S42, cyc_act, order=8)
        q^6 + 2*q^4 + q^2 + 2
        sage: CyclicSievingPolynomial([4,2])
        q^3 + 2*q^2 + q + 2

    TESTS:

    We check that :trac:`13997` is handled::

        sage: CyclicSievingPolynomial(S42, cyc_act, order=8, get_order=True)
        [q^6 + 2*q^4 + q^2 + 2, 8]
        sage: CyclicSievingPolynomial(S42, cyc_act, order=11)
        Traceback (most recent call last):
        ...
        ValueError: order is not a multiple of the order of the cyclic action
    """
    if cyc_act:
        orbits = orbit_decomposition(L, cyc_act)
    else:
        orbits = [list(range(k)) for k in L]

    R = PolynomialRing(ZZ, 'q')
    q = R.gen()
    p = R.zero()

    orbit_sizes = {}
    for orbit in orbits:
        length = len(orbit)
        if length in orbit_sizes:
            orbit_sizes[length] += 1
        else:
            orbit_sizes[length] = 1

    n = lcm(orbit_sizes)

    if order:
        if order.mod(n):
            raise ValueError("order is not a multiple of the order"
                             " of the cyclic action")
    else:
        order = n

    for i in range(n):
        if i == 0:
            j = sum(orbit_sizes.values())
        else:
            j = sum(orb for l, orb in orbit_sizes.items()
                    if not ZZ(i).mod(n // l))
        p += j * q**i

    p = p(q**(order // n))

    return [p, order] if get_order else p


def CyclicSievingCheck(L, cyc_act, f, order=None) -> bool:
    """
    Return whether the triple ``(L, cyc_act, f)`` exhibits
    the cyclic sieving phenomenon.

    If ``cyc_act`` is None, ``L`` is expected to contain the orbit lengths.

    INPUT:

    - ``L`` -- if ``cyc_act`` is ``None``: list of orbit sizes,
      otherwise list of objects

    - ``cyc_act`` -- (default:``None``) bijective function from ``L`` to ``L``

    - ``order`` -- (default:``None``) if set to an integer, this
        cyclic order of ``cyc_act`` is used (must be an integer
        multiple of the order of ``cyc_act``) otherwise, the order of
        ``cyc_action`` is used

    EXAMPLES::

        sage: from sage.combinat.cyclic_sieving_phenomenon import *
        sage: from sage.combinat.q_analogues import q_binomial
        sage: S42 = Subsets([1,2,3,4], 2)
        sage: def cyc_act(S): return Set(i.mod(4) + 1 for i in S)
        sage: cyc_act([1,3])
        {2, 4}
        sage: cyc_act([1,4])
        {1, 2}
        sage: p = q_binomial(4,2); p
        q^4 + q^3 + 2*q^2 + q + 1
        sage: CyclicSievingPolynomial( S42, cyc_act )
        q^3 + 2*q^2 + q + 2
        sage: CyclicSievingCheck( S42, cyc_act, p )
        True
    """
    p1, n = CyclicSievingPolynomial(L, cyc_act=cyc_act, order=order,
                                    get_order=True)
    R = p1.parent()
    q = R.gen()
    return p1 == R(f).mod(q**n - 1)


def orbit_decomposition(L, cyc_act) -> list[list]:
    """
    Return the orbit decomposition of ``L`` by the action of ``cyc_act``.

    INPUT:

    - ``L`` -- list

    - ``cyc_act`` -- bijective function from ``L`` to ``L``

    OUTPUT:

    - a list of lists, the orbits under the cyc_act acting on ``L``

    EXAMPLES::

        sage: from sage.combinat.cyclic_sieving_phenomenon import *
        sage: S42 = Subsets([1,2,3,4], 2); S42
        Subsets of {1, 2, 3, 4} of size 2
        sage: def cyc_act(S): return Set(i.mod(4) + 1 for i in S)
        sage: cyc_act([1,3])
        {2, 4}
        sage: cyc_act([1,4])
        {1, 2}
        sage: orbits = orbit_decomposition(S42, cyc_act)
        sage: sorted([sorted(orb, key=sorted) for orb in orbits], key=len)
        [[{1, 3}, {2, 4}], [{1, 2}, {1, 4}, {2, 3}, {3, 4}]]
    """
    orbits = []
    L_prime = set(L)
    while L_prime:
        obj = L_prime.pop()
        orbit = [obj]
        obj = cyc_act(obj)
        while obj in L_prime:
            orbit.append(obj)
            L_prime.remove(obj)
            obj = cyc_act(obj)
        orbits.append(orbit)
    return orbits
