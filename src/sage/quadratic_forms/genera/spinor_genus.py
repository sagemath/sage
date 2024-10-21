# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.groups
r"""
Spinor genus computations.

This file defines the group of spinor operators used
for the computation of spinor genera.
It is meant for internal use only.

EXAMPLES::

    sage: from sage.quadratic_forms.genera.spinor_genus import SpinorOperators
    sage: A = SpinorOperators((2, 3, 7))

AUTHORS:

- Simon Brandhorst (2020-11-1): initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.groups.abelian_gps.abelian_group_gap import (AbelianGroupGap,
                                                       AbelianGroupElement_gap)
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class SpinorOperator(AbelianGroupElement_gap):
    r"""
    A spinor operator seen as a tuple of square classes.

    For `2` the square class is represented as one of `1,3,5,7` and for
    `p` odd it is `1` for a `p`-adic unit square and `-1` for a non-square.

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.spinor_genus import *
        sage: A = SpinorOperators((2, 3, 7))
        sage: A.an_element()
        [2:7, 3:-1, 7:-1]
    """

    def _repr_(self) -> str:
        r"""
        Return the print representation.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.spinor_genus import *
            sage: A = SpinorOperators((2, 3, 7))
            sage: A.an_element()
            [2:7, 3:-1, 7:-1]
        """
        e = self.exponents()
        p = self.parent()._primes
        s = "[2:"
        if e[0] == 0 == e[1]:
            s += "1"
        elif e[0] == 1 and e[1] == 0:
            s += "3"
        elif e[0] == 0 and e[1] == 1:
            s += "5"
        elif e[0] == 1 == e[1]:
            s += "7"
        for k in range(1, len(p)):
            s += f", {p[k]}:{(-1)**e[k + 1]}"
        s += "]"
        return s


class SpinorOperators(AbelianGroupGap):
    r"""
    The group of spinor operators of a genus.

    It is a product of `p`-adic unit square classes
    used for spinor genus computations.

    INPUT:

    - a tuple of primes `(p_1=2,\dots, p_n`)

    EXAMPLES::

        sage: from sage.quadratic_forms.genera.spinor_genus import *
        sage: SpinorOperators((2, 3, 7))
        Group of SpinorOperators at primes (2, 3, 7)
    """
    def __init__(self, primes):
        r"""
        Initialize the group of spinor operators.

        TESTS::

            sage: from sage.quadratic_forms.genera.spinor_genus import *
            sage: S = SpinorOperators((2, 3, 7))
            sage: TestSuite(S).run()
        """
        if primes[0] != 2:
            raise ValueError("first prime must be 2")
        self._primes = tuple(ZZ(p) for p in primes)
        orders = len(self._primes) * [2] + [2]
        # 3, 5, unit_p1, unit_p2,...
        orders = tuple(orders)
        AbelianGroupGap.__init__(self, orders)

    def __reduce__(self):
        r"""
        Implement pickling.

        OUTPUT:

        a tuple ``f`` such that this element is ``f[0](*f[1])``

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.spinor_genus import SpinorOperators
            sage: S = SpinorOperators((2, 3, 7))
            sage: S == loads(dumps(S))
            True
        """
        return SpinorOperators, (self._primes,)

    Element = SpinorOperator

    def _repr_(self) -> str:
        r"""
        Return the print representation of ``self``.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.spinor_genus import SpinorOperators
            sage: SpinorOperators((2, 3, 7))
            Group of SpinorOperators at primes (2, 3, 7)
        """
        return f"Group of SpinorOperators at primes {self._primes}"

    def to_square_class(self, x, p):
        r"""
        Return `(1,...,1,x,1,...,1)` with the square class of `x` at position `p`.

        INPUT:

        - ``p`` -- a prime

        - ``x`` -- nonzero rational number

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.spinor_genus import SpinorOperators
            sage: AS = SpinorOperators((2, 3, 7))
            sage: AS.to_square_class(5, 7)
            [2:1, 3:1, 7:-1]
            sage: AS.to_square_class(5, 2)
            [2:5, 3:1, 7:1]
            sage: AS.to_square_class(-5, 2)
            [2:3, 3:1, 7:1]
            sage: AS.to_square_class(7, 2)
            [2:7, 3:1, 7:1]
        """
        x = QQ(x)
        if x == 0:
            raise ValueError("x must be nonzero")
        if p not in self._primes:
            raise ValueError("not a coordinate prime")
        v, u = x.val_unit(p)
        v = v % 2
        if v != 0:
            raise ValueError(f"x(={x}) must be a p-adic unit")
        y = self.one()
        if p == 2:
            u = u % 8
            if u == 3:
                y *= self.gens()[0]
            if u == 5:
                y *= self.gens()[1]
            if u == 7:
                y *= self.gens()[0] * self.gens()[1]
            return y
        i = 1 + self._primes.index(p)
        if not u.is_padic_square(p):
            y *= self.gens()[i]
        return y

    def delta(self, r, prime=None):
        r"""
        Diagonal embedding of rational square classes.

        INPUT:

        - ``r`` -- a nonzero integer; if ``prime`` is ``None``, ``r`` must not
          be divisible by the defining primes of ``self``

        - ``prime`` -- (default: ``None``) a prime or `-1`

        OUTPUT:

        If a prime `p` is given, the method returns
        `\Delta_p(r)`
        otherwise returns `\Delta(r)`
        where both are as defined by Conway-Sloane in
        Chapter 15 9.3 of [CS1988]_.

        EXAMPLES::

            sage: from sage.quadratic_forms.genera.spinor_genus import SpinorOperators
            sage: AS = SpinorOperators((2, 3, 7))
            sage: AS.delta(5)
            [2:5, 3:-1, 7:-1]
            sage: AS.delta(2, prime=3)
            [2:1, 3:-1, 7:1]
            sage: AS.delta(11)
            [2:3, 3:-1, 7:1]
            sage: AS.delta(3, prime=7)
            [2:1, 3:1, 7:-1]
        """
        r = ZZ(r)
        if prime is None:
            if any(p.divides(r) for p in self._primes):
                raise ValueError(f"r must not be divisible by {self._primes}")
            return self.prod([self.to_square_class(r, p)
                              for p in self._primes])
        prime = ZZ(prime)
        if prime == -1:
            r = r.sign()
            return self.prod([self.to_square_class(r, p)
                              for p in self._primes])
        if prime not in self._primes:
            raise ValueError("prime must be among %s" % self._primes)
        v, u = r.val_unit(prime)
        pv = prime**v
        y = self.prod([self.to_square_class(pv, q)
                       for q in self._primes if q != prime])
        if prime in self._primes:
            y *= self.to_square_class(u, p=prime)
        return y
