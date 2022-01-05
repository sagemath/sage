r"""
Linear feedback shift register (LFSR) sequence commands

Stream ciphers have been used for a long time as a source of pseudo-random
number generators.

S. Golomb [Go1967]_ gives a list of three statistical properties that a sequence of
numbers `{\bf a}=\{a_n\}_{n=1}^\infty`, `a_n\in \{0,1\}` should display to be
considered "random". Define the autocorrelation of `{\bf a}` to be

.. MATH::

     C(k)=C(k,{\bf a})=\lim_{N\rightarrow \infty} \frac{1}{N}\sum_{n=1}^N (-1)^{a_n+a_{n+k}}.

In the case where `{\bf a}` is periodic with period `P`, then this reduces to

.. MATH::

     C(k)=\frac{1}{P}\sum_{n=1}^P (-1)^{a_n+a_{n+k}}.

Assume `{\bf a}` is periodic with period `P`.

-  balance: `|\sum_{n=1}^P(-1)^{a_n}|\leq 1`.

-  low autocorrelation:

   .. MATH::

      C(k)= \left\{ \begin{array}{cc} 1,& k=0,\\ \epsilon, & k\not= 0. \end{array} \right.

   (For sequences satisfying these first two properties, it is known that
   `\epsilon=-1/P` must hold.)

-  proportional runs property: In each period, half the runs have length `1`,
   one-fourth have length `2`, etc.  Moreover, there are as many runs of `1`'s as
   there are of `0`'s.

A general feedback shift register is a map `f:{\bf F}_q^d\rightarrow {\bf
F}_q^d` of the form

.. MATH::

     \begin{array}{c} f(x_0,...,x_{n-1})=(x_1,x_2,...,x_n),\\ x_n=C(x_0,...,x_{n-1}), \end{array}


where `C:{\bf F}_q^d\rightarrow {\bf F}_q` is a given function. When `C` is of
the form

.. MATH::

     C(x_0,...,x_{n-1})=a_0x_0+...+a_{n-1}x_{n-1},

for some given constants `a_i\in {\bf F}_q`, the map is called a linear
feedback shift register (LFSR).

Example of an LFSR: Let

.. MATH::

     f(x)=a_{{0}}+a_{{1}}x+...+a_{{n}}{x}^n+...,

.. MATH::

     g(x)=b_{{0}}+b_{{1}}x+...+b_{{n}}{x}^n+...,

be given polynomials in `{\bf F}_2[x]` and let

.. MATH::

     h(x)=\frac{f(x)}{g(x)}=c_0+c_1x+...+c_nx^n+... \ .

We can compute a recursion formula which allows us to rapidly compute the
coefficients of `h(x)` (take `f(x)=1`):

.. MATH::

     c_{n}=\sum_{i=1}^n {\frac{-b_i}{b_0}c_{n-i}}.

The coefficients of `h(x)` can, under certain conditions on `f(x)` and `g(x)`,
be considered "random" from certain statistical points of view.

Example: For instance, if

.. MATH::

     f(x)=1,\ \ \ \ g(x)=x^4+x+1,

then

.. MATH::

     h(x)=1+x+x^2+x^3+x^5+x^7+x^8+...\ .

The coefficients of `h` are

.. MATH::

     \begin{array}{c} 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, \\
     1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, ...\ . \end{array}

The sequence of `0,1`'s is periodic with period `P=2^4-1=15` and satisfies
Golomb's three randomness conditions. However, this sequence of period 15 can
be "cracked" (i.e., a procedure to reproduce `g(x)`) by knowing only 8 terms!
This is the function of the Berlekamp-Massey algorithm [Mas1969]_, implemented
in ``berlekamp_massey.py``.

AUTHORS:

- David Joyner (2005-11-24): initial creation.

- Timothy Brock (2005-11): added ``lfsr_sequence`` with code modified from
  Python Cookbook, http://aspn.activestate.com/ASPN/Python/Cookbook/

- Timothy Brock (2006-04-17): added ``lfsr_autocorrelation`` and
  ``lfsr_connection_polynomial``.

"""

###########################################################################
#  Copyright (C) 2006 Timothy Brock
#  and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
###########################################################################

import copy

from sage.structure.all import Sequence
from sage.rings.all import Integer, PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField


def lfsr_sequence(key, fill, n):
    r"""
    Create an LFSR sequence.

    INPUT:

    - ``key`` -- a list of finite field elements, `[c_0, c_1,\dots, c_k]`

    - ``fill`` -- the list of the initial terms of the LFSR sequence, `[x_0,x_1,\dots,x_k]`

    - ``n`` -- number of terms of the sequence that the function returns

    OUTPUT:

    The LFSR sequence defined by `x_{n+1} = c_kx_n+...+c_0x_{n-k}` for `n \geq k`.

    EXAMPLES::

        sage: F = GF(2); l = F(1); o = F(0)
        sage: F = GF(2); S = LaurentSeriesRing(F,'x'); x = S.gen()
        sage: fill = [l,l,o,l]; key = [1,o,o,l]; n = 20
        sage: L = lfsr_sequence(key,fill,20); L
        [1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0]
        sage: from sage.matrix.berlekamp_massey import berlekamp_massey
        sage: g = berlekamp_massey(L); g
        x^4 + x^3 + 1
        sage: (1)/(g.reverse()+O(x^20))
        1 + x + x^2 + x^3 + x^5 + x^7 + x^8 + x^11 + x^15 + x^16 + x^17 + x^18 + O(x^20)
        sage: (1+x^2)/(g.reverse()+O(x^20))
        1 + x + x^4 + x^8 + x^9 + x^10 + x^11 + x^13 + x^15 + x^16 + x^19 + O(x^20)
        sage: (1+x^2+x^3)/(g.reverse()+O(x^20))
        1 + x + x^3 + x^5 + x^6 + x^9 + x^13 + x^14 + x^15 + x^16 + x^18 + O(x^20)
        sage: fill = [l,l,o,l]; key = [l,o,o,o]; n = 20
        sage: L = lfsr_sequence(key,fill,20); L
        [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1]
        sage: g = berlekamp_massey(L); g
        x^4 + 1
        sage: (1+x)/(g.reverse()+O(x^20))
        1 + x + x^4 + x^5 + x^8 + x^9 + x^12 + x^13 + x^16 + x^17 + O(x^20)
        sage: (1+x+x^3)/(g.reverse()+O(x^20))
        1 + x + x^3 + x^4 + x^5 + x^7 + x^8 + x^9 + x^11 + x^12 + x^13 + x^15 + x^16 + x^17 + x^19 + O(x^20)

    """
    if not isinstance(key, list):
        raise TypeError("key must be a list")
    key = Sequence(key)
    F = key.universe()
    if not is_FiniteField(F):
        raise TypeError("universe of sequence must be a finite field")

    s = fill
    k = len(fill)
    L = []
    for i in range(n):
        s0 = copy.copy(s)
        L.append(s[0])
        s = s[1:k]
        s.append(sum([key[i] * s0[i] for i in range(k)]))
    return L


def lfsr_autocorrelation(L, p, k):
    """
    INPUT:

    - ``L`` -- a periodic sequence of elements of ZZ or GF(2); must have length `p`

    - ``p`` -- the period of `L`

    - ``k`` -- an integer between `0` and `p`

    OUTPUT: autocorrelation sequence of `L`

    EXAMPLES::

        sage: F = GF(2)
        sage: o = F(0)
        sage: l = F(1)
        sage: key = [l,o,o,l]; fill = [l,l,o,l]; n = 20
        sage: s = lfsr_sequence(key,fill,n)
        sage: lfsr_autocorrelation(s,15,7)
        4/15
        sage: lfsr_autocorrelation(s,int(15),7)
        4/15
    """
    if not isinstance(L, list):
        raise TypeError("L (=%s) must be a list" % L)
    p = Integer(p)
    _p = int(p)
    k = int(k)
    L0 = L[:_p]     # slices makes a copy
    L0 = L0 + L0[:k]
    return sum([int(L0[i]) * int(L0[i + k]) / p for i in range(_p)])


def lfsr_connection_polynomial(s):
    """
    INPUT:

    - ``s`` -- a sequence of elements of a finite field of even length

    OUTPUT:

    - ``C(x)`` -- the connection polynomial of the minimal LFSR.

    This implements the algorithm in section 3 of J. L. Massey's article
    [Mas1969]_.

    EXAMPLES::

        sage: F = GF(2)
        sage: F
        Finite Field of size 2
        sage: o = F(0); l = F(1)
        sage: key = [l,o,o,l]; fill = [l,l,o,l]; n = 20
        sage: s = lfsr_sequence(key,fill,n); s
        [1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0]
        sage: lfsr_connection_polynomial(s)
        x^4 + x + 1
        sage: from sage.matrix.berlekamp_massey import berlekamp_massey
        sage: berlekamp_massey(s)
        x^4 + x^3 + 1

    Notice that ``berlekamp_massey`` returns the reverse of the connection
    polynomial (and is potentially must faster than this implementation).
    """
    # Initialization:
    FF = s[0].base_ring()
    R = PolynomialRing(FF, "x")
    x = R.gen()
    C = R.one()
    B = R.one()
    m = 1
    b = FF.one()
    L = 0
    N = 0

    while N < len(s):
        if L > 0:
            r = min(L + 1, C.degree() + 1)
            d = s[N] + sum([(C.list())[i] * s[N - i] for i in range(1, r)])
        if L == 0:
            d = s[N]
        if d == 0:
            m += 1
            N += 1
        if d > 0:
            if 2 * L > N:
                C = C - d*b**(-1)*x**m*B
                m += 1
                N += 1
            else:
                T = C
                C = C - d*b**(-1)*x**m*B
                L = N + 1 - L
                m = 1
                b = d
                B = T
                N += 1
    return C
