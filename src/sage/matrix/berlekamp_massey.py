"""
Minimal Polynomials of Linear Recurrence Sequences

AUTHORS:

- William Stein
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.rational_field


def berlekamp_massey(a):
    r"""
    Use the Berlekamp-Massey algorithm to find the minimal polynomial
    of a linear recurrence sequence `a`.

    The minimal polynomial of a linear recurrence `\{a_r\}` is
    by definition the unique monic polynomial `g`, such that if
    `\{a_r\}` satisfies a linear recurrence
    `a_{j+k} + b_{j-1} a_{j-1+k} + \cdots + b_0 a_k=0`
    (for all `k\geq 0`), then `g` divides the
    polynomial `x^j + \sum_{i=0}^{j-1} b_i x^i`.

    INPUT:

    - ``a`` -- list of even length of elements of a field (or domain)

    OUTPUT:

    the minimal polynomial of the sequence, as a polynomial over the
    field in which the entries of `a` live

    .. WARNING::

         The result is only guaranteed to be correct on the full
         sequence if there exists a linear recurrence of length less
         than half the length of `a`.

    EXAMPLES::

        sage: from sage.matrix.berlekamp_massey import berlekamp_massey
        sage: berlekamp_massey([1,2,1,2,1,2])
        x^2 - 1
        sage: berlekamp_massey([GF(7)(1), 19, 1, 19])
        x^2 + 6
        sage: berlekamp_massey([2,2,1,2,1,191,393,132])
        x^4 - 36727/11711*x^3 + 34213/5019*x^2 + 7024942/35133*x - 335813/1673
        sage: berlekamp_massey(prime_range(2, 38))                                      # needs sage.libs.pari
        x^6 - 14/9*x^5 - 7/9*x^4 + 157/54*x^3 - 25/27*x^2 - 73/18*x + 37/9

    TESTS::

        sage: berlekamp_massey("banana")
        Traceback (most recent call last):
        ...
        TypeError: argument must be a list or tuple
        sage: berlekamp_massey([1,2,5])
        Traceback (most recent call last):
        ...
        ValueError: argument must have an even number of terms
    """
    if not isinstance(a, (list, tuple)):
        raise TypeError("argument must be a list or tuple")
    if len(a) % 2:
        raise ValueError("argument must have an even number of terms")

    M = len(a) // 2

    try:
        K = a[0].parent().fraction_field()
    except AttributeError:
        K = sage.rings.rational_field.RationalField()

    R, x = K['x'].objgen()
    f0, f1 = R(a), x**(2 * M)
    s0, s1 = 1, 0
    while f1.degree() >= M:
        f0, (q, f1) = f1, f0.quo_rem(f1)
        s0, s1 = s1, s0 - q * s1
    return s1.reverse().monic()
