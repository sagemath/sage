r"""
Miscellaneous Functions

This file contains several miscellaneous functions used by `p`-adics.

- ``gauss_sum`` -- compute Gauss sums using the Gross-Koblitz formula.
- ``min`` -- a version of ``min`` that returns `\infty` on empty input.
- ``max`` -- a version of ``max`` that returns `-\infty` on empty input.

AUTHORS:

- David Roe
- Adriana Salerno
- Ander Steele
- Kiran Kedlaya (modified gauss_sum 2017/09)
"""
# ****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import infinity

python_min = min
python_max = max


def gauss_sum(a, p, f, prec=20, factored=False, algorithm='pari', parent=None):
    r"""
    Return the Gauss sum `g_q(a)` as a `p`-adic number.

    The Gauss sum `g_q(a)` is defined by

    .. MATH::

        g_q(a)= \sum_{u\in F_q^*} \omega(u)^{-a} \zeta_q^u,

    where `q = p^f`, `\omega` is the Teichmüller character and
    `\zeta_q` is some arbitrary choice of primitive `q`-th root of
    unity. The computation is adapted from the main theorem in Alain
    Robert's paper *The Gross-Koblitz formula revisited*,
    Rend. Sem. Mat. Univ. Padova 105 (2001), 157--170.

    Let `p` be a prime, `f` a positive integer, `q=p^f`, and `\pi` be
    the unique root of `f(x) = x^{p-1}+p` congruent to `\zeta_p - 1` modulo
    `(\zeta_p - 1)^2`. Let `0\leq a < q-1`. Then the
    Gross-Koblitz formula gives us the value of the Gauss sum `g_q(a)`
    as a product of `p`-adic Gamma functions as follows:

    .. MATH::

        g_q(a) = -\pi^s \prod_{0\leq i < f} \Gamma_p(a^{(i)}/(q-1)),

    where `s` is the sum of the digits of `a` in base `p` and the
    `a^{(i)}` have `p`-adic expansions obtained from cyclic
    permutations of that of `a`.

    INPUT:

    - ``a`` -- integer

    - ``p`` -- prime

    - ``f`` -- positive integer

    - ``prec`` -- positive integer (default: 20)

    - ``factored`` -- boolean (default: ``False``)

    - ``algorithm`` -- flag passed to `p`-adic Gamma function (default: ``'pari'``)

    OUTPUT:

    If ``factored`` is ``False``, returns a `p`-adic number in an Eisenstein extension of `\QQ_p`.
    This number has the form `pi^e * z` where `pi` is as above, `e` is some nonnegative
    integer, and `z` is an element of `\ZZ_p`; if ``factored`` is ``True``, the pair `(e,z)`
    is returned instead, and the Eisenstein extension is not formed.

    .. NOTE::

        This is based on GP code written by Adriana Salerno.

    EXAMPLES:

    In this example, we verify that `g_3(0) = -1`::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: -gauss_sum(0, 3, 1)                                                       # needs sage.libs.ntl sage.rings.padics
        1 + O(pi^40)

    Next, we verify that `g_5(a) g_5(-a) = 5 (-1)^a`::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: gauss_sum(2,5,1)^2 - 5                                                    # needs sage.libs.ntl
        O(pi^84)
        sage: gauss_sum(1,5,1)*gauss_sum(3,5,1) + 5                                     # needs sage.libs.ntl
        O(pi^84)

    Finally, we compute a non-trivial value::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: gauss_sum(2,13,2)                                                         # needs sage.libs.ntl
        6*pi^2 + 7*pi^14 + 11*pi^26 + 3*pi^62 + 6*pi^74 + 3*pi^86 + 5*pi^98 +
        pi^110 + 7*pi^134 + 9*pi^146 + 4*pi^158 + 6*pi^170 + 4*pi^194 +
        pi^206 + 6*pi^218 + 9*pi^230 + O(pi^242)
        sage: gauss_sum(2,13,2, prec=5, factored=True)                                  # needs sage.rings.padics
        (2, 6 + 6*13 + 10*13^2 + O(13^5))

    .. SEEALSO::

        - :func:`sage.arith.misc.gauss_sum` for general finite fields
        - :meth:`sage.modular.dirichlet.DirichletCharacter.gauss_sum`
          for prime finite fields
        - :meth:`sage.modular.dirichlet.DirichletCharacter.gauss_sum_numerical`
          for prime finite fields
    """
    from sage.rings.padics.factory import Zp
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    q = p**f
    a = a % (q-1)
    if parent is None:
        R = Zp(p, prec)
    else:
        R = parent
    out = -R.one()
    if a != 0:
        t = R(1/(q-1))
        for i in range(f):
            out *= (a*t).gamma(algorithm)
            a = (a*p) % (q-1)
    s = sum(a.digits(base=p))
    if factored:
        return s, out
    X = PolynomialRing(R, name='X').gen()
    pi = R.ext(X**(p - 1) + p, names='pi').gen()
    out *= pi**s
    return out


def min(*L):
    r"""
    Return the minimum of the inputs, where the minimum of the empty
    list is `\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import min
        sage: min()
        +Infinity
        sage: min(2,3)
        2
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity


def max(*L):
    r"""
    Return the maximum of the inputs, where the maximum of the empty
    list is `-\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import max
        sage: max()
        -Infinity
        sage: max(2,3)
        3
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity


def precprint(prec_type, prec_cap, p):
    """
    String describing the precision mode on a `p`-adic ring or field.

    EXAMPLES::

        sage: from sage.rings.padics.misc import precprint
        sage: precprint('capped-rel', 12, 2)
        'with capped relative precision 12'
        sage: precprint('capped-abs', 11, 3)
        'with capped absolute precision 11'
        sage: precprint('floating-point', 1234, 5)
        'with floating precision 1234'
        sage: precprint('fixed-mod', 1, 17)
        'of fixed modulus 17^1'
    """
    precD = {'capped-rel':'with capped relative precision %s' % prec_cap,
             'capped-abs':'with capped absolute precision %s' % prec_cap,
             'floating-point':'with floating precision %s' % prec_cap,
             'fixed-mod':'of fixed modulus %s^%s' % (p, prec_cap),
             'lattice-cap':'with lattice-cap precision',
             'lattice-float':'with lattice-float precision',
             'relaxed':'handled with relaxed arithmetics'}
    return precD[prec_type]


def trim_zeros(L):
    r"""
    Strip trailing zeros/empty lists from a list.

    EXAMPLES::

        sage: from sage.rings.padics.misc import trim_zeros
        sage: trim_zeros([1,0,1,0])
        [1, 0, 1]
        sage: trim_zeros([[1],[],[2],[],[]])
        [[1], [], [2]]
        sage: trim_zeros([[],[]])
        []
        sage: trim_zeros([])
        []

    Zeros are also trimmed from nested lists (one deep)::

        sage: trim_zeros([[1,0]])
        [[1]]
        sage: trim_zeros([[0],[1]])
        [[], [1]]
    """
    strip_trailing = True
    n = len(L)
    for i, c in zip(reversed(range(len(L))), reversed(L)):
        if strip_trailing and (c == 0 or c == []):
            n = i
        elif isinstance(c, list):
            strip_trailing = False
            m = len(c)
            # strip trailing zeros from the sublists
            for j, d in zip(reversed(range(len(c))), reversed(c)):
                if d == 0:
                    m = j
                else:
                    break
            L[i] = c[:m]
        else:
            break
    return L[:n]
