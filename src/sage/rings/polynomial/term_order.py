r"""
Term orders

Sage supports the following term orders:

Lexicographic (lex)
    `x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
    This term order is called 'lp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
        sage: x > y
        True
        sage: x > y^2
        True
        sage: x > 1
        True
        sage: x^1*y^2 > y^3*z^4
        True
        sage: x^3*y^2*z^4 < x^3*y^2*z^1
        False

Degree reverse lexicographic (degrevlex)
    Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
    `x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
    there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
    This term order is called 'dp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='degrevlex')
        sage: x > y
        True
        sage: x > y^2*z
        False
        sage: x > 1
        True
        sage: x^1*y^5*z^2 > x^4*y^1*z^3
        True
        sage: x^2*y*z^2 > x*y^3*z
        False

Degree lexicographic (deglex)
    Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
    `x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
    there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
    This term order is called 'Dp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='deglex')
        sage: x > y
        True
        sage: x > y^2*z
        False
        sage: x > 1
        True
        sage: x^1*y^2*z^3 > x^3*y^2*z^0
        True
        sage: x^2*y*z^2 > x*y^3*z
        True

Inverse lexicographic (invlex)
    `x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i < b_i`.
    This order is called 'rp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='invlex')
        sage: x > y
        False
        sage: y > x^2
        True
        sage: x > 1
        True
        sage: x*y > z
        False

    This term order only makes sense in a non-commutative setting
    because if P is the ring `k[x_1, \dots, x_n]` and term
    order 'invlex' then it is equivalent to the ring
    `k[x_n, \dots, x_1]` with term order 'lex'.

Negative lexicographic (neglex)
    `x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i > b_i`.
    This term order is called 'ls' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='neglex')
        sage: x > y
        False
        sage: x > 1
        False
        sage: x^1*y^2 > y^3*z^4
        False
        sage: x^3*y^2*z^4 < x^3*y^2*z^1
        True
        sage: x*y > z
        False

Negative degree reverse lexicographic (negdegrevlex)
    Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
    `x^a < x^b` if and only if `\deg(x^a) > \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
    there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
    This term order is called 'ds' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='negdegrevlex')
        sage: x > y
        True
        sage: x > x^2
        True
        sage: x > 1
        False
        sage: x^1*y^2 > y^3*z^4
        True
        sage: x^2*y*z^2 > x*y^3*z
        False

Negative degree lexicographic (negdeglex)
    Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
    `x^a < x^b` if and only if `\deg(x^a) > \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
    there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
    This term order is called 'Ds' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='negdeglex')
        sage: x > y
        True
        sage: x > x^2
        True
        sage: x > 1
        False
        sage: x^1*y^2 > y^3*z^4
        True
        sage: x^2*y*z^2 > x*y^3*z
        True

Weighted degree reverse lexicographic (wdegrevlex), positive integral weights
    Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
    `x^a < x^b` if and only if `\deg_w(x^a) < \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
    there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
    This term order is called 'wp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order=TermOrder('wdegrevlex',(1,2,3)))
        sage: x > y
        False
        sage: x > x^2
        False
        sage: x > 1
        True
        sage: x^1*y^2 > x^2*z
        True
        sage: y*z > x^3*y
        False

Weighted degree lexicographic (wdeglex), positive integral weights
    Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
    `x^a < x^b` if and only if `\deg_w(x^a) < \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
    there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
    This term order is called 'Wp' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order=TermOrder('wdeglex',(1,2,3)))
        sage: x > y
        False
        sage: x > x^2
        False
        sage: x > 1
        True
        sage: x^1*y^2 > x^2*z
        False
        sage: y*z > x^3*y
        False

Negative weighted degree reverse lexicographic (negwdegrevlex), positive integral weights
    Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
    `x^a < x^b` if and only if `\deg_w(x^a) > \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
    there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
    This term order is called 'ws' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order=TermOrder('negwdegrevlex',(1,2,3)))
        sage: x > y
        True
        sage: x > x^2
        True
        sage: x > 1
        False
        sage: x^1*y^2 > x^2*z
        True
        sage: y*z > x^3*y
        False

Degree negative lexicographic (degneglex)
    Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
    `x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
    there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i > b_i`.
    This term order is called 'dp_asc' in PolyBoRi.
    Singular has the extra weight vector ordering ``(a(1:n),ls)`` for this
    purpose.

    EXAMPLES:

    ::

        sage: t = TermOrder('degneglex')
        sage: P.<x,y,z> = PolynomialRing(QQ, order=t)
        sage: x*y > y*z # indirect doctest
        False
        sage: x*y > x
        True

Negative weighted degree lexicographic (negwdeglex), positive integral weights
    Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
    `x^a < x^b` if and only if `\deg_w(x^a) > \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
    there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
    This term order is called 'Ws' in Singular.

    EXAMPLES:

    ::

        sage: P.<x,y,z> = PolynomialRing(QQ, 3, order=TermOrder('negwdeglex',(1,2,3)))
        sage: x > y
        True
        sage: x > x^2
        True
        sage: x > 1
        False
        sage: x^1*y^2 > x^2*z
        False
        sage: y*z > x^3*y
        False

Of these, only 'degrevlex', 'deglex', 'degneglex', 'wdegrevlex', 'wdeglex',
'invlex' and 'lex' are global orders.

Sage also supports matrix term order. Given a square matrix `A`,

    `x^a <_A x^b` if and only if `Aa < Ab`

where `<` is the lexicographic term order.

EXAMPLES::

    sage: # needs sage.modules
    sage: m = matrix(2, [2,3,0,1]); m
    [2 3]
    [0 1]
    sage: T = TermOrder(m); T
    Matrix term order with matrix
    [2 3]
    [0 1]
    sage: P.<a,b> = PolynomialRing(QQ, 2, order=T)
    sage: P
    Multivariate Polynomial Ring in a, b over Rational Field
    sage: a > b
    False
    sage: a^3 < b^2
    True
    sage: S = TermOrder('M(2,3,0,1)')
    sage: T == S
    True

Additionally all these monomial orders may be combined to product or block
orders, defined as:

Let `x = (x_1, x_2, \dots, x_n)` and `y = (y_1, y_2, \dots, y_m)` be two
ordered sets of variables, `<_1` a monomial order on `k[x]` and `<_2` a
monomial order on `k[y]`.

The product order (or block order) `<` `:=` `(<_1,<_2)` on `k[x,y]` is defined
as: `x^a y^b < x^A y^B` if and only if `x^a <_1 x^A` or (`x^a =x^A` and `y^b
<_2 y^B`).

These block orders are constructed in Sage by giving a comma separated list of
monomial orders with the length of each block attached to them.

EXAMPLES:

As an example, consider constructing a block order where the first four
variables are compared using the degree reverse lexicographical order while the
last two variables in the second block are compared using negative
lexicographical order.

::

    sage: P.<a,b,c,d,e,f> = PolynomialRing(QQ, 6, order='degrevlex(4),neglex(2)')
    sage: a > c^4
    False
    sage: a > e^4
    True
    sage: e > f^2
    False

The same result can be achieved by::

    sage: T1 = TermOrder('degrevlex',4)
    sage: T2 = TermOrder('neglex',2)
    sage: T = T1 + T2
    sage: P.<a,b,c,d,e,f> = PolynomialRing(QQ, 6, order=T)
    sage: a > c^4
    False
    sage: a > e^4
    True

If any other unsupported term order is given the provided string can be forced
to be passed through as is to Singular, Macaulay2, and Magma.  This ensures
that it is for example possible to calculate a Groebner basis with respect to
some term order Singular supports but Sage doesn't::

    sage: T = TermOrder("royalorder")
    Traceback (most recent call last):
    ...
    ValueError: unknown term order 'royalorder'
    sage: T = TermOrder("royalorder", force=True)
    sage: T
    royalorder term order
    sage: T.singular_str()
    'royalorder'

AUTHORS:

- David Joyner and William Stein: initial version of multi_polynomial_ring

- Kiran S. Kedlaya: added macaulay2 interface

- Martin Albrecht: implemented native term orders, refactoring

- Kwankyu Lee: implemented matrix and weighted degree term orders

- Simon King (2011-06-06): added termorder_from_singular
"""
#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import re
from sage.structure.sage_object import SageObject

print_name_mapping = {
    'lex'           : 'Lexicographic',
    'invlex'        : 'Inverse lexicographic',
    'degrevlex'     : 'Degree reverse lexicographic',
    'deglex'        : 'Degree lexicographic',
    'neglex'        : 'Negative lexicographic',
    'negdegrevlex'  : 'Negative degree reverse lexicographic',
    'negdeglex'     : 'Negative degree lexicographic',
    'degneglex'     : 'Degree negative lexicographic',
    'wdegrevlex'    : 'Weighted degree reverse lexicographic',
    'wdeglex'       : 'Weighted degree lexicographic',
    'negwdegrevlex' : 'Negative weighted degree reverse lexicographic',
    'negwdeglex'    : 'Negative weighted degree lexicographic',
}

singular_name_mapping = {
    'lex'           : 'lp',
    'invlex'        : 'ip',
    'degrevlex'     : 'dp',
    'deglex'        : 'Dp',
    'neglex'        : 'ls',
    'negdegrevlex'  : 'ds',
    'negdeglex'     : 'Ds',
    'degneglex'     : '(a(1:%(ngens)i),ls(%(ngens)i))',
    'wdegrevlex'    : 'wp',
    'wdeglex'       : 'Wp',
    'negwdegrevlex' : 'ws',
    'negwdeglex'    : 'Ws',
}

inv_singular_name_mapping = dict(zip(singular_name_mapping.values(),singular_name_mapping.keys()))

macaulay2_name_mapping = {
    'lex'           : 'Lex',
    'revlex'        : 'RevLex, Global=>false',
    'degrevlex'     : 'GRevLex',
    'deglex'        : 'GLex',
}

inv_macaulay2_name_mapping = dict(zip(macaulay2_name_mapping.values(),macaulay2_name_mapping.keys()))

magma_name_mapping = {
    'lex'           : '"lex"',
    'degrevlex'     : '"grevlex"',
    'deglex'        : '"glex"',
}

inv_magma_name_mapping = dict(zip(magma_name_mapping.values(),magma_name_mapping.keys()))

lex_description = r"""
Lexicographic (lex) term order.

`x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
"""

invlex_description = r"""
Inverse lexicographic (invlex) term order.

`x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i < b_i`.
"""

degrevlex_description = r"""
Degree reverse lexicographic (degrevlex) term order.

Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
`x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
"""

deglex_description = r"""
Degree lexicographic (deglex) term order.

Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
`x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
"""

neglex_description = r"""
Negative lexicographic (neglex) term order.

`x^a < x^b` if and only if there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i > b_i`.
"""

negdegrevlex_description = r"""
Negative degree reverse lexicographic (negdegrevlex) term order.

Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
`x^a < x^b` if and only if `\deg(x^a) > \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
"""

negdeglex_description = r"""
Negative degree lexicographic (negdeglex) term order.

Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
`x^a < x^b` if and only if `\deg(x^a) > \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
"""

degneglex_description = r"""
Degree negative lexicographic (degneglex) term order.

Let `\deg(x^a) = a_1 + a_2 + \dots + a_n`, then
`x^a < x^b` if and only if `\deg(x^a) < \deg(x^b)` or `\deg(x^a) = \deg(x^b)` and
there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i > b_i`.
"""

wdegrevlex_description = r"""
Weighted degree reverse lexicographic (wdegrevlex) term order.

Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
`x^a < x^b` if and only if `\deg_w(x^a) < \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
"""

wdeglex_description = r"""
Weighted degree lexicographic (wdeglex) term order.

Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
`x^a < x^b` if and only if `\deg_w(x^a) < \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
"""

negwdegrevlex_description = r"""
Negative weighted degree reverse lexicographic (negwdegrevlex) term order.

Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
`x^a < x^b` if and only if `\deg_w(x^a) > \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
there exists `1 \le i \le n` such that `a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i > b_i`.
"""

negwdeglex_description = r"""
Negative weighted degree lexicographic (negwdeglex) term order.

Let `\deg_w(x^a) = a_1w_1 + a_2w_2 + \dots + a_nw_n` with weights `w`, then
`x^a < x^b` if and only if `\deg_w(x^a) > \deg_w(x^b)` or `\deg_w(x^a) = \deg_w(x^b)` and
there exists `1 \le i \le n` such that `a_1 = b_1, \dots, a_{i-1} = b_{i-1}, a_i < b_i`.
"""

matrix_description = """
Matrix term order defined by a matrix A.

`x^a < x^b` if and only if `x^{Aa} < x^{Ab}` where `<` is the lexicographic term order.
"""

block_description = r"""
Block term order defined by term orders `<_1, <_2, \dots, <_n`.

`x^a < x^b` if and only if `a = b` with respect to the first `n-1` term orders and `a <_n b`
with respect to the `n`-th term order `<_n`.
"""

description_mapping = {
    'lex'           : lex_description,
    'invlex'        : invlex_description,
    'degrevlex'     : degrevlex_description,
    'deglex'        : deglex_description,
    'neglex'        : neglex_description,
    'negdegrevlex'  : negdegrevlex_description,
    'negdeglex'     : negdeglex_description,
    'degneglex'     : degneglex_description,
    'wdeglex'       : wdeglex_description,
    'wdegrevlex'    : wdegrevlex_description,
    'negwdegrevlex' : negwdegrevlex_description,
    'negwdeglex'    : negwdeglex_description,
    'matrix'        : matrix_description,
    'block'         : block_description,
}

class TermOrder(SageObject):
    """
    A term order.

    See ``sage.rings.polynomial.term_order`` for details on supported
    term orders.
    """
    def __setstate__(self, dict):
        """
        Translate old pickled TermOrder objects.

        See Issue :issue:`11316`.

        EXAMPLES::

            sage: t = TermOrder('lex')
            sage: t2 = loads(dumps(t))
            sage: t2._weights is None
            True
        """
        if '_weights' not in dict:
            name = dict['_TermOrder__name']
            n = dict['_TermOrder__length']
            t = TermOrder(name,n)
            self.__dict__.update(t.__dict__)
        else:
            self.__dict__.update(dict)

    def __init__(self, name='lex', n=0, force=False):
        """
        Construct a new term order object.

        INPUT:

        - ``name`` -- name of the term order (default: ``'lex'``)

        - ``n`` -- number of variables (default: `0`) weights for
          weighted degree orders. The weights are converted to
          integers and must be positive.

        - ``force`` -- ignore unknown term orders

        See the ``sage.rings.polynomial.term_order`` module
        for help which names and orders are available.

        EXAMPLES::

            sage: t = TermOrder('lex')
            sage: t
            Lexicographic term order
            sage: loads(dumps(t)) == t
            True

        We can construct block orders directly as::

            sage: TermOrder('degrevlex(3),neglex(2)')
            Block term order with blocks:
            (Degree reverse lexicographic term order of length 3,
             Negative lexicographic term order of length 2)

        or by adding together the blocks::

            sage: t1 = TermOrder('degrevlex',3)
            sage: t2 = TermOrder('neglex',2)
            sage: t1 + t2
            Block term order with blocks:
            (Degree reverse lexicographic term order of length 3,
             Negative lexicographic term order of length 2)
            sage: t3 = TermOrder('wdeglex',(1,2))
            sage: t1 + t3
            Block term order with blocks:
            (Degree reverse lexicographic term order of length 3,
             Weighted degree lexicographic term order with weights (1, 2))
            sage: t4 = TermOrder('degneglex')
            sage: t4
            Degree negative lexicographic term order

        We allow blocks of length 0, these are simply ignored::

            sage: TermOrder('lex(0),degrevlex(5),deglex(0),deglex(2)')
            Block term order with blocks:
            (Degree reverse lexicographic term order of length 5,
             Degree lexicographic term order of length 2)

        .. NOTE::

           The optional `n` parameter is not necessary if only
           non-block orders like `deglex` are
           constructed. However, it is useful if block orders are
           to be constructed from this ``TermOrder`` object later.

        TESTS:

        We demonstrate that nonpositive weights are refused and non-integral weights
        are converted to integers (and potentially rounded)::

            sage: N.<a,b,c> = PolynomialRing(QQ, 3, order=TermOrder('wdeglex',[-1,2,-3]))
            Traceback (most recent call last):
            ...
            ValueError: the degree weights must be positive integers
            sage: N.<a,b,c> = PolynomialRing(QQ, 3, order=TermOrder('wdeglex',[1.1,2,3]))
            sage: a.degree()
            1

        We enforce consistency when calling the copy constructor (cf.
        :issue:`12748`)::

            sage: T = TermOrder('degrevlex', 6) + TermOrder('degrevlex',10)
            sage: R.<x0,y0,z0,x1,y1,z1,a0,a1,a2,a3,a4,a5,a6,a7,a8> = PolynomialRing(QQ, order=T)
            Traceback (most recent call last):
            ...
            ValueError: the length of the given term order (16) differs from the number of variables (15)

        If ``_singular_ringorder_column`` attribute is set, then it is regarded
        as a different term order::

            sage: T = TermOrder('degrevlex')
            sage: R.<x,y,z> = PolynomialRing(QQ, order=T)
            sage: R._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            // coefficients: QQ
            // number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C
            sage: T2 = copy(T)
            sage: T2 == T
            True
            sage: T2._singular_ringorder_column = 0
            sage: T2 == T
            False
            sage: S = R.change_ring(order=T2)
            sage: S == T
            False
            sage: S._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            // coefficients: QQ
            // number of vars : 3
            //        block   1 : ordering C
            //        block   2 : ordering dp
            //                  : names    x y z

        Check that :issue:`29635` is fixed::

            sage: T = PolynomialRing(GF(101^5), 'u,v,w',                                # needs sage.rings.finite_rings
            ....:                    order=TermOrder('degneglex')).term_order()
            sage: T.singular_str()                                                      # needs sage.rings.finite_rings
            '(a(1:3),ls(3))'
            sage: (T + T).singular_str()                                                # needs sage.rings.finite_rings
            '(a(1:3),ls(3),a(1:3),ls(3))'
        """
        if isinstance(name, TermOrder):
            self.__copy(name)
            if n:
                if not name.is_block_order() and not name.is_weighted_degree_order():
                    self._length = n
                    if self._length != 0:
                        self._singular_str = (self._singular_str
                                              % dict(ngens=self._length))
                elif self._length != n:
                    raise ValueError("the length of the given term order ({}) differs from the number of variables ({})"
                            .format(self._length, n))
            return

        if isinstance(name, str):
            name = name.lower()
        else:
            try:
                if not isinstance(name, (tuple, list)):
                    name = name.list() # name may be a matrix
                name = tuple(name)
            except Exception:
                raise ValueError("{!r} is not a valid term order".format(name))

        self._blocks = tuple()
        self._weights = None
        self._matrix = None
        self._singular_moreblocks = 0

        # Singular has special C block in ring order, that compares columns of
        # monomials. This does not affect ordering of monomials of a ring, but
        # does ordering of monomials of modules over the ring.
        #
        # If the attribute below is set to an integer i, then i divided by 2
        # gives the position of the C block in the ring order; negative i means
        # the default position at the end. The remainder of i divided by 2
        # decides the ascending order if 0, or descending order if 1, of the
        # column generators. Note that the ascending and descending orders of
        # generators are written as 'C'  and 'c', respectively, in Singular
        # ring order string. See singular_str() method of this class.
        #
        # This is to be used for internal code that constructs Singular modules
        # over the ring.
        self._singular_ringorder_column = None

        if name == "block": # block term order with blocks in a list
            length = 0
            blocks = []
            singular_str = []
            macaulay2_str = []

            for t in n:
                if not isinstance(t, TermOrder):
                    t = TermOrder(t, force=True)
                if t.name() == 'block':
                    blocks = blocks + list(t.blocks())
                    singular_str.append("%s" % (t.singular_str()[1:-1],))  # [1:-1] is needed to remove parenthesis
                    macaulay2_str.append("%s" % (t.macaulay2_str()[1:-1],))
                else:
                    if len(t) == 0:
                        raise ArithmeticError("Can only concatenate term orders with length attribute.")
                    blocks.append(t)
                    if t.is_weighted_degree_order(): # true if t is a matrix order as well
                        singular_str.append("%s" % (t.singular_str(),))
                    elif t.name() == 'degneglex':
                        singular_str.append("%s" % (t.singular_str()[1:-1],))  # [1:-1] to remove (,)
                    else:
                        singular_str.append("%s(%d)" % (t.singular_str(), len(t)))
                    macaulay2_str.append("%s => %d" % (t.macaulay2_str(), len(t)))
                length += len(t)
                self._singular_moreblocks += t.singular_moreblocks()

            self._length = length
            self._name = "block"
            self._singular_str = "(" + ",".join(singular_str) + ")"
            self._macaulay2_str = "{" + ",".join(macaulay2_str) + "}"
            self._magma_str = "" # Magma does not support block order
            self._blocks = tuple(blocks)
        elif isinstance(name, str) and not isinstance(n, (tuple, list)):
            # string representation of simple or block orders
            if force:
                self._length = n
                self._name = name
                self._singular_str = singular_name_mapping.get(name,name)
                self._macaulay2_str = macaulay2_name_mapping.get(name,name)
                self._magma_str = magma_name_mapping.get(name,name)
            else:
                split_pattern = r"([^(),]+(?:\([^()]*\)[^(),]*)*)" # split by outermost commas
                block_names = re.findall(split_pattern,name)

                if len(block_names) == 0:
                    raise ValueError("no term order specified")
                elif len(block_names) == 1:
                    name = block_names[0]
                    match = re.match(r'm\(([-+0-9,]+)\)$',name)
                    if match: # matrix term order
                        m = [int(_) for _ in match.groups()[0].split(',')] # replace match.groups()[0]  with match.group(1) later
                        self.__copy(TermOrder(m))
                    else: # simple order
                        if name not in print_name_mapping.keys() and name not in singular_name_mapping.values():
                            raise ValueError("unknown term order {!r}".format(name))
                        self._length = n
                        self._name = name
                        self._singular_str = singular_name_mapping.get(name,name)
                        self._macaulay2_str = macaulay2_name_mapping.get(name,name)
                        self._magma_str = magma_name_mapping.get(name,name)
                else: # len(block_names) > 1, and hence block order represented by a string
                    length = 0
                    blocks = []
                    singular_str = []
                    macaulay2_str = []

                    length_pattern = re.compile(r"\(([0-9]+)\)$") # match with parenthesized block length at end
                    for block in block_names:
                        try:
                            block_name, block_length, _ = re.split(length_pattern,block.strip())
                            block_length = int(block_length)
                            if block_length > 0:  # ignore blocks with length 0
                                blocks.append( TermOrder(block_name, block_length, force=force) )
                                singular_str.append("%s(%d)" % (singular_name_mapping.get(block_name, block_name), block_length))
                                macaulay2_str.append("%s => %d" % (macaulay2_name_mapping.get(block_name, block_name), block_length))
                                length += block_length
                        except ValueError:
                            block_name = block.strip()
                            if block_name.lower() != "c":
                                raise ValueError("{!r} is not a valid term order (wrong part: {!r})".format(name, block))

                    if n and length != n:
                        raise ValueError("term order length does not match the number of generators")
                    self.__copy(TermOrder('block', blocks))
        elif isinstance(name, str) and isinstance(n, (tuple, list)):
            # weighted degree term orders
            if name not in print_name_mapping.keys() and name not in singular_name_mapping.values() and not force:
                raise ValueError("unknown term order {!r}".format(name))
            weights = tuple(int(w) for w in n) # n is a tuple of weights
            if any(w <= 0 for w in weights):
                raise ValueError("the degree weights must be positive integers")

            self._length = len(weights)
            self._name = name
            self._singular_str = singular_name_mapping.get(name,name) + '(' + ','.join(str(w) for w in weights) + ')'
            self._macaulay2_str = ""
            self._magma_str = ""
            self._weights = weights # defined only for weighted degree orders
        elif isinstance(name, tuple): # name represents a matrix
            if not n:
                from math import sqrt
                n = int(sqrt(len(name)))
            if n*n != len(name):
                raise ValueError("{} does not specify a square matrix".format(name))

            int_str = ','.join(str(int(e)) for e in name)

            self._length = n
            self._name = "matrix"
            self._singular_str = "M(%s)" % (int_str,)
            self._macaulay2_str = "" # Macaulay2 does not support matrix term order directly
            self._magma_str = '"weight",[%s]' % (int_str,)

            from sage.matrix.constructor import matrix
            self._matrix = matrix(n,name)  # defined only for matrix term order
            self._matrix.set_immutable()
            self._weights = name[:n] # the first row of the matrix gives weights
        else:
            raise ValueError("{!r} is not a valid term order".format(name))

        if self._length != 0:
            self._singular_str = self._singular_str % dict(ngens=self._length)
        if self._name == 'degneglex':
            self._singular_moreblocks += 1

        self.__doc__ = description_mapping.get(self._name, "No description available")

    def __hash__(self):
        r"""
        A hash function.

        EXAMPLES::

            sage: _=hash(TermOrder('lex'))
        """
        return hash((self._name, self._blocks, self._weights, self._matrix))

    def __copy(self, other):
        """
        Copy ``other`` term order to ``self``.

        EXAMPLES::

            sage: t = TermOrder('lex')
            sage: s = TermOrder(t)
            sage: t == s            # indirect test
            True
        """
        self.__dict__ = other.__dict__.copy()

    @property
    def sortkey(self):
        """
        The default ``sortkey`` method for this term order.

        EXAMPLES::

            sage: O = TermOrder()
            sage: O.sortkey.__func__ is O.sortkey_lex.__func__
            True
            sage: O = TermOrder('deglex')
            sage: O.sortkey.__func__ is O.sortkey_deglex.__func__
            True
        """

        return getattr(self, 'sortkey_' + self._name)

    @property
    def greater_tuple(self):
        """
        The default ``greater_tuple`` method for this term order.

        EXAMPLES::

            sage: O = TermOrder()
            sage: O.greater_tuple.__func__ is O.greater_tuple_lex.__func__
            True
            sage: O = TermOrder('deglex')
            sage: O.greater_tuple.__func__ is O.greater_tuple_deglex.__func__
            True
        """

        return getattr(self, 'greater_tuple_' + self._name)

    def sortkey_matrix(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the matrix
        term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='m(1,3,1,0)')                # needs sage.rings.number_field
            sage: y > x^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: y > x^3                                                               # needs sage.rings.number_field
            False
        """
        return tuple(sum(l * r for l, r in zip(row, f))
                     for row in self._matrix)

    def sortkey_lex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='lex')                       # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: x > 1                                                                 # needs sage.rings.number_field
            True
        """
        return f

    def sortkey_invlex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the inversed
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='invlex')                    # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x > 1                                                                 # needs sage.rings.number_field
            True
        """
        return f.reversed()

    def sortkey_deglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the degree
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='deglex')                    # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x > 1                                                                 # needs sage.rings.number_field
            True
        """
        return (sum(f.nonzero_values(sort=False)), f)

    def sortkey_degrevlex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        degree reversed lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='degrevlex')                 # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x > 1                                                                 # needs sage.rings.number_field
            True
        """
        return (sum(f.nonzero_values(sort=False)),
                f.reversed().emul(-1))
                # tuple(-v for v in f.reversed()))

    def sortkey_neglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the negative
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='neglex')                    # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x > 1                                                                 # needs sage.rings.number_field
            False
        """
        return tuple(-v for v in f)

    def sortkey_negdegrevlex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        negative degree reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='negdegrevlex')              # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: x > 1                                                                 # needs sage.rings.number_field
            False
        """
        return (-sum(f.nonzero_values(sort=False)),
                tuple(-v for v in f.reversed()))

    def sortkey_negdeglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        negative degree lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='negdeglex')                 # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: x > 1                                                                 # needs sage.rings.number_field
            False
        """
        return (-sum(f.nonzero_values(sort=False)), f)

    def sortkey_degneglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        degree negative lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='degneglex')               # needs sage.rings.number_field
            sage: x*y > y*z  # indirect doctest                                         # needs sage.rings.number_field
            False
            sage: x*y > x                                                               # needs sage.rings.number_field
            True
        """
        return (sum(f.nonzero_values(sort=False)), tuple(-v for v in f))

    def sortkey_wdegrevlex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        weighted degree reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('wdegrevlex',(3,2))
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order=t)                           # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x^2 > y^3                                                             # needs sage.rings.number_field
            True
        """
        return (sum(l * r for (l, r) in zip(f, self._weights)),
                tuple(-v for v in f.reversed()))

    def sortkey_wdeglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        weighted degree lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('wdeglex',(3,2))
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order=t)                           # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: x > y                                                                 # needs sage.rings.number_field
            True
        """
        return (sum(l * r for (l, r) in zip(f, self._weights)), f)

    def sortkey_negwdeglex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        negative weighted degree lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('negwdeglex',(3,2))
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order=t)                           # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: x^2 > y^3                                                             # needs sage.rings.number_field
            True
        """
        return (-sum(l * r for (l, r) in zip(f, self._weights)), f)

    def sortkey_negwdegrevlex(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        negative weighted degree reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('negwdegrevlex',(3,2))
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order=t)                           # needs sage.rings.number_field
            sage: x > y^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: x^2 > y^3                                                             # needs sage.rings.number_field
            True
        """
        return (-sum(l * r for (l, r) in zip(f, self._weights)),
                tuple(-v for v in f.reversed()))

    def sortkey_block(self, f):
        """
        Return the sortkey of an exponent tuple with respect to the
        block order as specified when constructing this element.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: P.<a,b,c,d,e,f> = PolynomialRing(QQbar, 6,                            # needs sage.rings.number_field
            ....:                                  order='degrevlex(3),degrevlex(3)')
            sage: a > c^4  # indirect doctest                                           # needs sage.rings.number_field
            False
            sage: a > e^4                                                               # needs sage.rings.number_field
            True

        TESTS:

        Check that the issue in :issue:`27139` is fixed::

            sage: R.<x,y,z,t> = PolynomialRing(AA, order='lex(2),lex(2)')               # needs sage.rings.number_field
            sage: x > y                                                                 # needs sage.rings.number_field
            True
        """
        key = tuple()
        n = 0
        for block in self:
            r = getattr(block, "sortkey_" + block.name())(f[n:n + len(block)])
            key += tuple(r)
            n += len(block)
        return key

    def greater_tuple_matrix(self,f,g):
        """
        Return the greater exponent tuple with respect to the matrix
        term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='m(1,3,1,0)')                # needs sage.rings.number_field
            sage: y > x^2  # indirect doctest                                           # needs sage.rings.number_field
            True
            sage: y > x^3                                                               # needs sage.rings.number_field
            False
        """
        for row in self._matrix:
            sf = sum(l*r for (l,r) in zip(row,f))
            sg = sum(l*r for (l,r) in zip(row,g))

            if sf > sg:
                return f
            elif sf < sg:
                return g
        return g

    def greater_tuple_lex(self,f,g):
        """
        Return the greater exponent tuple with respect to the
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='lex')                     # needs sage.rings.number_field
            sage: f = x + y^2; f.lm()  # indirect doctest                               # needs sage.rings.number_field
            x

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        return f > g and f or g

    def greater_tuple_invlex(self,f,g):
        """
        Return the greater exponent tuple with respect to the inversed
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='invlex')                  # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            y
            sage: f = y + x^2; f.lm()                                                   # needs sage.rings.number_field
            y

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        return f.reversed() > g.reversed() and f or g

    def greater_tuple_deglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the total degree
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='deglex')                  # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            x
            sage: f = x + y^2*z; f.lm()                                                 # needs sage.rings.number_field
            y^2*z

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        return ( sf > sg or ( sf == sg and f > g )) and f or g

    def greater_tuple_degrevlex(self,f,g):
        """
        Return the greater exponent tuple with respect to the total degree
        reversed lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='degrevlex')               # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            x
            sage: f = x + y^2*z; f.lm()                                                 # needs sage.rings.number_field
            y^2*z

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        return ( sf > sg or ( sf == sg and f.reversed() < g.reversed() )) and f or g

    def greater_tuple_negdegrevlex(self,f,g):
        """
        Return the greater exponent tuple with respect to the negative
        degree reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='negdegrevlex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^2*y*z^2 + x*y^3*z; f.lm()
            x*y^3*z

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        return ( sf < sg or ( sf == sg and f.reversed() < g.reversed() )) and f or g

    def greater_tuple_negdeglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the negative
        degree lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='negdeglex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^2*y*z^2 + x*y^3*z; f.lm()
            x^2*y*z^2

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        return ( sf < sg or ( sf == sg and f > g )) and f or g

    def greater_tuple_degneglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the degree negative
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='degneglex')               # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            y
            sage: f = x + y^2*z; f.lm()                                                 # needs sage.rings.number_field
            y^2*z

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        return ( sf > sg or ( sf == sg and f < g )) and f or g

    def greater_tuple_neglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the negative
        lexicographical term order.

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<a,b,c,d,e,f> = PolynomialRing(QQbar, 6,                            # needs sage.rings.number_field
            ....:                                  order='degrevlex(3),degrevlex(3)')
            sage: f = a + c^4; f.lm()  # indirect doctest                               # needs sage.rings.number_field
            c^4
            sage: g = a + e^4; g.lm()                                                   # needs sage.rings.number_field
            a
        """
        return (f < g) and f or g

    def greater_tuple_wdeglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the weighted degree
        lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('wdeglex',(1,2,3))
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order=t)                         # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            y
            sage: f = x*y + z; f.lm()                                                   # needs sage.rings.number_field
            x*y

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(l*r for (l,r) in zip(f,self._weights))
        sg = sum(l*r for (l,r) in zip(g,self._weights))
        return (sf > sg or ( sf == sg and f > g )) and f or g

    def greater_tuple_wdegrevlex(self,f,g):
        """
        Return the greater exponent tuple with respect to the weighted degree
        reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('wdegrevlex',(1,2,3))
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order=t)                         # needs sage.rings.number_field
            sage: f = x + y; f.lm()  # indirect doctest                                 # needs sage.rings.number_field
            y
            sage: f = x + y^2*z; f.lm()                                                 # needs sage.rings.number_field
            y^2*z

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(l*r for (l,r) in zip(f,self._weights))
        sg = sum(l*r for (l,r) in zip(g,self._weights))
        return (sf > sg or ( sf == sg and f.reversed() < g.reversed())) and f or g

    def greater_tuple_negwdeglex(self,f,g):
        """
        Return the greater exponent tuple with respect to the negative
        weighted degree lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: t = TermOrder('negwdeglex',(1,2,3))
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order=t)
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^3 + z; f.lm()
            x^3

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(l*r for (l,r) in zip(f,self._weights))
        sg = sum(l*r for (l,r) in zip(g,self._weights))
        return (sf < sg or ( sf == sg and f > g )) and f or g

    def greater_tuple_negwdegrevlex(self,f,g):
        """
        Return the greater exponent tuple with respect to the negative
        weighted degree reverse lexicographical term order.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: t = TermOrder('negwdegrevlex',(1,2,3))
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order=t)
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^3 + z; f.lm()
            x^3

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.
        """
        sf = sum(l*r for (l,r) in zip(f,self._weights))
        sg = sum(l*r for (l,r) in zip(g,self._weights))
        return (sf < sg or ( sf == sg and f.reversed() < g.reversed() )) and f or g

    def greater_tuple_block(self, f,g):
        """
        Return the greater exponent tuple with respect to the block
        order as specified when constructing this element.

        This method is called by the lm/lc/lt methods of
        ``MPolynomial_polydict``.

        INPUT:

        - ``f`` -- exponent tuple

        - ``g`` -- exponent tuple

        EXAMPLES::

            sage: P.<a,b,c,d,e,f> = PolynomialRing(QQbar, 6,                            # needs sage.rings.number_field
            ....:                                  order='degrevlex(3),degrevlex(3)')
            sage: f = a + c^4; f.lm()  # indirect doctest                               # needs sage.rings.number_field
            c^4
            sage: g = a + e^4; g.lm()                                                   # needs sage.rings.number_field
            a
        """
        n = 0
        for block in self:
            keyfn = getattr(block, "sortkey_" + block.name())
            f_key = keyfn(f[n:n + len(block)])
            g_key = keyfn(g[n:n + len(block)])
            if f_key != g_key:
                if f_key < g_key:
                    return g
                else:
                    return f
            n += len(block)
        return f

    def tuple_weight(self, f):
        """
        Return the weight of tuple f.

        INPUT:

        - ``f`` -- exponent tuple

        EXAMPLES::

            sage: t = TermOrder('wdeglex',(1,2,3))
            sage: P.<a,b,c> = PolynomialRing(QQbar, order=t)                            # needs sage.rings.number_field
            sage: P.term_order().tuple_weight([3,2,1])                                  # needs sage.rings.number_field
            10
        """
        return sum(l*r for (l,r) in zip(f,self._weights))

    def name(self):
        """
        EXAMPLES::

            sage: TermOrder('lex').name()
            'lex'
        """
        return self._name

    def _repr_(self):
        """
        EXAMPLES::

            sage: TermOrder('lex') # indirect doctest
            Lexicographic term order
        """
        if self._name == 'matrix':
            return 'Matrix term order with matrix\n%s' % (self._matrix,)
        elif self._name == 'block':
            s = []
            for t in self._blocks:
                if not t.is_weighted_degree_order():
                    s.append('%s of length %d' % (t,len(t)))
                else: # includes matrix order
                    s.append('%s' % (t,))
            return 'Block term order with blocks:\n(%s)' % (',\n '.join(s),)
        else:
            s = print_name_mapping.get(self._name,self._name) + ' term order'
            if self.is_weighted_degree_order():
                s = s + ' with weights %s' % (self._weights,)
            return s

    def singular_str(self):
        """
        Return a SINGULAR representation of ``self``.

        Used to convert polynomial rings to their SINGULAR representation.

        EXAMPLES::

            sage: P = PolynomialRing(GF(127), 10, names='x',
            ....:                    order='lex(3),deglex(5),lex(2)')
            sage: T = P.term_order()
            sage: T.singular_str()
            '(lp(3),Dp(5),lp(2))'
            sage: P._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: ZZ/127
            //   number of vars : 10
            //        block   1 : ordering lp
            //                  : names    x0 x1 x2
            //        block   2 : ordering Dp
            //                  : names    x3 x4 x5 x6 x7
            //        block   3 : ordering lp
            //                  : names    x8 x9
            //        block   4 : ordering C

        The ``degneglex`` ordering is somehow special, it looks like a block
        ordering in SINGULAR::

            sage: T = TermOrder("degneglex", 2)
            sage: P = PolynomialRing(QQ,2, names='x', order=T)
            sage: T = P.term_order()
            sage: T.singular_str()
            '(a(1:2),ls(2))'

            sage: T = TermOrder("degneglex", 2) + TermOrder("degneglex", 2)
            sage: P = PolynomialRing(QQ,4, names='x', order=T)
            sage: T = P.term_order()
            sage: T.singular_str()
            '(a(1:2),ls(2),a(1:2),ls(2))'
            sage: P._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 4
            //        block   1 : ordering a
            //                  : names    x0 x1
            //                  : weights   1  1
            //        block   2 : ordering ls
            //                  : names    x0 x1
            //        block   3 : ordering a
            //                  : names    x2 x3
            //                  : weights   1  1
            //        block   4 : ordering ls
            //                  : names    x2 x3
            //        block   5 : ordering C

        The position of the ``ordering C`` block can be controlled by setting
        ``_singular_ringorder_column`` attribute to an integer::

            sage: T = TermOrder("degneglex", 2) + TermOrder("degneglex", 2)
            sage: T._singular_ringorder_column = 0
            sage: P = PolynomialRing(QQ, 4, names='x', order=T)
            sage: P._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            // coefficients: QQ
            // number of vars : 4
            //        block   1 : ordering C
            //        block   2 : ordering a
            //                  : names    x0 x1
            //                  : weights   1  1
            //        block   3 : ordering ls
            //                  : names    x0 x1
            //        block   4 : ordering a
            //                  : names    x2 x3
            //                  : weights   1  1
            //        block   5 : ordering ls
            //                  : names    x2 x3

            sage: T._singular_ringorder_column = 1
            sage: P = PolynomialRing(QQ, 4, names='y', order=T)
            sage: P._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            // coefficients: QQ
            // number of vars : 4
            //        block   1 : ordering c
            //        block   2 : ordering a
            //                  : names    y0 y1
            //                  : weights   1  1
            //        block   3 : ordering ls
            //                  : names    y0 y1
            //        block   4 : ordering a
            //                  : names    y2 y3
            //                  : weights   1  1
            //        block   5 : ordering ls
            //                  : names    y2 y3

            sage: T._singular_ringorder_column = 2
            sage: P = PolynomialRing(QQ, 4, names='z', order=T)
            sage: P._singular_()                                                        # needs sage.libs.singular
            polynomial ring, over a field, global ordering
            // coefficients: QQ
            // number of vars : 4
            //        block   1 : ordering a
            //                  : names    z0 z1
            //                  : weights   1  1
            //        block   2 : ordering C
            //        block   3 : ordering ls
            //                  : names    z0 z1
            //        block   4 : ordering a
            //                  : names    z2 z3
            //                  : weights   1  1
            //        block   5 : ordering ls
            //                  : names    z2 z3
        """
        if self._singular_ringorder_column is not None:
            singular_str = self._singular_str
            if singular_str.startswith('('):
                singular_str = singular_str[1:-1]  # remove parenthesis
            split_pattern = r"([^(),]+(?:\([^()]*\)[^(),]*)*)"  # regex to split by outermost commas
            singular_str_blocks = re.findall(split_pattern, singular_str)
            if (self._singular_ringorder_column < 0 or
                self._singular_ringorder_column >= 2*len(singular_str_blocks)+2):
                singular_str_blocks.append("C")
            else:
                singular_str_blocks.insert(self._singular_ringorder_column // 2,
                       "C" if self._singular_ringorder_column % 2 == 0 else "c")
            return "(" + ",".join(singular_str_blocks) + ")"
        else:
            return self._singular_str

    def singular_moreblocks(self):
        """
        Return a the number of additional blocks SINGULAR needs to allocate
        for handling non-native orderings like ``degneglex``.

        EXAMPLES::

            sage: P = PolynomialRing(GF(127), 10, names='x',
            ....:                    order='lex(3),deglex(5),lex(2)')
            sage: T = P.term_order()
            sage: T.singular_moreblocks()
            0
            sage: P = PolynomialRing(GF(127), 10, names='x',
            ....:                    order='lex(3),degneglex(5),lex(2)')
            sage: T = P.term_order()
            sage: T.singular_moreblocks()
            1
            sage: P = PolynomialRing(GF(127), 10, names='x',
            ....:                    order='degneglex(5),degneglex(5)')
            sage: T = P.term_order()
            sage: T.singular_moreblocks()
            2

        TESTS:

        The 'degneglex' ordering is somehow special: SINGULAR handles it
        using an extra weight vector block. ::

            sage: T = TermOrder("degneglex", 2)
            sage: P = PolynomialRing(QQ,2, names='x', order=T)
            sage: T = P.term_order()
            sage: T.singular_moreblocks()
            1
            sage: T = TermOrder("degneglex", 2) + TermOrder("degneglex", 2)
            sage: P = PolynomialRing(QQ,4, names='x', order=T)
            sage: T = P.term_order()
            sage: T.singular_moreblocks()
            2
        """
        return self._singular_moreblocks

    def macaulay2_str(self):
        """
        Return a Macaulay2 representation of ``self``.

        Used to convert polynomial rings to their Macaulay2
        representation.

        EXAMPLES::

            sage: P = PolynomialRing(GF(127), 8, names='x', order='degrevlex(3),lex(5)')
            sage: T = P.term_order()
            sage: T.macaulay2_str()
            '{GRevLex => 3,Lex => 5}'
            sage: P._macaulay2_().options()['MonomialOrder']    # optional - macaulay2
            {MonomialSize => 16  }
            {GRevLex => {1, 1, 1}}
            {Lex => 5            }
            {Position => Up      }
        """
        return self._macaulay2_str

    def magma_str(self):
        """
        Return a MAGMA representation of ``self``.

        Used to convert polynomial rings to their MAGMA representation.

        EXAMPLES::

            sage: P = PolynomialRing(GF(127), 10, names='x', order='degrevlex')
            sage: magma(P)                                      # optional - magma
            Polynomial ring of rank 10 over GF(127)
            Order: Graded Reverse Lexicographical
            Variables: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9

        ::

            sage: T = P.term_order()
            sage: T.magma_str()
            '"grevlex"'
        """
        return self._magma_str

    def blocks(self):
        """
        Return the term order blocks of ``self``.

        NOTE:

        This method has been added in :issue:`11316`. There used
        to be an *attribute* of the same name and the same content.
        So, it is a backward incompatible syntax change.

        EXAMPLES::

            sage: t = TermOrder('deglex',2) + TermOrder('lex',2)
            sage: t.blocks()
            (Degree lexicographic term order, Lexicographic term order)
        """
        if self._blocks: # self is a block order
            return self._blocks
        else:
            return [self]

    def matrix(self):
        """
        Return the matrix defining matrix term order.

        EXAMPLES::

            sage: t = TermOrder("M(1,2,0,1)")                                           # needs sage.modules
            sage: t.matrix()                                                            # needs sage.modules
            [1 2]
            [0 1]
        """
        return self._matrix

    def weights(self):
        """
        Return the weights for weighted term orders.

        EXAMPLES::

            sage: t = TermOrder('wdeglex',(2,3))
            sage: t.weights()
            (2, 3)
        """
        return self._weights

    def __eq__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` are equal.

        EXAMPLES::

            sage: TermOrder('lex') == TermOrder('lex',3)
            True

        ::

            sage: TermOrder('degrevlex') == TermOrder('lex')
            False

        ::

            sage: T1 = TermOrder('lex',2) + TermOrder('lex',3)
            sage: T2 = TermOrder('lex',3) + TermOrder('lex',2)
            sage: T1 == T2
            False

        ::

            sage: T1 = TermOrder('lex',2) + TermOrder('neglex',3)
            sage: T2 = TermOrder('lex',2) + TermOrder('neglex',3)
            sage: T1 == T2
            True

        TESTS:

        We assert that comparisons take into account the block size of
        orderings (cf. :issue:`24981`)::

            sage: R = PolynomialRing(QQ, 6, 'x', order="lex(1),degrevlex(5)")
            sage: S = R.change_ring(order="lex(2),degrevlex(4)")
            sage: R == S
            False
            sage: S.term_order() == R.term_order()
            False
            sage: S.term_order() == TermOrder('lex', 2) + TermOrder('degrevlex', 4)
            True
        """
        if not isinstance(other, TermOrder):
            try:
                other = TermOrder(other, force=True)
            except Exception:
                return False

        return (self._name == other._name
            and self._blocks == other._blocks
            and (not self.is_block_order()
                or all(len(t1) == len(t2) for (t1, t2) in zip(self._blocks, other._blocks)))
            and self._weights == other._weights
            and self._matrix == other._matrix
            and self._singular_ringorder_column == other._singular_ringorder_column)

    def __ne__(self, other):
        """
        Return ``True`` if ``self`` and ``other`` are not equal.

        EXAMPLES::

            sage: T1 = TermOrder('lex',2) + TermOrder('lex',3)
            sage: T2 = TermOrder('lex',3) + TermOrder('lex',2)
            sage: T1 != T2
            True
        """
        return not self == other

    def __add__(self, other):
        """
        Construct a block order combining ``self`` and ``other``.

        INPUT:

        - ``other`` -- a term order

        OUTPUT: a block order

        EXAMPLES::

            sage: from sage.rings.polynomial.term_order import TermOrder
            sage: TermOrder('deglex',2) + TermOrder('degrevlex(3),neglex(3)')
            Block term order with blocks:
            (Degree lexicographic term order of length 2,
             Degree reverse lexicographic term order of length 3,
             Negative lexicographic term order of length 3)
        """
        if isinstance(other, TermOrder):
            return TermOrder('block',[self,other])
        else:
            return self

    def __len__(self):
        """
        Return the length of this term order, i.e. the number of
        variables it covers. This may be zero for indefinitely many
        variables.

        EXAMPLES::

            sage: T = TermOrder('lex')
            sage: len(T)
            0
            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: len(T)
            5
        """
        return self._length

    def __getitem__(self, i):
        r"""
        Return the `i`-th block of this term order.

        INPUT:

        - ``i`` -- index

        EXAMPLES::

            sage: T = TermOrder('lex')
            sage: T[0]
            Lexicographic term order

        ::

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: T[1]
            Degree reverse lexicographic term order

        Note that ``len(self)`` does not count blocks but
        variables.

        ::

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: T[len(T)-1]
            Traceback (most recent call last):
            \dots
            IndexError: tuple index out of range
        """
        return self.blocks()[i]

    def __iter__(self):
        r"""
        Iterate over the blocks of this term order.

        EXAMPLES::

            sage: T = TermOrder('lex')
            sage: list(T) # indirect doctest
            [Lexicographic term order]

        ::

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: list(T)
            [Lexicographic term order, Degree reverse lexicographic term order]

        Note that ``len(self)`` and
        ``len(list(self))`` are not the same. The former counts
        the number of variables in ``self`` while the latter
        counts the number of blocks.
        """
        return iter(self.blocks())

    def is_global(self):
        r"""
        Return ``True`` if this term order is definitely
        global. Return false otherwise, which includes
        unknown term orders.

        EXAMPLES::

            sage: T = TermOrder('lex')
            sage: T.is_global()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('degrevlex', 3)
            sage: T.is_global()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_global()
            False
            sage: T = TermOrder('degneglex', 3)
            sage: T.is_global()
            True
            sage: T = TermOrder('invlex', 3)
            sage: T.is_global()
            True
        """
        if self.name() in ('lex', 'degrevlex', 'deglex', 'degneglex',
                           'wdegrevlex', 'wdeglex', 'invlex'):
            return True
        elif self.name() == 'block':
            return all(t.is_global() for t in self.blocks())
        else:
            return False

    def is_local(self):
        r"""
        Return ``True`` if this term order is definitely
        local. Return false otherwise, which includes
        unknown term orders.

        EXAMPLES::

            sage: T = TermOrder('lex')
            sage: T.is_local()
            False
            sage: T = TermOrder('negdeglex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_local()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_local()
            False
        """
        if (self.name() in ('neglex', 'negdegrevlex', 'negdeglex',
                            'negwdegrevlex', 'negwdeglex') or
            self.singular_str() in ('ls', 'ds', 'Ds', 'ws', 'Ws')):
            return True
        elif self.name() == 'block':
            return all(t.is_local() for t in self.blocks())
        else:
            return False

    def is_block_order(self):
        """
        Return ``True`` if ``self`` is a block term order.

        EXAMPLES::

            sage: t = TermOrder('deglex',2) + TermOrder('lex',2)
            sage: t.is_block_order()
            True
        """
        return self._name == 'block'

    def is_weighted_degree_order(self):
        """
        Return ``True`` if ``self`` is a weighted degree term order.

        EXAMPLES::

            sage: t = TermOrder('wdeglex',(2,3))
            sage: t.is_weighted_degree_order()
            True
        """
        return self._weights is not None


def termorder_from_singular(S):
    """
    Return the Sage term order of the basering in the given Singular interface.

    INPUT:

    An instance of the Singular interface.

    EXAMPLES::

        sage: from sage.rings.polynomial.term_order import termorder_from_singular
        sage: singular.eval('ring r1 = (9,x),(a,b,c,d,e,f),(M((1,2,3,0)),wp(2,3),lp)')  # needs sage.libs.singular
        ''
        sage: termorder_from_singular(singular)                                         # needs sage.libs.singular
        Block term order with blocks:
        (Matrix term order with matrix
           [1 2]
           [3 0],
         Weighted degree reverse lexicographic term order with weights (2, 3),
         Lexicographic term order of length 2)

    A term order in Singular also involves information on orders for modules.
    This information is reflected in ``_singular_ringorder_column`` attribute of
    the term order. ::

        sage: # needs sage.libs.singular
        sage: singular.ring(0, '(x,y,z,w)', '(C,dp(2),lp(2))')
        polynomial ring, over a field, global ordering
        // coefficients: QQ
        // number of vars : 4
        //        block   1 : ordering C
        //        block   2 : ordering dp
        //                  : names    x y
        //        block   3 : ordering lp
        //                  : names    z w
        sage: T = termorder_from_singular(singular)
        sage: T
        Block term order with blocks:
        (Degree reverse lexicographic term order of length 2,
         Lexicographic term order of length 2)
        sage: T._singular_ringorder_column
        0

        sage: # needs sage.libs.singular
        sage: singular.ring(0, '(x,y,z,w)', '(c,dp(2),lp(2))')
        polynomial ring, over a field, global ordering
        // coefficients: QQ
        // number of vars : 4
        //        block   1 : ordering c
        //        block   2 : ordering dp
        //                  : names    x y
        //        block   3 : ordering lp
        //                  : names    z w
        sage: T = termorder_from_singular(singular)
        sage: T
        Block term order with blocks:
        (Degree reverse lexicographic term order of length 2,
         Lexicographic term order of length 2)
        sage: T._singular_ringorder_column
        1

    TESTS:

    Check that ``degneglex`` term orders are converted correctly
    (:issue:`29635`)::

        sage: # needs sage.libs.singular
        sage: _ = singular.ring(0, '(x,y,z,w)', '(a(1:4),ls(4))')
        sage: termorder_from_singular(singular).singular_str()
        '(a(1:4),ls(4))'
        sage: _ = singular.ring(0, '(x,y,z,w)', '(a(1:2),ls(2),a(1:2),ls(2))')
        sage: termorder_from_singular(singular).singular_str()
        '(a(1:2),ls(2),a(1:2),ls(2))'
        sage: _ = singular.ring(0, '(x,y,z,w)', '(a(1:2),ls(2),C,a(1:2),ls(2))')
        sage: termorder_from_singular(singular).singular_str()
        '(a(1:2),ls(2),C,a(1:2),ls(2))'
        sage: PolynomialRing(QQ, 'x,y', order='degneglex')('x^2')._singular_().sage()
        x^2
    """
    from sage.rings.integer_ring import ZZ
    singular = S
    T = singular('ringlist(basering)[3]')
    order = []
    ringorder_column = None
    weights_one_block = False
    for idx, block in enumerate(T):
        blocktype = singular.eval('%s[1]' % block.name())
        if blocktype in ['a']:
            weights = list(block[2].sage())
            weights_one_block = all(w == 1 for w in weights)
            continue
        elif blocktype == 'c':
            ringorder_column = 2*idx + 1
        elif blocktype == 'C':
            if idx < len(T) - 1:  # skip Singular default
                ringorder_column = 2*idx
        elif blocktype == 'M':
            from sage.matrix.constructor import matrix
            coefs = list(block[2].sage())
            n = ZZ(len(coefs)).sqrt()
            order.append(TermOrder(matrix(n,coefs)))
        elif weights_one_block and blocktype == 'ls':
            # 'degneglex' is encoded as '(a(1:n),ls(n))'
            n = ZZ(singular.eval("size(%s[2])" % block.name()))
            order.append(TermOrder('degneglex', n))
        elif blocktype[0] in ['w','W']:
            order.append(TermOrder(inv_singular_name_mapping[blocktype], list(block[2].sage())))
        else:
            order.append(TermOrder(inv_singular_name_mapping[blocktype], ZZ(singular.eval("size(%s[2])" % block.name()))))
        weights_one_block = False

    if not order:
        raise ValueError("invalid term order in Singular")
    out = order.pop(0)
    while order:
        out = out + order.pop(0)

    if ringorder_column is not None:
        out._singular_ringorder_column = ringorder_column

    return out
