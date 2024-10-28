"""
Fast Arithmetic Functions
"""

# ****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_check

from sage.libs.gmp.mpz cimport mpz_lcm, mpz_set_ui
from sage.rings.integer cimport Integer
from sage.structure.coerce cimport coercion_model


def lcm(a, b=None):
    r"""
    The least common multiple of a and b, or if a is a list and b is
    omitted the least common multiple of all elements of a.

    Note that LCM is an alias for lcm.

    INPUT:

    - ``a``, ``b`` -- two elements of a ring with lcm or

    - ``a`` -- list; tuple or iterable of elements of a ring with lcm

    OUTPUT:

    First, the given elements are coerced into a common parent. Then,
    their least common multiple *in that parent* is returned.

    EXAMPLES::

        sage: lcm(97, 100)
        9700
        sage: LCM(97, 100)
        9700
        sage: LCM(0, 2)
        0
        sage: LCM(-3, -5)
        15
        sage: LCM([1,2,3,4,5])
        60
        sage: v = LCM(range(1, 10000))   # *very* fast!
        sage: len(str(v))
        4349

    TESTS:

    The following tests against a bug that was fixed in :issue:`10771`::

        sage: lcm(4/1,2)
        4

    The following shows that indeed coercion takes place before
    computing the least common multiple::

        sage: R.<x> = QQ[]
        sage: S.<x> = ZZ[]
        sage: p = S.random_element(degree=(0,5))
        sage: q = R.random_element(degree=(0,5))
        sage: parent(lcm([1/p,q]))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field

    Make sure we try `\QQ` and not merely `\ZZ` (:issue:`13014`)::

        sage: bool(lcm(2/5, 3/7) == lcm(SR(2/5), SR(3/7)))                              # needs sage.symbolic
        True

    Make sure that the lcm of Expressions stays symbolic::

        sage: parent(lcm(2, 4))
        Integer Ring
        sage: parent(lcm(SR(2), 4))                                                     # needs sage.symbolic
        Symbolic Ring
        sage: parent(lcm(2, SR(4)))                                                     # needs sage.symbolic
        Symbolic Ring
        sage: parent(lcm(SR(2), SR(4)))                                                 # needs sage.symbolic
        Symbolic Ring

    Verify that objects without lcm methods but which can't be
    coerced to `\ZZ` or `\QQ` raise an error::

        sage: F.<x,y> = FreeMonoid(2)                                                   # needs sage.groups
        sage: lcm(x,y)                                                                  # needs sage.groups
        Traceback (most recent call last):
        ...
        TypeError: unable to find lcm of x and y

    Check rational and integers (:issue:`17852`)::

        sage: lcm(1/2, 4)
        4
        sage: lcm(4, 1/2)
        4

    Check that we do not mutate the list (:issue:`22630`)::

        sage: L = [int(1), int(2)]
        sage: lcm(L)
        2
        sage: [type(x).__name__ for x in L]
        ['int', 'int']
    """
    if b is None:
        return LCM_list(a)

    try:
        return a.lcm(b)
    except (AttributeError, TypeError):
        pass
    try:
        return Integer(a).lcm(Integer(b))
    except TypeError:
        pass
    raise TypeError(f"unable to find lcm of {a!r} and {b!r}")


cpdef LCM_list(v):
    """
    Return the LCM of an iterable ``v``.

    Elements of ``v`` are converted to Sage objects if they aren't
    already.

    This function is used, e.g., by :func:`~sage.arith.functions.lcm`.

    INPUT:

    - ``v`` -- an iterable

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.arith.functions import LCM_list
        sage: w = LCM_list([3,9,30]); w
        90
        sage: type(w)
        <class 'sage.rings.integer.Integer'>

    The inputs are converted to Sage integers::

        sage: w = LCM_list([int(3), int(9), int(30)]); w
        90
        sage: type(w)
        <class 'sage.rings.integer.Integer'>

    TESTS::

        sage: from sage.structure.sequence import Sequence
        sage: from sage.arith.functions import LCM_list
        sage: l = Sequence(())
        sage: LCM_list(l)
        1
        sage: LCM_list([])
        1

    This is because ``lcm(0,x) = 0`` for all ``x`` (by convention)::

        sage: LCM_list(Sequence(srange(100)))
        0
        sage: LCM_list(range(100))
        0

    So for the lcm of all integers up to 10 you must do this::

        sage: LCM_list(Sequence(srange(1,100)))
        69720375229712477164533808935312303556800

    Note that the following example did not work in QQ[] as of 2.11,
    but does in 3.1.4; the answer is different, though equivalent::

        sage: R.<X> = ZZ[]
        sage: LCM_list(Sequence((2*X+4,2*X^2,2)))
        2*X^3 + 4*X^2
        sage: R.<X> = QQ[]
        sage: LCM_list(Sequence((2*X+4,2*X^2,2)))
        X^3 + 2*X^2
    """
    cdef Integer x
    cdef Integer z = <Integer>(Integer.__new__(Integer))
    mpz_set_ui(z.value, 1)

    itr = iter(v)
    for elt in itr:
        sig_check()
        if isinstance(elt, Integer):
            x = <Integer>elt
        elif isinstance(elt, int):
            x = Integer(elt)
        else:
            # The result is no longer an Integer, pass to generic code
            a, b = coercion_model.canonical_coercion(z, elt)
            return LCM_generic(itr, a.lcm(b))
        mpz_lcm(z.value, z.value, x.value)

    return z


cdef LCM_generic(itr, ret):
    """
    Return the least common multiple of the element ``ret`` and the
    elements in the iterable ``itr``.
    """
    for vi in itr:
        sig_check()
        a, b = coercion_model.canonical_coercion(vi, ret)
        ret = a.lcm(b)
        if not ret:
            return ret
    return ret
