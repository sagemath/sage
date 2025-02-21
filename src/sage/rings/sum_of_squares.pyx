# distutils: libraries = m
r"""
Fast decomposition of small integers into sums of squares

Implement fast version of decomposition of (small) integers into sum of squares
by direct method not relying on factorisation.

AUTHORS:

- Vincent Delecroix (2014): first implementation (:issue:`16374`)
"""
#*****************************************************************************
#       Copyright (C) 2014 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.math cimport sqrt
from cysignals.signals cimport sig_on, sig_off

cimport sage.rings.integer as integer
from sage.rings import integer

cdef int two_squares_c(uint_fast32_t n, uint_fast32_t res[2]) noexcept:
    r"""
    Return ``1`` if ``n`` is a sum of two squares and ``0`` otherwise.

    If ``1`` is returned, the value of ``res[0]`` and ``res[1]`` are set to the
    lexicographically smallest solution of `a^2 + b^2 = n`.
    """
    cdef uint_fast32_t fac,i,ii,j,jj,nn

    if n == 0:
        res[0] = res[1] = 0
        return 1

    # if n = 0 mod 4 then i and j must be even
    # hence, we first remove the maximum power of 4 from n and will then
    # multiply by the corresponding power of 2 the solution
    fac = 0
    while n%4 == 0:
        n >>= 2
        fac += 1

    # now, n is congruent to 1,2 or 3 mod 4.
    # As a square is congruent to 0,1 mod 4, a sum of square is congruent to
    # 0,1,2 mod 4.
    if n%4 == 3:
        return 0

    # if n=1 mod 4 then exactly one of i or j must be even
    # if n=2 mod 4 then i and j must be odd
    if n%4 == 1:
        i = ii = 0
        j = <uint_fast32_t> sqrt(<double> n)
        jj = j*j
        while ii <= jj:
            nn = n - ii
            while jj > nn:
                j -= 1
                # strangely enough, the 1-by-1 decreasing above is much faster
                # than integer Newton iteration:
                # j = (j+nn/j)/2
                jj = j*j
            if jj == nn:
                res[0] = i<<fac
                res[1] = j<<fac
                return 1
            i += 1
            ii = i*i
    else: # n mod 4 = 2
        i = ii = 1
        j = <uint_fast32_t> sqrt(<double> n)
        j += 1 - j%2
        jj = j*j
        while ii <= jj:
            nn = n - ii
            while jj > nn:
                j -= 2
                # strangely enough, the 2-by-2 decreasing above is much faster
                # than integer Newton iteration:
                # j = (j+nn/j)/2
                jj = j*j
            if jj == nn:
                res[0] = i<<fac
                res[1] = j<<fac
                return 1
            i += 2
            ii = i*i

    return 0


cdef int three_squares_c(uint_fast32_t n, uint_fast32_t res[3]) noexcept:
    r"""
    Return `1` if `n` is a sum of three squares and `0` otherwise.

    If `1` is returned, then the values of ``res[0]``, ``res[1]`` and
    ``res[2]`` are set to a solution of `a^2 + b^2 + c^2 = n` such
    that `a \leq b \leq c`.
    """
    cdef uint_fast32_t fac,i

    if n == 0:
        res[0] = res[1] = res[2] = 0
        return 1

    # if n == 0 mod 4 then i,j,k must be even
    # hence we remove from n the maximum power of 4 and at the very end we
    # multiply each term of the solution by the appropriate power of 2
    fac = 0
    while n%4 == 0:
        n >>= 2
        fac += 1

    # Legendre's three-square theorem: n is a sum of three squares if and only
    # if it is not of the form 4^a(8b + 7)
    if n%8 == 7:
        return 0

    i = <uint_fast32_t> sqrt(<double> n)
    while not two_squares_c(n-i*i, res):
        i -= 1
    res[0] <<= fac
    res[1] <<= fac
    res[2] = i<<fac

    return 1


def two_squares_pyx(uint32_t n):
    r"""
    Return a pair of nonnegative integers ``(i,j)`` such that `i^2 + j^2 = n`.

    If ``n`` is not a sum of two squares, a :exc:`ValueError` is raised. The input
    must be less than `2^{32}=4294967296`, otherwise an :exc:`OverflowError` is
    raised.

    .. SEEALSO::

        :func:`~sage.arith.misc.two_squares` is much more suited for large inputs

    EXAMPLES::

        sage: from sage.rings.sum_of_squares import two_squares_pyx
        sage: two_squares_pyx(0)
        (0, 0)
        sage: two_squares_pyx(1)
        (0, 1)
        sage: two_squares_pyx(2)
        (1, 1)
        sage: two_squares_pyx(3)
        Traceback (most recent call last):
        ...
        ValueError: 3 is not a sum of 2 squares
        sage: two_squares_pyx(106)
        (5, 9)

        sage: two_squares_pyx(2**32)
        Traceback (most recent call last):
        ...
        OverflowError: ...

    TESTS::

        sage: s = lambda t: sum(i^2 for i in t)
        sage: for ij in Subsets(Subsets(45000, 15).random_element(), 2):
        ....:     if s(two_squares_pyx(s(ij))) != s(ij):
        ....:         print("hey")

        sage: for n in range(1,65536):
        ....:     if two_squares_pyx(n^2) != (0, n):
        ....:         print("hey")
        ....:     if two_squares_pyx(n^2 + 1) != (1, n):
        ....:         print("ho")
    """
    cdef uint_fast32_t i[2]

    sig_on()
    if two_squares_c(n, i):
        sig_off()
        return (integer.smallInteger(i[0]), integer.smallInteger(i[1]))
    sig_off()

    raise ValueError("%d is not a sum of 2 squares" % n)


def is_sum_of_two_squares_pyx(uint32_t n):
    r"""
    Return ``True`` if ``n`` is a sum of two squares and ``False`` otherwise.

    The input must be smaller than `2^{32} = 4294967296`, otherwise an
    :exc:`OverflowError` is raised.

    EXAMPLES::

        sage: from sage.rings.sum_of_squares import is_sum_of_two_squares_pyx
        sage: [x for x in range(30) if is_sum_of_two_squares_pyx(x)]
        [0, 1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20, 25, 26, 29]

        sage: is_sum_of_two_squares_pyx(2**32)
        Traceback (most recent call last):
        ...
        OverflowError: ...
    """
    cdef uint_fast32_t i[2]

    sig_on()
    if two_squares_c(n, i):
        sig_off()
        return True
    else:
        sig_off()
        return False


def three_squares_pyx(uint32_t n):
    r"""
    If ``n`` is a sum of three squares return a 3-tuple ``(i,j,k)`` of Sage
    integers such that `i^2 + j^2 + k^2 = n` and `i \leq j \leq k`. Otherwise
    raise a :exc:`ValueError`.

    The input must be less than `2^{32}=4294967296`, otherwise an
    :exc:`OverflowError` is raised.

    EXAMPLES::

        sage: from sage.rings.sum_of_squares import three_squares_pyx
        sage: three_squares_pyx(0)
        (0, 0, 0)
        sage: three_squares_pyx(1)
        (0, 0, 1)
        sage: three_squares_pyx(2)
        (0, 1, 1)
        sage: three_squares_pyx(3)
        (1, 1, 1)
        sage: three_squares_pyx(4)
        (0, 0, 2)
        sage: three_squares_pyx(5)
        (0, 1, 2)
        sage: three_squares_pyx(6)
        (1, 1, 2)
        sage: three_squares_pyx(7)
        Traceback (most recent call last):
        ...
        ValueError: 7 is not a sum of 3 squares
        sage: three_squares_pyx(107)
        (1, 5, 9)

        sage: three_squares_pyx(2**32)
        Traceback (most recent call last):
        ...
        OverflowError: ...

    TESTS::

        sage: s = lambda t: sum(i^2 for i in t)
        sage: for ijk in Subsets(Subsets(35000,15).random_element(),3):
        ....:     if s(three_squares_pyx(s(ijk))) != s(ijk):
        ....:         print("hey")
    """
    cdef uint_fast32_t i[3]

    sig_on()
    if three_squares_c(n, i):
        sig_off()
        return (integer.smallInteger(i[0]), integer.smallInteger(i[1]), integer.smallInteger(i[2]))
    sig_off()

    raise ValueError("%d is not a sum of 3 squares" % n)


def four_squares_pyx(uint32_t n):
    r"""
    Return a 4-tuple of nonnegative integers ``(i,j,k,l)`` such that `i^2 + j^2
    + k^2 + l^2 = n` and `i \leq j \leq k \leq l`.

    The input must be less than `2^{32}=4294967296`, otherwise an
    :exc:`OverflowError` is raised.

    .. SEEALSO::

        :func:`~sage.arith.misc.four_squares` is much more suited for large input

    EXAMPLES::

        sage: from sage.rings.sum_of_squares import four_squares_pyx
        sage: four_squares_pyx(15447)
        (2, 5, 17, 123)
        sage: 2^2 + 5^2 + 17^2 + 123^2
        15447

        sage: four_squares_pyx(523439)
        (3, 5, 26, 723)
        sage: 3^2 + 5^2 + 26^2 + 723^2
        523439

        sage: four_squares_pyx(2**32)
        Traceback (most recent call last):
        ...
        OverflowError: ...

    TESTS::

        sage: four_squares_pyx(0)
        (0, 0, 0, 0)

        sage: s = lambda t: sum(i^2 for i in t)
        sage: all(s(four_squares_pyx(n)) == n for n in range(5000,10000))
        True
    """
    cdef uint_fast32_t fac, j
    cdef uint_fast32_t i[3]

    if n == 0:
        return (integer.smallInteger(0),)*4

    # division by power of 4
    fac = 0
    while n % 4 == 0:
        n >>= 2
        fac += 1

    sig_on()
    # we pick the largest square we can for j
    j = <uint_fast32_t> sqrt(<double> n)
    while not three_squares_c(n-j*j, i):
        j -= 1
    sig_off()

    return (integer.smallInteger((i[0])<<fac), integer.smallInteger((i[1])<<fac),
            integer.smallInteger((i[2])<<fac), integer.smallInteger(j<<fac))
