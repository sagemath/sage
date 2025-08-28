# distutils: extra_compile_args = -D_XPG6
"""
FLINT Arithmetic Functions
"""

#*****************************************************************************
#       Copyright (C) 2013 Fredrik Johansson <fredrik.johansson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpq cimport *
from sage.libs.flint.arith cimport *


from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational


def bell_number(unsigned long n):
    """
    Return the `n`-th Bell number.

    See :wikipedia:`Bell_number`.

    ALGORITHM:

    Uses :c:func:`arith_bell_number`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import bell_number
        sage: [bell_number(i) for i in range(10)]
        [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
        sage: bell_number(10)
        115975
        sage: bell_number(40)
        157450588391204931289324344702531067
        sage: bell_number(100)
        47585391276764833658790768841387207826363669686825611466616334637559114497892442622672724044217756306953557882560751
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans = Integer(0)

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()
    arith_bell_number(ans_fmpz, n)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    if n > 1000:
        sig_off()

    return ans


def bernoulli_number(unsigned long n):
    """
    Return the `n`-th Bernoulli number.

    See :wikipedia:`Bernoulli_number`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import bernoulli_number
        sage: [bernoulli_number(i) for i in range(10)]
        [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0]
        sage: bernoulli_number(10)
        5/66
        sage: bernoulli_number(40)
        -261082718496449122051/13530
        sage: bernoulli_number(100)
        -94598037819122125295227433069493721872702841533066936133385696204311395415197247711/33330
    """
    cdef fmpq_t ans_fmpq
    cdef Rational ans = <Rational>Rational.__new__(Rational)

    fmpq_init(ans_fmpq)
    sig_on()
    arith_bernoulli_number(ans_fmpq, n)
    sig_off()
    fmpq_get_mpq(ans.value, ans_fmpq)
    fmpq_clear(ans_fmpq)

    return ans


def euler_number(unsigned long n):
    """
    Return the Euler number of index `n`.

    See :wikipedia:`Euler_number`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import euler_number
        sage: [euler_number(i) for i in range(8)]
        [1, 0, -1, 0, 5, 0, -61, 0]
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans = Integer(0)

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()
    arith_euler_number(ans_fmpz, n)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    if n > 1000:
        sig_off()

    return ans


def stirling_number_1(long n, long k):
    """
    Return the unsigned Stirling number of the first kind.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import stirling_number_1
        sage: [stirling_number_1(8,i) for i in range(9)]
        [0, 5040, 13068, 13132, 6769, 1960, 322, 28, 1]
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans = Integer(0)

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()
    arith_stirling_number_1u(ans_fmpz, n, k)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    if n > 1000:
        sig_off()
    return ans


def stirling_number_2(long n, long k):
    """
    Return the Stirling number of the second kind.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import stirling_number_2
        sage: [stirling_number_2(8,i) for i in range(9)]
        [0, 1, 127, 966, 1701, 1050, 266, 28, 1]
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans = Integer(0)

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()
    arith_stirling_number_2(ans_fmpz, n, k)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    if n > 1000:
        sig_off()
    return ans


def number_of_partitions(unsigned long n):
    """
    Return the number of partitions of the integer `n`.

    See :wikipedia:`Partition_(number_theory)`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import number_of_partitions
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    TESTS::

        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0  # long time
        True
    """
    cdef fmpz_t ans_fmpz
    cdef Integer ans

    fmpz_init(ans_fmpz)

    if n > 1000:
        sig_on()

    arith_number_of_partitions(ans_fmpz, n)

    if n > 1000:
        sig_off()

    ans = Integer(0)
    fmpz_get_mpz(ans.value, ans_fmpz)
    fmpz_clear(ans_fmpz)
    return ans


def dedekind_sum(p, q):
    """
    Return the Dedekind sum `s(p, q)` where `p` and `q` are arbitrary integers.

    See :wikipedia:`Dedekind_sum`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import dedekind_sum
        sage: dedekind_sum(4, 5)
        -1/5
    """
    p = Integer(p)
    q = Integer(q)
    s = Rational(0)

    cdef fmpz_t p_fmpz, q_fmpz
    cdef fmpq_t s_fmpq

    fmpz_init(p_fmpz)
    fmpz_init(q_fmpz)
    fmpq_init(s_fmpq)

    fmpz_set_mpz(p_fmpz, (<Integer>p).value)
    fmpz_set_mpz(q_fmpz, (<Integer>q).value)

    arith_dedekind_sum(s_fmpq, p_fmpz, q_fmpz)

    fmpq_get_mpq((<Rational>s).value, s_fmpq)

    fmpz_clear(p_fmpz)
    fmpz_clear(q_fmpz)
    fmpq_clear(s_fmpq)

    return s


def harmonic_number(unsigned long n):
    """
    Return the harmonic number `H_n`.

    See :wikipedia:`Harmonic_number`.

    EXAMPLES::

        sage: from sage.libs.flint.arith_sage import harmonic_number
        sage: n = 500 + randint(0,500)
        sage: bool( sum(1/k for k in range(1,n+1)) == harmonic_number(n) )
        True
    """
    s = Rational(0)
    cdef fmpq_t s_fmpq

    fmpq_init(s_fmpq)

    sig_on()
    arith_harmonic_number(s_fmpq, n)

    fmpq_get_mpq((<Rational>s).value, s_fmpq)
    sig_off()

    fmpq_clear(s_fmpq)

    return s
