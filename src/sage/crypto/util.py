# sage.doctest: needs sage.combinat
"""
Utility Functions for Cryptography

Miscellaneous utility functions for cryptographic purposes.

AUTHORS:

- Minh Van Nguyen (2009-12): initial version with the following functions:
  ``ascii_integer``, ``ascii_to_bin``, ``bin_to_ascii``, ``has_blum_prime``,
  ``is_blum_prime``, ``least_significant_bits``, ``random_blum_prime``.
"""

#*****************************************************************************
#       Copyright (c) 2009, 2010 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.functions import lcm
from sage.arith.misc import is_prime, primes, random_prime
from sage.misc.lazy_import import lazy_import
from sage.rings.finite_rings.integer_mod import Mod as mod
from sage.rings.integer import Integer

lazy_import('sage.arith.misc', ('carmichael_lambda'), deprecation=34719)
lazy_import('sage.monoids.string_monoid', 'BinaryStrings')


def ascii_integer(B):
    r"""
    Return the ASCII integer corresponding to the binary string ``B``.

    INPUT:

    - ``B`` -- a non-empty binary string or a non-empty list of bits. The
      number of bits in ``B`` must be 8

    OUTPUT: the ASCII integer corresponding to the 8-bit block ``B``

    EXAMPLES:

    The ASCII integers of some binary strings::

        sage: from sage.crypto.util import ascii_integer
        sage: bin = BinaryStrings()
        sage: B = bin.encoding("A"); B
        01000001
        sage: ascii_integer(B)
        65
        sage: B = bin.encoding("C"); list(B)
        [0, 1, 0, 0, 0, 0, 1, 1]
        sage: ascii_integer(list(B))
        67
        sage: ascii_integer("01000100")
        68
        sage: ascii_integer([0, 1, 0, 0, 0, 1, 0, 1])
        69

    TESTS:

    The input ``B`` must be a non-empty string or a non-empty list of bits::

        sage: from sage.crypto.util import ascii_integer
        sage: ascii_integer("")
        Traceback (most recent call last):
        ...
        ValueError: B must consist of 8 bits.
        sage: ascii_integer([])
        Traceback (most recent call last):
        ...
        ValueError: B must consist of 8 bits.

    The input ``B`` must be an 8-bit string or a list of 8 bits::

        sage: from sage.crypto.util import ascii_integer
        sage: ascii_integer("101")
        Traceback (most recent call last):
        ...
        ValueError: B must consist of 8 bits.
        sage: ascii_integer([1, 0, 1, 1])
        Traceback (most recent call last):
        ...
        ValueError: B must consist of 8 bits.
    """
    if len(B) != 8:
        raise ValueError("B must consist of 8 bits.")
    L = [int(str(x)) for x in list(B)]
    return sum([L[7], L[6]*2, L[5]*4, L[4]*8,
                L[3]*16, L[2]*32, L[1]*64, L[0]*128])

def ascii_to_bin(A):
    r"""
    Return the binary representation of the ASCII string ``A``.

    INPUT:

    - ``A`` -- string or list of ASCII characters

    OUTPUT: the binary representation of ``A``

    ALGORITHM:

    Let `A = a_0 a_1 \cdots a_{n-1}` be an ASCII string, where each `a_i` is
    an ASCII character. Let `c_i` be the ASCII integer corresponding to `a_i`
    and let `b_i` be the binary representation of `c_i`. The binary
    representation `B` of `A` is `B = b_0 b_1 \cdots b_{n-1}`.

    EXAMPLES:

    The binary representation of some ASCII strings::

        sage: from sage.crypto.util import ascii_to_bin
        sage: ascii_to_bin("A")
        01000001
        sage: ascii_to_bin("Abc123")
        010000010110001001100011001100010011001000110011

    The empty string is different from the string with one space character.
    For the empty string and the empty list, this function returns the same
    result::

        sage: from sage.crypto.util import ascii_to_bin
        sage: ascii_to_bin("")
        <BLANKLINE>
        sage: ascii_to_bin(" ")
        00100000
        sage: ascii_to_bin([])
        <BLANKLINE>

    This function also accepts a list of ASCII characters. You can also pass
    in a list of strings::

        sage: from sage.crypto.util import ascii_to_bin
        sage: ascii_to_bin(["A", "b", "c", "1", "2", "3"])
        010000010110001001100011001100010011001000110011
        sage: ascii_to_bin(["A", "bc", "1", "23"])
        010000010110001001100011001100010011001000110011

    TESTS:

    For a list of ASCII characters or strings, do not mix characters or
    strings with integers::

        sage: from sage.crypto.util import ascii_to_bin
        sage: ascii_to_bin(["A", "b", "c", 1, 2, 3])
        Traceback (most recent call last):
        ...
        TypeError: sequence item 3: expected str..., sage.rings.integer.Integer found
        sage: ascii_to_bin(["Abc", 1, 2, 3])
        Traceback (most recent call last):
        ...
        TypeError: sequence item 1: expected str..., sage.rings.integer.Integer found
    """
    bin = BinaryStrings()
    return bin.encoding("".join(list(A)))

def bin_to_ascii(B):
    r"""
    Return the ASCII representation of the binary string ``B``.

    INPUT:

    - ``B`` -- a non-empty binary string or a non-empty list of bits. The
      number of bits in ``B`` must be a multiple of 8

    OUTPUT: the ASCII string corresponding to ``B``

    ALGORITHM:

    Consider a block of bits `B = b_0 b_1 \cdots b_{n-1}` where each
    sub-block `b_i` is a binary string of length 8. Then the total number
    of bits is a multiple of 8 and is given by `8n`. Let `c_i` be the
    integer representation of `b_i`. We can consider `c_i` as the integer
    representation of an ASCII character. Then the ASCII representation
    `A` of `B` is `A = a_0 a_1 \cdots a_{n-1}`.

    EXAMPLES:

    Convert some ASCII strings to their binary representations and recover
    the ASCII strings from the binary representations::

        sage: from sage.crypto.util import ascii_to_bin
        sage: from sage.crypto.util import bin_to_ascii
        sage: A = "Abc"
        sage: B = ascii_to_bin(A); B
        010000010110001001100011
        sage: bin_to_ascii(B)
        'Abc'
        sage: bin_to_ascii(B) == A
        True

    ::

        sage: A = "123 \" #"
        sage: B = ascii_to_bin(A); B
        00110001001100100011001100100000001000100010000000100011
        sage: bin_to_ascii(B)
        '123 " #'
        sage: bin_to_ascii(B) == A
        True

    This function also accepts strings and lists of bits::

        sage: from sage.crypto.util import bin_to_ascii
        sage: bin_to_ascii("010000010110001001100011")
        'Abc'
        sage: bin_to_ascii([0, 1, 0, 0, 0, 0, 0, 1])
        'A'

    TESTS:

    The number of bits in ``B`` must be non-empty and a multiple of 8::

        sage: from sage.crypto.util import bin_to_ascii
        sage: bin_to_ascii("")
        Traceback (most recent call last):
        ...
        ValueError: B must be a non-empty binary string.
        sage: bin_to_ascii([])
        Traceback (most recent call last):
        ...
        ValueError: B must be a non-empty binary string.
        sage: bin_to_ascii(" ")
        Traceback (most recent call last):
        ...
        ValueError: The number of bits in B must be a multiple of 8.
        sage: bin_to_ascii("101")
        Traceback (most recent call last):
        ...
        ValueError: The number of bits in B must be a multiple of 8.
        sage: bin_to_ascii([1, 0, 1])
        Traceback (most recent call last):
        ...
        ValueError: The number of bits in B must be a multiple of 8.
    """
    # sanity checks
    n = len(B)
    if n == 0:
        raise ValueError("B must be a non-empty binary string.")
    if mod(n, 8) != 0:
        raise ValueError("The number of bits in B must be a multiple of 8.")
    # perform conversion to ASCII string
    b = [int(str(x)) for x in list(B)]
    A = []
    # the number of 8-bit blocks
    k = n // 8
    for i in range(k):
        # Convert from 8-bit string to ASCII integer. Then convert the
        # ASCII integer to the corresponding ASCII character.
        A.append(chr(ascii_integer(b[8*i: 8*(i+1)])))
    return "".join(A)


def has_blum_prime(lbound, ubound):
    r"""
    Determine whether or not there is a Blum prime within the specified closed
    interval.

    INPUT:

    - ``lbound`` -- positive integer; the lower bound on how small a
      Blum prime can be. The lower bound must be distinct from the upper
      bound.

    - ``ubound`` -- positive integer; the upper bound on how large a
      Blum prime can be. The lower bound must be distinct from the upper
      bound.

    OUTPUT:

    - ``True`` if there is a Blum prime ``p`` such that
      ``lbound <= p <= ubound``. ``False`` otherwise.

    ALGORITHM:

    Let `L` and `U` be distinct positive integers. Let `P` be the set of all
    odd primes `p` such that `L \leq p \leq U`. Our main focus is on Blum
    primes, i.e. odd primes that are congruent to 3 modulo 4, so we assume
    that the lower bound `L > 2`. The closed interval `[L, U]` has a Blum
    prime if and only if the set `P` has a Blum prime.

    EXAMPLES:

    Testing for the presence of Blum primes within some closed intervals.
    The interval `[4, 100]` has a Blum prime, the smallest such prime being
    7. The interval `[24, 28]` has no primes, hence no Blum primes. ::

        sage: from sage.crypto.util import has_blum_prime
        sage: from sage.crypto.util import is_blum_prime
        sage: has_blum_prime(4, 100)                                                    # needs sage.libs.pari
        True
        sage: for n in range(4, 100):
        ....:     if is_blum_prime(n):
        ....:         print(n)
        ....:         break
        7
        sage: has_blum_prime(24, 28)                                                    # needs sage.libs.pari
        False

    TESTS:

    Both the lower and upper bounds must be greater than 2::

        sage: from sage.crypto.util import has_blum_prime
        sage: has_blum_prime(2, 3)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.
        sage: has_blum_prime(3, 2)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.
        sage: has_blum_prime(2, 2)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.

    The lower and upper bounds must be distinct from each other::

        sage: has_blum_prime(3, 3)
        Traceback (most recent call last):
        ...
        ValueError: The lower and upper bounds must be distinct.

    The lower bound must be less than the upper bound::

        sage: has_blum_prime(4, 3)
        Traceback (most recent call last):
        ...
        ValueError: The lower bound must be less than the upper bound.
    """
    # sanity checks
    if (lbound < 3) or (ubound < 3):
        raise ValueError("Both the lower and upper bounds must be > 2.")
    if lbound == ubound:
        raise ValueError("The lower and upper bounds must be distinct.")
    if lbound > ubound:
        raise ValueError("The lower bound must be less than the upper bound.")
    # now test for presence of a Blum prime
    for p in primes(lbound, ubound + 1):
        if mod(p, 4).lift() == 3:
            return True
    return False

def is_blum_prime(n):
    r"""
    Determine whether or not ``n`` is a Blum prime.

    INPUT:

    - ``n`` -- a positive prime

    OUTPUT: ``True`` if ``n`` is a Blum prime; ``False`` otherwise

    Let `n` be a positive prime. Then `n` is a Blum prime if `n` is
    congruent to 3 modulo 4, i.e. `n \equiv 3 \pmod{4}`.

    EXAMPLES:

    Testing some integers to see if they are Blum primes::

        sage: from sage.crypto.util import is_blum_prime
        sage: from sage.crypto.util import random_blum_prime
        sage: is_blum_prime(101)
        False
        sage: is_blum_prime(7)
        True
        sage: p = random_blum_prime(10**3, 10**5)                                       # needs sage.libs.pari
        sage: is_blum_prime(p)                                                          # needs sage.libs.pari
        True
    """
    if n < 0:
        return False
    if is_prime(n):
        if mod(n, 4).lift() == 3:
            return True
        else:
            return False
    else:
        return False

def least_significant_bits(n, k):
    r"""
    Return the ``k`` least significant bits of ``n``.

    INPUT:

    - ``n`` -- integer

    - ``k`` -- positive integer

    OUTPUT:

    - The ``k`` least significant bits of the integer ``n``. If ``k=1``,
      then return the parity bit of the integer ``n``. Let `b` be the
      binary representation of ``n``, where `m` is the length of the
      binary string `b`. If `k \geq m`, then return the binary
      representation of ``n``.

    EXAMPLES:

    Obtain the parity bits of some integers::

        sage: from sage.crypto.util import least_significant_bits
        sage: least_significant_bits(0, 1)
        [0]
        sage: least_significant_bits(2, 1)
        [0]
        sage: least_significant_bits(3, 1)
        [1]
        sage: least_significant_bits(-2, 1)
        [0]
        sage: least_significant_bits(-3, 1)
        [1]

    Obtain the 4 least significant bits of some integers::

        sage: least_significant_bits(101, 4)
        [0, 1, 0, 1]
        sage: least_significant_bits(-101, 4)
        [0, 1, 0, 1]
        sage: least_significant_bits(124, 4)
        [1, 1, 0, 0]
        sage: least_significant_bits(-124, 4)
        [1, 1, 0, 0]

    The binary representation of 123::

        sage: n = 123; b = n.binary(); b
        '1111011'
        sage: least_significant_bits(n, len(b))
        [1, 1, 1, 1, 0, 1, 1]
    """
    return [int(_) for _ in list(n.binary()[-k:])]

def random_blum_prime(lbound, ubound, ntries=100):
    r"""
    A random Blum prime within the specified bounds.

    Let `p` be a positive prime. Then `p` is a Blum prime if `p` is
    congruent to 3 modulo 4, i.e. `p \equiv 3 \pmod{4}`.

    INPUT:

    - ``lbound`` -- positive integer; the lower bound on how small a
      random Blum prime `p` can be. So we have
      ``0 < lbound <= p <= ubound``. The lower bound must be distinct from
      the upper bound.

    - ``ubound`` -- positive integer; the upper bound on how large a
      random Blum prime `p` can be. So we have
      ``0 < lbound <= p <= ubound``. The lower bound must be distinct
      from the upper bound.

    - ``ntries`` -- (default: ``100``) the number of attempts to generate
      a random Blum prime. If ``ntries`` is a positive integer, then
      perform that many attempts at generating a random Blum prime. This
      might or might not result in a Blum prime.

    OUTPUT: a random Blum prime within the specified lower and upper bounds

    .. NOTE::

        Beware that there might not be any primes between the lower and
        upper bounds. So make sure that these two bounds are
        "sufficiently" far apart from each other for there to be primes
        congruent to 3 modulo 4. In particular, there should be at least
        two distinct Blum primes within the specified bounds.

    EXAMPLES:

    Choose a random prime and check that it is a Blum prime::

        sage: from sage.crypto.util import random_blum_prime
        sage: p = random_blum_prime(10**4, 10**5)                                       # needs sage.libs.pari
        sage: is_prime(p)                                                               # needs sage.libs.pari
        True
        sage: mod(p, 4) == 3                                                            # needs sage.libs.pari
        True

    TESTS:

    Make sure that there is at least one Blum prime between the lower and
    upper bounds. In the following example, we have ``lbound=24`` and
    ``ubound=30`` with 29 being the only prime within those bounds. But 29
    is not a Blum prime. ::

        sage: from sage.crypto.util import random_blum_prime
        sage: random_blum_prime(24, 30, ntries=10)                                      # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: No Blum primes within the specified closed interval.
        sage: random_blum_prime(24, 28)                                                 # needs sage.libs.pari
        Traceback (most recent call last):
        ...
        ValueError: No Blum primes within the specified closed interval.
    """
    # sanity check
    if not has_blum_prime(lbound, ubound):
        raise ValueError("No Blum primes within the specified closed interval.")
    # Now we know that there is a Blum prime within the closed interval
    # [lbound, ubound]. Pick one such prime at random.
    p = random_prime(ubound, lbound=lbound, proof=True)
    n = 1
    while mod(p, 4) != 3:
        p = random_prime(ubound, lbound=lbound, proof=True)
        n += 1
        if n > ntries:
            raise ValueError("Maximum number of attempts exceeded.")
    return p
