# sage.doctest: needs sage.combinat sage.rings.finite_rings
"""
Stream Cryptosystems
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.misc import gcd, power_mod
from sage.crypto.cryptosystem import SymmetricKeyCryptosystem
from sage.crypto.stream_cipher import LFSRCipher, ShrinkingGeneratorCipher
from sage.crypto.util import random_blum_prime
from sage.monoids.string_monoid import BinaryStrings
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.finite_rings.integer_mod_ring import IntegerModFactory
from sage.rings.polynomial.polynomial_element import Polynomial


IntegerModRing = IntegerModFactory("IntegerModRing")

class LFSRCryptosystem(SymmetricKeyCryptosystem):
    """
    Linear feedback shift register cryptosystem class
    """
    def __init__(self, field=None):
        """
        Create a linear feedback shift cryptosystem.

        INPUT:

        - ``field`` -- (default: ``None``) string monoid over a binary alphabet

        EXAMPLES::

            sage: E = LFSRCryptosystem(FiniteField(2))
            sage: E
            LFSR cryptosystem over Finite Field of size 2

        TESTS::

            sage: E = LFSRCryptosystem(FiniteField(2))
            sage: E == loads(dumps(E))
            True

        TODO: Implement LFSR cryptosystem for arbitrary rings. The current
        implementation is limited to the finite field of 2 elements only
        because of the dependence on binary strings.
        """
        if field is None:
            field = FiniteField(2)
        if field.cardinality() != 2:
            raise NotImplementedError("Not yet implemented.")
        S = BinaryStrings()
        SymmetricKeyCryptosystem.__init__(self, S, S, None)
        self._field = field

    def __eq__(self, right):
        return type(self) is type(right) and self._field == right._field

    def __call__(self, key):
        """
        Create a LFSR cipher.

        INPUT:

        - ``key`` -- a polynomial and initial state of the LFSR
        """
        if not isinstance(key, (list, tuple)) and len(key) == 2:
            raise TypeError("Argument key (= %s) must be a list of tuple of length 2" % key)
        poly, IS = key
        if not isinstance(poly, Polynomial):
            raise TypeError("poly (= %s) must be a polynomial." % poly)
        if not isinstance(IS, (list, tuple)):
            raise TypeError("IS (= %s) must be an initial in the key space." % IS)
        if len(IS) != poly.degree():
            raise TypeError("The length of IS (= %s) must equal the degree of poly (= %s)" % (IS, poly))
        return LFSRCipher(self, poly, IS)

    def _repr_(self):
        r"""
        Return the string representation of this LFSR cryptosystem.

        EXAMPLES::

            sage: LFSRCryptosystem(FiniteField(2))
            LFSR cryptosystem over Finite Field of size 2
        """
        return "LFSR cryptosystem over %s" % self._field

    def encoding(self,M):
        S = self.cipher_domain()
        try:
            return S.encoding(M)
        except Exception:
            raise TypeError("Argument M = %s does not encode in the cipher domain" % M)

class ShrinkingGeneratorCryptosystem(SymmetricKeyCryptosystem):
    """
    Shrinking generator cryptosystem class
    """
    def __init__(self, field=None):
        """
        Create a shrinking generator cryptosystem.

        INPUT:

        - ``field`` -- (default: ``None``) string monoid over a binary alphabet

        EXAMPLES::

            sage: E = ShrinkingGeneratorCryptosystem()
            sage: E
            Shrinking generator cryptosystem over Finite Field of size 2
        """
        if field is None:
            field = FiniteField(2)
        if field.cardinality() != 2:
            raise NotImplementedError("Not yet implemented.")
        S = BinaryStrings()
        SymmetricKeyCryptosystem.__init__(self, S, S, None)
        self._field = field

    def __call__(self, key):
        """
        Create a Shrinking generator cipher.

        INPUT:

        - ``key`` -- list or tuple consisting of two LFSR ciphers (e1,e2)

        OUTPUT: the shrinking generator cipher with key stream generator e1
        and decimating cipher e2
        """
        if not isinstance(key, (list, tuple)) and len(key) == 2:
            raise TypeError("Argument key (= %s) must be a list of tuple of length 2" % key)
        e1 = key[0]
        e2 = key[1]
        if not isinstance(e1, LFSRCipher) or not isinstance(e2, LFSRCipher):
            raise TypeError("The key (= (%s,%s)) must be a tuple of two LFSR ciphers." % key)
        return ShrinkingGeneratorCipher(self, e1, e2)

    def _repr_(self):
        r"""
        Return the string representation of this shrinking generator
        cryptosystem.

        EXAMPLES::

            sage: ShrinkingGeneratorCryptosystem()
            Shrinking generator cryptosystem over Finite Field of size 2
        """
        return "Shrinking generator cryptosystem over %s" % self._field

    def encoding(self,M):
        S = self.cipher_domain()
        try:
            return S.encoding(M)
        except Exception:
            raise TypeError("Argument M = %s does not encode in the cipher domain" % M)

def blum_blum_shub(length, seed=None, p=None, q=None,
                   lbound=None, ubound=None, ntries=100):
    r"""
    The Blum-Blum-Shub (BBS) pseudorandom bit generator.

    See the original paper by Blum, Blum and Shub [BBS1986]_. The
    BBS algorithm is also discussed in section 5.5.2 of [MvOV1996]_.

    INPUT:

    - ``length`` -- positive integer; the number of bits in the output
      pseudorandom bit sequence

    - ``seed`` -- (default: ``None``) if `p` and `q` are Blum primes, then
      ``seed`` is a quadratic residue in the multiplicative group
      `(\ZZ/n\ZZ)^{\ast}` where `n = pq`. If ``seed=None``, then the function
      would generate its own random quadratic residue in `(\ZZ/n\ZZ)^{\ast}`.
      If you provide a value for ``seed``, then it is your responsibility to
      ensure that the seed is a quadratic residue in the multiplicative group
      `(\ZZ/n\ZZ)^{\ast}`.

    - ``p`` -- (default: ``None``) a large positive prime congruent to 3
      modulo 4. Both ``p`` and ``q`` must be distinct. If ``p=None``, then
      a value for ``p`` will be generated, where
      ``0 < lower_bound <= p <= upper_bound``.

    - ``q`` -- (default: ``None``) a large positive prime congruence to 3
      modulo 4. Both ``p`` and ``q`` must be distinct. If ``q=None``, then
      a value for ``q`` will be generated, where
      ``0 < lower_bound <= q <= upper_bound``.

    - ``lbound`` -- (positive integer, default: ``None``) the lower
      bound on how small each random primes `p` and `q` can be. So we
      have ``0 < lbound <= p, q <= ubound``. The lower bound must be
      distinct from the upper bound.

    - ``ubound`` -- (positive integer, default: ``None``) the upper
      bound on how large each random primes `p` and `q` can be. So we have
      ``0 < lbound <= p, q <= ubound``. The lower bound must be distinct
      from the upper bound.

    - ``ntries`` -- (default: ``100``) the number of attempts to generate
      a random Blum prime. If ``ntries`` is a positive integer, then
      perform that many attempts at generating a random Blum prime. This
      might or might not result in a Blum prime.

    OUTPUT: a pseudorandom bit sequence whose length is specified by ``length``

    Here is a common use case for this function. If you want this
    function to use pre-computed values for `p` and `q`, you should pass
    those pre-computed values to this function. In that case, you only need
    to specify values for ``length``, ``p`` and ``q``, and you do not need
    to worry about doing anything with the parameters ``lbound`` and
    ``ubound``. The pre-computed values `p` and `q` must be Blum primes.
    It is your responsibility to check that both `p` and `q` are Blum primes.

    Here is another common use case. If you want the function to generate
    its own values for `p` and `q`, you must specify the lower and upper
    bounds within which these two primes must lie. In that case, you must
    specify values for ``length``, ``lbound`` and ``ubound``, and you do
    not need to worry about values for the parameters ``p`` and ``q``. The
    parameter ``ntries`` is only relevant when you want this function to
    generate ``p`` and ``q``.

    .. NOTE::

        Beware that there might not be any primes between the lower and
        upper bounds. So make sure that these two bounds are
        "sufficiently" far apart from each other for there to be primes
        congruent to 3 modulo 4. In particular, there should be at least
        two distinct primes within these bounds, each prime being congruent
        to 3 modulo 4. This function uses the function
        :func:`random_blum_prime() <sage.crypto.util.random_blum_prime>` to
        generate random primes that are congruent to 3 modulo 4.

    ALGORITHM:

    The BBS algorithm as described below is adapted from the presentation
    in Algorithm 5.40, page 186 of [MvOV1996]_.

    #. Let `L` be the desired number of bits in the output bit sequence.
       That is, `L` is the desired length of the bit string.
    #. Let `p` and `q` be two large distinct primes, each congruent to 3
       modulo 4.
    #. Let `n = pq` be the product of `p` and `q`.
    #. Select a random seed value `s \in (\ZZ/n\ZZ)^{\ast}`, where
       `(\ZZ/n\ZZ)^{\ast}` is the multiplicative group of `\ZZ/n\ZZ`.
    #. Let `x_0 = s^2 \bmod n`.
    #. For `i` from 1 to `L`, do

       #. Let `x_i = x_{i-1}^2 \bmod n`.
       #. Let `z_i` be the least significant bit of `x_i`.

    #. The output pseudorandom bit sequence is `z_1, z_2, \dots, z_L`.

    EXAMPLES:

    A BBS pseudorandom bit sequence with a specified seed::

        sage: from sage.crypto.stream import blum_blum_shub
        sage: blum_blum_shub(length=6, seed=3, p=11, q=19)
        110000

    You could specify the length of the bit string, with given values for
    ``p`` and ``q``::

        sage: blum_blum_shub(length=6, p=11, q=19)  # random
        001011

    Or you could specify the length of the bit string, with given values for
    the lower and upper bounds::

        sage: blum_blum_shub(length=6, lbound=10**4, ubound=10**5)  # random
        110111

    Under some reasonable hypotheses, Blum-Blum-Shub [BBS1982]_
    sketch a proof that the period of the BBS stream cipher is equal to
    `\lambda(\lambda(n))`, where `\lambda(n)` is the Carmichael function of
    `n`. This is verified below in a few examples by using the function
    :func:`lfsr_connection_polynomial() <sage.crypto.lfsr.lfsr_connection_polynomial>`
    (written by Tim Brock) which computes the connection polynomial of a
    linear feedback shift register sequence. The degree of that polynomial
    is the period. ::

        sage: from sage.crypto.stream import blum_blum_shub
        sage: from sage.arith.misc import carmichael_lambda
        sage: carmichael_lambda(carmichael_lambda(7*11))
        4
        sage: s = [GF(2)(int(str(x))) for x in blum_blum_shub(60, p=7, q=11, seed=13)]
        sage: lfsr_connection_polynomial(s)
        x^3 + x^2 + x + 1
        sage: carmichael_lambda(carmichael_lambda(11*23))
        20
        sage: s = [GF(2)(int(str(x))) for x in blum_blum_shub(60, p=11, q=23, seed=13)]
        sage: lfsr_connection_polynomial(s)
        x^19 + x^18 + x^17 + x^16 + x^15 + x^14 + x^13 + x^12 + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1

    TESTS:

    Make sure that there is at least one Blum prime between the lower and
    upper bounds. In the following example, we have ``lbound=24`` and
    ``ubound=30`` with 29 being the only prime within those bounds. But 29
    is not a Blum prime. ::

        sage: from sage.crypto.stream import blum_blum_shub
        sage: blum_blum_shub(6, lbound=24, ubound=30, ntries=10)
        Traceback (most recent call last):
        ...
        ValueError: No Blum primes within the specified closed interval.

    Both the lower and upper bounds must be greater than 2::

        sage: blum_blum_shub(6, lbound=2, ubound=3)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.
        sage: blum_blum_shub(6, lbound=3, ubound=2)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.
        sage: blum_blum_shub(6, lbound=2, ubound=2)
        Traceback (most recent call last):
        ...
        ValueError: Both the lower and upper bounds must be > 2.

    The lower and upper bounds must be distinct from each other::

        sage: blum_blum_shub(6, lbound=3, ubound=3)
        Traceback (most recent call last):
        ...
        ValueError: The lower and upper bounds must be distinct.

    The lower bound must be less than the upper bound::

        sage: blum_blum_shub(6, lbound=4, ubound=3)
        Traceback (most recent call last):
        ...
        ValueError: The lower bound must be less than the upper bound.
    """
    # sanity checks
    if length < 0:
        raise ValueError("The length of the bit string must be positive.")
    if (p is None) and (p == q == lbound == ubound):
        raise ValueError("Either specify values for p and q, or specify values for the lower and upper bounds.")
    # Use pre-computed Blum primes. Both the parameters p and q are
    # assumed to be Blum primes. No attempts are made to ensure that they
    # are indeed Blum primes.
    randp = 0
    randq = 0
    if (p is not None) and (q is not None):
        randp = p
        randq = q
    # generate random Blum primes within specified bounds
    elif (lbound is not None) and (ubound is not None):
        randp = random_blum_prime(lbound, ubound, ntries=ntries)
        randq = random_blum_prime(lbound, ubound, ntries=ntries)
        while randp == randq:
            randq = random_blum_prime(lbound, ubound, ntries=ntries)
    # no pre-computed primes given, and no appropriate bounds given
    else:
        raise ValueError("Either specify values for p and q, or specify values for the lower and upper bounds.")
    # By now, we should have two distinct Blum primes.
    n = randp * randq
    # If no seed is provided, select a random seed.
    x0 = seed
    if seed is None:
        zmod = IntegerModRing(n)
        s = zmod.random_element().lift()
        while gcd(s, n) != 1:
            s = zmod.random_element().lift()
        x0 = power_mod(s, 2, n)
    # start generating pseudorandom bits
    z = []
    for i in range(length):
        x1 = power_mod(x0, 2, n)
        z.append(x1 % 2)
        x0 = x1
    bin = BinaryStrings()
    return bin(z)
