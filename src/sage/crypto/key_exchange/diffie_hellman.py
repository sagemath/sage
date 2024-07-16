r"""
Diffie-Hellman Public Key Exchange Scheme

This module contains a toy implementation of the Diffie-Hellman public key exchange
scheme.

AUTHORS:

- Vincent Macri (2024-07-15): initial version
"""
# ****************************************************************************
#       Copyright (C) 2024 Vincent Macri <vincent.macri@ucalgary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.superseded import experimental

from sage.crypto.key_exchange.pke import KeyExchangeScheme

from sage.arith.misc import is_prime
from sage.misc.prandom import randint
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_rings.integer_mod import IntegerMod_int
from sage.structure.sage_object import SageObject

from typing import Union

class DiffieHellman(KeyExchangeScheme):

    @experimental(37305)
    def __init__(self, p: Integer, g: Union[Integer, IntegerMod_int]):
        """
        Create an instance of the Diffie-Hellman key exchange scheme using the
        given prime ``p`` and base ``g``.

        INPUT:

        - ``p`` -- prime integer defining the field `\\GF{p}`` that the key exchanges
          will be performed over, must be at least 5

        - ``g`` -- base for the key exchange, (coerceable to) an element of `\\GF{p}` from `2` to `p - 2`

        .. WARNING::

            This is a toy implementation for educational use only! Do not use this
            implementation, or any cryptographic features of Sage, in any setting where
            security is needed!

        REFERENCES:

        For more information, see [PP2010]_, section 8.1.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(13, 2)
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/37305 for details.

        TESTS::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(3, 2)
            Traceback (most recent call last):
            ...
            ValueError: p must be at least 5

            sage: DH = DiffieHellman(5, 0)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)

            sage: DH = DiffieHellman(5, 1)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)

            sage: DH = DiffieHellman(5, 4)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)
        """

        if p < 5:
            raise ValueError('p must be at least 5')

        # The modn implementation takes care of checking that ``p`` is prime
        self._field = GF(p, impl='modn')

        self._p = p
        self._g = self._field(g)

        # While these values won't cause mathematical problems, they do completely
        # break the security of the Diffie-Hellman scheme.
        # g = 0 makes every secret key and shared secret 0
        # g = 1 makes every secret key and shared secret 1
        # g = -1 makes every secret key and shared secret 1 or -1
        if self._g == 0 or self._g == 1 or self._g == p - 1:
            raise ValueError('g cannot be 0, 1, or p - 1 (mod p)')

    def field(self) -> FiniteField_prime_modn:
        """
        Return the field this Diffie-Hellman instance is working over.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(5, 2)
            sage: DH.field()
            Finite Field of size 5
        """
        return self._field

    def prime(self) -> Integer:
        """
        Return the prime ``p`` for this Diffie-Hellman instance.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(7, 3)
            sage: DH.prime()
            7
        """
        return self._p

    def generator(self) -> IntegerMod_int:
        """
        Return the generator ``g`` for this Diffie-Hellman instance.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(7, 3)
            sage: DH.generator()
            3
        """
        return self._g

    def parameters(self) -> tuple[Integer, IntegerMod_int]:
        """
        Output the parameters ``(p, g)`` for this Diffie-Hellman instance.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(7, 3)
            sage: DH.parameters()
            (7, 3)
        """
        return (self._p, self._g)

    def generate_secret_key(self) -> Integer:
        """
        Generate a random Diffie-Hellman secret key.

        TESTS:

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(7, 2)
            sage: keys = [DH.generate_secret_key() for i in range(10)]
            sage: all(2 <= i <= 5 for i in keys)
            True
        """
        return randint(2, self._p - 2)

    def generate_public_key(self, secret_key: Integer) -> IntegerMod_int:
        """
        Generate a Diffie-Hellman public key using the given secret key.

        INPUT:

        - ``secret_key`` -- the secret key to generate the public key with

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(13, 2)
            sage: DH.generate_public_key(4)
            3
        """
        return self._g**secret_key

    def compute_shared_secret(self, alice_pk: IntegerMod_int, bob_sk: Integer) -> IntegerMod_int:
        """
        Compute Alice and Bob's shared secret using Alice's public key and
        Bob's secret key.

        INPUT:

        - ``alice_pk`` -- Alice's public key

        - ``bob_sk`` -- Bob's secret key

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(17, 3)
            sage: DH.compute_shared_secret(13, 11)
            4
        """
        return self._field(alice_pk**bob_sk)

    def subgroup_size(self) -> Integer:
        """
        Calculates the size of the subgroup of `\\GF{p}` generated by ``self.generator()``.

        EXAMPLES:

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(47, 2)
            sage: DH.subgroup_size()
            23

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH = DiffieHellman(47, 5)
            sage: DH.subgroup_size()
            46
        """
        return self._g.multiplicative_order()

    def __len__(self):
        return int(self.subgroup_size())

    def __eq__(self, other):
        """
        Check if two Diffie-Hellman instances have the same parameter set.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH1 = DiffieHellman(5, 2)
            sage: DH2 = DiffieHellman(5, 2)
            sage: DH1 == DH2
            True
        """
        return self.parameters() == other.parameters()

    def __hash__(self):
        """
        Compute the hash value of a Diffie-Hellman instance.

        EXAMPLES::

            sage: from sage.crypto.key_exchange.diffie_hellman import DiffieHellman
            sage: DH1 = DiffieHellman(7, 3)
            sage: DH2 = DiffieHellman(7, 3)
            sage: s = set([DH1, DH2])
            sage: len(s)
            1
        """
        return hash((self._p, self._g))

    def _repr_(self):
        return f'Diffie-Hellman key exchange over {self._field} with generator {self._g}'

    def _latex_(self):
        return f'\\text{{Diffie-Hellman key exchange over }}{self._field}\\text{{ with generator }}{self._g}'
