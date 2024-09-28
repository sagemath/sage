r"""
Diffie-Hellman Key Exchange Scheme

This module contains a toy implementation of the Diffie-Hellman key exchange
scheme.

AUTHORS:

- Vincent Macri (2024-07-30): initial version
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

from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme

from sage.arith.misc import is_prime
from sage.misc.prandom import randint
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.finite_field_prime_modn import \
    FiniteField_prime_modn
from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
from sage.structure.proof.proof import WithProof

from typing import Union

class DiffieHellman(KeyExchangeScheme):

    @experimental(37305)
    def __init__(self, p: Integer, g: Union[Integer, IntegerMod_abstract],
                 proof: bool = True) -> None:
        r"""
        Create an instance of the Diffie-Hellman key exchange scheme using the
        given prime ``p`` and base ``g``.

        INPUT:

        - ``p`` -- prime integer defining the field `\GF{p}` that the key
          exchanges will be performed over, must be at least 5

        - ``g`` -- base for the key exchange, (coerceable to) an element of
          `\GF{p}` from `2` to `p - 2`

        - ``proof`` -- (default: ``True``) whether to require a proof that
          ``p`` is prime. If ``False``, a probabilistic test can be used for
          checking that ``p`` is prime. This should be set to ``False``
          when using large (cryptographic size) primes, otherwise checking
          primality will take too long.

        .. WARNING::

            This is a toy implementation for educational use only! Do not use
            this implementation, or any cryptographic features of Sage, in any
            setting where security is needed!

        REFERENCES:

        For more information, see Section 8.1 of [PP2010]_.

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(13, 2)
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/37305 for details.

        This is an example of a full key exchange using a cryptographically
        large prime. This is the prime from the 8192-bit MODP group in RFC 3526
        (see [KK2003]_)::

            sage: p = 2^8192 - 2^8128 - 1 + 2^64 * (round(2^8062 * pi) + 4743158)
            sage: DH = key_exchange.DiffieHellman(p, 2, proof=False)
            sage: alice_sk = DH.generate_secret_key()
            sage: alice_pk = DH.generate_public_key(alice_sk)
            sage: bob_sk = DH.generate_secret_key()
            sage: bob_pk = DH.generate_public_key(bob_sk)
            sage: alice_shared_secret = DH.compute_shared_secret(bob_pk, alice_sk)
            sage: bob_shared_secret = DH.compute_shared_secret(alice_pk, bob_sk)
            sage: alice_shared_secret == bob_shared_secret
            True

        TESTS::

            sage: DH = key_exchange.DiffieHellman(3, 2)
            Traceback (most recent call last):
            ...
            ValueError: p must be at least 5

            sage: DH = key_exchange.DiffieHellman(5, 0)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)

            sage: DH = key_exchange.DiffieHellman(5, 1)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)

            sage: DH = key_exchange.DiffieHellman(5, 4)
            Traceback (most recent call last):
            ...
            ValueError: g cannot be 0, 1, or p - 1 (mod p)
        """

        if p < 5:
            raise ValueError('p must be at least 5')

        if proof:
            # The modn implementation checks that ``p`` is prime
            self._field = GF(p, impl='modn')
        else:
            with WithProof('arithmetic', False):
                self._field = GF(p, impl='modn')

        self._p = p
        self._g = self._field(g)

        # While these values won't cause mathematical problems, they do
        # completely break the security of the Diffie-Hellman scheme.
        # g = 0 makes every secret key and shared secret 0
        # g = 1 makes every secret key and shared secret 1
        # g = -1 makes every secret key and shared secret 1 or -1
        if self._g == 0 or self._g == 1 or self._g == p - 1:
            raise ValueError('g cannot be 0, 1, or p - 1 (mod p)')

    def field(self) -> FiniteField_prime_modn:
        """
        Return the field this ``DiffieHellman`` instance is working over.

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(5, 2)
            sage: DH.field()
            Finite Field of size 5
        """
        return self._field

    def prime(self) -> Integer:
        """
        Return the prime ``p`` for this ``DiffieHellman`` instance.

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(7, 3)
            sage: DH.prime()
            7
        """
        return self._p

    def generator(self) -> IntegerMod_abstract:
        """
        Return the generator ``g`` for this ``DiffieHellman`` instance.

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(7, 3)
            sage: DH.generator()
            3
        """
        return self._g

    def parameters(self) -> tuple[Integer, IntegerMod_abstract]:
        """
        Get the parameters ``(p, g)`` for this ``DiffieHellman`` instance.

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(7, 3)
            sage: DH.parameters()
            (7, 3)
        """
        return (self._p, self._g)

    def generate_secret_key(self) -> Integer:
        """
        Generate a random Diffie-Hellman secret key.

        TESTS:

            sage: DH = key_exchange.DiffieHellman(7, 2)
            sage: keys = [DH.generate_secret_key() for i in range(10)]
            sage: all(2 <= i <= 5 for i in keys)
            True
        """
        return randint(2, self._p - 2)

    def generate_public_key(self, secret_key: Integer) -> IntegerMod_abstract:
        """
        Generate a Diffie-Hellman public key using the given secret key.

        INPUT:

        - ``secret_key`` -- the secret key to generate the public key with

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(13, 2)
            sage: DH.generate_public_key(4)
            3
        """
        return self._g**secret_key

    def compute_shared_secret(self, pk: IntegerMod_abstract,
                              sk: Integer) -> IntegerMod_abstract:
        """
        Compute the shared secret using the given public key and secret keys.

        INPUT:

        - ``pk`` -- public key

        - ``sk`` -- secret key

        EXAMPLES::

            sage: DH = key_exchange.DiffieHellman(17, 3)
            sage: DH.compute_shared_secret(13, 11)
            4
        """
        return self._field(pk**sk)

    def subgroup_size(self) -> Integer:
        """
        Calculates the size of the subgroup of `\\GF{p}` generated by
        ``self.generator()``.

        EXAMPLES:

        This is an example of a ``DiffieHellman`` instance where the subgroup
        size is `(p - 1) / 2`::

            sage: DH = key_exchange.DiffieHellman(47, 2)
            sage: DH.subgroup_size()
            23

        This is an example of a ``DiffieHellman`` instance where the subgroup
        size is `p - 1`::

            sage: DH = key_exchange.DiffieHellman(47, 5)
            sage: DH.subgroup_size()
            46
        """
        return self._g.multiplicative_order()

    def __len__(self) -> int:
        """
        Calculates the size of the subgroup of `\\GF{p}` generated by
        ``self.generator()``. This is a wrapper around `subgroup_size`.

        TESTS::

            sage: DH = key_exchange.DiffieHellman(53, 9)
            sage: len(DH)
            26
        """
        return int(self.subgroup_size())

    def __eq__(self, other) -> bool:
        """
        Check if two ``DiffieHellman`` instances have the same parameter set.

        TESTS::

            sage: DH1 = key_exchange.DiffieHellman(5, 2)
            sage: DH2 = key_exchange.DiffieHellman(5, 2)
            sage: DH1 == DH2
            True
            sage: DH1 == 5
            False
        """
        if isinstance(other, DiffieHellman):
            return self.parameters() == other.parameters()
        return False

    def __hash__(self) -> int:
        """
        Compute the hash value of a ``DiffieHellman`` instance.

        TESTS::

            sage: DH1 = key_exchange.DiffieHellman(7, 3)
            sage: DH2 = key_exchange.DiffieHellman(7, 3)
            sage: s = set([DH1, DH2])
            sage: len(s)
            1
        """
        return hash((self._p, self._g))

    def _repr_(self) -> str:
        """
        Get the string representation of the ``DiffieHellman`` instance.

        TESTS::

            sage: DH = key_exchange.DiffieHellman(7, 3)
            sage: DH
            Diffie-Hellman key exchange over Finite Field of size 7 with generator 3
        """
        return ('Diffie-Hellman key exchange over '
                f'{self._field} '
                'with generator '
                f'{self._g}')

    def _latex_(self) -> str:
        r"""
        Get the LaTeX representation of the ``DiffieHellman`` instance.

        TESTS::

            sage: DH = key_exchange.DiffieHellman(7, 3)
            sage: latex(DH)
            \text{Diffie-Hellman key exchange over }\Bold{F}_{7}\text{ with generator }3
        """
        return ('\\text{Diffie-Hellman key exchange over }'
                f'{self._field._latex_()}'
                '\\text{ with generator }'
                f'{self._g}')
