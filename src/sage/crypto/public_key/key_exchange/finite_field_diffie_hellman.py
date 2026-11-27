r"""
Finite Field Diffie-Hellman

Toy implementation of Diffie-Hellman key exchange over finite fields `\Zmod{p}`.

.. WARNING::

    This is a toy implementation for educational use only! Do not use
    this implementation, or any cryptographic features of Sage, in any
    setting where security is needed!

AUTHORS:

- Vincent Macri (2024-07-30): initial version
- Brian Heckel (2025-12-06): converted to inherit from KeyExchangeBase class
"""

# ****************************************************************************
#       Copyright (C) 2024 Vincent Macri <vincent.macri@ucalgary.ca>
#       Copyright (C) 2025 Brian Heckel <heckelbri@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING

from sage.crypto.public_key.key_exchange.key_exchange_base import (
    CommutativeKeyExchangeBase,
)
from sage.misc.prandom import randint
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer

if TYPE_CHECKING:
    from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
    from sage.rings.finite_rings.integer_mod import IntegerMod_abstract


class FiniteFieldDiffieHellman(CommutativeKeyExchangeBase):
    def __init__(
        self,
        p: Integer | int,
        generator: Integer | IntegerMod_abstract | int,
        proof: bool | None = None,
    ) -> None:
        r"""
        Create an instance of the Diffie-Hellman key exchange scheme using the
        given prime ``p`` and base ``g``.

        INPUT:

        - ``p`` -- prime integer defining the field `\GF{p}` that the key
          exchanges will be performed over, must be at least 5

        - ``generator`` -- base for the key exchange, (coerceable to) an element of
          `\GF{p}` from `2` to `p - 2`

        - ``proof`` -- whether to require a proof that ``p`` is prime.
          If ``False``, a probabilistic test is used to check that ``p`` is
          prime. This should be set to ``False`` when using large
          (cryptographic size) primes, otherwise checking primality will take
          too long. If this is not specified, then the default behaviour is
          to use the current value of `proof.arithmetic()`.

        REFERENCES:

        For more information, see Section 8.1 of [PP2010]_.

        EXAMPLES::

            sage: DH = key_exchange.FiniteFieldDiffieHellman(13, 2)
            doctest:...: FutureWarning: SageMath's key exchange functionality is experimental and might change in the future.
                         See https://github.com/sagemath/sage/issues/41218 for details.

        This is an example of a full key exchange using a cryptographically
        large prime. This is the prime from the 8192-bit MODP group in RFC 3526
        (see [KK2003]_)::

            sage: p = 2^8192 - 2^8128 - 1 + 2^64 * (round(2^8062 * pi) + 4743158)
            sage: DH = key_exchange.FiniteFieldDiffieHellman(p, 2, proof=False)
            sage: alice_sk = DH.secret_key()
            sage: alice_pk = DH.public_key(alice_sk)
            sage: bob_sk = DH.secret_key()
            sage: bob_pk = DH.public_key(bob_sk)
            sage: alice_shared_secret = DH.compute_shared_secret(alice_sk, bob_pk)
            sage: bob_shared_secret = DH.compute_shared_secret(bob_sk, alice_pk)
            sage: alice_shared_secret == bob_shared_secret
            True

        TESTS::

            sage: DH = key_exchange.FiniteFieldDiffieHellman(5, 3)
            sage: TestSuite(DH).run()
        """
        if p < 5:
            raise ValueError('p must be at least 5')

        self._p = Integer(p)
        # The modn implementation checks that ``p`` is prime
        self._field = GF(self._p, impl='modn', proof=proof)

        self._generator = self._field(generator)

    def secret_key(self) -> IntegerMod_abstract:
        """
        Generate a random Diffie-Hellman secret key.

        TESTS:

            sage: DH = key_exchange.FiniteFieldDiffieHellman(7, 2)
            sage: keys = [DH.secret_key() for i in range(10)]
            sage: all(2 <= i <= 5 for i in keys)
            True
        """
        return self._field(randint(2, self._p - 2))

    def public_key(self, secret_key) -> IntegerMod_abstract:
        """
        Generate a Diffie-Hellman public key using the given secret key.

        INPUT:

        - ``secret_key`` -- the secret key to generate the public key with

        EXAMPLES::

            sage: DH = key_exchange.FiniteFieldDiffieHellman(13, 2)
            sage: DH.public_key(4)
            3
        """
        return self._generator**secret_key

    def compute_shared_secret(self, secret_key, public_key) -> IntegerMod_abstract:
        """
        Compute the shared secret using the given public key and secret keys.

        INPUT:

        - ``pk`` -- public key

        - ``sk`` -- secret key

        EXAMPLES::

            sage: DH = key_exchange.FiniteFieldDiffieHellman(17, 3)
            sage: DH.compute_shared_secret(11, 13)
            4
        """
        return self._field(public_key**secret_key)

    def parameters(self) -> tuple[FiniteField_prime_modn, IntegerMod_abstract]:
        """
        Returns a list of the prime ``p`` used, the base generator ``g`` and the field used.

        EXAMPLES::

            sage: DH = key_exchange.FiniteFieldDiffieHellman(17, 3)
            sage: DH.parameters()
            (Finite Field of size 17, 3)
        """
        return (self._field, self._generator)
