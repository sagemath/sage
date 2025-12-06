"""
Diffie-Hellman Key Exchange Scheme

This module contains a toy implementation of the Diffie-Hellman key exchange
scheme.

AUTHORS:

- Vincent Macri (2024-07-30): initial version
- Brian Heckel (2025-12-06): converted to inherit from KeyExchangeBase class
"""

# ****************************************************************************
#       Copyright (C) 2025 Vincent Macri <vincent.macri@ucalgary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.crypto.public_key.key_exchange.key_exchange_base import (
    CommutativeKeyExchangeBase,
)
from sage.misc.prandom import randint
from sage.misc.superseded import experimental
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer


class FiniteFieldDH(CommutativeKeyExchangeBase):
    @experimental(41218)
    def __init__(self, p: Integer, generator):
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

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: DH = FiniteFieldDH(13, 2)
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/41218 for details.

        This is an example of a full key exchange using a cryptographically
        large prime. This is the prime from the 8192-bit MODP group in RFC 3526
        (see [KK2003]_)::

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: p = 2^8192 - 2^8128 - 1 + 2^64 * (round(2^8062 * pi) + 4743158)
            sage: DH = FiniteFieldDH(13, 2)
            sage: alice_sk = DH.secret_key()
            sage: alice_pk = DH.public_key(alice_sk)
            sage: bob_sk = DH.secret_key()
            sage: bob_pk = DH.public_key(bob_sk)
            sage: alice_shared_secret = DH.compute_shared_secret(alice_sk, bob_pk)
            sage: bob_shared_secret = DH.compute_shared_secret(bob_sk, alice_pk)
            sage: alice_shared_secret == bob_shared_secret
            True

        TESTS::

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH

            sage: DH = FiniteFieldDH(5, 3)
        """
        self._field = GF(p, impl='modn')
        self.p = p
        self.generator = self._field(generator)

    def secret_key(self):
        """
        Generate a random Diffie-Hellman secret key.

        TESTS:

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: DH = FiniteFieldDH(7, 2)
            sage: keys = [DH.secret_key() for i in range(10)]
            sage: all(2 <= i <= 5 for i in keys)
            True
        """
        return self._field(randint(2, self.p - 2))

    def public_key(self, secret_key):
        """
        Generate a Diffie-Hellman public key using the given secret key.

        INPUT:

        - ``secret_key`` -- the secret key to generate the public key with

        EXAMPLES::

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: DH = FiniteFieldDH(13, 2)
            sage: DH.public_key(4)
            3
        """
        return self.generator ** secret_key

    def compute_shared_secret(self, secret_key, public_key):
        """
        Compute the shared secret using the given public key and secret keys.

        INPUT:

        - ``pk`` -- public key

        - ``sk`` -- secret key

        EXAMPLES::

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: DH = FiniteFieldDH(17, 3)
            sage: DH.compute_shared_secret(11, 13)
            4
        """
        return self._field(public_key ** secret_key)

    def parameters(self):
        """
        Returns a list of the prime ``p`` used, the base generator ``g`` and the field used.

        Examples::

            sage: from sage.crypto.public_key.key_exchange.finite_field_dh import FiniteFieldDH
            sage: DH = FiniteFieldDH(17, 3)
            sage: p = DH.parameters()
            sage: p[0]
            17
            sage: p[2]
            Finite Field of size 17
        """
        return [self.p, self.generator, self._field]
