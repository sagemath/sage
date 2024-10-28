r"""
Key Exchange Schemes

This module contains base classes for key exchange schemes. The classes defined
in this module should not be called directly. It is the responsibility of child
classes to implement specific key exchange schemes.


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
from sage.structure.sage_object import SageObject


class KeyExchangeScheme(SageObject):
    """
    Abstract base class for key exchange schemes.

    Currently experimental and subject to change.
    """

    @experimental(37305)
    def __init__(self):
        """
        Create a ``KeyExchangeScheme`` instance.

        TESTS::

            sage: from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme
            sage: K = KeyExchangeScheme()
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/37305 for details.
        """
        pass

    def generate_secret_key(self):
        """
        Generate a secret key.

        TESTS::

            sage: from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme
            sage: K = KeyExchangeScheme()
            sage: K.generate_secret_key()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def generate_public_key(self, secret_key):
        """
        Generate a public key using the given secret key.

        TESTS::

            sage: from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme
            sage: K = KeyExchangeScheme()
            sage: K.generate_public_key(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def compute_shared_secret(self, alice_pk, bob_sk):
        """
        Compute the shared secret using the given public key and secret keys.

        TESTS::

            sage: from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme
            sage: K = KeyExchangeScheme()
            sage: K.compute_shared_secret(None, None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def parameters(self):
        """
        Get the parameters for the ``KeyExchangeScheme`` instance.

        TESTS::

            sage: from sage.crypto.key_exchange.key_exchange_scheme import KeyExchangeScheme
            sage: K = KeyExchangeScheme()
            sage: K.parameters()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
