r"""
Base Classes for Key Exchange Schemes

This module contains base classes for key exchange schemes. The classes defined
in this module should not be initialized directly. It is the responsibility of
child classes to implement specific key exchange schemes.

A key exchange protocol establishes a shared secret value between two parties,
Alice and Bob. Either party is able to initiate the key exchange, in the sense
that either party can compute the shared secret using only their own private key
and the other party's public key.

AUTHORS:

- Brian Heckel (2025-11-26): initial version
- Vincent Macri (2025-12-18): add named_parameter_set method
"""

# ****************************************************************************
#       Copyright (C) 2025 Brian Heckel <heckelbri@gmail.com>
#       Copyright (C) 2025 Vincent Macri <vincent.macri@ucalgary.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from abc import abstractmethod
from typing import Any, Self

from sage.misc.superseded import experimental_warning
from sage.structure.sage_object import SageObject

experimental_warning(
    41218,
    "SageMath's key exchange functionality is experimental and might change in the future.",
)


class KeyExchangeBase(SageObject):
    r"""
    A base class for key exchange schemes.

    Implementers of this class must implement all abstract methods
    defined in :meth:`KeyExchangeBase`.

    If all ``alice`` methods are the same as ``bob`` methods,
    then the :class:`CommutativeKeyExchangeBase` might be easier to implement.
    """

    @abstractmethod
    def alice_secret_key(self):
        r"""
        Generate a valid secret key for Alice.
        """
        raise NotImplementedError

    @abstractmethod
    def alice_public_key(self, alice_secret_key):
        r"""
        Generate a valid public key for Alice.

        INPUT:

        - ``alice_secret_key`` -- Alice's secret key that will be used to generate
            the public key

        OUTPUT:

        A valid public key that will be sent to Bob
        """
        raise NotImplementedError

    @abstractmethod
    def bob_secret_key(self):
        r"""
        Generate a valid secret key for Bob.
        """
        raise NotImplementedError

    @abstractmethod
    def bob_public_key(self, bob_secret_key):
        r"""
        Generate a valid public key for Bob.

        INPUT:

        - ``bob_secret_key`` -- Bob's secret key that will be used to generate
            the public key

        OUTPUT:

        A valid public key that will be sent to Alice
        """
        raise NotImplementedError

    @abstractmethod
    def alice_compute_shared_secret(self, alice_secret_key, bob_public_key):
        r"""
        Compute Alice's shared secret.

        INPUT:

        - ``alice_secret_key`` -- Alice's secret key that is kept secret from all parties

        - ``bob_public_key`` -- Bob's public key that has been sent to Alice

        OUTPUT:

        A secret key that is shared between Alice and Bob
        """
        raise NotImplementedError

    @abstractmethod
    def bob_compute_shared_secret(self, bob_secret_key, alice_public_key) -> Any:
        r"""
        Compute Bob's shared secret.

        INPUT:

        - ``bob_secret_key`` -- Bob's secret key that is kept secret from all parties

        - ``alice_public_key`` -- Alice's public key that has been sent to Bob

        OUTPUT:

        The secret key that is shared between Alice and Bob
        """
        raise NotImplementedError

    @abstractmethod
    def parameters(self) -> tuple:
        """
        A tuple of the public parameters of the key exchange.

        :meth:`parameters` should a tuple of useful attributes of the instance
        which are sufficient to define the parameter set of the key exchange
        scheme that the instance represents. For some implementations this
        may simply be the parameters passed to ``__init__`` when the object was
        constructed. For some implementations we may wish to return additional
        information for convenience. For example, a key exchange scheme that
        works over an elliptic curve over a finite field may wish to return the
        characteristic of the finite field in addition to the elliptic curve, even
        though the finite field can be accessed via methods on elliptic curve objects.

        The default implementations of ``_eq_`` and ``__hash__`` for
        :class:`KeyExchangeBase` are implementing using :meth:`parameters`.
        Hence two key exchange instances ``a`` and ``b`` compare as equal
        if and only if ``a.parameters()`` == ``b.parameters()``. Similarly,
        a key exchange instance ``a`` is hashable if and only if ``a.parameters()``
        is hashable. This is a reasonable default that should work for most key
        exchange schemes, but user classes can override the ``_eq_`` and ``__hash__``
        methods if this is not desirable.

        OUTPUT:

        A tuple of public parameters used for the key exchange
        """
        raise NotImplementedError

    def alice_key_pair(self) -> tuple[Any, Any]:
        r"""
        Generate a valid (secret key, public key) key pair for Alice.

        OUTPUT:

        A 2-tuple (secret_key, public_key) which is Alice's
        secret and public keys
        """
        alice_sk = self.alice_secret_key()
        alice_pk = self.alice_public_key(alice_sk)
        return (alice_sk, alice_pk)

    def bob_key_pair(self) -> tuple[Any, Any]:
        r"""
        Generate a valid (secret key, public key) key pair for Bob.

        OUTPUT:

        A 2-tuple (secret_key, public_key) which is Bob's
        secret and public keys
        """
        bob_sk = self.bob_secret_key()
        bob_pk = self.bob_public_key(bob_sk)
        return (bob_sk, bob_pk)

    def do_key_exchange(self) -> tuple[Any, Any, Any, Any, Any]:
        r"""
        Do a full key exchange and returns all public keys, secret keys,
        and the computed shared secret between Alice and Bob. Raises
        an ``AssertException`` if Alice and Bob do not compute the same
        shared secret.

        OUTPUT:

        A 5-tuple ``(alice_secret_key, alice_public_key, bob_secret_key, bob_public_key, shared_secret)``
        """
        alice_sk, alice_pk = self.alice_key_pair()
        bob_sk, bob_pk = self.bob_key_pair()
        alice_shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
        if alice_shared_secret != self.bob_compute_shared_secret(bob_sk, alice_pk):
            raise RuntimeError('Alice and Bob did not arrive at the same shared secret')
        return (alice_sk, alice_pk, bob_sk, bob_pk, alice_shared_secret)

    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Convenience method to easily construct particular instances of a key exchange scheme
        for actual parameter sets that are used in practice and have names. Implementations
        may also wish to implement a parameter set named 'toy' of a size that is just large
        enough to be non-trivial but is nowhere near cryptographic size. Implementations may
        also wish to set a custom name on the key exchange instance before returning it.

        Sage library implementations of key exchange schemes should define a 'toy' implementation
        and use it for most tests to reduce testing time.
        """
        raise ValueError(f'Unknown parameter set name "{name}" for {cls}')

    def _repr_(self) -> str:
        return f'{type(self).__name__} with parameter set: {self.parameters()}'

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.parameters() == other.parameters()

    def __hash__(self) -> int:
        return hash(self.parameters())

    def _test_key_exchange(self, **options) -> None:
        r"""
        Test that the key exchange generates the same shared secret for both parties.
        """
        tester = self._tester(**options)
        alice_sk, alice_pk = self.alice_key_pair()
        bob_sk, bob_pk = self.bob_key_pair()
        alice_shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
        bob_shared_secret = self.bob_compute_shared_secret(bob_sk, alice_pk)
        tester.assertEqual(alice_shared_secret, bob_shared_secret)


class CommutativeKeyExchangeBase(KeyExchangeBase):
    r"""
    A base class for key exchange schemes such as Diffie-Hellman where Alice
    and Bob perform the same computations for generating public/secret keys and
    the shared secret key.

    Implementers of this class only need to implement the abstract methods
    defined in :class:`CommutativeKeyExchangeBase` and do not need to implement
    method defined in :class:`KeyExchangeBase`. This class is for convenience
    to reduce code duplication when implementing key exchange schemes where
    Alice and Bob perform the same calculations.
    """

    @abstractmethod
    def secret_key(self) -> Any:
        r"""
        Generate a secret key for the key exchange.
        """
        raise NotImplementedError

    @abstractmethod
    def public_key(self, secret_key) -> Any:
        r"""
        Generate a public key for the secret key that you have chosen.

        INPUT:

        - ``secret_key`` -- A secret key that has been chosen beforehand
        """
        raise NotImplementedError

    @abstractmethod
    def compute_shared_secret(self, secret_key, public_key) -> Any:
        r"""
        Generate the computed shared secret.

        INPUT:

        - ``secret_key`` -- A secret key that has been chosen beforehand
        - ``public_key`` -- A public key that has been sent to this party through
            an insecure channel

        OUTPUT:

        A shared secret key between the two parties
        """
        raise NotImplementedError

    def alice_secret_key(self) -> Any:
        r"""
        Alias of :meth:`secret_key` for compatibility with base class.

        :meta private:
        """
        return self.secret_key()

    def alice_public_key(self, alice_secret_key) -> Any:
        r"""
        Alias of :meth:`public_key` for compatibility with base class.

        :meta private:
        """
        return self.public_key(alice_secret_key)

    def bob_secret_key(self) -> Any:
        r"""
        Alias of :meth:`secret_key` for compatibility with base class.

        :meta private:
        """
        return self.secret_key()

    def bob_public_key(self, bob_secret_key) -> Any:
        r"""
        Alias of :meth:`public_key` for compatibility with base class.

        :meta private:
        """
        return self.public_key(bob_secret_key)

    def alice_compute_shared_secret(self, alice_sk, bob_pk) -> Any:
        """
        Alias of :meth:`compute_shared_secret` for compatibility with base class.

        :meta private:
        """
        return self.compute_shared_secret(alice_sk, bob_pk)

    def bob_compute_shared_secret(self, bob_sk, alice_pk) -> Any:
        r"""
        Alias of :meth:`compute_shared_secret` for compatibility with base class.

        :meta private:
        """
        return self.compute_shared_secret(bob_sk, alice_pk)
