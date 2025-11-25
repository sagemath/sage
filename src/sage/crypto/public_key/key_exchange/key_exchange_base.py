from abc import ABC, abstractmethod
from typing import Any

from sage.structure.sage_object import SageObject


class KeyExchangeBase(SageObject, ABC):
    """
    An Base class for Public Key Exchange Schemes.

    Implementors of this class must give implementations
    of the following methods
    - ``alice_secret_key(self)``
    - ``bob_secret_key(self)``
    - ``alice_public_key(self, alice_secret_key)``
    - ``bob_public_key(self, bob_secret_key)``
    - ``alice_compute_shared_secret(self, alice_secret_key, bob_public_key)``
    - ``bob_compute_shared_secret(self, bob_secret_key, alice_public_key)``
    If all `alice` methods are the same as `bob` methods,
    then the `CommutativeKeyExchange` might be eaiser to implement.
    """

    @abstractmethod
    def alice_secret_key(self):
        """
        Generate a valid secret key for Alice
        """
        raise NotImplementedError

    @abstractmethod
    def alice_public_key(self, alice_secret_key):
        """
        Generate a valid public key for Alice

        INPUT:

        - ``alice_secret_key``: Alice's secret key that will be used to generate
            the public key.

        OUTPUT:

        - A valid public key that will be sent to Bob
        """
        raise NotImplementedError

    @abstractmethod
    def bob_secret_key(self):
        """
        Generate a valid secret key for Bob
        """
        raise NotImplementedError

    @abstractmethod
    def bob_public_key(self, bob_secret_key):
        """
        Generate a valid public key for Bob

        INPUT:

        - ``bob_secret_key``: Bob's secret key that will be used to generate
            the public key.

        OUTPUT:

        - A valid public key that will be sent to Alice
        """
        raise NotImplementedError

    @abstractmethod
    def alice_compute_shared_secret(self, alice_secret_key, bob_public_key):
        """
        Compute Alice's shared secret

        INPUT:

        - ``alice_secret_key``: Alice's secret key that is kept secret from all parties.

        - ``bob_public_key``: Bob's public key that has been sent to Alice.

        OUTPUT:

        - A secret key that is shared between Alice and Bob

        """
        raise NotImplementedError

    @abstractmethod
    def bob_compute_shared_secret(self, bob_secret_key, alice_public_key) -> Any:
        """
        Compute Bob's shared secret

        INPUT:

        - ``bob_secret_key``: Bob's secret key that is kept secret from all parties.

        - ``alice_public_key``: Alice's public key that has been sent to Bob.

        OUTPUT:

        - The secret key that is shared between Alice and Bob

        """
        raise NotImplementedError

    def alice_key_generate(self) -> tuple[Any, Any]:
        """
        Generate a valid (secret key, public key) pair for Alice's
        key exchange

        OUTPUT:
            A two tuple (secret_key, public_key) which is Alice's
            secret and public keys
        """
        alice_sk = self.alice_secret_key()
        alice_pk = self.alice_public_key(alice_sk)
        return (alice_sk, alice_pk)

    def bob_key_generate(self) -> tuple[Any, Any]:
        """
        Generate a valid (secret key, public key) pair for Bob's
        key exchange

        OUTPUT:
            A 2-tuple (secret_key, public_key) which is Bob's
            secret and public keys
        """
        bob_sk = self.bob_secret_key()
        bob_pk = self.bob_public_key(bob_sk)
        return (bob_sk, bob_pk)

    def do_key_exchange(self) -> tuple[Any, Any, Any, Any, Any]:
        """
        Do a full key exchange and returns all public keys, secret keys,
        and the computed shared secret between Alice and Bob. Raises
        an AssertException if the computed shared secret between Alice
        and Bob are not the same.

        OUTPUT:
            A 5-tuple (alice_secret_key, alice_public_key, bob_secret_key, bob_public_key, shared_secret).
        """
        alice_sk, alice_pk = self.alice_key_generate()
        bob_sk, bob_pk = self.bob_key_generate()
        alice_shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
        assert alice_shared_secret == self.bob_compute_shared_secret(bob_sk, alice_pk)
        return (alice_sk, alice_pk, bob_sk, bob_pk, alice_shared_secret)

    def _test_key_exchange(self, **options):
        """
        Tests the key exchange generates the same shared secrets for both parties
        """
        tester = self._tester(**options)
        for _ in range(5):
            alice_sk, alice_pk = self.alice_key_generate()
            bob_sk, bob_pk = self.bob_key_generate()
            alice_shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
            bob_shared_secret = self.bob_compute_shared_secret(bob_sk, alice_pk)
            tester.assertEqual(alice_shared_secret, bob_shared_secret)


class CommutativeKeyExchangeBase(KeyExchangeBase):
    """
    A base class for Key Exchange schemes where alice and bob perform
    the same computations for generating public/secret keys and the
    shared secret key
    """

    @abstractmethod
    def secret_key(self) -> Any:
        """
        Generate a secret key for the key exchange
        """
        raise NotImplementedError

    @abstractmethod
    def public_key(self, secret_key) -> Any:
        """
        Generate a public key for the secret key that you have chosen

        INPUT:
            - ``secret_key``: A secret key that has been chosen beforehand
        """
        raise NotImplementedError

    @abstractmethod
    def compute_shared_secret(self, secret_key, public_key) -> Any:
        """
        Generate the computed shared secret

        INPUT:
            - ``secret_key``: A secret key that has been chosen beforehand
            - ``public_key``: A public key that has been sent to this party through
                an insecure channel
        OUTPUT:
            - A shared secret key between the two parties
        """
        raise NotImplementedError

    def alice_secret_key(self) -> Any:
        """
        Generate a valid secret key for Alice
        """
        return self.secret_key()

    def alice_public_key(self, alice_secret_key) -> Any:
        """
        Generate a valid public key for Alice

        INPUT:

        - ``alice_secret_key``: Alice's secret key that will be used to generate
            the public key.

        OUTPUT:

        - A valid public key that will be sent to Bob
        """
        return self.public_key(alice_secret_key)

    def bob_secret_key(self) -> Any:
        """
        Generate a valid secret key for Bob
        """
        return self.secret_key()

    def bob_public_key(self, bob_secret_key) -> Any:
        """
        Generate a valid public key for Bob

        INPUT:

        - ``bob_secret_key``: Bob's secret key that will be used to generate
            the public key.

        OUTPUT:

        - A valid public key that will be sent to Alice
        """
        return self.public_key(bob_secret_key)

    def alice_compute_shared_secret(self, alice_sk, bob_pk) -> Any:
        """
        Compute Alice's shared secret

        INPUT:

        - ``alice_secret_key``: Alice's secret key that is kept secret from all parties.

        - ``bob_public_key``: Bob's public key that has been sent to Alice.

        OUTPUT:

        - A secret key that is shared between Alice and Bob

        """
        return self.compute_shared_secret(alice_sk, bob_pk)

    def bob_compute_shared_secret(self, bob_sk, alice_pk) -> Any:
        """
        Compute Bob's shared secret

        INPUT:

        - ``bob_secret_key``: Bob's secret key that is kept secret from all parties.

        - ``alice_public_key``: Alice's public key that has been sent to Bob.

        OUTPUT:

        - The secret key that is shared between Alice and Bob
        """
        return self.compute_shared_secret(bob_sk, alice_pk)

    def __dir__(self):
        return ['secret_key', 'public_key', 'compute_shared_secret']
