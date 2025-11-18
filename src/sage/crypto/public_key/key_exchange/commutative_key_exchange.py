from sage.crypto.public_key.key_exchange import KeyExchange
from abc import ABC, abstractmethod

class CommutativeKeyExchange(KeyExchange, ABC):

    @abstractmethod
    def secret_key(self):
        raise NotImplementedError

    @abstractmethod
    def public_key(self, secret_key):
        raise NotImplementedError

    @abstractmethod
    def compute_shared_secret(self, sk, pk):
        raise NotImplementedError
    
    def alice_secret_key(self):
        return self.secret_key()

    def alice_public_key(self, alice_secret_key):
        return self.public_key(alice_secret_key)

    def bob_secret_key(self):
        return secret_key()

    def bob_public_key(self, bob_secret_key):
        return self.public_key(bob_secret_key)

    def alice_compute_shared_secret(self, alice_sk, bob_pk):
        return self.compute_shared_secret(alice_sk, bob_pk)

    def bob_compute_shared_secret(self, bob_sk, alice_pk):
        return self.compute_shared_secret(bob_sk, alice_pk)

