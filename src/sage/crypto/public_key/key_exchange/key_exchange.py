from abc import ABC, abstractmethod
from sage.structure.sage_object import SageObject

class KeyExchange(SageObject, ABC):

    @abstractmethod
    def alice_secret_key(self):
        raise NotImplementedError

    @abstractmethod
    def alice_public_key(self, alice_secret_key):
        raise NotImplementedError

    @abstractmethod
    def bob_secret_key(self):
        raise NotImplementedError

    @abstractmethod
    def bob_public_key(self, bob_secret_key):
        raise NotImplementedError

    @abstractmethod
    def alice_compute_shared_secret(self, alice_sk, bob_pk):
        raise NotImplementedError

    @abstractmethod
    def bob_compute_shared_secret(self, bob_sk, alice_pk):
        raise NotImplementedError

    def alice_keygen(self):
        alice_sk = self.alice_secret_key()
        alice_pk = self.alice_public_key(alice_sk)
        bob_sk = self.bob_secret_key()
        bob_pk = self.bob_public_key(bob_sk)
        shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
        return (alice_sk, alice_pk, shared_secret)

    def bob_keygen(self):
        alice_sk = self.alice_secret_key()
        alice_pk = self.alice_public_key(alice_sk)
        bob_sk = self.bob_secret_key()
        bob_pk = self.bob_public_key(bob_sk)
        shared_secret = self.alice_compute_shared_secret(bob_sk, alice_pk)
        return (bob_sk, bob_pk, shared_secret)

    def compute_all_keys(self):
        alice_sk = self.alice_secret_key()
        alice_pk = self.alice_public_key(alice_sk)
        bob_sk = self.bob_secret_key()
        bob_pk = self.bob_public_key(bob_sk)
        alice_shared_secret = self.alice_compute_shared_secret(alice_sk, bob_pk)
        bob_shared_secret = self.bob_compute_shared_secret(bob_sk, alice_pk)
        return (alice_sk, alice_pk, alice_shared_secret, bob_sk, bob_pk, bob_shared_secret)
       


