import random
from sage.crypto.public_key.key_exchange import KeyExchange

class AnshelAnshelGoldfeld(KeyExchange):
    def __init__(self, G, A, B):
        """
        - A is a list of elements in the group G
        - B is a list of elements in the group G
        """
        self._group = G
        # TODO: check that the group is
        # non-abelian and the subgroups
        # are subgroups of 
        self._alice_list = A
        self._bob_list = B
        self._alice_subgroup = G.subgroup(A)
        self._bob_subgroup = G.subgroup(B)
        self._alice_size = alice_size
        self._bob_size = bob_size

    def alice_secret_key(self):
        return random.choice(self._alice_subgroup)

    def alice_public_key(self, alice_secret_key):
        pub_keys = []
        sk_inverse = alice_secret_key ** (-1)
        sk = alice_secret_key
        for b in self._bob_list:
            conj_by_A = sk_inverse * b * sk
            pub_keys.append(conj_by_A)
        return pub_keys

    def bob_secret_key(self):
        return random.choice(self._bob_subgroup)

    def bob_public_key(self, bob_secret_key):
        pub_keys = []
        sk_inverse = bob_secret_key ** (-1)
        sk = bob_secret_key
        for b in self._alice_list:
            conj_by_A = sk_inverse * b * sk
            pub_keys.append(conj_by_A)
        return pub_keys

    def alice_compute_shared_secret(self, alice_sk, bob_pk):
        beta = 
        raise NotImplementedError

    def bob_compute_shared_secret(self, bob_sk, alice_pk):
        raise NotImplementedError
