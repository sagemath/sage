from sage.crypto.public_key.key_exchange.key_exchange import KeyExchange
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.misc.prandom import randint

class FiniteFieldDH(KeyExchange):
    def __init__(self, p, generator):
        self._field = GF(p, impl='modn')
        self.p = p
        self.generator = self._field(generator)

    def alice_secret_key(self):
        return randint(2, self.p - 2)

    def alice_public_key(self, sk):
        return self.generator ** sk

    def bob_secret_key(self):
        return randint(2, self.p - 2)

    def bob_public_key(self, sk):
        return self.generator ** sk

    def alice_compute_shared_secret(self, alice_sk, bob_pk):
        return bob_pk ** alice_sk

    def bob_compute_shared_secret(self, bob_sk, alice_pk):
        return alice_pk ** bob_sk

