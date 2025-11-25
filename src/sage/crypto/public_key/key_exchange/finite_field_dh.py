from sage.crypto.public_key.key_exchange.key_exchange import CommutativeKeyExchangeBase
from sage.misc.prandom import randint
from sage.rings.finite_rings.finite_field_constructor import GF


class FiniteFieldDH(CommutativeKeyExchangeBase):
    def __init__(self, p, generator):
        self._field = GF(p, impl='modn')
        self.p = p
        self.generator = self._field(generator)

    def secret_key(self):
        return self._field(randint(2, self.p - 2))

    def public_key(self, sk):
        return self.generator ** sk

    def compute_shared_secret(self, sk, pk):
        return pk ** sk


