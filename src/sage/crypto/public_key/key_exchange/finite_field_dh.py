r"""
Toy implementation of Diffie-Hellman key exchange over finite fields `\Zmod{p}`.

AUTHORS:

- Brian Heckel (2025-11-26): initial version
"""

# ****************************************************************************
#       Copyright (C) 2025 Brian Heckel <heckelbri@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.crypto.public_key.key_exchange.key_exchange_base import CommutativeKeyExchangeBase
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
