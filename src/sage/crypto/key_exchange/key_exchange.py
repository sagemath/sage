r"""
Public key exchange schemes.

This module contains base classes for public key exchange schemes. The classes
defined in this module should not be called directly. It is the responsibility
of child classes to implement specific cryptosystems.


AUTHORS:

- Vincent Macri (2024-07-15): initial version
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

    @experimental(37305)
    def __init__(self):
        pass

    def generate_secret_key(self):
        raise NotImplementedError

    def generate_public_key(self, secret_key):
        raise NotImplementedError

    def compute_shared_secret(self, alice_pk, bob_sk):
        raise NotImplementedError

    def parameters(self):
        raise NotImplementedError
