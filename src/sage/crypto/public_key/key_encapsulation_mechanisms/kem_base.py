r"""
KEM base class
"""

# Copyright (C) 2026 Vivek Kumar <viveksantoshkumardubey@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# https://www.gnu.org/licenses/

from abc import abstractmethod

from sage.misc.superseded import experimental_warning
from sage.structure.sage_object import SageObject

experimental_warning(
    41218,
    "SageMath's key encapsulation functionality is experimental and might change in the future.",
)


class KEMBase(SageObject):
    r"""
    Abstract base class for Key Encapsulation Mechanisms.

    Subclasses must implement:

    - ``keygen()`` -- return a (public key, secret key) pair
    - ``encaps(public_key)`` -- return a (ciphertext, shared secret) pair
    - ``decaps(secret_key, ciphertext)`` -- return a shared secret
    """

    @abstractmethod
    def keygen(self):
        """
        Generate a public and secret key pair.

        OUTPUT: tuple ``(public_key, secret_key)``
        """
        raise NotImplementedError("Subclasses must implement keygen()")

    @abstractmethod
    def encaps(self, public_key):
        """
        Encapsulate a shared secret using a public key.

        INPUT:

        - ``public_key`` -- the recipient's public key

        OUTPUT: tuple ``(ciphertext, shared_secret)``
        """
        raise NotImplementedError("Subclasses must implement encaps()")

    @abstractmethod
    def decaps(self, secret_key, ciphertext):
        """
        Decapsulate a ciphertext to recover the shared secret.

        INPUT:

        - ``secret_key`` -- the recipient's secret key
        - ``ciphertext`` -- the ciphertext from ``encaps()``

        OUTPUT: shared secret
        """
        raise NotImplementedError("Subclasses must implement decaps()")

    def _test_kem_correctness(self, **options):
        """
        Check that encaps and decaps produce the same shared secret.
        """
        tester = self._tester(**options)
        pk, sk = self.keygen()
        ct, ss1 = self.encaps(pk)
        ss2 = self.decaps(sk, ct)
        tester.assertEqual(ss1, ss2, "Shared secrets do not match")
