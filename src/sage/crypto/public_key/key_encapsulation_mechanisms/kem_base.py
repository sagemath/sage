"""
KEM base class

.. WARNING::
    This module is experimental. The API may change in future versions.
"""

from abc import ABC, abstractmethod

from sage.structure.sage_object import SageObject


class KEMBase(SageObject, ABC):
    """
    Abstract base class for Key Encapsulation Mechanisms.

    Subclasses must implement:
    - keygen() -> (public_key, secret_key)
    - encaps(public_key) -> (ciphertext, shared_secret)
    - decaps(secret_key, ciphertext) -> shared_secret
    """

    @abstractmethod
    def keygen(self):
        """
        Generate a public and secret key pair.

        OUTPUT: tuple (public_key, secret_key)
        """
        raise NotImplementedError("Subclasses must implement keygen()")

    @abstractmethod
    def encaps(self, public_key):
        """
        Encapsulate a shared secret using a public key.

        INPUT:
        - ``public_key`` -- the recipient's public key

        OUTPUT: tuple (ciphertext, shared_secret)
        """
        raise NotImplementedError("Subclasses must implement encaps()")

    @abstractmethod
    def decaps(self, secret_key, ciphertext):
        """
        Decapsulate a ciphertext to recover the shared secret.

        INPUT:
        - ``secret_key`` -- the recipient's secret key
        - ``ciphertext`` -- the ciphertext from encaps()

        OUTPUT: shared_secret
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
