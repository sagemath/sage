"""
KEM base class
"""

from sage.structure.sage_object import SageObject

class KEMBase(SageObject):
    """
    Abstract base class for Key Encapsulation Mechanisms.

    Subclasses must implement:
    - keygen() -> (public_key, secret_key)
    - encaps(public_key) -> (ciphertext, shared_secret)
    - decaps(secret_key, ciphertext) -> shared_secret
    """

    def keygen(self):
        raise NotImplementedError("Subclasses must implement keygen()")

    def encaps(self, public_key):
        raise NotImplementedError("Subclasses must implement encaps()")

    def decaps(self, secret_key, ciphertext):
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
