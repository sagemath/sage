
from typing import TYPE_CHECKING, Any

from sage.arith.misc import is_prime
from sage.misc.prandom import randint
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer

from .digital_signature_base import DigitalSignatureBase

if TYPE_CHECKING:
    from sage.rings.finite_rings.integer_mod import IntegerMod_abstract

class DigitalSignatureAlgorithm(DigitalSignatureBase):

    def __init__(self, p: Integer | int, q: Integer | int, generator: Integer | int | IntegerMod_abstract) -> None:
        r"""
        Create an instance of the Digital Signature Algorithm using
        primes ``p`` and ``q`` and generator ``g``

        INPUT:

        - ``q`` -- prime integer

        - ``p`` -- prime integer for which (p-1) is a multiple of ``q``

        - ``generator`` -- base non-one element of the digital signature `\frac{\mathbb{Z}}{q\mathbb{Z}}`.


        EXAMPLES::

            sage: DSA = digital_signature.DigitalSignatureAlgorithm(4007, 2003, 867)
            doctest:...: FutureWarning: SageMath's digital signature functionality is experimental and might change in the future.
                         See https://github.com/sagemath/sage/issues/41218 for details.
            sage: message = 300
            sage: public, secret, sig, message, result = DSA.do_signature(message)
            sage: result
            True
        """
        if (p-1) % q != 0:
            raise ValueError('p-1 must be a multiple of q')

        if not is_prime(q) or not is_prime(p):
            raise ValueError('p and q must be prime')
        if generator == 1:
            raise ValueError('generator must not be equal to one')
        self.q = q
        self.p = p
        self.pZmod = IntegerModRing(p)
        self.qZmod = IntegerModRing(q)
        self.generator = self.pZmod(generator)


    def generate_keys(self) -> tuple[IntegerMod_abstract, IntegerMod_abstract]:
        """
        Generates a keypair to be used for signatures

        OUTPUT:

        Returns a 2-tuple of (public_key, secret_key)
        """
        secret_key = randint(1, self.q-1)
        public_key = self.pZmod(self.generator)**secret_key
        return (public_key, secret_key)

    def sign(self, message, secret_key: IntegerMod_abstract) -> tuple[tuple[IntegerMod_abstract, IntegerMod_abstract], Any]:
        """
        Signs a message from the secret_key

        INPUT:

        - ``message`` -- Integer message value to be signed
        - ``secret_key`` -- The secret_key that only the signer has.

        OUTPUT:

        Returns a pair of (signature, message)
        """
        k = Integer(randint(1, self.q-1))
        r = self.qZmod(self.pZmod(self.generator**k))
        s = self.qZmod(k)**(-1) * (self.qZmod(message) + self.qZmod(secret_key) * r)
        return ((r, s), message)

    def verify(self, public_key: IntegerMod_abstract, signature: tuple[IntegerMod_abstract, IntegerMod_abstract], message) -> bool:
        """
        Verifies that the signature is valid

        INPUT:

        - ``public_key`` -- A public key that is used to verify that the signature is valid
        - ``signature`` -- The signature that the verifier is checking is valid
        - ``message`` -- The message being sent over, assumed to be an integer for the test suite

        OUTPUT:

        Returns True if it is a valid signature for the public_key and message and returns
        False if not.
        """
        r, s = signature
        s = self.qZmod(s)
        w = s**(-1)
        u1 = self.qZmod(message) * w
        u2 = r * w
        v = self.qZmod(
            self.generator**u1 * public_key**u2
        )
        return v == r

    def parameters(self) -> tuple[Any, Any, IntegerMod_abstract]:
        """
        Returns the public parameter set, which is of the form ``(p, q, generator)``

        OUTPUT:

        The public parameter set of the signature scheme ``(p, q, g)``.
        """
        return (self.p, self.q, self.generator)
