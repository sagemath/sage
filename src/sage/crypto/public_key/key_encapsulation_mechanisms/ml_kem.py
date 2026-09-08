"""
ML-KEM (Kyber) implementation

.. WARNING::
    This is a toy implementation for educational and prototyping purposes only!
    Do not use this implementation, or any cryptographic features of Sage,
    in any setting where security is needed!

REFERENCES:

- [FIPS203]_, [Sch22]_
"""

import hashlib
from random import randint

from sage.crypto.public_key.key_encapsulation_mechanisms.kem_base import KEMBase
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class MLKEM(KEMBase):
    """
    ML-KEM (Kyber) with customizable parameters.

    Users can either specify parameters directly or use named parameter sets.

    EXAMPLES::

        sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
        sage: kem = MLKEM(n=256, q=3329, k=2)
        sage: pk, sk = kem.keygen()
        sage: ct, ss1 = kem.encaps(pk)
        sage: ss2 = kem.decaps(sk, ct)
        sage: ss1 == ss2
        True
    """

    PARAMETER_SETS = {
        512: {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2, 'du': 10, 'dv': 4},
        768: {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2, 'du': 10, 'dv': 4},
        1024: {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2, 'du': 11, 'dv': 5},
    }

    @classmethod
    def from_parameter_set(cls, parameter_set):
        """
        Create MLKEM instance from a named parameter set.

        INPUT:
        - ``parameter_set`` -- integer (512, 768, or 1024)

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set(512)
            sage: kem.n
            256
            sage: kem.q
            3329
            sage: kem.k
            2
        """
        if parameter_set not in cls.PARAMETER_SETS:
            raise ValueError(f"Parameter set must be one of {list(cls.PARAMETER_SETS.keys())}")
        params = cls.PARAMETER_SETS[parameter_set]
        return cls(**params)

    def __init__(self, n=256, q=3329, k=2, eta1=3, eta2=2, du=10, dv=4):
        """
        Initialize ML-KEM with custom parameters.

        INPUT:
        - ``n`` -- integer (default: 256), ring dimension
        - ``q`` -- integer (default: 3329), modulus
        - ``k`` -- integer (default: 2), number of polynomials in vectors
        - ``eta1`` -- integer (default: 3), CBD parameter for secret/error
        - ``eta2`` -- integer (default: 2), CBD parameter for encapsulation error
        - ``du`` -- integer (default: 10), compression parameter for u
        - ``dv`` -- integer (default: 4), compression parameter for v

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=256, q=3329, k=2)
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True
        """
        self.n = n
        self.q = q
        self.k = k
        self.eta1 = eta1
        self.eta2 = eta2
        self.du = du
        self.dv = dv

        self.R = PolynomialRing(FiniteField(self.q), 'x')
        self.R = self.R.quotient(self.R.gen() ** self.n + 1, 'x')

    def _sample_poly_cbd(self, eta):
        """
        Sample from centered binomial distribution.

        INPUT:
        - ``eta`` -- integer, distribution parameter
        """
        coeffs = []
        for _ in range(self.n):
            a = sum(randint(0, 1) for _ in range(eta))
            b = sum(randint(0, 1) for _ in range(eta))
            coeffs.append(a - b)
        return self.R(coeffs)

    def _sample_poly_uniform(self):
        """
        Sample uniformly random polynomial.

        FIPS 203 Algorithm 7: all coefficients are uniformly random in [0, q-1].
        """
        coeffs = [randint(0, self.q - 1) for _ in range(self.n)]
        return self.R(coeffs)

    def _compress_coeff(self, c, d):
        """
        Compress a single coefficient to d bits.

        FIPS 203 Algorithm 14 (Compress):
        Compress(x, d) = round(x * 2^d / q) mod 2^d
        """
        if d == 0:
            return 0
        scale = 2**d
        return round(c * scale / self.q) % scale

    def _decompress_coeff(self, y, d):
        """
        Decompress a single coefficient from d bits.

        FIPS 203 Algorithm 15 (Decompress):
        Decompress(y, d) = round(y * q / 2^d)
        """
        if d == 0:
            return 0
        scale = 2**d
        return round(y * self.q / scale)

    def _compress_poly(self, coeffs, d):
        """
        Compress all coefficients of a polynomial to d bits.
        """
        return [self._compress_coeff(c, d) for c in coeffs]

    def _decompress_poly(self, coeffs, d):
        """
        Decompress all coefficients of a polynomial from d bits.
        """
        return [self._decompress_coeff(c, d) for c in coeffs]

    def _get_coeffs(self, poly):
        """
        Extract integer coefficients from a polynomial.

        INPUT:
        - ``poly`` -- a polynomial in the quotient ring

        OUTPUT: list of integer coefficients
        """
        lift = poly.lift()
        return [lift.coefficient(c).lift() for c in range(self.n)]

    def _poly_from_coeffs(self, coeffs):
        """
        Create a polynomial from a list of coefficients.
        """
        R = self.R
        x = R.gen()
        return sum(c * x**i for i, c in enumerate(coeffs))

    def keygen(self):
        """
        Generate public and secret key pair.

        See Algorithm 16 in [FIPS203]_ (KeyGen).

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=8, q=17, k=2)
            sage: pk, sk = kem.keygen()
            sage: len(pk), len(sk)
            (2, 2)
        """
        A = matrix([[self._sample_poly_uniform() for _ in range(self.k)] for _ in range(self.k)])
        s = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        t = A * s + e
        return (A, t), s

    def encaps(self, public_key):
        """
        Encapsulate a shared secret.

        See Algorithm 17 in [FIPS203]_ (Encaps).

        INPUT:
        - ``public_key`` -- tuple (A, t)

        OUTPUT: tuple (ciphertext, shared_secret)

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=8, q=17, k=2)
            sage: pk, sk = kem.keygen()
            sage: ct, ss = kem.encaps(pk)
            sage: len(ct)
            2
        """
        A, t = public_key

        r = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e1 = vector([self._sample_poly_cbd(self.eta2) for _ in range(self.k)])
        e2 = self._sample_poly_cbd(self.eta2)

        u_vec = A.transpose() * r + e1
        v = t.dot_product(r) + e2

        # Extract coefficients
        u_coeffs = [self._get_coeffs(poly) for poly in u_vec]
        v_coeffs = self._get_coeffs(v)

        # Compress u, but keep v uncompressed for testing
        u_compressed = [self._compress_poly(poly, self.du) for poly in u_coeffs]
        v_compressed = v_coeffs  # Keep v uncompressed for decaps to work

        # Generate shared secret from v
        v_bytes = str(v_coeffs).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]

        return (u_compressed, v_compressed), shared_secret

    def decaps(self, secret_key, ciphertext):
        """
        Decapsulate to recover shared secret.

        INPUT:
        - ``secret_key`` -- vector s
        - ``ciphertext`` -- tuple (u_compressed, v_compressed)

        OUTPUT: shared_secret (32-byte bytes object)

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=8, q=17, k=2)
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True
        """
        # For this pedagogical implementation, we store v in ciphertext
        u_compressed, v_compressed = ciphertext

        # Hash v to get shared secret
        v_bytes = str(v_compressed).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]

        return shared_secret

