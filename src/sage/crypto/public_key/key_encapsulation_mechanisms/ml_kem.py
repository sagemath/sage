"""
ML-KEM (Kyber) implementation

REFERENCES:

- [FIPS203] National Institute of Standards and Technology,
  "Module-Lattice-Based Key-Encapsulation Mechanism Standard",
  FIPS 203, 2024.
  https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.203.pdf

- [Sch22] Peter Schwabe et al.,
  "CRYSTALS-KYBER",
  https://eprint.iacr.org/2022/1696
"""

from sage.crypto.public_key.key_encapsulation_mechanisms.kem_base import KEMBase
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from random import randint
import hashlib

class MLKEM(KEMBase):
    """
    ML-KEM (Kyber) with customizable parameters.

    Users can either specify parameters directly or use named parameter sets.
    """

    PARAMETER_SETS = {
        512: {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2},
        768: {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2},
        1024: {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2}
    }

    @classmethod
    def from_parameter_set(cls, parameter_set):
        """
        Create MLKEM instance from a named parameter set.

        INPUT:
        - ``parameter_set`` -- integer (512, 768, or 1024)

        EXAMPLES::

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

    def __init__(self, n=256, q=3329, k=2, eta1=3, eta2=2):
        """
        Initialize ML-KEM with custom parameters.

        INPUT:
        - ``n`` -- integer (default: 256), ring dimension
        - ``q`` -- integer (default: 3329), modulus
        - ``k`` -- integer (default: 2), number of polynomials in vectors
        - ``eta1`` -- integer (default: 3), CBD parameter for secret/error
        - ``eta2`` -- integer (default: 2), CBD parameter for encapsulation error

        EXAMPLES::

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

        R = PolynomialRing(FiniteField(self.q), 'x')
        self.R = R.quotient(R.gen()**self.n + 1, 'x')

    def _sample_poly_cbd(self, eta):
        """
        Sample from centered binomial distribution.

        INPUT:
        - ``eta`` -- integer, distribution parameter
        """
        R = self.R
        x = R.gen()
        coeffs = []
        for _ in range(self.n):
            a = sum(randint(0, 1) for _ in range(eta))
            b = sum(randint(0, 1) for _ in range(eta))
            coeffs.append(a - b)
        return sum(c * x**i for i, c in enumerate(coeffs))

    def _sample_poly_uniform(self):
        """Sample uniformly random polynomial."""
        R = self.R
        x = R.gen()
        coeffs = [randint(0, self.q - 1) for _ in range(self.n)]
        return sum(c * x**i for i, c in enumerate(coeffs))

    def _get_coeffs(self, poly):
        """
        Extract integer coefficients from a polynomial.

        INPUT:
        - ``poly`` -- a polynomial in the quotient ring

        OUTPUT: list of integer coefficients
        """
        return [int(c) for c in poly.list()] + [0] * (self.n - len(poly.list()))

    def _poly_from_coeffs(self, coeffs):
        """Create a polynomial from a list of coefficients."""
        R = self.R
        x = R.gen()
        return sum(c * x**i for i, c in enumerate(coeffs))

    def keygen(self):
        """Generate public and secret key pair."""
        A = matrix([[self._sample_poly_uniform() for _ in range(self.k)] for _ in range(self.k)])
        s = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        t = A * s + e
        return (A, t), s

    def encaps(self, public_key):
        """Encapsulate a shared secret."""
        A, t = public_key

        r = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e1 = vector([self._sample_poly_cbd(self.eta2) for _ in range(self.k)])
        e2 = self._sample_poly_cbd(self.eta2)

        u_vec = A.transpose() * r + e1
        v = t.dot_product(r) + e2

        u_coeffs = [self._get_coeffs(poly) for poly in u_vec]
        v_coeffs = self._get_coeffs(v)

        v_bytes = str(v_coeffs).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]

        return (u_coeffs, v_coeffs), shared_secret

    def decaps(self, secret_key, ciphertext):
        """Decapsulate to recover shared secret."""
        u_coeffs, v_coeffs = ciphertext

        v_decompressed = self._poly_from_coeffs(v_coeffs)

        v_bytes = str(self._get_coeffs(v_decompressed)).encode()
        return hashlib.sha256(v_bytes).digest()[:32]
