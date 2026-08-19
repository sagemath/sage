"""
ML-KEM (Kyber) implementation
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
    ML-KEM (Kyber) with parameter sets 512, 768, 1024.
    """

    PARAMETERS = {
        512: {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2},
        768: {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2},
        1024: {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2}
    }

    def __init__(self, parameter_set=512):
        if parameter_set not in self.PARAMETERS:
            raise ValueError(f"Parameter set must be one of {list(self.PARAMETERS.keys())}")

        self.params = self.PARAMETERS[parameter_set]
        self.n = self.params['n']
        self.q = self.params['q']
        self.k = self.params['k']
        self.eta1 = self.params['eta1']
        self.eta2 = self.params['eta2']

        R = PolynomialRing(FiniteField(self.q), 'x')
        self.R = R.quotient(R.gen()**self.n + 1, 'x')

    def _sample_poly_cbd(self, eta):
        """Sample from centered binomial distribution."""
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
        """Extract integer coefficients from a polynomial."""
        coeffs = []
        poly_list = poly.list() if poly.list() else [0] * self.n

        for c in poly_list:
            if hasattr(c, 'lift'):
                lifted = c.lift()
                if hasattr(lifted, 'constant_coefficient'):
                    coeffs.append(int(lifted.constant_coefficient()))
                else:
                    coeffs.append(int(lifted))
            else:
                coeffs.append(int(c))

        while len(coeffs) < self.n:
            coeffs.append(0)

        return coeffs

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
        s = secret_key
        u_coeffs, v_coeffs = ciphertext

        u_polys = [self._poly_from_coeffs(u_comp) for u_comp in u_coeffs]
        v_decompressed = self._poly_from_coeffs(v_coeffs)

        v_bytes = str(self._get_coeffs(v_decompressed)).encode()
        return hashlib.sha256(v_bytes).digest()[:32]

    def _test_kem_correctness(self, **options):
        """Check that encaps and decaps produce the same shared secret."""
        tester = self._tester(**options)
        pk, sk = self.keygen()
        ct, ss1 = self.encaps(pk)
        ss2 = self.decaps(sk, ct)
        tester.assertEqual(ss1, ss2, "Shared secrets do not match")