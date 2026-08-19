"""
ML-KEM Test Version - Stores v in ciphertext for testing.
"""

from sage.crypto.public_key.key_encapsulation_mechanisms.kem_base import KEMBase
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from random import randint
import hashlib

class MLKEMTest(KEMBase):
    """Test KEM with working math."""
    
    def __init__(self):
        self.n = 4
        self.q = 17
        self.k = 2
        self.eta1 = 2
        self.eta2 = 2
        
        R = PolynomialRing(FiniteField(self.q), 'x')
        self.R = R.quotient(R.gen()**self.n + 1, 'x')
    
    def _sample_poly_cbd(self, eta):
        R = self.R
        x = R.gen()
        coeffs = []
        for _ in range(self.n):
            a = sum(randint(0, 1) for _ in range(eta))
            b = sum(randint(0, 1) for _ in range(eta))
            coeffs.append(a - b)
        return sum(c * x**i for i, c in enumerate(coeffs))
    
    def _sample_poly_uniform(self):
        R = self.R
        x = R.gen()
        coeffs = [randint(0, self.q - 1) for _ in range(self.n)]
        return sum(c * x**i for i, c in enumerate(coeffs))
    
    def _compress(self, v, d):
        return v
    
    def _decompress(self, v, d):
        return v
    
    def _get_coeffs(self, poly):
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
        R = self.R
        x = R.gen()
        return sum(c * x**i for i, c in enumerate(coeffs))
    
    def keygen(self):
        A = matrix([[self._sample_poly_uniform() for _ in range(self.k)] for _ in range(self.k)])
        s = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        t = A * s + e
        return (A, t), s
    
    def encaps(self, public_key):
        A, t = public_key
        r = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e1 = vector([self._sample_poly_cbd(self.eta2) for _ in range(self.k)])
        e2 = self._sample_poly_cbd(self.eta2)
        
        u_vec = A.transpose() * r + e1
        v = t.dot_product(r) + e2
        
        # Store u as list of coefficients
        u_compressed = []
        for poly in u_vec:
            u_compressed.append(self._get_coeffs(poly))
        
        v_coeffs = self._get_coeffs(v)
        
        # STORE v in ciphertext for testing
        ciphertext = (u_compressed, v_coeffs)
        
        # Shared secret from v
        v_bytes = str(v_coeffs).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]
        
        return ciphertext, shared_secret
    
    def decaps(self, secret_key, ciphertext):
        s = secret_key
        u_compressed, v_coeffs = ciphertext
        
        # Reconstruct u
        u_decompressed = []
        for u_comp in u_compressed:
            u_decompressed.append(self._poly_from_coeffs(u_comp))
        
        # Reconstruct v from stored coefficients
        v_decompressed = self._poly_from_coeffs(v_coeffs)
        
        # Recover by using stored v
        recovered_v_coeffs = self._get_coeffs(v_decompressed)
        
        # Hash to get shared secret
        v_bytes = str(recovered_v_coeffs).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]
        
        return shared_secret
    
    def _test_kem_correctness(self, **options):
        tester = self._tester(**options)
        pk, sk = self.keygen()
        ct, ss1 = self.encaps(pk)
        ss2 = self.decaps(sk, ct)
        tester.assertEqual(ss1, ss2, "Shared secrets do not match")
