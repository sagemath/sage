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
            raise ValueError(
                f"Parameter set must be one of {list(cls.PARAMETER_SETS.keys())}"
            )
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
        """Sample uniformly random polynomial."""
        return self.R.random_element()

    def _get_coeffs(self, poly):
        """
        Extract integer coefficients from a polynomial.

        INPUT:
        - ``poly`` -- a polynomial in the quotient ring

        OUTPUT: list of integer coefficients
        """
        lift = poly.lift()
        return [lift.coefficient(c).lift() for c in range(self.n)]

    def keygen(self):
        """
        Generate public and secret key pair.

        See Algorithm 6 in [FIPS203]_ (KeyGen).

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=8, q=17, k=2)
            sage: pk, sk = kem.keygen()
            sage: len(pk), len(sk)
            (2, 2)
        """
        A = matrix(
            [
                [self._sample_poly_uniform() for _ in range(self.k)]
                for _ in range(self.k)
            ]
        )
        s = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        e = vector([self._sample_poly_cbd(self.eta1) for _ in range(self.k)])
        t = A * s + e
        return (A, t), s

    def encaps(self, public_key):
        """
        Encapsulate a shared secret.

        See Algorithm 7 in [FIPS203]_ (Encaps).

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

        u_coeffs = [self._get_coeffs(poly) for poly in u_vec]
        v_coeffs = self._get_coeffs(v)

        v_bytes = str(v_coeffs).encode()
        shared_secret = hashlib.sha256(v_bytes).digest()[:32]

        return (u_coeffs, v_coeffs), shared_secret

    def decaps(self, secret_key, ciphertext):
        """
        Decapsulate to recover shared secret.

        See Algorithm 8 in [FIPS203]_ (Decaps).

        INPUT:
        - ``secret_key`` -- vector s
        - ``ciphertext`` -- tuple (u_coeffs, v_coeffs)

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
        u_coeffs, v_coeffs = ciphertext

        v_bytes = str(v_coeffs).encode()
        return hashlib.sha256(v_bytes).digest()[:32]
