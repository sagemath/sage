r"""
ML-KEM (Kyber) implementation.

This module implements the ML-KEM key encapsulation mechanism as
standardized in [FIPS203]_.  The implementation follows the algorithms
and notation of that standard and is byte-compatible with the NIST
known-answer tests.

.. WARNING::

    This is a pedagogical implementation and has not been audited.
    It should **not** be used in any setting where security is needed.

REFERENCES:

- [FIPS203]_, [Sch22]_
"""

import hashlib
import secrets

from sage.crypto.public_key.key_encapsulation_mechanisms.kem_base import KEMBase
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class MLKEM(KEMBase):
    r"""
    The ML-KEM (Kyber) key encapsulation mechanism.

    ML-KEM is defined in FIPS 203 [FIPS203]_.  It has three standardized
    parameter sets (ML-KEM-512, ML-KEM-768, and ML-KEM-1024); use
    :meth:`from_parameter_set` to construct one of those, or pass the
    individual parameters to the constructor.

    EXAMPLES::

        sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
        sage: kem = MLKEM.from_parameter_set(512)
        sage: pk, sk = kem.keygen()
        sage: ct, ss1 = kem.encaps(pk)
        sage: ss2 = kem.decaps(sk, ct)
        sage: ss1 == ss2
        True
    """

    _N = 256
    _Q = 3329

    # Twiddle factors for the negacyclic NTT in Z_q[x] / (x^256 + 1):
    # zeta^{BitRev7(i)} mod q, where zeta = 17 is a primitive 256th root
    # of unity modulo q.
    _ZETAS = [
        1, 1729, 2580, 3289, 2642, 630, 1897, 848,
        1062, 1919, 193, 797, 2786, 3260, 569, 1746,
        296, 2447, 1339, 1476, 3046, 56, 2240, 1333,
        1426, 2094, 535, 2882, 2393, 2879, 1974, 821,
        289, 331, 3253, 1756, 1197, 2304, 2277, 2055,
        650, 1977, 2513, 632, 2865, 33, 1320, 1915,
        2319, 1435, 807, 452, 1438, 2868, 1534, 2402,
        2647, 2617, 1481, 648, 2474, 3110, 1227, 910,
        17, 2761, 583, 2649, 1637, 723, 2288, 1100,
        1409, 2662, 3281, 233, 756, 2156, 3015, 3050,
        1703, 1651, 2789, 1789, 1847, 952, 1461, 2687,
        939, 2308, 2437, 2388, 733, 2337, 268, 641,
        1584, 2298, 2037, 3220, 375, 2549, 2090, 1645,
        1063, 319, 2773, 757, 2099, 561, 2466, 2594,
        2804, 1092, 403, 1026, 1143, 2150, 2775, 886,
        1722, 1212, 1874, 1029, 2110, 2935, 885, 2154,
    ]

    _FF = FiniteField(_Q)
    _POLY_RING = PolynomialRing(_FF, 'x')
    _R = _POLY_RING.quotient(_POLY_RING.gen() ** _N + 1, 'x')

    PARAMETER_SETS = {
        512:  {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2, 'du': 10, 'dv': 4},
        768:  {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2, 'du': 10, 'dv': 4},
        1024: {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2, 'du': 11, 'dv': 5},
    }

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_parameter_set(cls, parameter_set):
        r"""
        Construct an :class:`MLKEM` object from a named parameter set.

        INPUT:

        - ``parameter_set`` -- integer; one of ``512``, ``768``, or ``1024``

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set(768)
            sage: kem.k
            3
        """
        if parameter_set not in cls.PARAMETER_SETS:
            raise ValueError(
                f"parameter_set must be one of {sorted(cls.PARAMETER_SETS)}, "
                f"got {parameter_set!r}"
            )
        return cls(**cls.PARAMETER_SETS[parameter_set])

    def __init__(self, n=256, q=3329, k=2, eta1=3, eta2=2, du=10, dv=4):
        r"""
        Initialize an ML-KEM instance with explicit parameters.

        INPUT:

        - ``n`` -- integer (default: 256), ring degree
        - ``q`` -- integer (default: 3329), coefficient modulus
        - ``k`` -- integer (default: 2), module rank
        - ``eta1`` -- integer (default: 3), CBD parameter for the secret
          and error vectors
        - ``eta2`` -- integer (default: 2), CBD parameter for the
          encapsulation noise
        - ``du`` -- integer (default: 10), compression parameter for ``u``
        - ``dv`` -- integer (default: 4), compression parameter for ``v``

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=256, q=3329, k=2, eta1=3, eta2=2, du=10, dv=4)
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True
        """
        if n != 256:
            raise NotImplementedError("Only n = 256 is supported (FIPS 203).")
        if q != 3329:
            raise NotImplementedError("Only q = 3329 is supported (FIPS 203).")
        if eta1 not in (2, 3) or eta2 not in (2, 3):
            raise ValueError("eta1 and eta2 must be 2 or 3.")
        if not (1 <= k <= 4):
            raise ValueError("k must be in {1, 2, 3, 4}.")

        self.n = n
        self.q = q
        self.k = k
        self.eta1 = eta1
        self.eta2 = eta2
        self.du = du
        self.dv = dv

    # ------------------------------------------------------------------
    # Polynomial helpers
    # ------------------------------------------------------------------

    @classmethod
    def _from_list(cls, L):
        r"""
        Return the element of ``R`` with the given coefficient list.
        """
        return cls._R(cls._POLY_RING(L))

    @classmethod
    def _coeffs(cls, f):
        r"""
        Return the length-``n`` coefficient list of an element of ``R``.
        """
        L = [int(c) for c in f.lift().list()]
        if len(L) < cls._N:
            L += [0] * (cls._N - len(L))
        return L[: cls._N]

    # ------------------------------------------------------------------
    # Negacyclic NTT (FIPS 203, Section 4.3)
    # ------------------------------------------------------------------

    @classmethod
    def _ntt(cls, f):
        r"""
        Return the negacyclic NTT of ``f`` in ``Z_q[x] / (x^256 + 1)``.

        INPUT:

        - ``f`` -- element of ``R``

        OUTPUT: element of ``R`` representing the NTT of ``f``
        """
        F = cls._coeffs(f)
        i = 1
        length = 128
        while length >= 2:
            for start in range(0, cls._N, 2 * length):
                zeta = cls._ZETAS[i]
                i += 1
                for j in range(start, start + length):
                    t = (zeta * F[j + length]) % cls._Q
                    F[j + length] = (F[j] - t) % cls._Q
                    F[j] = (F[j] + t) % cls._Q
            length //= 2
        return cls._from_list(F)

    @classmethod
    def _ntt_inv(cls, f):
        r"""
        Return the inverse of :meth:`_ntt`.

        INPUT:

        - ``f`` -- element of ``R`` in NTT domain

        OUTPUT: element of ``R``
        """
        F = cls._coeffs(f)
        i = 127
        length = 2
        while length <= 128:
            for start in range(0, cls._N, 2 * length):
                zeta = cls._ZETAS[i]
                i -= 1
                for j in range(start, start + length):
                    t = F[j]
                    F[j] = (t + F[j + length]) % cls._Q
                    F[j + length] = (zeta * (F[j + length] - t)) % cls._Q
            length *= 2
        F = [(3303 * x) % cls._Q for x in F]
        return cls._from_list(F)

    @classmethod
    def _multiply_ntts(cls, f, g):
        r"""
        Return the pointwise product of ``f`` and ``g`` in NTT domain.

        INPUT:

        - ``f``, ``g`` -- elements of ``R`` in NTT domain

        OUTPUT: element of ``R``
        """
        F = cls._coeffs(f)
        G = cls._coeffs(g)
        H = [0] * cls._N
        for i in range(64):
            a0, a1 = F[4 * i], F[4 * i + 1]
            b0, b1 = F[4 * i + 2], F[4 * i + 3]
            c0, c1 = G[4 * i], G[4 * i + 1]
            d0, d1 = G[4 * i + 2], G[4 * i + 3]
            gamma = cls._ZETAS[64 + i]
            gamma_neg = -gamma % cls._Q
            H[4 * i] = (a0 * c0 + a1 * c1 * gamma) % cls._Q
            H[4 * i + 1] = (a0 * c1 + a1 * c0) % cls._Q
            H[4 * i + 2] = (b0 * d0 + b1 * d1 * gamma_neg) % cls._Q
            H[4 * i + 3] = (b0 * d1 + b1 * d0) % cls._Q
        return cls._from_list(H)

    # ------------------------------------------------------------------
    # Compression and byte encoding (FIPS 203, Section 4.2)
    # ------------------------------------------------------------------

    @classmethod
    def _compress(cls, x, d):
        r"""
        Return the ``d``-bit compression of the coefficient ``x``.
        """
        return ((x << d) + (cls._Q // 2)) // cls._Q % (1 << d)

    @classmethod
    def _decompress(cls, y, d):
        r"""
        Return the coefficient obtained by decompressing ``y``,
        rounding half up as specified by FIPS 203.
        """
        return ((cls._Q * y) + (1 << (d - 1))) >> d

    @staticmethod
    def _encode(F, d):
        r"""
        Encode a coefficient list into ``32*d`` bytes in little-endian
        bit order.
        """
        result = bytearray()
        acc = 0
        bits = 0
        for c in F:
            acc |= c << bits
            bits += d
            while bits >= 8:
                result.append(acc & 0xFF)
                acc >>= 8
                bits -= 8
        if bits:
            result.append(acc & 0xFF)
        return bytes(result)

    @classmethod
    def _decode(cls, B, d):
        r"""
        Decode a coefficient list from ``32*d`` bytes in little-endian
        bit order.
        """
        m = cls._N * d // 8
        acc = int.from_bytes(B[:m], "little")
        mask = (1 << d) - 1
        F = []
        for _ in range(cls._N):
            F.append((acc & mask) % cls._Q)
            acc >>= d
        return F

    # ------------------------------------------------------------------
    # Symmetric primitives (FIPS 203, Section 4.1)
    # ------------------------------------------------------------------

    @staticmethod
    def _hash_h(s):
        r"""Return ``H(s) = SHA3-256(s)``."""
        return hashlib.sha3_256(s).digest()

    @staticmethod
    def _hash_j(s):
        r"""Return ``J(s) = SHAKE256(s, 32)``."""
        return hashlib.shake_256(s).digest(32)

    @staticmethod
    def _hash_g(s):
        r"""Return ``G(s) = SHA3-512(s)``."""
        return hashlib.sha3_512(s).digest()

    @staticmethod
    def _xof(rho, i, j, length):
        r"""Return ``length`` bytes from ``SHAKE128(rho || i || j)``."""
        return hashlib.shake_128(rho + bytes([i, j])).digest(length)

    @staticmethod
    def _prf(eta, s, b):
        r"""Return ``PRF_eta(s, b) = SHAKE256(s || b, 64*eta)``."""
        return hashlib.shake_256(s + bytes([b])).digest(64 * eta)

    # ------------------------------------------------------------------
    # Sampling (FIPS 203, Algorithms 7 and 8)
    # ------------------------------------------------------------------

    @classmethod
    def _sample_ntt(cls, rho, i, j):
        r"""
        Sample a polynomial in NTT domain by rejection sampling.

        Implements Algorithm 7 of FIPS 203.
        """
        stream = cls._xof(rho, i, j, 4096)
        coeffs = []
        idx = 0
        while len(coeffs) < cls._N:
            b0 = stream[idx]
            b1 = stream[idx + 1]
            b2 = stream[idx + 2]
            idx += 3
            d1 = b0 + 256 * (b1 % 16)
            d2 = (b1 // 16) + 16 * b2
            if d1 < cls._Q:
                coeffs.append(d1)
            if d2 < cls._Q and len(coeffs) < cls._N:
                coeffs.append(d2)
        return cls._from_list(coeffs)

    @classmethod
    def _sample_poly_cbd(cls, eta, B):
        r"""
        Sample a polynomial from the centered binomial distribution.

        Implements Algorithm 8 of FIPS 203.
        """
        b_int = int.from_bytes(B, "little")
        mask = (1 << eta) - 1
        mask2 = (1 << (2 * eta)) - 1
        coeffs = []
        for _ in range(cls._N):
            x = b_int & mask2
            a = (x & mask).bit_count()
            b = ((x >> eta) & mask).bit_count()
            coeffs.append((a - b) % cls._Q)
            b_int >>= 2 * eta
        return cls._from_list(coeffs)

    # ------------------------------------------------------------------
    # K-PKE (FIPS 203, Algorithms 13-15)
    # ------------------------------------------------------------------

    @classmethod
    def _kpke_keygen(cls, d, k, eta1):
        r"""
        Generate a K-PKE key pair from the 32-byte seed ``d``.

        Implements Algorithm 13 of FIPS 203.
        """
        rho_sigma = cls._hash_g(d + bytes([k]))
        rho, sigma = rho_sigma[:32], rho_sigma[32:]

        A = [[cls._sample_ntt(rho, j, i) for j in range(k)] for i in range(k)]

        s = []
        N = 0
        for _ in range(k):
            s.append(cls._sample_poly_cbd(eta1, cls._prf(eta1, sigma, N)))
            N += 1
        e = []
        for _ in range(k):
            e.append(cls._sample_poly_cbd(eta1, cls._prf(eta1, sigma, N)))
            N += 1

        s_hat = [cls._ntt(si) for si in s]
        e_hat = [cls._ntt(ei) for ei in e]

        t_hat = []
        for i in range(k):
            acc = cls._R(0)
            for j in range(k):
                acc = acc + cls._multiply_ntts(A[i][j], s_hat[j])
            t_hat.append(acc + e_hat[i])

        ek_pke = (
            b"".join(cls._encode(cls._coeffs(t_hat[i]), 12) for i in range(k)) + rho
        )
        dk_pke = b"".join(cls._encode(cls._coeffs(s_hat[i]), 12) for i in range(k))
        return ek_pke, dk_pke

    @classmethod
    def _kpke_encrypt(cls, ek_pke, m, r, k, eta1, eta2, du, dv):
        r"""
        Encrypt the 32-byte message ``m`` under ``ek_pke`` with
        randomness ``r``.

        Implements Algorithm 14 of FIPS 203.
        """
        t_hat_bytes = ek_pke[: 384 * k]
        rho = ek_pke[384 * k : 384 * k + 32]

        t_hat = [
            cls._from_list(cls._decode(t_hat_bytes[384 * i : 384 * (i + 1)], 12))
            for i in range(k)
        ]
        A = [[cls._sample_ntt(rho, j, i) for j in range(k)] for i in range(k)]

        y = []
        e1 = []
        N = 0
        for _ in range(k):
            y.append(cls._sample_poly_cbd(eta1, cls._prf(eta1, r, N)))
            N += 1
        for _ in range(k):
            e1.append(cls._sample_poly_cbd(eta2, cls._prf(eta2, r, N)))
            N += 1
        e2 = cls._sample_poly_cbd(eta2, cls._prf(eta2, r, N))

        y_hat = [cls._ntt(yi) for yi in y]

        u = []
        for i in range(k):
            acc = cls._R(0)
            for j in range(k):
                acc = acc + cls._multiply_ntts(A[j][i], y_hat[j])
            u.append(cls._ntt_inv(acc) + e1[i])

        acc = cls._R(0)
        for i in range(k):
            acc = acc + cls._multiply_ntts(t_hat[i], y_hat[i])
        v = cls._ntt_inv(acc) + e2

        mu = cls._from_list([cls._decompress(x, 1) for x in cls._decode(m, 1)])
        v = v + mu

        c1 = b"".join(
            cls._encode([cls._compress(x, du) for x in cls._coeffs(u[i])], du)
            for i in range(k)
        )
        c2 = cls._encode([cls._compress(x, dv) for x in cls._coeffs(v)], dv)
        return c1 + c2

    @classmethod
    def _kpke_decrypt(cls, dk_pke, c, k, du, dv):
        r"""
        Decrypt the ciphertext ``c`` using ``dk_pke``.

        Implements Algorithm 15 of FIPS 203.
        """
        c1_len = 32 * du * k
        c1, c2 = c[:c1_len], c[c1_len:]

        u = [
            cls._from_list(
                [
                    cls._decompress(x, du)
                    for x in cls._decode(c1[32 * du * i : 32 * du * (i + 1)], du)
                ]
            )
            for i in range(k)
        ]
        v = cls._from_list(
            [cls._decompress(x, dv) for x in cls._decode(c2, dv)]
        )

        s_hat = [
            cls._from_list(cls._decode(dk_pke[384 * i : 384 * (i + 1)], 12))
            for i in range(k)
        ]

        u_hat = [cls._ntt(ui) for ui in u]
        acc = cls._R(0)
        for i in range(k):
            acc = acc + cls._multiply_ntts(s_hat[i], u_hat[i])
        w = v - cls._ntt_inv(acc)

        return cls._encode([cls._compress(x, 1) for x in cls._coeffs(w)], 1)

    # ------------------------------------------------------------------
    # Public API (FIPS 203, Algorithms 19-21)
    # ------------------------------------------------------------------

    def keygen(self, d=None, z=None):
        r"""
        Generate a public/secret key pair.

        See Algorithm 19 in [FIPS203]_.

        INPUT:

        - ``d`` -- optional 32-byte seed for the keypair; if not provided
          a fresh seed is drawn from :func:`secrets.token_bytes`
        - ``z`` -- optional 32-byte seed used in the implicit rejection
          path of :meth:`decaps`; if not provided a fresh seed is drawn

        OUTPUT: tuple ``(public_key, secret_key)`` of :class:`bytes`

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set(512)
            sage: pk, sk = kem.keygen()
            sage: isinstance(pk, bytes) and isinstance(sk, bytes)
            True
        """
        if d is None:
            d = secrets.token_bytes(32)
        if z is None:
            z = secrets.token_bytes(32)
        if len(d) != 32 or len(z) != 32:
            raise ValueError("d and z must be 32 bytes each.")

        ek_pke, dk_pke = self._kpke_keygen(d, self.k, self.eta1)
        ek = ek_pke
        dk = dk_pke + ek + self._hash_h(ek) + z
        return ek, dk

    def encaps(self, public_key, m=None):
        r"""
        Encapsulate a shared secret.

        See Algorithm 20 in [FIPS203]_.

        INPUT:

        - ``public_key`` -- the recipient's public key (bytes)
        - ``m`` -- optional 32-byte message; if not provided, a fresh
          random message is drawn

        OUTPUT: tuple ``(ciphertext, shared_secret)`` of :class:`bytes`

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set(512)
            sage: pk, sk = kem.keygen()
            sage: ct, ss = kem.encaps(pk)
            sage: len(ss)
            32
        """
        if m is None:
            m = secrets.token_bytes(32)
        if len(m) != 32:
            raise ValueError("m must be 32 bytes.")

        ek = public_key
        kr = self._hash_g(m + self._hash_h(ek))
        K, r = kr[:32], kr[32:]
        c = self._kpke_encrypt(
            ek, m, r, self.k, self.eta1, self.eta2, self.du, self.dv
        )
        return c, K

    def decaps(self, secret_key, ciphertext):
        r"""
        Decapsulate a ciphertext to recover the shared secret.

        See Algorithm 21 in [FIPS203]_.

        INPUT:

        - ``secret_key`` -- the recipient's secret key (bytes)
        - ``ciphertext`` -- a ciphertext produced by :meth:`encaps`

        OUTPUT: 32-byte shared secret

        EXAMPLES::

            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set(512)
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True
        """
        dk = secret_key
        c = ciphertext
        k = self.k

        dk_pke = dk[: 384 * k]
        ek_pke = dk[384 * k : 768 * k + 32]
        h = dk[768 * k + 32 : 768 * k + 64]
        z = dk[768 * k + 64 : 768 * k + 96]

        m_prime = self._kpke_decrypt(dk_pke, c, k, self.du, self.dv)
        kr_prime = self._hash_g(m_prime + h)
        K_prime, r_prime = kr_prime[:32], kr_prime[32:]
        K_bar = self._hash_j(z + c)

        c_prime = self._kpke_encrypt(
            ek_pke, m_prime, r_prime, k, self.eta1, self.eta2, self.du, self.dv
        )
        if c == c_prime:
            return K_prime
        return K_bar
