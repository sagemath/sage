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

# Copyright (C) 2026 Vivek Kumar <viveksantoshkumardubey@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# https://www.gnu.org/licenses/

import hashlib
import secrets

from sage.arith.misc import is_prime
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

    The implementation is primarily intended for prototyping and
    educational purposes.  The constructor accepts arbitrary ``n`` and
    ``q`` (subject to the arithmetic requirements of the scheme), so
    small toy parameters can be used to illustrate the algorithm.

    Keys and ciphertexts are represented as mathematical objects:
    elements of the ring `R = \GF(q)[x] / (x^{n} + 1)` and vectors
    thereof.  Byte-string serialization, as specified by FIPS 203, is
    kept separate and is only used as an input to the symmetric
    primitives `H`, `J`, `G` and for interoperability with other
    implementations.

    EXAMPLES::

        sage: import warnings; warnings.simplefilter('ignore')
        sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
        sage: kem = MLKEM.from_parameter_set("512")
        sage: pk, sk = kem.keygen()
        sage: ct, ss1 = kem.encaps(pk)
        sage: ss2 = kem.decaps(sk, ct)
        sage: ss1 == ss2
        True

    The constructor accepts arbitrary ``n`` and ``q`` (subject to the
    arithmetic requirements), so toy parameters can be used to
    illustrate the algorithm::

        sage: import warnings; warnings.simplefilter('ignore')
        sage: kem = MLKEM.from_parameter_set("toy")
        sage: pk, sk = kem.keygen()
        sage: ct, ss1 = kem.encaps(pk)
        sage: ss2 = kem.decaps(sk, ct)
        sage: ss1 == ss2
        True
    """

    PARAMETER_SETS = {
        "512":  {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2, 'du': 10, 'dv': 4},
        "768":  {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2, 'du': 10, 'dv': 4},
        "1024": {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2, 'du': 11, 'dv': 5},
        "toy":  {'n': 32, 'q': 3329, 'k': 2, 'eta1': 2, 'eta2': 2, 'du': 10, 'dv': 4},
    }

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_parameter_set(cls, parameter_set):
        r"""
        Construct an :class:`MLKEM` object from a named parameter set.

        INPUT:

        - ``parameter_set`` -- string; one of ``"512"``, ``"768"``,
          ``"1024"``, or ``"toy"``

        EXAMPLES::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set("768")
            sage: kem.k
            3
            sage: MLKEM.from_parameter_set("toy").n
            32
        """
        if parameter_set not in cls.PARAMETER_SETS:
            raise ValueError(
                f"parameter_set must be one of {list(cls.PARAMETER_SETS)}, "
                f"got {parameter_set!r}"
            )
        return cls(**cls.PARAMETER_SETS[parameter_set])

    def __init__(self, n=256, q=3329, k=2, eta1=3, eta2=2, du=10, dv=4):
        r"""
        Initialize an ML-KEM instance with explicit parameters.

        INPUT:

        - ``n`` -- integer (default: 256), ring degree; must be a power of 2
        - ``q`` -- integer (default: 3329), coefficient modulus; must be prime
        - ``k`` -- integer (default: 2), module rank
        - ``eta1`` -- integer (default: 3), CBD parameter for the secret
          and error vectors
        - ``eta2`` -- integer (default: 2), CBD parameter for the
          encapsulation noise
        - ``du`` -- integer (default: 10), compression parameter for ``u``
        - ``dv`` -- integer (default: 4), compression parameter for ``v``

        The parameters must satisfy the arithmetic requirements of the
        scheme: ``n`` must be a power of 2 and at least 8 (so the
        message length ``n / 8`` is positive), ``q`` must be prime and
        strictly less than 4096 (so that coefficients fit in 12 bits
        as required by :meth:`_encode_ek` and by the rejection sampler),
        and ``n`` must divide ``q - 1`` so that the negacyclic NTT is
        defined.

        EXAMPLES::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM(n=256, q=3329, k=2, eta1=3, eta2=2, du=10, dv=4)
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True

        Toy parameters are also accepted::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: kem = MLKEM.from_parameter_set("toy")
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True

        Invalid parameters are rejected::

            sage: MLKEM(n=255, q=3329)
            Traceback (most recent call last):
            ...
            ValueError: n must be a power of 2
            sage: MLKEM(n=256, q=3328)
            Traceback (most recent call last):
            ...
            ValueError: q must be prime
        """
        if n < 8:
            raise ValueError("n must be at least 8")
        if (n & (n - 1)) != 0:
            raise ValueError("n must be a power of 2")
        if not is_prime(q):
            raise ValueError("q must be prime")
        if q >= 4096:
            raise ValueError("q must be less than 4096 (12-bit coefficient encoding)")
        if (q - 1) % n != 0:
            raise ValueError("n must divide q - 1 (required for the NTT)")
        if eta1 not in (2, 3) or eta2 not in (2, 3):
            raise ValueError("eta1 and eta2 must be 2 or 3")
        if k < 1:
            raise ValueError("k must be positive")

        self.n = n
        self.q = q
        self.k = k
        self.eta1 = eta1
        self.eta2 = eta2
        self.du = du
        self.dv = dv

        self._FF = FiniteField(q)
        self._POLY_RING = PolynomialRing(self._FF, 'x')
        self._R = self._POLY_RING.quotient(
            self._POLY_RING.gen() ** n + 1, 'x'
        )
        self._ZETAS = self._compute_zetas()

    # ------------------------------------------------------------------
    # Twiddle factors
    # ------------------------------------------------------------------

    def _compute_zetas(self):
        r"""
        Compute the twiddle factors for the negacyclic NTT.

        Returns the list ``[zeta^BitRev_{log2(n)-1}(i) mod q for i in range(n//2)]``,
        where ``zeta`` is a primitive ``n``-th root of unity in ``GF(q)``.

        For the FIPS 203 parameters ``(n, q) = (256, 3329)`` this returns
        the list specified by the standard (with ``zeta = 17``).
        """
        n, q = self.n, self.q
        half = n // 2
        logn = n.bit_length() - 1

        if (n, q) == (256, 3329):
            zeta = 17
        else:
            FF = self._FF
            for candidate in range(2, q):
                z = FF(candidate)
                if z ** n == 1 and z ** half == -1:
                    zeta = candidate
                    break
            else:
                raise ValueError(
                    f"no primitive {n}-th root of unity found in GF({q})"
                )

        zetas = []
        for i in range(half):
            x = i
            br = 0
            for _ in range(logn - 1):
                br = (br << 1) | (x & 1)
                x >>= 1
            zetas.append(pow(zeta, br, q))
        return zetas

    # ------------------------------------------------------------------
    # Polynomial helpers
    # ------------------------------------------------------------------

    def _from_list(self, L):
        r"""
        Return the element of ``R`` with the given coefficient list.
        """
        return self._R(self._POLY_RING(L))

    def _coeffs(self, f):
        r"""
        Return the length-``n`` coefficient list of an element of ``R``.
        """
        L = [int(c) for c in f.lift().list()]
        if len(L) < self.n:
            L += [0] * (self.n - len(L))
        return L[: self.n]

    # ------------------------------------------------------------------
    # Negacyclic NTT (FIPS 203, Section 4.3)
    # ------------------------------------------------------------------

    def _ntt(self, f):
        r"""
        Return the negacyclic NTT of ``f``.

        INPUT:

        - ``f`` -- element of ``R``

        OUTPUT: element of ``R`` representing the NTT of ``f``
        """
        F = self._coeffs(f)
        i = 1
        length = self.n // 2
        while length >= 2:
            for start in range(0, self.n, 2 * length):
                zeta = self._ZETAS[i]
                i += 1
                for j in range(start, start + length):
                    t = (zeta * F[j + length]) % self.q
                    F[j + length] = (F[j] - t) % self.q
                    F[j] = (F[j] + t) % self.q
            length //= 2
        return self._from_list(F)

    def _ntt_inv(self, f):
        r"""
        Return the inverse of :meth:`_ntt`.
        """
        F = self._coeffs(f)
        i = self.n // 2 - 1
        length = 2
        while length <= self.n // 2:
            for start in range(0, self.n, 2 * length):
                zeta = self._ZETAS[i]
                i -= 1
                for j in range(start, start + length):
                    t = F[j]
                    F[j] = (t + F[j + length]) % self.q
                    F[j + length] = (
                        zeta * (F[j + length] - t)
                    ) % self.q
            length *= 2
        # Scale by (n/2)^{-1} mod q
        n_inv = pow(self.n // 2, -1, self.q)
        F = [(n_inv * x) % self.q for x in F]
        return self._from_list(F)

    # ------------------------------------------------------------------
    # Compression (FIPS 203, Section 4.2)
    # ------------------------------------------------------------------

    def _compress(self, x, d):
        r"""
        Return the ``d``-bit compression of the coefficient ``x``.
        """
        return ((x << d) + (self.q // 2)) // self.q % (1 << d)

    def _decompress(self, y, d):
        r"""
        Return the coefficient obtained by decompressing ``y``,
        rounding half up as specified by FIPS 203.
        """
        return ((self.q * y) + (1 << (d - 1))) >> d

    # ------------------------------------------------------------------
    # Byte encoding (FIPS 203, Section 4.2)
    # ------------------------------------------------------------------

    @staticmethod
    def _encode(F, d):
        r"""
        Encode a coefficient list into ``ceil(len(F) * d / 8)`` bytes in
        little-endian bit order.
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

    def _decode(self, B, d):
        r"""
        Decode ``n`` coefficients from ``ceil(n*d/8)`` bytes in
        little-endian bit order.
        """
        nbytes = (self.n * d + 7) // 8
        acc = int.from_bytes(B[:nbytes], "little")
        mask = (1 << d) - 1
        F = []
        for _ in range(self.n):
            F.append((acc & mask) % self.q)
            acc >>= d
        return F

    # ------------------------------------------------------------------
    # FIPS 203-specific serialization
    # ------------------------------------------------------------------

    def _encode_ek(self, rho, t):
        r"""
        Serialize the public key ``(rho, t)`` to FIPS 203 bytes.

        FIPS 203 serializes ``t`` in NTT form, so we apply the forward
        NTT before byte-encoding.
        """
        parts = []
        for ti in t:
            ti_hat = self._ntt(ti)
            parts.append(self._encode(self._coeffs(ti_hat), 12))
        return b"".join(parts) + rho

    def _encode_ct(self, u, v):
        r"""
        Serialize the ciphertext ``(u, v)`` to FIPS 203 bytes.
        """
        c1_parts = []
        for ui in u:
            coeffs = self._coeffs(ui)
            compressed = [self._compress(x, self.du) for x in coeffs]
            c1_parts.append(self._encode(compressed, self.du))
        v_coeffs = self._coeffs(v)
        v_compressed = [self._compress(x, self.dv) for x in v_coeffs]
        c2 = self._encode(v_compressed, self.dv)
        return b"".join(c1_parts) + c2

    def _compress_ct(self, u, v):
        r"""
        Apply Compress followed by Decompress to the ciphertext
        polynomials ``(u, v)``.

        FIPS 203 transmits the ciphertext in compressed form, i.e.
        ``ByteEncode_du(Compress_du(u)) || ByteEncode_dv(Compress_dv(v))``.
        We model this by applying ``Decompress(Compress(-, du))`` to
        each coefficient of ``u`` and ``Decompress(Compress(-, dv))``
        to each coefficient of ``v``, so that decapsulation operates
        on the same data as a real recipient would receive.  This
        exposes any decryption failure caused by the compression loss.
        """
        u_compressed = [
            self._from_list([
                self._decompress(self._compress(x, self.du), self.du)
                for x in self._coeffs(ui)
            ])
            for ui in u
        ]
        v_compressed = self._from_list([
            self._decompress(self._compress(x, self.dv), self.dv)
            for x in self._coeffs(v)
        ])
        return u_compressed, v_compressed

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

    def _prf(self, eta, s, b):
        r"""Return ``PRF_eta(s, b)``."""
        return hashlib.shake_256(s + bytes([b])).digest(self.n * eta // 4)

    # ------------------------------------------------------------------
    # Sampling (FIPS 203, Algorithms 7 and 8)
    # ------------------------------------------------------------------

    def _sample_ntt(self, rho, i, j):
        r"""
        Sample a polynomial in NTT domain by rejection sampling.

        Implements Algorithm 7 of FIPS 203.
        """
        length = max(1024, 3 * self.n * (4096 // self.q + 10))
        stream = self._xof(rho, i, j, length)
        coeffs = []
        idx = 0
        while len(coeffs) < self.n:
            if idx + 2 >= len(stream):
                length *= 2
                stream = self._xof(rho, i, j, length)
            b0 = stream[idx]
            b1 = stream[idx + 1]
            b2 = stream[idx + 2]
            idx += 3
            d1 = b0 + 256 * (b1 % 16)
            d2 = (b1 // 16) + 16 * b2
            if d1 < self.q:
                coeffs.append(d1)
            if d2 < self.q and len(coeffs) < self.n:
                coeffs.append(d2)
        return self._from_list(coeffs[: self.n])

    def _sample_poly_cbd(self, eta, B):
        r"""
        Sample a polynomial from the centered binomial distribution.

        Implements Algorithm 8 of FIPS 203.
        """
        b_int = int.from_bytes(B, "little")
        mask = (1 << eta) - 1
        mask2 = (1 << (2 * eta)) - 1
        coeffs = []
        for _ in range(self.n):
            x = b_int & mask2
            a = (x & mask).bit_count()
            b = ((x >> eta) & mask).bit_count()
            coeffs.append((a - b) % self.q)
            b_int >>= 2 * eta
        return self._from_list(coeffs)

    # ------------------------------------------------------------------
    # K-PKE (FIPS 203, Algorithms 13-15), plain polynomial arithmetic
    # ------------------------------------------------------------------

    def _kpke_keygen(self, d):
        r"""
        Generate a K-PKE key pair from the 32-byte seed ``d``.

        Implements Algorithm 13 of FIPS 203, using plain polynomial
        arithmetic in ``R``.

        OUTPUT: tuple ``(rho, t, s)``, where ``rho`` is 32 bytes and
        ``t``, ``s`` are lists of ``k`` elements of ``R``.
        """
        rho_sigma = self._hash_g(d + bytes([self.k]))
        rho, sigma = rho_sigma[:32], rho_sigma[32:]

        # A_hat[i][j] = SampleNTT(rho, j, i), then A = NTT^{-1}(A_hat)
        A_hat = [
            [self._sample_ntt(rho, j, i) for j in range(self.k)]
            for i in range(self.k)
        ]
        A = [[self._ntt_inv(A_hat[i][j]) for j in range(self.k)] for i in range(self.k)]

        s = [
            self._sample_poly_cbd(self.eta1, self._prf(self.eta1, sigma, i))
            for i in range(self.k)
        ]
        e = [
            self._sample_poly_cbd(
                self.eta1, self._prf(self.eta1, sigma, self.k + i)
            )
            for i in range(self.k)
        ]

        # t = A * s + e
        t = []
        for i in range(self.k):
            acc = self._R(0)
            for j in range(self.k):
                acc = acc + A[i][j] * s[j]
            t.append(acc + e[i])

        return rho, t, s

    def _kpke_encrypt(self, rho, t, m, r):
        r"""
        Encrypt the message ``m`` under the public key ``(rho, t)``.

        Implements Algorithm 14 of FIPS 203, using plain polynomial
        arithmetic in ``R``.

        OUTPUT: tuple ``(u, v)``, where ``u`` is a list of ``k``
        elements of ``R`` and ``v`` is an element of ``R``.  These are
        the *uncompressed* K-PKE ciphertext polynomials.
        """
        A_hat = [
            [self._sample_ntt(rho, j, i) for j in range(self.k)]
            for i in range(self.k)
        ]
        A = [[self._ntt_inv(A_hat[i][j]) for j in range(self.k)] for i in range(self.k)]

        y = [
            self._sample_poly_cbd(self.eta1, self._prf(self.eta1, r, i))
            for i in range(self.k)
        ]
        e1 = [
            self._sample_poly_cbd(self.eta2, self._prf(self.eta2, r, self.k + i))
            for i in range(self.k)
        ]
        e2 = self._sample_poly_cbd(self.eta2, self._prf(self.eta2, r, 2 * self.k))

        # u = A^T * y + e1
        u = []
        for i in range(self.k):
            acc = self._R(0)
            for j in range(self.k):
                acc = acc + A[j][i] * y[j]
            u.append(acc + e1[i])

        # v = t^T * y + e2 + mu
        acc = self._R(0)
        for i in range(self.k):
            acc = acc + t[i] * y[i]
        v = acc + e2

        m_coeffs = self._decode(m, 1)
        mu = self._from_list([self._decompress(x, 1) for x in m_coeffs])
        v = v + mu

        return u, v

    def _kpke_decrypt(self, s, u, v):
        r"""
        Decrypt the ciphertext ``(u, v)`` using the secret vector ``s``.

        Implements Algorithm 15 of FIPS 203, using plain polynomial
        arithmetic in ``R``.

        OUTPUT: the recovered message as bytes
        """
        # w = v - s^T * u
        acc = self._R(0)
        for i in range(self.k):
            acc = acc + s[i] * u[i]
        w = v - acc

        w_coeffs = self._coeffs(w)
        return self._encode([self._compress(x, 1) for x in w_coeffs], 1)

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

        OUTPUT: tuple ``(pk, sk)`` where

        - ``pk = (rho, t)``: ``rho`` is 32 bytes, and ``t`` is a list of
          ``k`` elements of ``R``.
        - ``sk = (s, rho, t, z)``: ``s`` is a list of ``k`` elements of
          ``R``, and ``z`` is 32 bytes.

        EXAMPLES::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set("512")
            sage: pk, sk = kem.keygen()
            sage: rho, t = pk
            sage: len(rho)
            32
            sage: len(t)
            2
        """
        if d is None:
            d = secrets.token_bytes(32)
        if z is None:
            z = secrets.token_bytes(32)
        if len(d) != 32 or len(z) != 32:
            raise ValueError("d and z must be 32 bytes each.")

        rho, t, s = self._kpke_keygen(d)
        return (rho, t), (s, rho, t, z)

    def encaps(self, public_key, m=None):
        r"""
        Encapsulate a shared secret.

        See Algorithm 20 in [FIPS203]_.

        INPUT:

        - ``public_key`` -- the recipient's public key, i.e. a tuple
          ``(rho, t)`` as returned by :meth:`keygen`
        - ``m`` -- optional message of ``n / 8`` bytes; if not provided,
          a fresh random message is drawn

        OUTPUT: tuple ``((u, v), K)``, where ``(u, v)`` is the
        ciphertext (a list of ``k`` elements of ``R`` followed by an
        element of ``R``) and ``K`` is the 32-byte shared secret.

        The ciphertext polynomials are returned in their compressed
        (and then decompressed) form, matching what a real recipient
        would receive.

        EXAMPLES::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set("512")
            sage: pk, sk = kem.keygen()
            sage: ct, ss = kem.encaps(pk)
            sage: len(ss)
            32
        """
        msg_len = self.n // 8
        if m is None:
            m = secrets.token_bytes(msg_len)
        if len(m) != msg_len:
            raise ValueError(f"m must be {msg_len} bytes.")

        rho, t = public_key
        ek = self._encode_ek(rho, t)
        kr = self._hash_g(m + self._hash_h(ek))
        K, r = kr[:32], kr[32:]
        u, v = self._kpke_encrypt(rho, t, m, r)
        u, v = self._compress_ct(u, v)
        return (u, v), K

    def decaps(self, secret_key, ciphertext):
        r"""
        Decapsulate a ciphertext to recover the shared secret.

        See Algorithm 21 in [FIPS203]_.

        INPUT:

        - ``secret_key`` -- the recipient's secret key, i.e. a tuple
          ``(s, rho, t, z)`` as returned by :meth:`keygen`
        - ``ciphertext`` -- a tuple ``(u, v)`` produced by :meth:`encaps`

        OUTPUT: 32-byte shared secret

        EXAMPLES::

            sage: import warnings; warnings.simplefilter('ignore')
            sage: from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM
            sage: kem = MLKEM.from_parameter_set("512")
            sage: pk, sk = kem.keygen()
            sage: ct, ss1 = kem.encaps(pk)
            sage: ss2 = kem.decaps(sk, ct)
            sage: ss1 == ss2
            True
        """
        s, rho, t, z = secret_key
        u, v = ciphertext

        ek = self._encode_ek(rho, t)
        h = self._hash_h(ek)

        m_prime = self._kpke_decrypt(s, u, v)
        kr_prime = self._hash_g(m_prime + h)
        K_prime, r_prime = kr_prime[:32], kr_prime[32:]

        c = self._encode_ct(u, v)
        u_prime, v_prime = self._kpke_encrypt(rho, t, m_prime, r_prime)
        if c == self._encode_ct(u_prime, v_prime):
            return K_prime
        return self._hash_j(z + c)
