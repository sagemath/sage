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

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

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


# ----------------------------------------------------------------------
# Negacyclic NTT (FIPS 203, Section 4.3)
# ----------------------------------------------------------------------

def _ntt(f):
    r"""
    Return the negacyclic NTT of ``f`` in ``Z_q[x] / (x^256 + 1)``.

    INPUT:

    - ``f`` -- list of 256 integers modulo ``q``

    OUTPUT: list of 256 integers modulo ``q``, in NTT domain
    """
    f = list(f)
    i = 1
    length = 128
    while length >= 2:
        for start in range(0, _N, 2 * length):
            zeta = _ZETAS[i]
            i += 1
            for j in range(start, start + length):
                t = (zeta * f[j + length]) % _Q
                f[j + length] = (f[j] - t) % _Q
                f[j] = (f[j] + t) % _Q
        length //= 2
    return f


def _ntt_inv(f):
    r"""
    Return the inverse of :func:`_ntt`.

    INPUT:

    - ``f`` -- list of 256 integers modulo ``q``, in NTT domain

    OUTPUT: list of 256 integers modulo ``q``
    """
    f = list(f)
    i = 127
    length = 2
    while length <= 128:
        for start in range(0, _N, 2 * length):
            zeta = _ZETAS[i]
            i -= 1
            for j in range(start, start + length):
                t = f[j]
                f[j] = (t + f[j + length]) % _Q
                f[j + length] = (zeta * (f[j + length] - t)) % _Q
        length *= 2
    # Scale by 128^{-1} mod q.
    return [(3303 * x) % _Q for x in f]


def _base_case_multiply(a0, a1, b0, b1, gamma):
    r"""
    Multiply ``a0 + a1*x`` and ``b0 + b1*x`` in ``Z_q[x] / (x^2 - gamma)``.

    INPUT:

    - ``a0``, ``a1``, ``b0``, ``b1`` -- integers modulo ``q``
    - ``gamma`` -- integer modulo ``q``

    OUTPUT: tuple ``(c0, c1)`` of integers modulo ``q``
    """
    c0 = (a0 * b0 + a1 * b1 * gamma) % _Q
    c1 = (a0 * b1 + a1 * b0) % _Q
    return c0, c1


def _multiply_ntts(f, g):
    r"""
    Return the product of ``f`` and ``g`` in NTT domain.

    INPUT:

    - ``f``, ``g`` -- lists of 256 integers modulo ``q``, in NTT domain

    OUTPUT: list of 256 integers modulo ``q``, in NTT domain
    """
    h = [0] * _N
    for i in range(64):
        a0, a1 = f[4 * i], f[4 * i + 1]
        b0, b1 = f[4 * i + 2], f[4 * i + 3]
        c0, c1 = g[4 * i], g[4 * i + 1]
        d0, d1 = g[4 * i + 2], g[4 * i + 3]

        gamma = _ZETAS[64 + i]

        h0, h1 = _base_case_multiply(a0, a1, c0, c1, gamma)
        h2, h3 = _base_case_multiply(b0, b1, d0, d1, -gamma % _Q)

        h[4 * i] = h0
        h[4 * i + 1] = h1
        h[4 * i + 2] = h2
        h[4 * i + 3] = h3
    return h


# ----------------------------------------------------------------------
# Compression and byte encoding (FIPS 203, Section 4.2)
# ----------------------------------------------------------------------

def _compress(x, d):
    r"""
    Return the ``d``-bit compression of the coefficient ``x``.

    INPUT:

    - ``x`` -- integer modulo ``q``
    - ``d`` -- positive integer

    OUTPUT: integer in ``range(2**d)``
    """
    return ((x << d) + (_Q // 2)) // _Q % (1 << d)


def _decompress(y, d):
    r"""
    Return the coefficient obtained by decompressing the ``d``-bit
    integer ``y``, rounding half up as specified by FIPS 203.

    INPUT:

    - ``y`` -- integer in ``range(2**d)``
    - ``d`` -- positive integer

    OUTPUT: integer modulo ``q``
    """
    return ((_Q * y) + (1 << (d - 1))) >> d


def _encode(F, d):
    r"""
    Encode the polynomial ``F`` into ``32*d`` bytes, in little-endian
    bit order.

    INPUT:

    - ``F`` -- list of 256 integers in ``range(2**d)``
    - ``d`` -- positive integer

    OUTPUT: bytes of length ``32*d``
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


def _decode(B, d):
    r"""
    Decode a polynomial from ``32*d`` bytes, in little-endian bit order.

    INPUT:

    - ``B`` -- bytes of length ``32*d``
    - ``d`` -- positive integer

    OUTPUT: list of 256 integers modulo ``q``
    """
    m = _N * d // 8
    acc = int.from_bytes(B[:m], "little")
    mask = (1 << d) - 1
    F = []
    for _ in range(_N):
        F.append((acc & mask) % _Q)
        acc >>= d
    return F


# ----------------------------------------------------------------------
# Symmetric primitives (FIPS 203, Section 4.1)
# ----------------------------------------------------------------------

def _hash_h(s):
    r"""
    Return ``H(s) = SHA3-256(s)``.
    """
    return hashlib.sha3_256(s).digest()


def _hash_j(s):
    r"""
    Return ``J(s) = SHAKE256(s, 32)``.
    """
    return hashlib.shake_256(s).digest(32)


def _hash_g(s):
    r"""
    Return ``G(s) = SHA3-512(s)``.
    """
    return hashlib.sha3_512(s).digest()


def _xof(rho, i, j, length):
    r"""
    Return ``length`` bytes from ``SHAKE128(rho || i || j)``.

    INPUT:

    - ``rho`` -- 32 bytes
    - ``i``, ``j`` -- integers in ``range(256)``
    - ``length`` -- positive integer
    """
    return hashlib.shake_128(rho + bytes([i, j])).digest(length)


def _prf(eta, s, b):
    r"""
    Return ``PRF_eta(s, b) = SHAKE256(s || b, 64*eta)``.

    INPUT:

    - ``eta`` -- 2 or 3
    - ``s`` -- 32 bytes
    - ``b`` -- integer in ``range(256)``
    """
    return hashlib.shake_256(s + bytes([b])).digest(64 * eta)


# ----------------------------------------------------------------------
# Sampling (FIPS 203, Algorithms 7 and 8)
# ----------------------------------------------------------------------

def _sample_ntt(rho, i, j):
    r"""
    Sample a polynomial in NTT domain by rejection sampling.

    Implements Algorithm 7 of FIPS 203.

    INPUT:

    - ``rho`` -- 32 bytes
    - ``i``, ``j`` -- integers in ``range(256)``

    OUTPUT: list of 256 integers modulo ``q``, in NTT domain
    """
    stream = _xof(rho, i, j, 4096)
    coeffs = []
    idx = 0
    while len(coeffs) < _N:
        b0 = stream[idx]
        b1 = stream[idx + 1]
        b2 = stream[idx + 2]
        idx += 3
        d1 = b0 + 256 * (b1 % 16)
        d2 = (b1 // 16) + 16 * b2
        if d1 < _Q:
            coeffs.append(d1)
        if d2 < _Q and len(coeffs) < _N:
            coeffs.append(d2)
    return coeffs


def _sample_poly_cbd(eta, B):
    r"""
    Sample a polynomial from the centered binomial distribution.

    Implements Algorithm 8 of FIPS 203.

    INPUT:

    - ``eta`` -- 2 or 3
    - ``B`` -- bytes of length ``64*eta``

    OUTPUT: list of 256 integers modulo ``q``
    """
    b_int = int.from_bytes(B, "little")
    mask = (1 << eta) - 1
    mask2 = (1 << (2 * eta)) - 1
    coeffs = []
    for _ in range(_N):
        x = b_int & mask2
        a = (x & mask).bit_count()
        b = ((x >> eta) & mask).bit_count()
        coeffs.append((a - b) % _Q)
        b_int >>= 2 * eta
    return coeffs


# ----------------------------------------------------------------------
# Polynomial arithmetic
# ----------------------------------------------------------------------

def _poly_add(f, g):
    r"""
    Return the coefficient-wise sum of ``f`` and ``g``, reduced mod ``q``.
    """
    return [(a + b) % _Q for a, b in zip(f, g)]


def _poly_sub(f, g):
    r"""
    Return the coefficient-wise difference of ``f`` and ``g``, reduced
    mod ``q``.
    """
    return [(a - b) % _Q for a, b in zip(f, g)]


def _matrix_a(rho, k):
    r"""
    Return the ``k x k`` matrix ``A`` in NTT domain, where
    ``A[i][j] = SampleNTT(rho, j, i)``.
    """
    return [[_sample_ntt(rho, j, i) for j in range(k)] for i in range(k)]


# ----------------------------------------------------------------------
# K-PKE (FIPS 203, Algorithms 13-15)
# ----------------------------------------------------------------------

def _kpke_keygen(d, k, eta1):
    r"""
    Generate a K-PKE key pair from the 32-byte seed ``d``.

    Implements Algorithm 13 of FIPS 203.

    OUTPUT: tuple ``(ek_pke, dk_pke)`` of bytes
    """
    rho_sigma = _hash_g(d + bytes([k]))
    rho, sigma = rho_sigma[:32], rho_sigma[32:]

    A = _matrix_a(rho, k)

    s = []
    e = []
    N = 0
    for _ in range(k):
        s.append(_sample_poly_cbd(eta1, _prf(eta1, sigma, N)))
        N += 1
    for _ in range(k):
        e.append(_sample_poly_cbd(eta1, _prf(eta1, sigma, N)))
        N += 1

    s_hat = [_ntt(si) for si in s]
    e_hat = [_ntt(ei) for ei in e]

    t_hat = []
    for i in range(k):
        acc = [0] * _N
        for j in range(k):
            acc = _poly_add(acc, _multiply_ntts(A[i][j], s_hat[j]))
        t_hat.append(_poly_add(acc, e_hat[i]))

    ek_pke = b"".join(_encode(t_hat[i], 12) for i in range(k)) + rho
    dk_pke = b"".join(_encode(s_hat[i], 12) for i in range(k))
    return ek_pke, dk_pke


def _kpke_encrypt(ek_pke, m, r, k, eta1, eta2, du, dv):
    r"""
    Encrypt the 32-byte message ``m`` under ``ek_pke`` with randomness
    ``r``.

    Implements Algorithm 14 of FIPS 203.

    OUTPUT: ciphertext as bytes
    """
    t_hat_bytes = ek_pke[: 384 * k]
    rho = ek_pke[384 * k : 384 * k + 32]

    t_hat = [_decode(t_hat_bytes[384 * i : 384 * (i + 1)], 12) for i in range(k)]
    A = _matrix_a(rho, k)

    y = []
    e1 = []
    N = 0
    for _ in range(k):
        y.append(_sample_poly_cbd(eta1, _prf(eta1, r, N)))
        N += 1
    for _ in range(k):
        e1.append(_sample_poly_cbd(eta2, _prf(eta2, r, N)))
        N += 1
    e2 = _sample_poly_cbd(eta2, _prf(eta2, r, N))

    y_hat = [_ntt(yi) for yi in y]

    u = []
    for i in range(k):
        acc = [0] * _N
        for j in range(k):
            acc = _poly_add(acc, _multiply_ntts(A[j][i], y_hat[j]))
        u.append(_poly_add(_ntt_inv(acc), e1[i]))

    acc = [0] * _N
    for i in range(k):
        acc = _poly_add(acc, _multiply_ntts(t_hat[i], y_hat[i]))
    v = _poly_add(_ntt_inv(acc), e2)

    mu = [_decompress(x, 1) for x in _decode(m, 1)]
    v = _poly_add(v, mu)

    c1 = b"".join(_encode([_compress(x, du) for x in u[i]], du) for i in range(k))
    c2 = _encode([_compress(x, dv) for x in v], dv)
    return c1 + c2


def _kpke_decrypt(dk_pke, c, k, du, dv):
    r"""
    Decrypt the ciphertext ``c`` using ``dk_pke``.

    Implements Algorithm 15 of FIPS 203.

    OUTPUT: 32-byte message
    """
    c1_len = 32 * du * k
    c1, c2 = c[:c1_len], c[c1_len:]

    u = [
        [_decompress(x, du) for x in _decode(c1[32 * du * i : 32 * du * (i + 1)], du)]
        for i in range(k)
    ]
    v = [_decompress(x, dv) for x in _decode(c2, dv)]

    s_hat = [_decode(dk_pke[384 * i : 384 * (i + 1)], 12) for i in range(k)]

    u_hat = [_ntt(ui) for ui in u]
    acc = [0] * _N
    for i in range(k):
        acc = _poly_add(acc, _multiply_ntts(s_hat[i], u_hat[i]))
    w = _poly_sub(v, _ntt_inv(acc))

    return _encode([_compress(x, 1) for x in w], 1)


# ----------------------------------------------------------------------
# ML-KEM (FIPS 203, Algorithms 19-21)
# ----------------------------------------------------------------------

def _mlkem_keygen(d, z, k, eta1, eta2, du, dv):
    r"""
    Generate an ML-KEM key pair from the seeds ``d`` and ``z``.

    Implements Algorithm 19 of FIPS 203.

    OUTPUT: tuple ``(ek, dk)`` of bytes
    """
    ek_pke, dk_pke = _kpke_keygen(d, k, eta1)
    ek = ek_pke
    dk = dk_pke + ek + _hash_h(ek) + z
    return ek, dk


def _mlkem_encaps(ek, m, k, eta1, eta2, du, dv):
    r"""
    Encapsulate a shared secret under the public key ``ek``.

    Implements Algorithm 20 of FIPS 203.

    OUTPUT: tuple ``(c, K)``, where ``c`` is the ciphertext and ``K`` is
    the 32-byte shared secret
    """
    kr = _hash_g(m + _hash_h(ek))
    K, r = kr[:32], kr[32:]
    c = _kpke_encrypt(ek, m, r, k, eta1, eta2, du, dv)
    return c, K


def _mlkem_decaps(dk, c, k, eta1, eta2, du, dv):
    r"""
    Decapsulate the ciphertext ``c`` using the secret key ``dk``.

    Implements Algorithm 21 of FIPS 203, including the implicit-rejection
    path for malformed ciphertexts.

    OUTPUT: 32-byte shared secret
    """
    dk_pke = dk[: 384 * k]
    ek_pke = dk[384 * k : 768 * k + 32]
    h = dk[768 * k + 32 : 768 * k + 64]
    z = dk[768 * k + 64 : 768 * k + 96]

    m_prime = _kpke_decrypt(dk_pke, c, k, du, dv)
    kr_prime = _hash_g(m_prime + h)
    K_prime, r_prime = kr_prime[:32], kr_prime[32:]
    K_bar = _hash_j(z + c)

    c_prime = _kpke_encrypt(ek_pke, m_prime, r_prime, k, eta1, eta2, du, dv)
    if c == c_prime:
        return K_prime
    return K_bar


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

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

    PARAMETER_SETS = {
        512:  {'n': 256, 'q': 3329, 'k': 2, 'eta1': 3, 'eta2': 2, 'du': 10, 'dv': 4},
        768:  {'n': 256, 'q': 3329, 'k': 3, 'eta1': 2, 'eta2': 2, 'du': 10, 'dv': 4},
        1024: {'n': 256, 'q': 3329, 'k': 4, 'eta1': 2, 'eta2': 2, 'du': 11, 'dv': 5},
    }

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
        return _mlkem_keygen(d, z, self.k, self.eta1, self.eta2, self.du, self.dv)

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
        return _mlkem_encaps(
            public_key, m, self.k, self.eta1, self.eta2, self.du, self.dv
        )

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
        return _mlkem_decaps(
            secret_key, ciphertext, self.k, self.eta1, self.eta2, self.du, self.dv
        )
