# sage.doctest: needs sage.rings.finite_rings
r"""
Ed25519 signatures

This module provides the protocol layer for Ed25519 as specified by RFC 8032.
The underlying curve group is the native twisted Edwards group provided by
:mod:`sage.schemes.elliptic_curves.ell_edwards`.  Keeping the two layers
separate is important: a generic twisted Edwards curve has an equation and a
group law, while Ed25519 additionally fixes the field, parameters, base point,
cofactor, hash function, scalar pruning, and byte encoding.

The public functions implement deterministic Ed25519, Ed25519ctx, and
Ed25519ph signing and verification using the 32-byte private-key seed format
from RFC 8032.  Ed25519ctx separates signatures by a caller-supplied context;
Ed25519ph hashes the message before signing and can also use a context.  The
three variants have different signature domains and are not interchangeable::

    sage: seed = bytes.fromhex("9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60")
    sage: public_key = ed25519_public_key(seed)
    sage: public_key.hex()
    'd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a'
    sage: signature = ed25519_sign(seed, b"")
    sage: signature.hex()
    'e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b'
    sage: ed25519_verify(public_key, b"", signature)
    True

For an object-oriented interface that bundles the variant and context into a
single scheme object with the ``secret_key`` / ``public_key`` / ``sign`` /
``verify`` operation shape of the digital-signature framework proposed in
:issue:`42856`, use :class:`Ed25519SignatureScheme`.

The curve constructor and base point are cached, so ordinary Sage additive
group operations can be used directly without a second hand-written scalar
multiplication implementation::

    sage: C = Ed25519()
    sage: B = Ed25519BasePoint()
    sage: C(B) == B
    True
    sage: ed25519_decode(ed25519_encode(B)) == B
    True

This module intentionally does not provide X25519 or an interchangeable
private-key object.  X25519 is a different protocol using the Montgomery form
and Ed25519's signing key is a seed, not the already-pruned scalar.

REFERENCES:

- [RFC8032]_

- [RFC8032]_

AUTHORS:

- SageMath developers (2026): initial Ed25519 protocol support
"""

# ****************************************************************************
#       Copyright (C) 2026 SageMath developers
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import hashlib

from sage.misc.cachefunc import cached_function
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.schemes.elliptic_curves.ell_edwards import (
    TwistedEdwardsCurve,
    TwistedEdwardsPoint,
)


# RFC 8032 section 5.1 fixes all of these values for Ed25519.
ED25519_FIELD_SIZE = 2**255 - 19
ED25519_SUBGROUP_ORDER = 2**252 + 27742317777372353535851937790883648493
ED25519_COFACTOR = 8
ED25519_A = -1
ED25519_D = (
    -121665 * pow(121666, ED25519_FIELD_SIZE - 2, ED25519_FIELD_SIZE)
) % ED25519_FIELD_SIZE
ED25519_BASE_X = (
    15112221349535400772501151409588531511454012693041857206046113283949847762202
)
ED25519_BASE_Y = (
    46316835694926478169428394003475163141307993866256225615783033603165251855960
)
ED25519_SQRT_M1 = pow(2, (ED25519_FIELD_SIZE - 1) // 4, ED25519_FIELD_SIZE)
_ED25519_Y_MASK = (1 << 255) - 1


@cached_function
def Ed25519():
    r"""Return the standard Ed25519 twisted Edwards curve.

    The returned parent is the curve

    .. MATH::

        -x^2 + y^2 = 1 + d x^2 y^2,
        \qquad d = -121665/121666

    over ``GF(2^255 - 19)``.  The protocol-specific encoding and signing
    functions are :func:`ed25519_encode`, :func:`ed25519_decode`,
    :func:`ed25519_sign`, and :func:`ed25519_verify`.

    EXAMPLES::

        sage: C = Ed25519()
        sage: C.base_ring().order() == ED25519_FIELD_SIZE
        True
        sage: C.a() == -1
        True
        sage: C.d() == C.base_ring()(-121665) / C.base_ring()(121666)
        True
    """
    F = FiniteField(ED25519_FIELD_SIZE)
    return TwistedEdwardsCurve(F, F(-1), F(-121665) / F(121666))


@cached_function
def Ed25519BasePoint():
    r"""Return the standard Ed25519 base point ``B``.

    EXAMPLES::

        sage: B = Ed25519BasePoint()
        sage: B.parent() is Ed25519()
        True
        sage: B.coordinates() == (ED25519_BASE_X, ED25519_BASE_Y)
        True
    """
    return Ed25519()(ED25519_BASE_X, ED25519_BASE_Y)


def _as_bytes(value, name):
    """Convert a bytes-like value or raise a useful ``TypeError``."""
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise TypeError("%s must be bytes-like" % name)
    return bytes(value)


def _require_length(value, length, name):
    """Check the length of a byte string."""
    if len(value) != length:
        raise ValueError("%s must be exactly %d bytes" % (name, length))


def ed25519_encode(point):
    r"""Encode an Ed25519 point as 32 little-endian bytes.

    The low 255 bits encode the affine ``y`` coordinate and the high bit
    encodes the least significant bit of ``x``.  Coordinates are required to
    be canonical field elements and ``point`` must belong to :func:`Ed25519`.

    EXAMPLES::

        sage: ed25519_encode(Ed25519BasePoint()).hex()
        '5866666666666666666666666666666666666666666666666666666666666666'
        sage: ed25519_decode(ed25519_encode(Ed25519BasePoint())) == Ed25519BasePoint()
        True
    """
    if not isinstance(point, TwistedEdwardsPoint) or point.parent() is not Ed25519():
        raise TypeError("point must belong to the standard Ed25519 curve")
    x, y = point.coordinates()
    x = int(x)
    y = int(y)
    if not (0 <= x < ED25519_FIELD_SIZE and 0 <= y < ED25519_FIELD_SIZE):
        raise ValueError("point coordinates must be canonical field elements")
    encoded = bytearray(y.to_bytes(32, "little"))
    encoded[31] |= (x & 1) << 7
    return bytes(encoded)


def ed25519_decode(data):
    r"""Decode a canonical 32-byte Ed25519 point encoding.

    The square root is recovered in ``GF(2^255 - 19)`` using the ``p = 5 mod
    8`` square-root method.  The returned point is a point on the full curve;
    callers that need a prime-order point must apply the protocol's cofactor
    handling.

    EXAMPLES::

        sage: ed25519_decode(bytes.fromhex('5866666666666666666666666666666666666666666666666666666666666666')) == Ed25519BasePoint()
        True
        sage: ed25519_decode(b'\xff' * 32)
        Traceback (most recent call last):
        ...
        ValueError: the encoded y-coordinate is not canonical
    """
    data = _as_bytes(data, "data")
    _require_length(data, 32, "encoded point")

    value = int.from_bytes(data, "little")
    sign = value >> 255
    y = value & _ED25519_Y_MASK
    if y >= ED25519_FIELD_SIZE:
        raise ValueError("the encoded y-coordinate is not canonical")

    yy = y * y % ED25519_FIELD_SIZE
    numerator = (yy - 1) % ED25519_FIELD_SIZE
    denominator = (ED25519_D * yy + 1) % ED25519_FIELD_SIZE
    if denominator == 0:
        raise ValueError("the encoded point has no affine x-coordinate")
    x_squared = numerator * pow(
        denominator, ED25519_FIELD_SIZE - 2, ED25519_FIELD_SIZE
    ) % ED25519_FIELD_SIZE

    x = pow(x_squared, (ED25519_FIELD_SIZE + 3) // 8, ED25519_FIELD_SIZE)
    if x * x % ED25519_FIELD_SIZE != x_squared:
        x = x * ED25519_SQRT_M1 % ED25519_FIELD_SIZE
    if x * x % ED25519_FIELD_SIZE != x_squared:
        raise ValueError("invalid Ed25519 point encoding")
    if x == 0 and sign:
        raise ValueError("negative zero is not a canonical point encoding")
    if (x & 1) != sign:
        x = ED25519_FIELD_SIZE - x
    return Ed25519()(x, y)


def _secret_scalar(seed):
    """Return the pruned scalar and nonce prefix derived from a seed."""
    seed = _as_bytes(seed, "seed")
    _require_length(seed, 32, "seed")
    digest = hashlib.sha512(seed).digest()
    scalar = bytearray(digest[:32])
    scalar[0] &= 248
    scalar[31] &= 63
    scalar[31] |= 64
    return int.from_bytes(scalar, "little"), digest[32:]


def _dom2(phflag, context):
    """Return the RFC 8032 Ed25519ctx/Ed25519ph domain separator."""
    context = _as_bytes(context, "context")
    if len(context) > 255:
        raise ValueError("context must be at most 255 bytes")
    return b"SigEd25519 no Ed25519 collisions" + bytes((phflag, len(context))) + context


def ed25519_public_key(seed):
    r"""Derive the 32-byte Ed25519 public key from a 32-byte seed.

    EXAMPLES::

        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: ed25519_public_key(seed).hex()
        'd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a'
    """
    scalar, _ = _secret_scalar(seed)
    return ed25519_encode(scalar * Ed25519BasePoint())


def _sign(seed, message, domain=b"", prehash=False):
    """Sign using the selected RFC 8032 Ed25519 message and domain."""
    message = _as_bytes(message, "message")
    if prehash:
        message = hashlib.sha512(message).digest()
    scalar, prefix = _secret_scalar(seed)
    public_key = ed25519_encode(scalar * Ed25519BasePoint())

    nonce = int.from_bytes(hashlib.sha512(domain + prefix + message).digest(), "little")
    nonce %= ED25519_SUBGROUP_ORDER
    R = nonce * Ed25519BasePoint()
    encoded_R = ed25519_encode(R)

    challenge = int.from_bytes(
        hashlib.sha512(domain + encoded_R + public_key + message).digest(), "little"
    )
    challenge %= ED25519_SUBGROUP_ORDER
    response = (nonce + challenge * scalar) % ED25519_SUBGROUP_ORDER
    return encoded_R + response.to_bytes(32, "little")


def ed25519_sign(seed, message):
    r"""Create a deterministic pure Ed25519 signature.

    ``seed`` is the 32-byte private-key seed and ``message`` is a bytes-like
    message.  The result is the 64-byte concatenation ``ENC(R) || ENC(S)``.

    EXAMPLES::

        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: ed25519_sign(seed, b'').hex()
        'e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b'
    """
    return _sign(seed, message)


def ed25519ctx_sign(seed, message, context):
    r"""Create an Ed25519ctx signature with a context of at most 255 bytes.

    An empty context is permitted, but still uses the Ed25519ctx domain and
    does not produce a pure Ed25519 signature.  The context is bytes-like and
    must match exactly when verifying.

    EXAMPLES::

        sage: seed = bytes.fromhex('0305334e381af78f141cb666f6199f57bc3495335a256a95bd2a55bf546663f6')
        sage: message = bytes.fromhex('f726936d19c800494e3fdaff20b276a8')
        sage: ed25519ctx_sign(seed, message, b'foo').hex()  # RFC 8032, Section 7.2
        '55a4cc2f70a54e04288c5f4cd1e45a7bb520b36292911876cada7323198dd87a8b36950b95130022907a7fb7c4e9b2d5f6cca685a587b4b21f4b888e4e7edb0d'
        sage: ed25519ctx_sign(seed, message, b'bar').hex()
        'fc60d5872fc46b3aa69f8b5b4351d5808f92bcc044606db097abab6dbcb1aee3216c48e8b3b66431b5b186d1d28f8ee15a5ca2df6668346291c2043d4eb3e90d'

    TESTS::

        sage: ed25519ctx_sign(seed, bytes.fromhex('508e9e6882b979fea900f62adceaca35'), b'foo').hex()
        '8b70c1cc8310e1de20ac53ce28ae6e7207f33c3295e03bb5c0732a1d20dc64908922a8b052cf99b7c4fe107a5abb5b2c4085ae75890d02df26269d8945f84b0b'
        sage: other_seed = bytes.fromhex('ab9c2853ce297ddab85c993b3ae14bcad39b2c682beabc27d6d4eb20711d6560')
        sage: ed25519ctx_sign(other_seed, message, b'foo').hex()
        '21655b5f1aa965996b3f97b3c849eafba922a0a62992f73b3d1b73106a84ad85e9b86a7b6005ea868337ff2d20a7f5fbd4cd10b0be49a68da2b2e0dc0ad8960f'
        sage: ed25519ctx_sign(seed, message, b'x' * 256)
        Traceback (most recent call last):
        ...
        ValueError: context must be at most 255 bytes
    """
    return _sign(seed, message, _dom2(0, context))


def ed25519ph_sign(seed, message, context=b""):
    r"""Create an Ed25519ph signature of a raw message.

    The message is hashed once with SHA-512 before the Ed25519ph nonce and
    challenge hashes.  Do not prehash it before calling this function.
    ``context`` is bytes-like, defaults to empty, and is limited to 255 bytes.

    EXAMPLES::

        sage: seed = bytes.fromhex('833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42')
        sage: ed25519ph_sign(seed, b'abc').hex()  # RFC 8032, Section 7.3
        '98a70222f0b8121aa9d30f813d683f809e462b469c7ff87639499bb94e6dae4131f85042463c2a355a2003d062adf5aaa10b8c61e636062aaad11c2a26083406'
        sage: import hashlib
        sage: ed25519ph_sign(seed, b'abc') != ed25519_sign(seed, hashlib.sha512(b'abc').digest())
        True
    """
    return _sign(seed, message, _dom2(1, context), prehash=True)


def _verify(public_key, message, signature, domain=b"", prehash=False):
    """Verify using the selected RFC 8032 Ed25519 message and domain."""
    try:
        public_key = _as_bytes(public_key, "public_key")
        message = _as_bytes(message, "message")
        signature = _as_bytes(signature, "signature")
        _require_length(public_key, 32, "public key")
        _require_length(signature, 64, "signature")
        if prehash:
            message = hashlib.sha512(message).digest()

        A = ed25519_decode(public_key)
        encoded_R = signature[:32]
        R = ed25519_decode(encoded_R)
        response = int.from_bytes(signature[32:], "little")
        if response >= ED25519_SUBGROUP_ORDER:
            return False

        challenge = int.from_bytes(
            hashlib.sha512(domain + encoded_R + public_key + message).digest(), "little"
        )
        challenge %= ED25519_SUBGROUP_ORDER

        B = Ed25519BasePoint()
        return ED25519_COFACTOR * (response * B) == ED25519_COFACTOR * (
            R + challenge * A
        )
    except (TypeError, ValueError, ZeroDivisionError):
        return False


def ed25519_verify(public_key, message, signature):
    r"""Verify an Ed25519 signature.

    The RFC 8032 verification equation is checked with the cofactor:

    .. MATH::

        [8][S]B = [8]R + [8][k]A,

    where ``k`` is the reduced SHA-512 challenge.  Malformed encodings and
    non-canonical ``S`` values return ``False``.

    EXAMPLES::

        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: public_key = ed25519_public_key(seed)
        sage: signature = ed25519_sign(seed, b'')
        sage: ed25519_verify(public_key, b'', signature)
        True
        sage: ed25519_verify(public_key, b'bad', signature)
        False
        sage: import operator
        sage: ed25519_verify(public_key, b'', signature[:-1] + bytes([operator.xor(signature[-1], 1)]))
        False
    """
    return _verify(public_key, message, signature)


def ed25519ctx_verify(public_key, message, signature, context):
    r"""Verify an Ed25519ctx signature under the given context.

    Return ``False`` for malformed input or a context longer than 255 bytes.

    EXAMPLES::

        sage: seed = bytes.fromhex('0305334e381af78f141cb666f6199f57bc3495335a256a95bd2a55bf546663f6')
        sage: public_key = ed25519_public_key(seed)
        sage: message = bytes.fromhex('f726936d19c800494e3fdaff20b276a8')
        sage: signature = ed25519ctx_sign(seed, message, b'foo')
        sage: ed25519ctx_verify(public_key, message, signature, b'foo')
        True
        sage: ed25519ctx_verify(public_key, message, signature, b'bar')
        False
        sage: ed25519_verify(public_key, message, signature)
        False
        sage: ed25519ctx_verify(public_key, message, signature, b'x' * 256)
        False
    """
    try:
        domain = _dom2(0, context)
    except (TypeError, ValueError):
        return False
    return _verify(public_key, message, signature, domain)


def ed25519ph_verify(public_key, message, signature, context=b""):
    r"""Verify an Ed25519ph signature of a raw message.

    ``context`` defaults to empty and must contain at most 255 bytes.
    Return ``False`` for malformed input.

    EXAMPLES::

        sage: seed = bytes.fromhex('833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42')
        sage: public_key = ed25519_public_key(seed)
        sage: signature = ed25519ph_sign(seed, b'abc')
        sage: ed25519ph_verify(public_key, b'abc', signature)
        True
        sage: ed25519ph_verify(public_key, b'abd', signature)
        False
        sage: ed25519ph_verify(public_key, b'abc', signature, b'foo')
        False
        sage: ed25519_verify(public_key, b'abc', signature)
        False

    TESTS::

        sage: signature = ed25519ph_sign(seed, b'abc', b'foo')
        sage: ed25519ph_verify(public_key, b'abc', signature, b'foo')
        True
        sage: ed25519ph_verify(public_key, b'abc', signature)
        False
    """
    try:
        domain = _dom2(1, context)
    except (TypeError, ValueError):
        return False
    return _verify(public_key, message, signature, domain, prehash=True)


# The digital-signature framework of gh-42856 (sage.crypto.public_key.
# digital_signature.DigitalSignatureBase) is still under review.  Import it
# opportunistically so that this class already provides the framework's four
# operations (secret_key, public_key, sign, verify) and automatically becomes
# a DigitalSignatureBase subclass once that framework is merged.
try:
    from sage.crypto.public_key.digital_signature.digital_signature_base import (
        DigitalSignatureBase,
    )
except ImportError:
    DigitalSignatureBase = object


class Ed25519SignatureScheme(DigitalSignatureBase):
    r"""Ed25519, Ed25519ctx and Ed25519ph as a digital-signature scheme object.

    This is the object-oriented counterpart of the module-level functions:
    one instance bundles the RFC 8032 variant and (for ctx/ph) the context
    string, which removes the risk of signing and verifying with mismatched
    domain parameters.  When the digital-signature framework of :issue:`42856`
    is available, this class is a
    :class:`~sage.crypto.public_key.digital_signature.digital_signature_base.DigitalSignatureBase`
    subclass and can be used uniformly with the other schemes of that
    framework.

    INPUT:

    - ``variant`` -- one of ``'ed25519'``, ``'ed25519ctx'``, ``'ed25519ph'``

    - ``context`` -- bytes-like application context (at most 255 bytes),
      only used by the ``'ed25519ctx'`` and ``'ed25519ph'`` variants;
      ignored for pure Ed25519

    Keys are the RFC 8032 32-byte seed (secret) and 32-byte encoded point
    (public).  Messages are bytes-like; Ed25519ph hashes them with SHA-512
    internally, so pass the raw message, not its prehash.

    EXAMPLES::

        sage: from sage.crypto.ed25519 import Ed25519SignatureScheme
        sage: ed = Ed25519SignatureScheme('ed25519ph', b'foo')
        sage: secret = ed.secret_key()
        sage: public = ed.public_key(secret)
        sage: signature, signed_message = ed.sign(b'attack at dawn', secret)
        sage: ed.verify(public, signature, signed_message)
        True

    A signature made under one context does not verify under another::

        sage: other = Ed25519SignatureScheme('ed25519ph', b'bar')
        sage: other.verify(public, signature, signed_message)
        False

    The pure variant reproduces the RFC 8032 section 7.1 test vector::

        sage: pure = Ed25519SignatureScheme('ed25519')
        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: signature, _ = pure.sign(b'', seed)
        sage: signature.hex()
        'e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b'

    TESTS::

        sage: Ed25519SignatureScheme('ed448')
        Traceback (most recent call last):
        ...
        ValueError: variant must be one of 'ed25519', 'ed25519ctx', 'ed25519ph'
        sage: Ed25519SignatureScheme('ed25519ctx', b'x' * 256)
        Traceback (most recent call last):
        ...
        ValueError: context must be at most 255 bytes
    """

    def __init__(self, variant='ed25519', context=b''):
        if variant not in ('ed25519', 'ed25519ctx', 'ed25519ph'):
            raise ValueError(
                "variant must be one of 'ed25519', 'ed25519ctx', 'ed25519ph'"
            )
        self._variant = variant
        self._context = _as_bytes(context, "context")
        self._domain = (
            b""
            if variant == 'ed25519'
            else _dom2(variant == 'ed25519ph', self._context)
        )
        self._prehash = variant == 'ed25519ph'

    def _repr_(self):
        if self._variant == 'ed25519':
            return "Ed25519 signature scheme"
        return "Ed25519%s signature scheme with context %r" % (
            self._variant[6:], bytes(self._context)
        )

    def variant(self):
        r"""Return the RFC 8032 variant of this scheme."""
        return self._variant

    def context(self):
        r"""Return the context string (empty for pure Ed25519)."""
        return bytes(self._context)

    def secret_key(self):
        r"""Return a fresh random 32-byte Ed25519 seed.

        EXAMPLES::

            sage: from sage.crypto.ed25519 import Ed25519SignatureScheme
            sage: ed = Ed25519SignatureScheme()
            sage: sk = ed.secret_key()
            sage: len(sk)
            32
        """
        import secrets
        return secrets.token_bytes(32)

    def public_key(self, secret_key):
        r"""Return the 32-byte public key of a 32-byte seed.

        EXAMPLES::

            sage: from sage.crypto.ed25519 import Ed25519SignatureScheme
            sage: ed = Ed25519SignatureScheme()
            sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
            sage: ed.public_key(seed).hex()
            'd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a'
        """
        return ed25519_public_key(secret_key)

    def sign(self, message, secret_key, nonce=None):
        r"""Sign ``message`` with a 32-byte seed and return ``(signature, message)``.

        Ed25519 nonces are derived deterministically from the secret key and
        the message, so an explicit ``nonce`` is not supported; the parameter
        is accepted and ignored to keep the signature-shape of the framework.

        EXAMPLES::

            sage: from sage.crypto.ed25519 import Ed25519SignatureScheme
            sage: ed = Ed25519SignatureScheme('ed25519ctx', b'foo')
            sage: seed = bytes.fromhex('0305334e381af78f141cb666f6199f57bc3495335a256a95bd2a55bf546663f6')
            sage: message = bytes.fromhex('f726936d19c800494e3fdaff20b276a8')
            sage: signature, signed = ed.sign(message, seed)
            sage: signature.hex() == ed25519ctx_sign(seed, message, b'foo').hex()
            True
        """
        return _sign(secret_key, message, self._domain, self._prehash), message

    def verify(self, public_key, signature, message):
        r"""Return whether ``signature`` is a valid signature of ``message``.

        EXAMPLES::

            sage: from sage.crypto.ed25519 import Ed25519SignatureScheme
            sage: ed = Ed25519SignatureScheme('ed25519ph')
            sage: seed = bytes.fromhex('833fe62409237b9d62ec77587520911e9a759cec1d19755b7da901b96dca3d42')
            sage: public = ed.public_key(seed)
            sage: signature, message = ed.sign(b'abc', seed)
            sage: ed.verify(public, signature, message)
            True
            sage: ed.verify(public, signature, b'abd')
            False
        """
        return _verify(
            public_key, message, signature, self._domain, self._prehash
        )


__all__ = [
    "ED25519_A",
    "ED25519_BASE_X",
    "ED25519_BASE_Y",
    "ED25519_COFACTOR",
    "ED25519_D",
    "ED25519_FIELD_SIZE",
    "ED25519_SQRT_M1",
    "ED25519_SUBGROUP_ORDER",
    "Ed25519",
    "Ed25519BasePoint",
    "Ed25519SignatureScheme",
    "ed25519_decode",
    "ed25519_encode",
    "ed25519_public_key",
    "ed25519_sign",
    "ed25519_verify",
    "ed25519ctx_sign",
    "ed25519ctx_verify",
    "ed25519ph_sign",
    "ed25519ph_verify",
]
