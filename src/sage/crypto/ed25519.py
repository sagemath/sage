# sage.doctest: needs sage.rings.finite_rings
r"""
Ed25519 signatures

This module provides the protocol layer for Ed25519 as specified by RFC 8032.
The underlying curve group is the native twisted Edwards group provided by
:mod:`sage.schemes.elliptic_curves.ell_edwards`.  Keeping the two layers
separate is important: a generic twisted Edwards curve has an equation and a
group law, while Ed25519 additionally fixes the field, parameters, base point,
cofactor, hash function, scalar pruning, and byte encoding.

The public functions implement deterministic Ed25519 signing and verification
using the 32-byte private-key seed format from RFC 8032::

    sage: seed = bytes.fromhex("9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60")
    sage: public_key = ed25519_public_key(seed)
    sage: public_key.hex()
    'd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a'
    sage: signature = ed25519_sign(seed, b"")
    sage: signature.hex()
    'e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b'
    sage: ed25519_verify(public_key, b"", signature)
    True

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

.. [RFC8032] S. Josefsson and I. Liusvaara, *Edwards-Curve Digital Signature
   Algorithm (EdDSA)*, RFC 8032, January 2017.

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
        sage: ed25519_decode(b'\\xff' * 32)
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


def ed25519_public_key(seed):
    r"""Derive the 32-byte Ed25519 public key from a 32-byte seed.

    EXAMPLES::

        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: ed25519_public_key(seed).hex()
        'd75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a'
    """
    scalar, _ = _secret_scalar(seed)
    return ed25519_encode(scalar * Ed25519BasePoint())


def ed25519_sign(seed, message):
    r"""Create a deterministic Ed25519 signature.

    ``seed`` is the 32-byte private-key seed and ``message`` is a bytes-like
    message.  The result is the 64-byte concatenation ``ENC(R) || ENC(S)``.

    EXAMPLES::

        sage: seed = bytes.fromhex('9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60')
        sage: ed25519_sign(seed, b'').hex()
        'e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555fb8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b'
    """
    message = _as_bytes(message, "message")
    scalar, prefix = _secret_scalar(seed)
    public_key = ed25519_encode(scalar * Ed25519BasePoint())

    nonce = int.from_bytes(hashlib.sha512(prefix + message).digest(), "little")
    nonce %= ED25519_SUBGROUP_ORDER
    R = nonce * Ed25519BasePoint()
    encoded_R = ed25519_encode(R)

    challenge = int.from_bytes(
        hashlib.sha512(encoded_R + public_key + message).digest(), "little"
    )
    challenge %= ED25519_SUBGROUP_ORDER
    response = (nonce + challenge * scalar) % ED25519_SUBGROUP_ORDER
    return encoded_R + response.to_bytes(32, "little")


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
        sage: ed25519_verify(public_key, b'', signature[:-1] + bytes([signature[-1] ^ 1]))
        False
    """
    try:
        public_key = _as_bytes(public_key, "public_key")
        message = _as_bytes(message, "message")
        signature = _as_bytes(signature, "signature")
        _require_length(public_key, 32, "public key")
        _require_length(signature, 64, "signature")

        A = ed25519_decode(public_key)
        encoded_R = signature[:32]
        R = ed25519_decode(encoded_R)
        response = int.from_bytes(signature[32:], "little")
        if response >= ED25519_SUBGROUP_ORDER:
            return False

        challenge = int.from_bytes(
            hashlib.sha512(encoded_R + public_key + message).digest(), "little"
        )
        challenge %= ED25519_SUBGROUP_ORDER

        B = Ed25519BasePoint()
        return ED25519_COFACTOR * (response * B) == ED25519_COFACTOR * (
            R + challenge * A
        )
    except (TypeError, ValueError, ZeroDivisionError):
        return False


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
    "ed25519_decode",
    "ed25519_encode",
    "ed25519_public_key",
    "ed25519_sign",
    "ed25519_verify",
]
