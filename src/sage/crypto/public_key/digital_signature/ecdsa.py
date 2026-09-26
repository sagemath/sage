r"""
Elliptic Curve Digital Signature Algorithm

Toy implementation of the Elliptic Curve Digital Signature Algorithm (ECDSA)
over an elliptic curve defined over a finite field of prime cardinality.

.. WARNING::

    This is a toy implementation for educational use only! Do not use
    this implementation, or any cryptographic features of Sage, in any
    setting where security is needed!

AUTHORS:

- Brian Heckel (2026-09-24): initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Brian Heckel <heckelbri@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Self

from sage.misc.prandom import randint
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from .digital_signature_base import DigitalSignatureBase

if TYPE_CHECKING:
    from sage.rings.finite_rings.integer_mod import IntegerMod_abstract
    from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field


class ECDSA(DigitalSignatureBase):
    r"""
    Elliptic Curve Digital Signature Algorithm.

    Create an instance of ECDSA given elliptic curve ``E`` and base point ``G``.

    INPUT:

    - ``E`` -- elliptic curve over a finite field of prime cardinality
    - ``G`` -- base point on the elliptic curve of prime order
    - ``order`` -- (default: ``None``) the order `n` of ``G``. If ``None``, the
      order is computed, which is expensive for cryptographic size curves. When
      given, it is checked against ``G`` rather than trusted blindly.

    REFERENCES:

    For more information, see Chapter 11 of [MvOV1996]_.

    EXAMPLES:

    A toy ECDSA instance over `\GF{17}`, using a base point of prime
    order 19::

        sage: E = EllipticCurve(GF(17), [2, 2])
        sage: G = E(5, 1)
        sage: G.order()
        19
        sage: ecdsa = digital_signature.ECDSA(E, G)
        doctest:...: FutureWarning: SageMath's digital signature functionality is experimental and might change in the future.
                     See https://github.com/sagemath/sage/issues/41218 for details.
        sage: ecdsa
        ECDSA with parameter set: (Elliptic Curve defined by y^2 = x^3 + 2*x + 2 over Finite Field of size 17, (5 : 1 : 1), 19)

    Signing a message and verifying the resulting signature::

        sage: public_key, secret_key = ecdsa.key_generation()
        sage: signature, message = ecdsa.sign(13, secret_key)
        sage: ecdsa.verify(public_key, signature, message)
        True

    The same toy instance is available under the name ``'toy'``, and
    cryptographic size parameter sets are available by name as well::

        sage: digital_signature.ECDSA.named_parameter_set('toy')
        ecdsa-toy
        sage: digital_signature.ECDSA.named_parameter_set('P-256')
        ecdsa-P-256

    This is a known answer test against the NIST P-256 example of [RFC6979]_,
    Appendix A.2.5, using SHA-256 and the message ``b'sample'``. Passing
    ``nonce`` reproduces the published signature exactly::

        sage: import hashlib
        sage: p256 = digital_signature.ECDSA.named_parameter_set('P-256')
        sage: d = 0xC9AFA9D845BA75166B5C215767B1D6934E50C3DB36E89B127B8A622B120F6721
        sage: Q = p256.public_key(d)
        sage: ZZ(Q[0]) == 0x60FED4BA255A9D31C961EB74C6356D68C049B8923B61FA6CE669622E60F29FB6
        True
        sage: ZZ(Q[1]) == 0x7903FE1008B8BC99A41AE9E95628BC64F2F1B20C2D7E9F5177A3C294D4462299
        True
        sage: z = ZZ(int.from_bytes(hashlib.sha256(b'sample').digest(), 'big'))
        sage: k = 0xA6E3C57DD01ABE90086538398355DD4C3B17AA873382B0F24D6129493D8AAD60
        sage: (r, s), _ = p256.sign(z, d, nonce=k)
        sage: ZZ(r) == 0xEFD48B2AACB6A8FD1140DD9CD45E81D69D2C877B56AAF991C34D0EA84EAF3716
        True
        sage: ZZ(s) == 0xF7CB1C942D657C41D436C7A1B6E29F65F3E900DBB9AFF4064DC4AB2F843ACDA8
        True
        sage: p256.verify(Q, (r, s), z)
        True

    Verification fails for a different message::

        sage: z2 = ZZ(int.from_bytes(hashlib.sha256(b'test').digest(), 'big'))
        sage: p256.verify(Q, (r, s), z2)
        False

    TESTS::

        sage: TestSuite(ecdsa).run()
        sage: TestSuite(digital_signature.ECDSA.named_parameter_set('toy')).run()
        sage: ecdsa == digital_signature.ECDSA.named_parameter_set('toy')
        True

    The base point must have prime order. Otherwise `\Zmod{n}` is not a
    field, and signing fails with an obscure ``ZeroDivisionError`` instead::

        sage: E23 = EllipticCurve(GF(23), [1, 1])
        sage: E23(0, 1).order()
        28
        sage: digital_signature.ECDSA(E23, E23(0, 1))
        Traceback (most recent call last):
        ...
        ValueError: (0 : 1 : 1) must have prime order
        sage: digital_signature.ECDSA(E23, E23(0))
        Traceback (most recent call last):
        ...
        ValueError: (0 : 1 : 0) must have prime order

    The base field must be a finite field of prime cardinality, since ECDSA
    reduces the `x`-coordinate of a curve point modulo `n`::

        sage: K = GF(11^2, 'a')
        sage: Eext = EllipticCurve(K, [1, 1])
        sage: digital_signature.ECDSA(Eext, Eext(0, 1))
        Traceback (most recent call last):
        ...
        ValueError: E must be defined over a finite field of prime cardinality

    A supplied ``order`` is checked against ``G`` rather than trusted::

        sage: digital_signature.ECDSA(E, G, order=23)
        Traceback (most recent call last):
        ...
        ValueError: (5 : 1 : 1) does not have order 23
        sage: digital_signature.ECDSA(E, G, order=18)
        Traceback (most recent call last):
        ...
        ValueError: (5 : 1 : 1) must have prime order
        sage: digital_signature.ECDSA(E, G, order=19) == ecdsa
        True

    An unknown parameter set name is rejected::

        sage: digital_signature.ECDSA.named_parameter_set('P-999')
        Traceback (most recent call last):
        ...
        ValueError: Unknown parameter set name "P-999" for <class 'sage.crypto.public_key.digital_signature.ecdsa.ECDSA'>
    """

    def __init__(
        self,
        E: EllipticCurve_finite_field,
        G: EllipticCurvePoint_finite_field,
        order: Integer | int | None = None,
    ) -> None:
        K = E.base_field()
        if not K.is_finite() or not K.is_prime_field():
            raise ValueError('E must be defined over a finite field of prime cardinality')

        self._E = E

        if G not in self._E:
            raise ValueError(f'{G} is not on {self._E}')

        self._G = G

        if G.is_zero():
            raise ValueError(f'{G} must have prime order')

        self._n = Integer(order) if order is not None else G.order()

        if not self._n.is_prime():
            raise ValueError(f'{G} must have prime order')

        if order is not None and not (self._n * G).is_zero():
            raise ValueError(f'{G} does not have order {self._n}')

        self._R = IntegerModRing(self._n)

    def secret_key(self) -> IntegerMod_abstract:
        r"""
        Generate a random ECDSA secret key.

        OUTPUT:

        A uniformly random element of `\Zmod{n}` in the range `1` to `n - 1`

        TESTS::

            sage: ecdsa = digital_signature.ECDSA.named_parameter_set('toy')
            sage: keys = [ecdsa.secret_key() for _ in range(20)]
            sage: all(1 <= ZZ(d) <= 18 for d in keys)
            True
        """
        return self._R(randint(1, self._n - 1))

    def public_key(self, secret_key) -> EllipticCurvePoint_finite_field:
        """
        Return the ECDSA public key for ``secret_key``.

        INPUT:

        - ``secret_key`` -- the secret key, an integer mod `n`

        OUTPUT:

        The curve point ``secret_key * G``

        EXAMPLES::

            sage: ecdsa = digital_signature.ECDSA.named_parameter_set('toy')
            sage: ecdsa.public_key(7)
            (0 : 6 : 1)
        """
        return Integer(secret_key) * self._G

    def sign(self, message, secret_key, nonce=None) -> tuple[tuple[Any, Any], Any]:
        """
        Return the ECDSA (signature, message) pair signed by ``secret_key``.

        INPUT:

        - ``message`` -- the message to be signed; an integer mod `n`
        - ``secret_key`` -- the secret key of the ECDSA instance
        - ``nonce`` -- (default: ``None``) the per-signature nonce `k`. If
          ``None``, a fresh random nonce is used. An explicit nonce must never be
          reused across two messages, as that reveals the secret key.

        OUTPUT:

        A pair ``((r, s), message)``

        EXAMPLES:

        An explicit nonce makes signing deterministic::

            sage: ecdsa = digital_signature.ECDSA.named_parameter_set('toy')
            sage: ecdsa.sign(13, 7, nonce=11)
            ((13, 6), 13)

        TESTS:

        A nonce that would produce a degenerate signature is rejected rather
        than silently replaced::

            sage: ecdsa.sign(13, 7, nonce=0)
            Traceback (most recent call last):
            ...
            ValueError: nonce must not be zero
        """
        z = self._R(message)
        d = self._R(secret_key)

        while True:
            if nonce is None:
                k = self._R(randint(1, self._n - 1))
            else:
                k = self._R(nonce)
                if k.is_zero():
                    raise ValueError('nonce must not be zero')

            P = Integer(k) * self._G
            r = self._R(Integer(P.x()))
            s = (z + r * d) / k

            if not r.is_zero() and not s.is_zero():
                return ((r, s), message)

            if nonce is not None:
                raise ValueError('the given nonce produces a degenerate signature')

    def verify(self, public_key, signature, message) -> bool:
        """
        Return whether ``signature`` is a valid ECDSA signature on ``message``.

        INPUT:

        - ``public_key`` -- the public key of the ECDSA signature
        - ``signature`` -- a candidate ECDSA signature, a pair ``(r, s)``
        - ``message`` -- the message that was signed

        OUTPUT:

        ``True`` if the signature is valid and ``False`` otherwise

        TESTS:

        A malformed signature is rejected rather than raising, and the result
        is a ``bool``::

            sage: ecdsa = digital_signature.ECDSA.named_parameter_set('toy')
            sage: public_key, secret_key = ecdsa.key_generation()
            sage: ecdsa.verify(public_key, (0, 0), 13)
            False
            sage: ecdsa.verify(public_key, (1, 19), 13)
            False
            sage: ecdsa.verify(public_key, (-1, 1), 13)
            False
            sage: signature, message = ecdsa.sign(13, secret_key)
            sage: ecdsa.verify(public_key, signature, message)
            True
        """
        r, s = signature

        # A verifier must never raise on a malformed signature, and components
        # outside the range 1 to n-1 are invalid. See FIPS 186-4, 6.4.2.
        r = Integer(r)
        s = Integer(s)
        if not (0 < r < self._n and 0 < s < self._n):
            return False

        r = self._R(r)
        s = self._R(s)
        z = self._R(message)

        u_1 = z / s
        u_2 = r / s
        Q = Integer(u_1) * self._G + Integer(u_2) * public_key

        if Q.is_zero():
            return False

        return r == self._R(Integer(Q.x()))

    def parameters(self) -> tuple[EllipticCurve_finite_field, EllipticCurvePoint_finite_field, Integer]:
        """
        Return the elliptic curve ``E``, the base point ``G``, and the order ``n`` of ``G``.

        EXAMPLES::

            sage: digital_signature.ECDSA.named_parameter_set('toy').parameters()
            (Elliptic Curve defined by y^2 = x^3 + 2*x + 2 over Finite Field of size 17,
             (5 : 1 : 1),
             19)
        """
        return (self._E, self._G, self._n)

    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Return an ECDSA instance corresponding to a named parameter set.

        INPUT:

        - ``name`` -- one of the following:

            - ``"toy"``: a toy curve over `\GF{17}` with a base point of
              order 19, small enough to compute with by hand.
            - ``"P-256"``: the NIST curve P-256, also known as secp256r1
              and prime256v1.
            - ``"secp256k1"``: the curve used by Bitcoin.

        EXAMPLES::

            sage: p256 = digital_signature.ECDSA.named_parameter_set('P-256')
            sage: E, G, n = p256.parameters()
            sage: E.base_field().cardinality() == 2^256 - 2^224 + 2^192 + 2^96 - 1
            True
            sage: n == 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
            True
            sage: ZZ(G[0]) == 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
            True

        TESTS::

            sage: k1 = digital_signature.ECDSA.named_parameter_set('secp256k1')
            sage: E, G, n = k1.parameters()
            sage: E.base_field().cardinality() == 2^256 - 2^32 - 977
            True
            sage: E.a4(), E.a6()
            (0, 7)
            sage: n == 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
            True
            sage: public_key, secret_key = k1.key_generation()
            sage: signature, message = k1.sign(12345, secret_key)
            sage: k1.verify(public_key, signature, message)
            True
        """
        if name == 'toy':
            p = Integer(17)
            a, b = Integer(2), Integer(2)
            Gx, Gy = Integer(5), Integer(1)
            n = Integer(19)
        elif name == 'P-256':
            p = Integer(2)**256 - Integer(2)**224 + Integer(2)**192 + Integer(2)**96 - 1
            a = Integer(-3)
            b = Integer(0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B)
            Gx = Integer(0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296)
            Gy = Integer(0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5)
            n = Integer(0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551)
        elif name == 'secp256k1':
            p = Integer(2)**256 - Integer(2)**32 - 977
            a, b = Integer(0), Integer(7)
            Gx = Integer(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798)
            Gy = Integer(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
            n = Integer(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141)
        else:
            return super().named_parameter_set(name)

        # These primes are small enough that proving primality is instant, so we
        # use the default GF(p). That also keeps the cached field object, so an
        # instance from here compares equal to an equivalent hand-built one.
        E = EllipticCurve(GF(p), [a, b])
        ecdsa = cls(E, E(Gx, Gy), order=n)
        ecdsa.rename(f'ecdsa-{name}')
        return ecdsa
