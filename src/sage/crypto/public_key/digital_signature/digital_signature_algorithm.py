r"""
Digital Signature Algorithm

Toy implementation of the Digital Signature Algorithm (DSA) over the
multiplicative group of `\Zmod{p}`, using a subgroup of prime order `q`.

.. WARNING::

    This is a toy implementation for educational use only! Do not use
    this implementation, or any cryptographic features of Sage, in any
    setting where security is needed!

AUTHORS:

- Brian Heckel (2026-06-11): initial version
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
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer

from .digital_signature_base import DigitalSignatureBase

if TYPE_CHECKING:
    from sage.rings.finite_rings.integer_mod import IntegerMod_abstract


class DigitalSignatureAlgorithm(DigitalSignatureBase):
    r"""
    Digital Signature Algorithm.

    Create an instance of the Digital Signature Algorithm using
    primes ``p`` and ``q`` and generator ``g``.

    INPUT:

    - ``q`` -- prime integer

    - ``p`` -- prime integer for which (p-1) is a multiple of ``q``

    - ``generator`` -- element of `\Zmod{p}` of multiplicative order ``q``

    - ``proof`` -- whether to require a proof that ``p`` and ``q`` are prime.
      If ``False``, a probabilistic test is used. This should be set to
      ``False`` when using large (cryptographic size) primes, otherwise
      checking primality will take too long. If this is not specified, then the
      default behaviour is to use the current value of `proof.arithmetic()`.

    REFERENCES:

    For more information, see Chapter 11 of [MvOV1996]_.

    EXAMPLES::

        sage: DSA = digital_signature.DigitalSignatureAlgorithm(4007, 2003, 867)
        doctest:...: FutureWarning: SageMath's digital signature functionality is experimental and might change in the future.
                     See https://github.com/sagemath/sage/issues/41218 for details.
        sage: DSA
        DigitalSignatureAlgorithm with parameter set: (4007, 2003, 867)
        sage: message = 300
        sage: public, secret, sig, message, result = DSA.do_signature(message)
        sage: result
        True

    The same toy parameters are available by name::

        sage: digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
        dsa-toy

    This is a known answer test against the 512-bit example of [FIPS186-2]_,
    Appendix 5, which signs the SHA-1 digest of ``b'abc'``. Passing ``nonce``
    reproduces the published signature exactly::

        sage: import hashlib
        sage: DSA2 = digital_signature.DigitalSignatureAlgorithm.named_parameter_set('fips186-2')
        sage: x = 0x2070b3223dba372fde1c0ffc7b2e3b498b260614
        sage: y = DSA2.public_key(x)
        sage: ZZ(y).hex() == '19131871d75b1612a819f29d78d1b0d7346f7aa77bb62a859bfd6c5675da9d21'\
        ....:               '2d3a36ef1672ef660b8c7c255cc0ec74858fba33f44c06699630a76b030ee333'
        True
        sage: z = ZZ(int.from_bytes(hashlib.sha1(b'abc').digest(), 'big'))
        sage: k = 0x358dad571462710f50e254cf1a376b2bdeaadfbf
        sage: (r, s), _ = DSA2.sign(z, x, nonce=k)
        sage: ZZ(r) == 0x8bac1ab66410435cb7181f95b16ab97c92b341c0
        True
        sage: ZZ(s) == 0x41e2345f1f56df2458f426d155b4ba2db6dcd8c8
        True
        sage: DSA2.verify(y, (r, s), z)
        True

    Verification fails for a different message::

        sage: z2 = ZZ(int.from_bytes(hashlib.sha1(b'abd').digest(), 'big'))
        sage: DSA2.verify(y, (r, s), z2)
        False

    TESTS::

        sage: TestSuite(DSA).run()
        sage: TestSuite(digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')).run()
        sage: DSA == digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
        True

    The generator must have order exactly ``q``. Without this check, signing
    and verification silently disagree for roughly half of all signatures::

        sage: Zmod(4007).multiplicative_generator().multiplicative_order()
        4006
        sage: digital_signature.DigitalSignatureAlgorithm(4007, 2003, 5)
        Traceback (most recent call last):
        ...
        ValueError: generator must have order q
        sage: digital_signature.DigitalSignatureAlgorithm(4007, 2003, 0)
        Traceback (most recent call last):
        ...
        ValueError: generator must have order q
        sage: digital_signature.DigitalSignatureAlgorithm(4007, 2003, 1)
        Traceback (most recent call last):
        ...
        ValueError: generator must not be equal to one
    """

    def __init__(
        self,
        p: Integer | int,
        q: Integer | int,
        generator: Integer | int | IntegerMod_abstract,
        proof: bool | None = None,
    ) -> None:
        self._p = Integer(p)
        self._q = Integer(q)

        if (self._p - 1) % self._q != 0:
            raise ValueError('p-1 must be a multiple of q')

        # proof=False uses a probabilistic primality test, which is necessary
        # for cryptographic size primes where a proof would take far too long.
        if not self._q.is_prime(proof=proof) or not self._p.is_prime(proof=proof):
            raise ValueError('p and q must be prime')

        self._Zp = IntegerModRing(self._p)
        self._Zq = IntegerModRing(self._q)
        self._generator = self._Zp(generator)

        if self._generator == 1:
            raise ValueError('generator must not be equal to one')

        if self._generator**self._q != 1:
            raise ValueError('generator must have order q')

    def secret_key(self) -> IntegerMod_abstract:
        r"""
        Generate a random DSA secret key.

        OUTPUT:

        A uniformly random element of `\Zmod{q}` in the range `1` to `q - 1`

        TESTS::

            sage: DSA = digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
            sage: keys = [DSA.secret_key() for _ in range(20)]
            sage: all(1 <= ZZ(x) <= 2002 for x in keys)
            True
        """
        return self._Zq(randint(1, self._q - 1))

    def public_key(self, secret_key) -> IntegerMod_abstract:
        """
        Return the DSA public key for ``secret_key``.

        INPUT:

        - ``secret_key`` -- the secret key, an integer mod `q`

        OUTPUT:

        The group element ``generator^secret_key`` in `\\Zmod{p}`

        EXAMPLES::

            sage: DSA = digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
            sage: DSA.public_key(1234)
            1340
        """
        return self._generator**Integer(secret_key)

    def sign(self, message, secret_key, nonce=None) -> tuple[tuple[IntegerMod_abstract, IntegerMod_abstract], Any]:
        """
        Return the DSA (signature, message) pair signed by ``secret_key``.

        INPUT:

        - ``message`` -- integer message value to be signed
        - ``secret_key`` -- the secret key that only the signer has
        - ``nonce`` -- (default: ``None``) the per-signature nonce `k`. If
          ``None``, a fresh random nonce is used. An explicit nonce must never be
          reused across two messages, as that reveals the secret key.

        OUTPUT:

        A pair ``((r, s), message)``

        EXAMPLES:

        An explicit nonce makes signing deterministic::

            sage: DSA = digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
            sage: DSA.sign(300, 1234, nonce=567)
            ((1707, 1649), 300)

        TESTS::

            sage: DSA.sign(300, 1234, nonce=0)
            Traceback (most recent call last):
            ...
            ValueError: nonce must not be zero
        """
        z = self._Zq(message)
        x = self._Zq(secret_key)

        while True:
            if nonce is None:
                k = self._Zq(randint(1, self._q - 1))
            else:
                k = self._Zq(nonce)
                if k.is_zero():
                    raise ValueError('nonce must not be zero')

            r = self._Zq(Integer(self._generator**Integer(k)))
            s = k**(-1) * (z + x * r)

            if not r.is_zero() and not s.is_zero():
                return ((r, s), message)

            if nonce is not None:
                raise ValueError('the given nonce produces a degenerate signature')

    def verify(self, public_key: IntegerMod_abstract, signature: tuple[IntegerMod_abstract, IntegerMod_abstract], message) -> bool:
        """
        Return whether ``signature`` is a valid DSA signature on ``message``.

        INPUT:

        - ``public_key`` -- the public key that the signature is checked against
        - ``signature`` -- a candidate DSA signature, a pair ``(r, s)``
        - ``message`` -- the message that was signed; assumed to be an integer
          for the test suite

        OUTPUT:

        ``True`` if the signature is valid and ``False`` otherwise

        TESTS:

        A malformed signature is rejected rather than raising an exception. In
        particular ``s`` congruent to 0 has no inverse modulo ``q``::

            sage: DSA = digital_signature.DigitalSignatureAlgorithm(4007, 2003, 867)
            sage: public_key, secret_key = DSA.key_generation()
            sage: DSA.verify(public_key, (1, 0), 300)
            False
            sage: DSA.verify(public_key, (0, 1), 300)
            False
            sage: DSA.verify(public_key, (1, 5000), 300)
            False
            sage: DSA.verify(public_key, (-1, 1), 300)
            False
        """
        r, s = signature

        r = Integer(r)
        s = Integer(s)
        if not (0 < r < self._q and 0 < s < self._q):
            return False

        r = self._Zq(r)
        s = self._Zq(s)
        w = s**(-1)
        u1 = self._Zq(message) * w
        u2 = r * w
        v = self._Zq(
            Integer(self._generator**Integer(u1) * self._Zp(public_key)**Integer(u2))
        )
        return v == r

    def parameters(self) -> tuple[Integer, Integer, IntegerMod_abstract]:
        """
        Return the public parameter set, which is of the form ``(p, q, generator)``.

        OUTPUT:

        The public parameter set of the signature scheme ``(p, q, g)``

        EXAMPLES::

            sage: digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy').parameters()
            (4007, 2003, 867)
        """
        return (self._p, self._q, self._generator)

    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Return a DSA instance corresponding to a named parameter set.

        INPUT:

        - ``name`` -- one of the following:

            - ``"toy"``: a toy parameter set with `p = 4007` and `q = 2003`,
              nowhere near cryptographic size.
            - ``"fips186-2"``: the 512-bit example parameter set from
              [FIPS186-2]_, Appendix 5.

        EXAMPLES::

            sage: digital_signature.DigitalSignatureAlgorithm.named_parameter_set('toy')
            dsa-toy
            sage: DSA = digital_signature.DigitalSignatureAlgorithm.named_parameter_set('fips186-2')
            sage: p, q, g = DSA.parameters()
            sage: p.nbits(), q.nbits()
            (512, 160)

        TESTS::

            sage: digital_signature.DigitalSignatureAlgorithm.named_parameter_set('nope')
            Traceback (most recent call last):
            ...
            ValueError: Unknown parameter set name "nope" for <class 'sage.crypto.public_key.digital_signature.digital_signature_algorithm.DigitalSignatureAlgorithm'>
        """
        if name == 'toy':
            p = Integer(4007)
            q = Integer(2003)
            g = Integer(867)
        elif name == 'fips186-2':
            p = Integer('8df2a494492276aa3d25759bb06869cbeac0d83afb8d0cf7cbb8324f0d7882e5'
                        'd0762fc5b7210eafc2e9adac32ab7aac49693dfbf83724c2ec0736ee31c80291', 16)
            q = Integer('c773218c737ec8ee993b4f2ded30f48edace915f', 16)
            g = Integer('626d027839ea0a13413163a55b4cb500299d5522956cefcb3bff10f399ce2c2e'
                        '71cb9de5fa24babf58e5b79521925c9cc42e9f6f464b088cc572af53e6d78802', 16)
        else:
            return super().named_parameter_set(name)

        dsa = cls(p, q, g, proof=False)
        dsa.rename(f'dsa-{name}')
        return dsa
