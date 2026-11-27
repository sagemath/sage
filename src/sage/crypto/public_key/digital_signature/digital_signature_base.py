r"""
Base Class for Digital Signature Schemes

This module contains base classes for digital signature schemes. The class defined
in this module should not be initialized directly. It is the responsibility of
child classes to implement specific signature schemes.

A digital signature scheme produces a signature for a message using a public key
and a secret key.

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


from abc import abstractmethod
from typing import Any, Self

from sage.misc.prandom import randint
from sage.misc.superseded import experimental_warning
from sage.structure.sage_object import SageObject

experimental_warning(
    41218,
    "SageMath's digital signature functionality is experimental and might change in the future.",
)


class DigitalSignatureBase(SageObject):
    r"""
    A base class for digital signature schemes.

    Implementers of this class must implement all abstract methods
    defined in :class:`DigitalSignatureBase`.

    .. NOTE::

        Digital signatures typically sign arbitrary bytes as messages. However,
        for teaching purposes most schemes will apply a cryptographic hash to the
        message to obtain an integer. To make these classes easier to use, we
        therefore suggest that messages are represented as some 32-byte integer.
        An implementer may choose otherwise, but must then override the
        `_test_signature` method so that the test suite uses a
        representation the scheme accepts.
    """

    @abstractmethod
    def secret_key(self):
        """
        Generate a valid secret key for the signer.

        OUTPUT:

        A secret key, which must be kept secret from all other parties
        """
        raise NotImplementedError

    @abstractmethod
    def public_key(self, secret_key):
        """
        Return the public key corresponding to ``secret_key``.

        INPUT:

        - ``secret_key`` -- the signer's secret key

        OUTPUT:

        The public key that verifiers use to check signatures made with
        ``secret_key``
        """
        raise NotImplementedError

    @abstractmethod
    def sign(self, message, secret_key, nonce=None) -> tuple[Any, Any]:
        """
        Sign a message with the secret key.

        INPUT:

        - ``message`` -- the message that will be signed; assumed to be an
          integer for the test suite
        - ``secret_key`` -- the secret key that only the signer has
        - ``nonce`` -- (default: ``None``) the per-signature nonce to use. If
          ``None``, a nonce is chosen uniformly at random, which is what callers
          should normally do. Passing an explicit nonce makes signing
          deterministic, which is needed to reproduce published test vectors.

          .. WARNING::

              Reusing a nonce across two different messages reveals the secret
              key. Only pass an explicit nonce when reproducing a known answer.

        OUTPUT:

        Returns a pair of (signature, message)
        """
        raise NotImplementedError

    @abstractmethod
    def verify(self, public_key, signature, message) -> bool:
        """
        Verify that a signature is valid.

        This must return ``False`` for a malformed signature rather than raising
        an exception, since a verifier acts on data supplied by another party.

        INPUT:

        - ``public_key`` -- the public key that the signature is checked against
        - ``signature`` -- the candidate signature that the verifier is checking
        - ``message`` -- the message that was signed; assumed to be an integer
          for the test suite

        OUTPUT:

        ``True`` if ``signature`` is a valid signature on ``message`` under
        ``public_key``, and ``False`` otherwise
        """
        raise NotImplementedError

    @abstractmethod
    def parameters(self) -> tuple:
        """
        A tuple of the public parameters of the digital signature scheme.

        :meth:`parameters` should return a tuple of useful attributes of the
        instance which are sufficient to define the parameter set of the digital
        signature scheme that the instance represents. For some implementations this
        may simply be the parameters passed to ``__init__`` when the object was
        constructed. For some implementations we may wish to return additional
        information for convenience. For example, a digital signature scheme that
        works over an elliptic curve over a finite field may wish to return the
        characteristic of the finite field in addition to the elliptic curve, even
        though the finite field can be accessed via methods on elliptic curve objects.

        The default implementations of ``__eq__`` and ``__hash__`` for
        :class:`DigitalSignatureBase` are implemented using :meth:`parameters`.
        Hence two digital signature instances ``a`` and ``b`` compare as equal
        if and only if ``a.parameters() == b.parameters()``. Similarly,
        a digital signature instance ``a`` is hashable if and only if
        ``a.parameters()`` is hashable. This is a reasonable default that should
        work for most digital signature schemes, but user classes can override the
        ``__eq__`` and ``__hash__`` methods if this is not desirable.

        OUTPUT:

        A tuple of public parameters used for the digital signature
        """
        raise NotImplementedError

    def key_generation(self) -> tuple[Any, Any]:
        """
        Generate a keypair to be used for signatures.

        OUTPUT:

        A 2-tuple ``(public_key, secret_key)``
        """
        secret_key = self.secret_key()
        return (self.public_key(secret_key), secret_key)

    def do_signature(self, message) -> tuple[Any, Any, Any, Any, bool]:
        """
        Run the digital signature protocol for one message, and output
        all values computed.

        INPUT:

        - ``message`` -- the message to sign

        OUTPUT:

        A 5-tuple ``(public_key, secret_key, signature, message, result_of_verification)``
        """
        public_key, secret_key = self.key_generation()
        signature, message = self.sign(message, secret_key)
        result_of_verification = self.verify(public_key, signature, message)
        return public_key, secret_key, signature, message, result_of_verification

    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Convenience method to easily construct particular instances of a digital
        signature scheme for actual parameter sets that are used in practice and
        have names. Implementations may also wish to implement a parameter set
        named 'toy' of a size that is just large enough to be non-trivial but is
        nowhere near cryptographic size. Implementations may also wish to set a
        custom name on the digital signature instance before returning it.

        Sage library implementations of digital signature schemes should define a
        'toy' implementation and use it for most tests to reduce testing time.

        INPUT:

        - ``name`` -- the name of the parameter set to construct
        """
        raise ValueError(f'Unknown parameter set name "{name}" for {cls}')

    def _repr_(self) -> str:
        return f'{type(self).__name__} with parameter set: {self.parameters()}'

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.parameters() == other.parameters()

    def __hash__(self) -> int:
        return hash(self.parameters())

    def _test_signature(self, **options) -> None:
        """
        Test that the signature scheme verifies a correct signature for a random
        integer message.
        """
        tester = self._tester(**options)
        message = randint(2, 2000)
        *_, result = self.do_signature(message)
        tester.assertTrue(result)
