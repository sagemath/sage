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
from typing import Any

from sage.misc.superseded import experimental_warning
from sage.structure.sage_object import SageObject

experimental_warning(
    41218,
    "SageMath's digital signature functionality is experimental and might change in the future.",
)


class DigitalSignatureBase(SageObject):
    r"""
    A base class for Digial Signature Schemes

    Implementers of this class must implement all abstract methods
    defined in :meth:`DigitalSignatureBase`.

    NOTE:

    Typically Digital Signatures Sign arbitrary bytes as messages, however for teaching purposes
    most schemes will use a cryptographic hash on the message to generate an integer.
    Thus to make it easier for users we suggest that messages are represented as some 32-byte
    integer. This can be changed by the implementor however if so then the implementing class
    must override the ``_test_sign()`` method to be able to test with different messages.
    """

    @abstractmethod
    def generate_keys(self) -> tuple[Any, Any]:
        """
        Generates a keypair to be used for signatures

        OUTPUT:

        Returns a 2-tuple of (public_key, secret_key)
        """
        raise NotImplementedError

    @abstractmethod
    def sign(self, message, secret_key) -> tuple[Any, Any]:
        """
        Signs a message from the secret_key

        INPUT:

        - ``message`` -- The message that will be signed, these are assumed to be integers
            for the test suite
        - ``secret_key`` -- The secret_key that only the signer has.

        OUTPUT:

        Returns a pair of (signature, message)
        """
        raise NotImplementedError

    @abstractmethod
    def verify(self, public_key, signature, message) -> bool:
        """
        Verifies that the signature is valid

        INPUT:

        - ``public_key`` -- A public key that is used to verify that the signature is valid
        - ``signature`` -- The signature that the verifier is checking is valid
        - ``message`` -- The message being sent over, assumed to be an integer for the test suite

        OUTPUT:

        Returns True if it is a valid signature for the public_key and message and returns
        False if not.
        """
        raise NotImplementedError

    @abstractmethod
    def parameters(self):
        """
        Returns a tuple of the public parameter set for the Digital Signature

        :meth:`parameters` should a tuple of useful attributes of the instance
        which are sufficient to define the parameter set of the digital signature
        scheme that the instance represents. For some implementations this
        may simply be the parameters passed to ``__init__`` when the object was
        constructed. For some implementations we may wish to return additional
        information for convenience. For example, a digital signature scheme that
        works over an elliptic curve over a finite field may wish to return the
        characteristic of the finite field in addition to the elliptic curve, even
        though the finite field can be accessed via methods on elliptic curve objects.

        The default implementations of ``_eq_`` and ``__hash__`` for
        :class:`DigitalSignatureBase` are implementing using :meth:`parameters`.
        Hence two key exchange instances ``a`` and ``b`` compare as equal
        if and only if ``a.parameters()`` == ``b.parameters()``. Similarly,
        a key exchange instance ``a`` is hashable if and only if ``a.parameters()``
        is hashable. This is a reasonable default that should work for most key
        exchange schemes, but user classes can override the ``_eq_`` and ``__hash__``
        methods if this is not desirable.

        OUTPUT:

        A tuple of public parameters used for the digital signature
        """
        raise NotImplementedError

    def do_signature(self, message) -> tuple[Any, Any, Any, Any, bool]:
        """
        Runs the digital signature protocol for one message, and outputs
        all values computed.

        A 5-tuple ``(public_key, secret_key, signature, message, result_of_verification)``
        """
        public_key, secret_key = self.generate_keys()
        signature, message = self.sign(message, secret_key)
        result_of_verification = self.verify(public_key, signature, message)
        return public_key, secret_key, signature, message, result_of_verification

    
    @classmethod
    def named_parameter_set(cls, name: str) -> Self:
        r"""
        Convenience method to easily construct particular instances of a key exchange scheme
        for actual parameter sets that are used in practice and have names. Implementations
        may also wish to implement a parameter set named 'toy' of a size that is just large
        enough to be non-trivial but is nowhere near cryptographic size. Implementations may
        also wish to set a custom name on the key exchange instance before returning it.

        Sage library implementations of key exchange schemes should define a 'toy' implementation
        and use it for most tests to reduce testing time.
        """
        raise ValueError(f'Unknown parameter set name "{name}" for {cls}')

    def _repr_(self) -> str:
        return f'{type(self).__name__} with parameter set: {self.parameters()}'

    def __eq__(self, other) -> bool:
        return isinstance(other, type(self)) and self.parameters() == other.parameters()

    def __hash__(self) -> int:
        return hash(self.parameters())

    def _test_signature(self, **options):
        """
        Tests that the signature scheme verifies a correct signature for a random
        integer message.
        """
        tester = self._tester(**options)
        message = randint(2, 2000)
        public_key, secret_key, signature, message, result = self.do_signature(message)
        tester.assertTrue(result)



