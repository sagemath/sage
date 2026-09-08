"""
ML-KEM Tests

This module tests the ML-KEM implementation using pytest.

Tests include:
1. Consistency tests: encapsulation and decapsulation produce matching shared secrets
"""

import pytest

from sage.crypto.public_key.key_encapsulation_mechanisms.ml_kem import MLKEM


@pytest.mark.parametrize("params", [512, 768, 1024])
def test_consistency(params):
    """Test that encapsulation and decapsulation produce matching shared secrets."""
    kem = MLKEM.from_parameter_set(params)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, f"ML-KEM-{params}: shared secrets do not match"

