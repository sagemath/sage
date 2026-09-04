"""
ML-KEM Tests

This module tests the ML-KEM implementation using pytest.

Tests include:
1. Consistency tests: encapsulation and decapsulation produce matching shared secrets
2. Known Answer Tests (KATs) against the Kyber reference implementation

KAT source: https://gist.github.com/itzmeanjan/c8f5bc9640d0f0bdd2437dfe364d7710
"""

import pytest

from sage.crypto.public_key.key_encapsulation_mechanisms.ml_kem import MLKEM

# ---------- Full KAT Values (first test case from each file) ----------

# ML-KEM-512 KAT
KAT_512 = {'ss': 'ea636ce31b73f40229572146b97e590f1605fdadd1c3781861530effcf2b1e18'}

# ML-KEM-768 KAT
KAT_768 = {'ss': '8bdb8b7da6af99a68647983d18ef82d0278ba1edb9647e3bb15d30fec2ee826c'}

# ML-KEM-1024 KAT
KAT_1024 = {'ss': '7133aff8fc9e3b14e476971d9651976a1b41a289b54fa6040dcc820c96d55500'}


# ---------- Consistency Tests ----------


@pytest.mark.parametrize("params", [512, 768, 1024])
def test_consistency(params):
    """Test that encapsulation and decapsulation produce matching shared secrets."""
    kem = MLKEM.from_parameter_set(params)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, f"ML-KEM-{params}: shared secrets do not match"


# ---------- KAT Tests ----------


def test_kat_512():
    """Test ML-KEM-512 against known answer test."""
    kem = MLKEM.from_parameter_set(512)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)

    expected_ss = bytes.fromhex(KAT_512['ss'])
    assert ss1 == expected_ss, "ML-KEM-512 KAT failed"
    assert ss2 == expected_ss, "ML-KEM-512 decaps KAT failed"


def test_kat_768():
    """Test ML-KEM-768 against known answer test."""
    kem = MLKEM.from_parameter_set(768)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)

    expected_ss = bytes.fromhex(KAT_768['ss'])
    assert ss1 == expected_ss, "ML-KEM-768 KAT failed"
    assert ss2 == expected_ss, "ML-KEM-768 decaps KAT failed"


def test_kat_1024():
    """Test ML-KEM-1024 against known answer test."""
    kem = MLKEM.from_parameter_set(1024)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)

    expected_ss = bytes.fromhex(KAT_1024['ss'])
    assert ss1 == expected_ss, "ML-KEM-1024 KAT failed"
    assert ss2 == expected_ss, "ML-KEM-1024 decaps KAT failed"
