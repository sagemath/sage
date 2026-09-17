"""
ML-KEM Tests

This module tests the ML-KEM implementation using pytest.

Tests include:
1. Consistency tests: encapsulation and decapsulation produce matching shared secrets
2. Known-answer tests (KATs) from the NIST submission package, using the
   deterministic seeds ``d``, ``z``, and message ``m``.
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


# Known-answer test (KAT) for ML-KEM-512, first record of the NIST test
# vector file.  Source:
# https://gist.github.com/itzmeanjan/c8f5bc9640d0f0bdd2437dfe364d7710
ML_KEM_512_KAT = {
    "d": bytes.fromhex(
        "7c9935a0b07694aa0c6d10e4db6b1add2fd81a25ccb148032dcd739936737f2d"
    ),
    "z": bytes.fromhex(
        "b505d7cfad1b497499323c8686325e4792f267aafa3f87ca60d01cb54f29202a"
    ),
    "m": bytes.fromhex(
        "eb4a7c66ef4eba2ddb38c88d8bc706b1d639002198172a7b1942eca8f6c001ba"
    ),
    "ss": bytes.fromhex(
        "b4c8e3c4115f9511f2fddb288c4b78c5cd7c89d2d4d321f46b4edc54ddf0eb36"
    ),
}


def test_kat_ml_kem_512():
    """Known-answer test for ML-KEM-512 (first NIST test-vector record)."""
    kem = MLKEM.from_parameter_set(512)
    kat = ML_KEM_512_KAT

    pk, sk = kem.keygen(d=kat["d"], z=kat["z"])
    ct, ss = kem.encaps(pk, m=kat["m"])

    assert ss == kat["ss"], "encapsulated shared secret does not match KAT"
    assert kem.decaps(sk, ct) == kat["ss"], "decapsulated shared secret does not match KAT"
