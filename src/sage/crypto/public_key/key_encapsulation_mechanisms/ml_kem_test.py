"""
ML-KEM Tests

This module tests the ML-KEM implementation using pytest.

Tests include consistency tests (encapsulation and decapsulation produce
matching shared secrets) and known-answer tests (KATs) from the NIST
submission package, using the deterministic seeds ``d``, ``z``, and
message ``m``.
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


# Known-answer tests (KATs) from the NIST submission package.  The
# expected values are the first record of the official test-vector files
# for each parameter set; see the KAT files at
# https://gist.github.com/itzmeanjan/c8f5bc9640d0f0bdd2437dfe364d7710
# The first record uses the same deterministic seeds across all three
# parameter sets; only the expected shared secret differs.
_ML_KEM_KAT_D = bytes.fromhex(
    "7c9935a0b07694aa0c6d10e4db6b1add2fd81a25ccb148032dcd739936737f2d"
)
_ML_KEM_KAT_Z = bytes.fromhex(
    "b505d7cfad1b497499323c8686325e4792f267aafa3f87ca60d01cb54f29202a"
)
_ML_KEM_KAT_M = bytes.fromhex(
    "eb4a7c66ef4eba2ddb38c88d8bc706b1d639002198172a7b1942eca8f6c001ba"
)

_ML_KEM_KAT_SS = {
    512: bytes.fromhex(
        "b4c8e3c4115f9511f2fddb288c4b78c5cd7c89d2d4d321f46b4edc54ddf0eb36"
    ),
    768: bytes.fromhex(
        "ac865f839fef1bf3d528dd7504bed2f64b5502b0fa81d1c32763658e4aac5037"
    ),
    1024: bytes.fromhex(
        "ea636ce31b73f40229572146b97e590f1605fdadd1c3781861530effcf2b1e18"
    ),
}


@pytest.mark.parametrize("params", [512, 768, 1024])
def test_kat(params):
    """Known-answer test for the first record of each NIST KAT file."""
    kem = MLKEM.from_parameter_set(params)

    pk, sk = kem.keygen(d=_ML_KEM_KAT_D, z=_ML_KEM_KAT_Z)
    ct, ss = kem.encaps(pk, m=_ML_KEM_KAT_M)

    assert ss == _ML_KEM_KAT_SS[params], (
        f"ML-KEM-{params}: encapsulated shared secret does not match KAT"
    )
    assert kem.decaps(sk, ct) == _ML_KEM_KAT_SS[params], (
        f"ML-KEM-{params}: decapsulated shared secret does not match KAT"
    )
