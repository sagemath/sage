"""
Consistency Tests for ML-KEM

These tests verify that encapsulation and decapsulation produce
matching shared secrets. They are not official NIST KATs.
"""

from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM


def test_ml_kem_512():
    """
    Test ML-KEM-512 consistency.
    """
    kem = MLKEM.from_parameter_set(512)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-512 consistency test failed"
    print("ML-KEM-512 consistency test passed")


def test_ml_kem_768():
    kem = MLKEM.from_parameter_set(768)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-768 consistency test failed"
    print("ML-KEM-768 consistency test passed")


def test_ml_kem_1024():
    kem = MLKEM.from_parameter_set(1024)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-1024 consistency test failed"
    print("ML-KEM-1024 consistency test passed")


if __name__ == "__main__":
    test_ml_kem_512()
    test_ml_kem_768()
    test_ml_kem_1024()
    print("All consistency tests passed!")
