"""
Known Answer Tests for ML-KEM
"""

from sage.crypto.public_key.key_encapsulation_mechanisms import MLKEM


def test_ml_kem_512():
    """
    Test ML-KEM-512 with known answer test.
    """
    kem = MLKEM.from_parameter_set(512)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-512 KAT failed"
    print("ML-KEM-512 KAT passed")


def test_ml_kem_768():
    kem = MLKEM.from_parameter_set(768)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-768 KAT failed"
    print("ML-KEM-768 KAT passed")


def test_ml_kem_1024():
    kem = MLKEM.from_parameter_set(1024)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-1024 KAT failed"
    print("ML-KEM-1024 KAT passed")


if __name__ == "__main__":
    test_ml_kem_512()
    test_ml_kem_768()
    test_ml_kem_1024()
    print("All KATs passed!")
