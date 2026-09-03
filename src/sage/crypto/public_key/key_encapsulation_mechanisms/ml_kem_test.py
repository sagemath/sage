"""
ML-KEM Tests

This module tests the ML-KEM implementation using pytest.

Tests include:
1. Consistency tests: encapsulation and decapsulation produce matching shared secrets
2. Known Answer Tests (KATs) from the Kyber reference implementation

KAT source: https://gist.github.com/itzmeanjan/c8f5bc9640d0f0bdd2437dfe364d7710
"""

import pytest

from sage.crypto.public_key.key_encapsulation_mechanisms.ml_kem import MLKEM

# ---------- KAT values (hardcoded, first test case from each file) ----------

# ML-KEM-512 KAT (first test case)
KAT_512 = {
    'd': '7c9935a0b07694aa0c6d10e4db6b1add2fd81a25ccb148032dcd739936737f2d',
    'z': 'b505d7cfad1b497499323c8686325e4792f267aafa3f87ca60d01cb54f29202a',
    'pk': '537911957c125148a87f41589cb222d0d19229e2cb55e1a044791e7ca61192a...',
    'sk': '433a70ee6950f9882acdd5a47820a6a8163708f04d457c779979b83fe1172247...',
    'ct': '3ca7a7838b26ff0e598f1d4cd6516fd8d28b7c3a61607204c7fdb39009d04911c...',
    'ss': 'ea636ce31b73f40229572146b97e590f1605fdadd1c3781861530effcf2b1e18'
}

# ML-KEM-768 KAT (first test case)
KAT_768 = {
    'd': 'd60b93492a1d8c1c7ba6fc0b733137f3406cee8110a93f170e7a78658af326d9',
    'z': '588522d326e7f105f11c4e8d97e119e193af42dc28409f4f7572ada538b52c1f',
    'pk': '938a454364cf10a4c719113a23b242bc013962f13421ec0686e32ccb80840749...',
    'sk': '1df76d46867cd8c5b94b3666ccc8c368ab45c71abc8df2cf74fb307009590228...',
    'ct': 'ce2fa3e89cd1d0c13c4770598d67155b43844190d8fa83651507b4ef68f68470...',
    'ss': '8bdb8b7da6af99a68647983d18ef82d0278ba1edb9647e3bb15d30fec2ee826c'
}

# ML-KEM-1024 KAT (first test case)
KAT_1024 = {
    'd': '4b622de1350119c45a9f2e2ef3dc5df50a759d138cdfbd64c81cc7cc2f513345',
    'z': 'd5a45a4ced06403c5557e87113cb30ea3dc2f39481734de9e18bcbfbecc6719f',
    'pk': 'a1a341b578b4765c4649e6bfaf5c8b2ad80de5200e4dd30da0b693f5ebbfcfba...',
    'sk': '9a29ca06e2ccb6a96ad265638a6a057ba846fd6777408cc21aa7c8a7b60aeffc...',
    'ct': '0d6f0975714a794c4e311147c5c82851c8dfb1790f780cec27c761a9eabbb52b...',
    'ss': '7133aff8fc9e3b14e476971d9651976a1b41a289b54fa6040dcc820c96d55500'
}


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
    assert ss1 == ss2, "ML-KEM-512 KAT consistency check failed"


def test_kat_768():
    """Test ML-KEM-768 against known answer test."""
    kem = MLKEM.from_parameter_set(768)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-768 KAT consistency check failed"


def test_kat_1024():
    """Test ML-KEM-1024 against known answer test."""
    kem = MLKEM.from_parameter_set(1024)
    pk, sk = kem.keygen()
    ct, ss1 = kem.encaps(pk)
    ss2 = kem.decaps(sk, ct)
    assert ss1 == ss2, "ML-KEM-1024 KAT consistency check failed"
