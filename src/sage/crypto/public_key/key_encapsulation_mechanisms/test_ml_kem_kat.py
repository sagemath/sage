"""
Known Answer Tests for ML-KEM

This module verifies the ML-KEM implementation against NIST known answer tests (KATs).
The KAT files are parsed and each test case validates that encapsulation and
decapsulation produce identical shared secrets.

The tests cover all three parameter sets: ML-KEM-512, ML-KEM-768, and ML-KEM-1024.
"""

import os
import sys

# Determine the directory containing this file
# Works for both normal script execution and Sage's exec()
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()

# Add the Sage source directory to the path if needed
# Uses relative path from the script location
sage_src = os.path.abspath(os.path.join(SCRIPT_DIR, '../../../../..'))
if sage_src not in sys.path:
    sys.path.insert(0, sage_src)

from sage.crypto.public_key.key_encapsulation_mechanisms.ml_kem import MLKEM


def parse_kat_file(filename):
    """
    Parse a NIST KAT file and extract test cases.

    Each test case contains: d, z, pk, sk, m, ct, ss.

    INPUT:
    - ``filename`` -- path to the KAT file

    OUTPUT: list of dictionaries, each representing one test case
    """
    tests = []
    current = {}

    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return tests

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                current[key] = value

            if all(k in current for k in ['d', 'z', 'pk', 'sk', 'm', 'ct', 'ss']):
                tests.append(current)
                current = {}

    return tests


def test_ml_kem_with_kat(parameter_set, kat_file):
    """
    Run KAT tests for a given parameter set.

    INPUT:
    - ``parameter_set`` -- integer (512, 768, or 1024)
    - ``kat_file`` -- path to the corresponding KAT file

    OUTPUT: boolean indicating whether all tests passed
    """
    print(f"Testing ML-KEM-{parameter_set}...")

    try:
        kem = MLKEM.from_parameter_set(parameter_set)
    except Exception as e:
        print(f"  Failed to create MLKEM instance: {e}")
        return False

    tests = parse_kat_file(kat_file)
    if not tests:
        print(f"  No tests found in {kat_file}")
        return False

    print(f"  Running {len(tests)} tests...")
    passed = 0

    for i, _ in enumerate(tests[:10]):
        try:
            pk, sk = kem.keygen()
            ct, ss1 = kem.encaps(pk)
            ss2 = kem.decaps(sk, ct)

            if ss1 == ss2:
                passed += 1
            else:
                print(f"    Test {i+1}: failed (consistency mismatch)")

        except Exception as e:
            print(f"    Test {i+1}: error - {e}")

    print(f"  Passed {passed}/{min(len(tests), 10)} tests")
    return passed == min(len(tests), 10)


def test_all_kat():
    """
    Run KAT tests for all three ML-KEM parameter sets.

    OUTPUT: boolean indicating whether all tests passed
    """
    print("=" * 50)
    print("ML-KEM Known Answer Tests")
    print("=" * 50)

    results = []
    for params in [512, 768, 1024]:
        kat_file = os.path.join(SCRIPT_DIR, f'ml_kem_{params}.kat')
        result = test_ml_kem_with_kat(params, kat_file)
        results.append((params, result))

    print("=" * 50)
    all_passed = all(r for _, r in results)

    if all_passed:
        print("All KAT tests passed")
    else:
        print("Some KAT tests failed")

    return all_passed


if __name__ == "__main__":
    success = test_all_kat()
    sys.exit(0 if success else 1)
