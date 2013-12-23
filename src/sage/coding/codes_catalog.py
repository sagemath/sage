r"""
Index of Codes

The ``codes`` object may be used to access the codes that Sage can build.

- :func:`codes.BCHCode <sage.coding.code_constructions.BCHCode>`
- :func:`codes.BinaryGolayCode <sage.coding.code_constructions.BinaryGolayCode>`
- :func:`codes.BinaryReedMullerCode <sage.coding.guava.BinaryReedMullerCode>`
- :func:`codes.CyclicCode <sage.coding.code_constructions.CyclicCode>`
- :func:`codes.CyclicCodeFromGeneratingPolynomial <sage.coding.code_constructions.CyclicCodeFromGeneratingPolynomial>`
- :func:`codes.CyclicCodeFromCheckPolynomial <sage.coding.code_constructions.CyclicCodeFromCheckPolynomial>`
- :func:`codes.DuadicCodeEvenPair <sage.coding.code_constructions.DuadicCodeEvenPair>`
- :func:`codes.DuadicCodeOddPair <sage.coding.code_constructions.DuadicCodeOddPair>`
- :func:`codes.ExtendedBinaryGolayCode <sage.coding.code_constructions.ExtendedBinaryGolayCode>`
- :func:`codes.ExtendedQuadraticResidueCode <sage.coding.code_constructions.ExtendedQuadraticResidueCode>`
- :func:`codes.ExtendedTernaryGolayCode <sage.coding.code_constructions.ExtendedTernaryGolayCode>`
- :func:`codes.HammingCode <sage.coding.code_constructions.HammingCode>`
- :func:`codes.LinearCodeFromCheckMatrix <sage.coding.code_constructions.LinearCodeFromCheckMatrix>`
- :func:`codes.QuadraticResidueCode <sage.coding.code_constructions.QuadraticResidueCode>`
- :func:`codes.QuadraticResidueCodeEvenPair <sage.coding.code_constructions.QuadraticResidueCodeEvenPair>`
- :func:`codes.QuadraticResidueCodeOddPair <sage.coding.code_constructions.QuadraticResidueCodeOddPair>`
- :func:`codes.QuasiQuadraticResidueCode <sage.coding.guava.QuasiQuadraticResidueCode>`
- :func:`codes.RandomLinearCode <sage.coding.code_constructions.RandomLinearCode>`
- :func:`codes.RandomLinearCodeGuava <sage.coding.guava.RandomLinearCodeGuava>`
- :func:`codes.ReedSolomonCode <sage.coding.code_constructions.ReedSolomonCode>`
- :func:`codes.TernaryGolayCode <sage.coding.code_constructions.TernaryGolayCode>`
- :func:`codes.ToricCode <sage.coding.code_constructions.ToricCode>`
- :func:`codes.TrivialCode <sage.coding.code_constructions.TrivialCode>`
- :func:`codes.WalshCode <sage.coding.code_constructions.WalshCode>`

"""

# Implementation note:
#
# This module is imported as "codes" in all.py so that codes.<tab> is available
# in the global namespace.

from code_constructions import (BCHCode, BinaryGolayCode, CyclicCodeFromGeneratingPolynomial,
                                CyclicCode, CyclicCodeFromCheckPolynomial, DuadicCodeEvenPair,
                                DuadicCodeOddPair, ExtendedBinaryGolayCode,
                                ExtendedQuadraticResidueCode, ExtendedTernaryGolayCode,
                                HammingCode, LinearCodeFromCheckMatrix,
                                QuadraticResidueCode, QuadraticResidueCodeEvenPair,
                                QuadraticResidueCodeOddPair, RandomLinearCode,
                                ReedSolomonCode, TernaryGolayCode,
                                ToricCode, TrivialCode, WalshCode)

from guava import BinaryReedMullerCode, QuasiQuadraticResidueCode, RandomLinearCodeGuava
