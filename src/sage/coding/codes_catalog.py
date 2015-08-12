r"""
Index of Codes

The ``codes`` object may be used to access the codes that Sage can build.

{INDEX_OF_FUNCTIONS}

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.codes_catalog import *

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

from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__], only_local_functions=False))
