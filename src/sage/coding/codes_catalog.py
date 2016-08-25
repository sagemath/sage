r"""
Index of codes

The ``codes`` object may be used to access the codes that Sage can build.

{INDEX_OF_FUNCTIONS}

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.codes_catalog import *

"""
from __future__ import absolute_import

# Implementation note:
#
# This module is imported as "codes" in all.py so that codes.<tab> is available
# in the global namespace.

from sage.misc.lazy_import import lazy_import as _lazy_import
from .code_constructions import (BCHCode, BinaryGolayCode, CyclicCodeFromGeneratingPolynomial,
                                CyclicCode, CyclicCodeFromCheckPolynomial, DuadicCodeEvenPair,
                                DuadicCodeOddPair, ExtendedBinaryGolayCode,
                                ExtendedQuadraticResidueCode, ExtendedTernaryGolayCode,
                                LinearCode, LinearCodeFromCheckMatrix,
                                QuadraticResidueCode, QuadraticResidueCodeEvenPair,
                                QuadraticResidueCodeOddPair, RandomLinearCode,
                                ReedSolomonCode, TernaryGolayCode,
                                ToricCode, TrivialCode, WalshCode)

from .grs import GeneralizedReedSolomonCode
from .reed_muller_code import ReedMullerCode, BinaryReedMullerCode
from .extended_code import ExtendedCode
from .subfield_subcode import SubfieldSubcode

from .guava import QuasiQuadraticResidueCode, RandomLinearCodeGuava
_lazy_import('sage.coding.punctured_code', 'PuncturedCode')
from .hamming_code import HammingCode
from . import decoders_catalog as decoders
from . import encoders_catalog as encoders
from . import bounds_catalog as bounds
from sage.misc.rest_index_of_methods import gen_rest_table_index as _gen_rest_table_index
import sys as _sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=_gen_rest_table_index(_sys.modules[__name__], only_local_functions=False))
