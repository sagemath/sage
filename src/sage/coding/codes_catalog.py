r"""
Index of code constructions

The ``codes`` object may be used to access the codes that Sage can build.

{INDEX_OF_FUNCTIONS}

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.codes_catalog import *

TESTS::

    sage: import sage.coding.codes_catalog
    sage: 'absolute_import' in dir(sage.coding.codes_catalog)
    False
"""
#*****************************************************************************
#       Copyright (C) 2009 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


# Implementation note:
#
# This module is imported as "codes" in all.py so that codes.<tab> is available
# in the global namespace.

from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .linear_code import LinearCode
from sage.coding.linear_code import LinearCode

lazy_import('sage.coding.code_constructions',
        ['DuadicCodeEvenPair', 'DuadicCodeOddPair',
         'ExtendedQuadraticResidueCode', 'from_parity_check_matrix',
         'QuadraticResidueCode', 'QuadraticResidueCodeEvenPair',
         'QuadraticResidueCodeOddPair', 'random_linear_code',
         'ToricCode', 'WalshCode'])

lazy_import('sage.coding.bch', 'BCHCode')
lazy_import('sage.coding.cyclic_code', 'CyclicCode')
lazy_import('sage.coding.extended_code', 'ExtendedCode')
lazy_import('sage.coding.golay_code', 'GolayCode')
lazy_import('sage.coding.grs', ['GeneralizedReedSolomonCode', 'ReedSolomonCode'])
lazy_import('sage.coding.guava', ['QuasiQuadraticResidueCode',
                                    'RandomLinearCodeGuava'])
lazy_import('sage.coding.hamming_code', 'HammingCode')
lazy_import('sage.coding.parity_check_code', 'ParityCheckCode')
lazy_import('sage.coding.punctured_code', 'PuncturedCode')
lazy_import('sage.coding.reed_muller_code', ['BinaryReedMullerCode',
                                              'ReedMullerCode'])
lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcode')

from . import decoders_catalog as decoders
from . import encoders_catalog as encoders
from . import bounds_catalog as bounds

lazy_import('sage.coding','databases')

from sage.misc.rest_index_of_methods import gen_rest_table_index
import sys
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(sys.modules[__name__], only_local_functions=False))

# We don't want this to appear in tab completion
del absolute_import, lazy_import, sys, gen_rest_table_index
