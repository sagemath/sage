r"""
Index of code constructions

The ``codes`` object may be used to access the codes that Sage can build.

Families of Codes (Rich representation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: @

    :meth:`~sage.coding.parity_check_code.ParityCheckCode` @ Parity check codes
    :meth:`~sage.coding.cyclic_code.CyclicCode` @ Cyclic codes
    :meth:`~sage.coding.bch_code.BCHCode` @ BCH Codes
    :meth:`~sage.coding.grs_code.GeneralizedReedSolomonCode` @ Generalized Reed-Solomon codes
    :meth:`~sage.coding.grs_code.ReedSolomonCode` @ Reed-Solomon codes
    :meth:`~sage.coding.reed_muller_code.BinaryReedMullerCode` @ Binary Reed-Muller codes
    :meth:`~sage.coding.reed_muller_code.ReedMullerCode` @ q-ary Reed-Muller codes
    :meth:`~sage.coding.hamming_code.HammingCode` @ Hamming codes
    :meth:`~sage.coding.golay_code.GolayCode` @ Golay codes
    :meth:`~sage.coding.goppa_code.GoppaCode` @ Goppa codes
    :meth:`~sage.coding.kasami_codes.KasamiCode` @ Kasami codes

Families of Codes (Generator matrix representation)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: @

    :meth:`~sage.coding.code_constructions.DuadicCodeEvenPair` @ Duadic codes, even pair
    :meth:`~sage.coding.code_constructions.DuadicCodeOddPair` @ Duadic codes, odd pair
    :meth:`~sage.coding.code_constructions.QuadraticResidueCode` @  Quadratic residue codes
    :meth:`~sage.coding.code_constructions.ExtendedQuadraticResidueCode` @ Extended quadratic residue codes
    :meth:`~sage.coding.code_constructions.QuadraticResidueCodeEvenPair` @ Even-like quadratic residue codes
    :meth:`~sage.coding.code_constructions.QuadraticResidueCodeOddPair` @ Odd-like quadratic residue codes
    :meth:`~sage.coding.guava.QuasiQuadraticResidueCode` @ Quasi quadratic residue codes (Requires GAP/Guava)
    :meth:`~sage.coding.code_constructions.ToricCode` @ Toric codes
    :meth:`~sage.coding.code_constructions.WalshCode` @ Walsh codes
    :meth:`~sage.coding.code_constructions.from_parity_check_matrix` @ Construct a code from a parity check matrix
    :meth:`~sage.coding.code_constructions.random_linear_code` @ Construct a random linear code
    :meth:`~sage.coding.guava.RandomLinearCodeGuava` @ Construct a random linear code through Guava (Requires GAP/Guava)


Derived Codes
^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: @

    :meth:`~sage.coding.subfield_subcode.SubfieldSubcode` @ Subfield subcodes
    :meth:`~sage.coding.extended_code.ExtendedCode` @ Extended codes
    :meth:`~sage.coding.punctured_code.PuncturedCode` @ Puncturedcodes

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.codes_catalog import *

"""
#*****************************************************************************
#       Copyright (C) 2009 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

# This module is imported as "codes" in all.py so that codes.<tab> is
# available in the global namespace.

from sage.misc.lazy_import import lazy_import as _lazy_import

from .linear_code import LinearCode

_lazy_import('sage.coding.code_constructions',
        ['DuadicCodeEvenPair', 'DuadicCodeOddPair',
         'ExtendedQuadraticResidueCode', 'from_parity_check_matrix',
         'QuadraticResidueCode', 'QuadraticResidueCodeEvenPair',
         'QuadraticResidueCodeOddPair', 'random_linear_code',
         'ToricCode', 'WalshCode'])

_lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcode')
_lazy_import('sage.coding.extended_code', 'ExtendedCode')
_lazy_import('sage.coding.punctured_code', 'PuncturedCode')

_lazy_import('sage.coding.parity_check_code', 'ParityCheckCode')
_lazy_import('sage.coding.cyclic_code', 'CyclicCode')
_lazy_import('sage.coding.bch_code', 'BCHCode')
_lazy_import('sage.coding.grs_code', ['GeneralizedReedSolomonCode', 'ReedSolomonCode'])
_lazy_import('sage.coding.reed_muller_code', ['BinaryReedMullerCode', 'ReedMullerCode'])
_lazy_import('sage.coding.hamming_code', 'HammingCode')
_lazy_import('sage.coding.golay_code', 'GolayCode')
_lazy_import('sage.coding.goppa_code', 'GoppaCode')
_lazy_import('sage.coding.kasami_codes', 'KasamiCode')
_lazy_import('sage.coding.linear_rank_metric', 'LinearRankMetricCode')
_lazy_import('sage.coding.gabidulin_code', 'GabidulinCode')
_lazy_import('sage.coding.ag_code', ['EvaluationAGCode', 'DifferentialAGCode', 'CartierCode'])

_lazy_import('sage.coding.guava', ['QuasiQuadraticResidueCode', 'RandomLinearCodeGuava'])

from . import decoders_catalog as decoders
from . import encoders_catalog as encoders
from . import bounds_catalog as bounds

_lazy_import('sage.coding','databases')
