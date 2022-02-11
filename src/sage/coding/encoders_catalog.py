r"""
Index of encoders

The ``codes.encoders`` object may be used to access the encoders that Sage can build.

**Cyclic code encoders**

- :class:`cyclic_code.CyclicCodePolynomialEncoder <sage.coding.cyclic_code.CyclicCodePolynomialEncoder>`
- :class:`cyclic_code.CyclicCodeVectorEncoder <sage.coding.cyclic_code.CyclicCodeVectorEncoder>`

**Extended code encoder**

- :class:`extended_code.ExtendedCodeExtendedMatrixEncoder <sage.coding.extended_code.ExtendedCodeExtendedMatrixEncoder>`

**Generic encoders**

- :class:`linear_code.LinearCodeGeneratorMatrixEncoder <sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder>`
- :class:`linear_code_no_metric.LinearCodeSystematicEncoder <sage.coding.linear_code_no_metric.LinearCodeSystematicEncoder>`

**Generalized Reed-Solomon code encoders**

- :class:`grs_code.GRSEvaluationVectorEncoder <sage.coding.grs_code.GRSEvaluationVectorEncoder>`
- :class:`grs_code.GRSEvaluationPolynomialEncoder <sage.coding.grs_code.GRSEvaluationPolynomialEncoder>`

**Punctured code encoder**

- :class:`punctured_code.PuncturedCodePuncturedMatrixEncoder <sage.coding.punctured_code.PuncturedCodePuncturedMatrixEncoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.encoders_catalog import *
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#                     2015 David Lucas <david.lucas@inria.fr>
#                     2016 Tania Richmond <tania.richmond@univ-tln.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import('sage.coding.cyclic_code', ['CyclicCodePolynomialEncoder',
                                         'CyclicCodeVectorEncoder'])
_lazy_import('sage.coding.extended_code', 'ExtendedCodeExtendedMatrixEncoder')
_lazy_import('sage.coding.grs_code', ['GRSEvaluationVectorEncoder',
                                      'GRSEvaluationPolynomialEncoder'])
_lazy_import('sage.coding.linear_code', 'LinearCodeGeneratorMatrixEncoder')
_lazy_import('sage.coding.linear_code_no_metric', 'LinearCodeSystematicEncoder')
_lazy_import('sage.coding.punctured_code', 'PuncturedCodePuncturedMatrixEncoder')
_lazy_import('sage.coding.reed_muller_code', ['ReedMullerVectorEncoder',
                                              'ReedMullerPolynomialEncoder'])
_lazy_import('sage.coding.subfield_subcode', 'SubfieldSubcodeParityCheckEncoder')
_lazy_import('sage.coding.parity_check_code', ['ParityCheckCodeGeneratorMatrixEncoder',
                                               'ParityCheckCodeStraightforwardEncoder'])
_lazy_import('sage.coding.goppa_code', ['GoppaCodeEncoder'])
_lazy_import('sage.coding.gabidulin_code', ['GabidulinVectorEvaluationEncoder',
                                            'GabidulinPolynomialEvaluationEncoder'])
_lazy_import('sage.coding.ag_code_decoders', ['EvaluationAGCodeEncoder',
                                              'DifferentialAGCodeEncoder'])
