# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Bijection between rigged configurations and KR tableaux

Functions which are big switch statements to create the bijection class of the
correct type.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version
- Travis Scrimshaw (2012-12-21): Added all non-exceptional bijection types
- Travis Scrimshaw (2014-09-10): Added type `D_4^{(3)}`
"""

#*****************************************************************************
#       Copyright (C) 2011-2015 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA

from sage.combinat.rigged_configurations.bij_type_B import KRTToRCBijectionTypeB
from sage.combinat.rigged_configurations.bij_type_B import RCToKRTBijectionTypeB

from sage.combinat.rigged_configurations.bij_type_C import KRTToRCBijectionTypeC
from sage.combinat.rigged_configurations.bij_type_C import RCToKRTBijectionTypeC

from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD

from sage.combinat.rigged_configurations.bij_type_D_twisted import KRTToRCBijectionTypeDTwisted
from sage.combinat.rigged_configurations.bij_type_D_twisted import RCToKRTBijectionTypeDTwisted

from sage.combinat.rigged_configurations.bij_type_A2_even import KRTToRCBijectionTypeA2Even
from sage.combinat.rigged_configurations.bij_type_A2_even import RCToKRTBijectionTypeA2Even

from sage.combinat.rigged_configurations.bij_type_A2_dual import KRTToRCBijectionTypeA2Dual
from sage.combinat.rigged_configurations.bij_type_A2_dual import RCToKRTBijectionTypeA2Dual

from sage.combinat.rigged_configurations.bij_type_A2_odd import KRTToRCBijectionTypeA2Odd
from sage.combinat.rigged_configurations.bij_type_A2_odd import RCToKRTBijectionTypeA2Odd

from sage.combinat.rigged_configurations.bij_type_D_tri import KRTToRCBijectionTypeDTri
from sage.combinat.rigged_configurations.bij_type_D_tri import RCToKRTBijectionTypeDTri

from sage.combinat.rigged_configurations.bij_type_E67 import KRTToRCBijectionTypeE67
from sage.combinat.rigged_configurations.bij_type_E67 import RCToKRTBijectionTypeE67


def KRTToRCBijection(tp_krt):
    r"""
    Return the correct KR tableaux to rigged configuration bijection helper class.

    TESTS::

        sage: KRT = crystals.TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
        sage: from sage.combinat.rigged_configurations.bijection import KRTToRCBijection
        sage: bijection = KRTToRCBijection(KRT(pathlist=[[5,2]]))
    """
    ct = tp_krt.cartan_type()
    typ = ct.type()
    if ct.is_untwisted_affine():
        if typ == 'A':
            return KRTToRCBijectionTypeA(tp_krt)
        if typ == 'B':
            return KRTToRCBijectionTypeB(tp_krt)
        if typ == 'C':
            return KRTToRCBijectionTypeC(tp_krt)
        if typ == 'D':
            return KRTToRCBijectionTypeD(tp_krt)
        if typ == 'E':
            if ct.classical().rank() < 8:
                return KRTToRCBijectionTypeE67(tp_krt)
        #if typ == 'F':
        #if typ == 'G':
    else:
        if typ == 'BC': # A_{2n}^{(2)}
            return KRTToRCBijectionTypeA2Even(tp_krt)
        typ = ct.dual().type()
        if typ == 'BC': # A_{2n}^{(2)\dagger}
            return KRTToRCBijectionTypeA2Dual(tp_krt)
        if typ == 'B': # A_{2n-1}^{(2)}
            return KRTToRCBijectionTypeA2Odd(tp_krt)
        if typ == 'C': # D_{n+1}^{(2)}
            return KRTToRCBijectionTypeDTwisted(tp_krt)
        #if typ == 'F': # E_6^{(2)}
        if typ == 'G': # D_4^{(3)}
            return KRTToRCBijectionTypeDTri(tp_krt)
    raise NotImplementedError


def RCToKRTBijection(rigged_configuration_elt):
    r"""
    Return the correct rigged configuration to KR tableaux bijection helper class.

    TESTS::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
        sage: from sage.combinat.rigged_configurations.bijection import RCToKRTBijection
        sage: bijection = RCToKRTBijection(RC(partition_list=[[1],[1],[1],[1]]))
    """
    ct = rigged_configuration_elt.cartan_type()
    typ = ct.type()
    if not ct.is_affine() or ct.is_untwisted_affine():
        if typ == 'A':
            return RCToKRTBijectionTypeA(rigged_configuration_elt)
        if typ == 'B':
            return RCToKRTBijectionTypeB(rigged_configuration_elt)
        if typ == 'C':
            return RCToKRTBijectionTypeC(rigged_configuration_elt)
        if typ == 'D':
            return RCToKRTBijectionTypeD(rigged_configuration_elt)
        if typ == 'E':
            if ct.classical().rank() < 8:
                return RCToKRTBijectionTypeE67(rigged_configuration_elt)
        #if typ == 'F':
        #if typ == 'G':
    else:
        if typ == 'BC': # A_{2n}^{(2)}
            return RCToKRTBijectionTypeA2Even(rigged_configuration_elt)
        typ = ct.dual().type()
        if typ == 'BC': # A_{2n}^{(2)\dagger}
            return RCToKRTBijectionTypeA2Dual(rigged_configuration_elt)
        if typ == 'B': # A_{2n-1}^{(2)}
            return RCToKRTBijectionTypeA2Odd(rigged_configuration_elt)
        if typ == 'C': # D_{n+1}^{(2)}
            return RCToKRTBijectionTypeDTwisted(rigged_configuration_elt)
        #if typ == 'F': # E_6^{(2)}
        if typ == 'G': # D_4^{(3)}
            return RCToKRTBijectionTypeDTri(rigged_configuration_elt)
    raise NotImplementedError
