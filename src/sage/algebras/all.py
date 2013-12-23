"""
Algebras
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
from sage.misc.lazy_import import lazy_import

from quatalg.all import *

# Algebra base classes
from algebra import Algebra, is_Algebra

# Ring element base classes
from algebra_element import AlgebraElement, is_AlgebraElement


from free_algebra import FreeAlgebra, is_FreeAlgebra
from free_algebra_quotient import FreeAlgebraQuotient

from steenrod.all import *

from group_algebra_new import GroupAlgebra

from iwahori_hecke_algebra import IwahoriHeckeAlgebra, IwahoriHeckeAlgebraT
from affine_nil_temperley_lieb import AffineNilTemperleyLiebTypeA
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra')

lazy_import('sage.algebras.hall_algebra', 'HallAlgebra')

lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra')
