"""
Lie Algebras
"""

# *****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.lie_algebras.lie_algebra', 'LieAlgebra')
# from kac_moody import KacMoodyAlgebra
lazy_import('sage.algebras.lie_algebras', 'examples', 'lie_algebras')
del lazy_import
