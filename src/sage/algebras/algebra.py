"""
Abstract base class for algebras
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

from sage.rings.ring import Algebra
from sage.categories.algebras import Algebras

def is_Algebra(x):
    r"""
    Return True if x is an Algebra.

    EXAMPLES::

        sage: from sage.algebras.algebra import is_Algebra
        sage: R.<x,y> = FreeAlgebra(QQ,2)
        sage: is_Algebra(R)
        True
    """
    try:
        return isinstance(x, Algebra) or x in Algebras(x.base_ring())
    except Exception:
        return False

