"""
Abstract base class for commutative algebras
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

from sage.rings.ring import CommutativeAlgebra

def is_CommutativeAlgebra(x):
    """
    Check to see if ``x`` is a :class:`CommutativeAlgebra`.

    EXAMPLES::

        sage: from sage.rings.commutative_algebra import is_CommutativeAlgebra
        sage: from sage.rings.ring import CommutativeAlgebra
        sage: is_CommutativeAlgebra(CommutativeAlgebra(ZZ))
        True
    """
    return isinstance(x, CommutativeAlgebra)
