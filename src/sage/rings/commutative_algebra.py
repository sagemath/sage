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
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.commutative_algebras import CommutativeAlgebras


def is_CommutativeAlgebra(x):
    """
    Check to see if ``x`` is in the category of ``CommutativeAlgebras``.

    EXAMPLES::

        sage: from sage.rings.commutative_algebra import is_CommutativeAlgebra
        sage: is_CommutativeAlgebra(QQ['x'])
        doctest:warning...
        DeprecationWarning: the function is_CommutativeAlgebra is deprecated; use '... in Algebras(base_ring).Commutative()' instead
        See https://github.com/sagemath/sage/issues/35999 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(35999, "the function is_CommutativeAlgebra is deprecated; use '... in Algebras(base_ring).Commutative()' instead")
    return x in CommutativeAlgebras(x.base_ring())
