# sage.doctest: needs sage.combinat sage.modules
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

from sage.categories.algebras import Algebras


def is_Algebra(x):
    r"""
    Return ``True`` if `x` is an Algebra.

    EXAMPLES::

        sage: from sage.algebras.algebra import is_Algebra
        sage: R.<x,y> = FreeAlgebra(QQ,2)
        sage: is_Algebra(R)
        doctest:warning...
        DeprecationWarning: the function is_Algebra is deprecated; use '... in Algebras(base_ring)' instead
        See https://github.com/sagemath/sage/issues/35253 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(35253, "the function is_Algebra is deprecated; use '... in Algebras(base_ring)' instead")
    return x in Algebras(x.base_ring())
