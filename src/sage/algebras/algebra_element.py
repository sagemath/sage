"""
Base class for algebra elements
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

from sage.misc.superseded import deprecation
deprecation(21141, "the module sage.algebras.algebra_element is deprecated, import from sage.structure.element instead")

from sage.structure.element import AlgebraElement

def is_AlgebraElement(x):
    r"""
    Return true if x is an AlgebraElement

    EXAMPLES:

    """
    return isinstance(x, AlgebraElement)
