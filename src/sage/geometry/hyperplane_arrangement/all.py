# sage_setup: distribution = sagemath-polyhedra
"""
Hyperplane arrangements
"""

# *****************************************************************************
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <was@math.harvard.edu>
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
# *****************************************************************************

from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.hyperplane_arrangement.arrangement', 'HyperplaneArrangements')

lazy_import('sage.geometry.hyperplane_arrangement.ordered_arrangement', 'OrderedHyperplaneArrangements')

lazy_import('sage.geometry.hyperplane_arrangement.library', 'hyperplane_arrangements')
