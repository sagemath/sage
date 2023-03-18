"""
Matrix Group Homsets

AUTHORS:

- William Stein (2006-05-07): initial version

- Volker Braun (2013-1) port to new Parent, libGAP
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.misc.lazy_import import lazy_import
from sage.misc.superseded import deprecation

def is_MatrixGroupHomset(x):
    r"""
    Test whether ``x`` is a matrix group homset.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.homset import is_MatrixGroupHomset
        sage: is_MatrixGroupHomset(4)
        doctest:...: DeprecationWarning:
        Importing MatrixGroupHomset from here is deprecated; please use
        "from sage.groups.libgap_morphism import GroupHomset_libgap as MatrixGroupHomset" instead.
        See https://github.com/sagemath/sage/issues/25444 for details.
        False

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens)
        sage: from sage.groups.matrix_gps.homset import MatrixGroupHomset
        sage: M = MatrixGroupHomset(G, G)
        sage: is_MatrixGroupHomset(M)
        True
    """
    deprecation(25444, "MatrixGroupHomset is deprecated. "
                "Use GroupHomset_libgap instead.")
    return isinstance(x, MatrixGroupHomset)


lazy_import('sage.groups.libgap_morphism', 'GroupHomset_libgap',
            'MatrixGroupHomset', deprecation=25444)
