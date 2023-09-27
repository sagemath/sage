"""
Homomorphisms Between Matrix Groups

Deprecated May, 2018; use :class:`sage.groups.libgap_morphism` instead.
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_import import lazy_import

def to_libgap(x):
    """
    Helper to convert ``x`` to a LibGAP matrix or matrix group
    element.

    Deprecated; use the ``x.gap()`` method or ``libgap(x)`` instead.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.morphism import to_libgap
        sage: to_libgap(GL(2,3).gen(0))
        doctest:...: DeprecationWarning: this function is deprecated.
         Use x.gap() or libgap(x) instead.
        See https://github.com/sagemath/sage/issues/25444 for details.
        [ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
        sage: to_libgap(matrix(QQ, [[1,2],[3,4]]))
        [ [ 1, 2 ], [ 3, 4 ] ]
    """
    from sage.misc.superseded import deprecation
    deprecation(25444, "this function is deprecated."
                " Use x.gap() or libgap(x) instead.")
    try:
        return x.gap()
    except AttributeError:
        from sage.libs.gap.libgap import libgap
        return libgap(x)


lazy_import('sage.groups.libgap_morphism', 'GroupMorphism_libgap',
            'MatrixGroupMorphism_im_gens', deprecation=25444)

lazy_import('sage.categories.morphism', 'Morphism',
            'MatrixGroupMap', deprecation=25444)

lazy_import('sage.categories.morphism', 'Morphism',
            'MatrixGroupMorphism', deprecation=25444)
