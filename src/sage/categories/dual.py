"""
Dual functorial construction

AUTHORS:

 - Nicolas M. Thiery (2009-2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# could do SelfDualCategory

from sage.categories.covariant_functorial_construction import CovariantFunctorialConstruction, CovariantConstructionCategory

class DualFunctor(CovariantFunctorialConstruction):
    """
    A singleton class for the dual functor
    """
    _functor_name = "dual"
    _functor_category = "DualObjects"
    symbol = "^*"

class DualObjectsCategory(CovariantConstructionCategory):

    _functor_category = "DualObjects"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).DualObjects() # indirect doctest
            Category of duals of vector spaces over Rational Field
        """
        # Just to remove the `objects`
        return "duals of %s"%(self.base_category()._repr_object_names())
