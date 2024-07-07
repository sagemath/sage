# sage_setup: distribution = sagemath-objects
"""
Subquotient Functorial Construction

AUTHORS:

 - Nicolas M. Thiery (2010): initial revision
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class SubquotientsCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Subquotients"
