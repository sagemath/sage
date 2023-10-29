r"""
A subset of sage.categories.all with just the basic categories needed
for sage startup (i.e. to define ZZ, QQ, ...).
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.posets import Posets

# For backward compatibility; will be deprecated at some point
PartiallyOrderedSets = Posets
OrderedSets = Posets


from sage.categories.partially_ordered_monoids import PartiallyOrderedMonoids

# For backward compatibility; might be deprecated at some point
OrderedMonoids = PartiallyOrderedMonoids



