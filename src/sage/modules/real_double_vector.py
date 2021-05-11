"""
Pickling for the old RDF vector class

AUTHORS:

- Jason Grout

TESTS::

    sage: v = vector(RDF, [1,2,3,4])
    sage: loads(dumps(v)) == v
    True
"""
###############################################################################
#       Copyright (C) 2008 Jason Grout <jason-sage@creativetrax.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################
from .vector_real_double_dense import Vector_real_double_dense

RealDoubleVectorSpaceElement = Vector_real_double_dense
