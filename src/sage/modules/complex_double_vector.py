r"""
Pickling for the old CDF vector class

AUTHORS:

- Jason Grout

TESTS::

    sage: v = vector(CDF,[(1,-1), (2,pi), (3,5)]); v                                    # needs sage.symbolic
    (1.0 - 1.0*I, 2.0 + 3.141592653589793*I, 3.0 + 5.0*I)
    sage: loads(dumps(v)) == v                                                          # needs sage.symbolic
    True
"""

###############################################################################
#       Copyright (C) 2008 Jason Grout <jason-sage@creativetrax.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################

from sage.modules.vector_complex_double_dense import Vector_complex_double_dense

ComplexDoubleVectorSpaceElement = Vector_complex_double_dense
