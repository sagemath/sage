"""
T-test using R

TESTS::

    sage: import rpy2                                        # optional - rpy2
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Mike Hansen   <mhansen@gmail.com>
#                     2008 Harald Schilly <harald.schilly@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.interfaces.r import R

# my own copy of an R interface
myR = R()

def ttest(x,y,conf_level = 0.95, **kw):
   """
   T-Test using R

   Arguments:

   - x, y -- vectors of same length
   - conf_level -- confidence level of the interval, [0,1) in percent

   Result:

      Tuple: (p-value, R return object)

   Example::

      sage: a, b = ttest([1,2,3,4,5],[1,2,3,3.5,5.121]); a # abs tol 1e-12  # optional - rpy2
      0.9410263720274274
   """
   if len(x) != len(y):
      raise AttributeError("vectors x and y must be of same length")

   test = myR.t_test(x,y,conf_level = conf_level, **kw)._sage_()
   t = test.get('DATA').get('p.value')
   return t, test
