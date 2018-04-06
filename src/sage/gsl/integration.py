"""
Numerical Integration

AUTHORS:

- Josh Kantor (2007-02): first version

- William Stein (2007-02): rewrite of docs, conventions, etc.

- Robert Bradshaw (2008-08): fast float integration

- Jeroen Demeyer (2011-11-23): :trac:`12047`: return 0 when the
  integration interval is a point; reformat documentation and add to
  the reference manual.
"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006,2007 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.integration has moved"
            " to sage.calculus.integration")

from sage.calculus.integration import *
