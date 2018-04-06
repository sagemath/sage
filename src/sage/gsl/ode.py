r"""
Solving ODE numerically by GSL

AUTHORS:

- Joshua Kantor (2004-2006)

- Robert Marik (2010 - fixed docstrings)

"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006 Joshua Kantor <kantor.jm@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.ode has moved to sage.calculus.ode")

from sage.calculus.ode import *
