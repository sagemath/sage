"""
Fast Fourier Transforms Using GSL

AUTHORS:

- William Stein (2006-9): initial file (radix2)
- D. Joyner (2006-10): Minor modifications (from radix2 to general case\
                        and some documentation).
- M. Hansen (2013-3): Fix radix2 backwards transformation
- L.F. Tabera Alonso (2013-3): Documentation
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.fft has moved"
            " to sage.calculus.transforms.fft")

from sage.calculus.transforms.fft import *
