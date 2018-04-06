"""
Discrete Wavelet Transform

Wraps ``GSL's gsl_wavelet_transform_forward()``,
and ``gsl_wavelet_transform_inverse()`` and creates plot methods.

AUTHOR:

- Josh Kantor (2006-10-07)  - initial version
- David Joyner (2006-10-09) - minor changes to docstrings and examples.

"""

#*****************************************************************************
#       Copyright (C) 2006 Joshua Kantor <jkantor@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(9084, "the module sage.gsl.dwt has moved to sage.calculus.transforms.dwt")

from sage.calculus.transforms.dwt import *
