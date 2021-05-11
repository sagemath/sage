"""
All imports for the Sage (IPython) command line

This is all.py (load all sage functions) plus set-up for the Sage commandline.
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

sage_mode = 'cmdline'

from sage.all import *
from sage.calculus.predefined import x

sage.misc.session.init()

