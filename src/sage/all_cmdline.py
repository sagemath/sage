"""
All imports for the Sage (IPython) command line

This is all.py (load all sage functions) plus set-up for the Sage commandline.
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

sage_mode = 'cmdline'

from sage.all import *
from sage.calculus.predefined import x

sage.misc.session.init()

