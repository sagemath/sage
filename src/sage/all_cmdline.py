"""
All imports for the Sage (IPython) command line

This is all.py (load all sage functions) plus set-up for the Sage commandline.
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
sage_mode = 'cmdline'

from sage.all import *

try:
    from sage.calculus.predefined import x
except ImportError:
    pass

from sage.misc.lazy_import import lazy_import

for pkg in ['axiom', 'fricas', 'gap', 'gap3', 'giac', 'gp',
            'gnuplot', 'kash', 'magma', 'macaulay2', 'maple',
            'mathematica', 'mathics', 'matlab',
            'mupad', 'mwrank', 'octave', 'qepcad', 'singular',
            'sage0', 'lie', 'r']:
    lazy_import(f'sage.interfaces.{pkg}', f'{pkg}_console')
del pkg

lazy_import('sage.interfaces.maxima_abstract', 'maxima_console')

sage.misc.session.init()
