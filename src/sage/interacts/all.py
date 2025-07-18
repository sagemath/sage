"""
Interacts included with sage

AUTHORS:

- Harald Schilly (2011-01-16): initial version (#9623) partially based on work by Lauri Ruotsalainen
"""

# ****************************************************************************
#       Copyright (C) 2011 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.lazy_import import lazy_import

from sage.interacts import calculus
from sage.interacts import geometry
from sage.interacts import statistics
from sage.interacts import fractals
from sage.interacts import algebra
lazy_import('sage.interacts.library', 'demo')
del lazy_import
