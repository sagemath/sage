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

from . import calculus
from . import geometry
from . import statistics
from . import fractals
from . import algebra
lazy_import('sage.interacts.library', 'demo')
