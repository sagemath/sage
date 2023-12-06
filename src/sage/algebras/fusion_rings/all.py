# sage_setup: distribution = sagemath-singular
"""
Fusion Rings
"""
# ****************************************************************************
#  Copyright (C) 2022 Guillermo Aboumrad <gh_willieab>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import

lazy_import('sage.algebras.fusion_rings.fusion_ring', ['FusionRing'])
lazy_import('sage.algebras.fusion_rings.fusion_double', ['FusionDouble'])
del lazy_import
