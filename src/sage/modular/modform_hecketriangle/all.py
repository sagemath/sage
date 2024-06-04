r"""
AUTHORS:

- Jonas Jermann (2013): initial version
"""
# ****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup

from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor

from sage.modular.modform_hecketriangle.graded_ring import (QuasiMeromorphicModularFormsRing,
                                                            QuasiWeakModularFormsRing, QuasiModularFormsRing,
                                                            QuasiCuspFormsRing, MeromorphicModularFormsRing,
                                                            WeakModularFormsRing,
                                                            ModularFormsRing, CuspFormsRing)

from sage.modular.modform_hecketriangle.space import (QuasiMeromorphicModularForms, QuasiWeakModularForms,
                                                      QuasiModularForms, QuasiCuspForms,
                                                      MeromorphicModularForms, WeakModularForms, ModularForms,
                                                      CuspForms, ZeroForm)

from sage.modular.modform_hecketriangle.subspace import ModularFormsSubSpace
