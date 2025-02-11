r"""
Homsets of finitely generated free graded left modules

For an overview, see the :mod:`free graded modules documentation
<sage.modules.fp_graded.free_module>`.

EXAMPLES::

    sage: from sage.modules.fp_graded.free_module import FreeGradedModule
    sage: A = SteenrodAlgebra(2)
    sage: F1 = FreeGradedModule(A, (1,3), names='g')
    sage: F2 = FreeGradedModule(A, (2,3), names='h')
    sage: homset = Hom(F1, F2)
    sage: homset
    Set of Morphisms from Free graded left module on 2 generators ...
    sage: homset([F2((Sq(1), 1)), F2((0, Sq(2)))])
    Module morphism:
      From: Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      To:   Free graded left module on 2 generators over mod 2 Steenrod algebra, milnor basis
      Defn: g[1] |--> Sq(1)*h[2] + h[3]
            g[3] |--> Sq(2)*h[3]
    sage: TestSuite(homset).run()

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# *****************************************************************************
#       Copyright (C) 2021 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.fp_graded.free_morphism import FreeGradedModuleMorphism
from sage.modules.fp_graded.homspace import FPModuleHomspace


class FreeGradedModuleHomspace(FPModuleHomspace):
    """
    Homspace between two free graded modules.
    """
    Element = FreeGradedModuleMorphism
