r"""
Homsets of finitely presented graded modules over the Steenrod algebra

EXAMPLES::

    sage: from sage.modules.fp_graded.steenrod.module import SteenrodFPModule
    sage: from sage.misc.sage_unittest import TestSuite
    sage: A = SteenrodAlgebra(2, profile=(3,2,1))
    sage: F = SteenrodFPModule(A, [1,3], names='c')
    sage: L = SteenrodFPModule(A, [2,3], [[Sq(2),Sq(1)], [0,Sq(2)]], names='h')
    sage: homset = Hom(F, L); homset
    Set of Morphisms from Free graded left module on 2 generators ...
    sage: homset.an_element()
    Module morphism:
      From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: c[1] |--> 0
            c[3] |--> Sq(1)*h[2]
    sage: f = homset([L((Sq(1), 1)), L((0, Sq(2)))]); f
    Module morphism:
      From: Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Finitely presented left module on 2 generators and 2 relations over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: c[1] |--> Sq(1)*h[2] + h[3]
            c[3] |--> Sq(2)*h[3]
    sage: f.kernel_inclusion()
    Module morphism:
      From: Finitely presented left module on 2 generators and 1 relation over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      To:   Free graded left module on 2 generators over sub-Hopf algebra of mod 2 Steenrod algebra, milnor basis, profile function [3, 2, 1]
      Defn: g[3] |--> c[3]
            g[4] |--> Sq(0,1)*c[1]
    sage: TestSuite(homset).run()

AUTHORS:

- Robert R. Bruner, Michael J. Catanzaro (2012): Initial version.
- Sverre Lunoee--Nielsen and Koen van Woerden (2019-11-29): Updated the
  original software to Sage version 8.9.
- Sverre Lunoee--Nielsen (2020-07-01): Refactored the code and added
  new documentation and tests.
"""

# ****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modules.fp_graded.homspace import FPModuleHomspace
from sage.modules.fp_graded.free_homspace import FreeGradedModuleHomspace
from .morphism import SteenrodFPModuleMorphism, SteenrodFreeModuleMorphism


class SteenrodFPModuleHomspace(FPModuleHomspace):
    Element = SteenrodFPModuleMorphism


class SteenrodFreeModuleHomspace(SteenrodFPModuleHomspace, FreeGradedModuleHomspace):
    Element = SteenrodFreeModuleMorphism
