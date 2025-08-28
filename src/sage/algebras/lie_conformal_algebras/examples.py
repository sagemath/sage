r"""
Examples of Lie Conformal Algebras

We implement the following examples of Lie conformal algebras:

- :mod:`Abelian Lie conformal algebra<.abelian_lie_conformal_algebra>`
- :mod:`Affine Lie conformal algebra<.affine_lie_conformal_algebra>`
- :mod:`Bosonic Ghosts<.bosonic_ghosts_lie_conformal_algebra>`
- :mod:`Fermionic Ghosts<.fermionic_ghosts_lie_conformal_algebra>`
- :mod:`Free Bosons<.free_bosons_lie_conformal_algebra>`
- :mod:`Free Fermions<.free_fermions_lie_conformal_algebra>`
- :mod:`N=2 super Lie Conformal algebra<.n2_lie_conformal_algebra>`
- :mod:`Neveu-Schwarz super Lie conformal
  algebra<.neveu_schwarz_lie_conformal_algebra>`
- :mod:`Virasoro Lie conformal algebra<.virasoro_lie_conformal_algebra>`
- :mod:`Weyl Lie conformal algebra<.weyl_lie_conformal_algebra>`

AUTHORS:

- Reimundo Heluani (2020-06-15): Initial implementation.
"""

# *****************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .abelian_lie_conformal_algebra import AbelianLieConformalAlgebra as Abelian
from .affine_lie_conformal_algebra import AffineLieConformalAlgebra as Affine
from .bosonic_ghosts_lie_conformal_algebra import BosonicGhostsLieConformalAlgebra as BosonicGhosts
from .fermionic_ghosts_lie_conformal_algebra import FermionicGhostsLieConformalAlgebra as FermionicGhosts
from .free_bosons_lie_conformal_algebra import FreeBosonsLieConformalAlgebra as FreeBosons
from .free_fermions_lie_conformal_algebra import FreeFermionsLieConformalAlgebra as FreeFermions
from .n2_lie_conformal_algebra import N2LieConformalAlgebra as N2
from .neveu_schwarz_lie_conformal_algebra import NeveuSchwarzLieConformalAlgebra as NeveuSchwarz
from .virasoro_lie_conformal_algebra import VirasoroLieConformalAlgebra as Virasoro
from .weyl_lie_conformal_algebra import WeylLieConformalAlgebra as Weyl
