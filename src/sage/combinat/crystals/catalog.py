r"""
Catalog Of Crystals

Let `I` be an index set and let `(A,\Pi,\Pi^\vee,P,P^\vee)` be a Cartan datum
associated with generalized Cartan matrix `A = (a_{ij})_{i,j\in I}`.  An
*abstract crystal* associated to this Cartan datum is a set `B` together with
maps

.. MATH::

    e_i,f_i \colon B \to B \cup \{0\}, \qquad
    \varepsilon_i,\varphi_i\colon B \to \ZZ \cup \{-\infty\}, \qquad
    \mathrm{wt}\colon B \to P,

subject to the following conditions:

    1. `\varphi_i(b) = \varepsilon_i(b) + \langle h_i, \mathrm{wt}(b) \rangle` for all `b \in B` and `i \in I`;
    2. `\mathrm{wt}(e_ib) = \mathrm{wt}(b) + \alpha_i` if `e_ib \in B`;
    3. `\mathrm{wt}(f_ib) = \mathrm{wt}(b) - \alpha_i` if `f_ib \in B`;
    4. `\varepsilon_i(e_ib) = \varepsilon_i(b) - 1`, `\varphi_i(e_ib) = \varphi_i(b) + 1` if `e_ib \in B`;
    5. `\varepsilon_i(f_ib) = \varepsilon_i(b) + 1`, `\varphi_i(f_ib) = \varphi_i(b) - 1` if `f_ib \in B`;
    6. `f_ib = b'` if and only if `b = e_ib'` for `b,b' \in B` and `i\in I`;
    7. if `\varphi_i(b) = -\infty` for `b\in B`, then `e_ib = f_ib = 0`.

.. SEEALSO::

    - :mod:`sage.categories.crystals`
    - :mod:`sage.combinat.crystals.crystals`

Catalog
-------

This is a catalog of crystals that are currently implemented in Sage:

* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassical`
* :class:`~sage.combinat.crystals.affine.AffineCrystalFromClassicalAndPromotion`
* :class:`AffineFactorization <sage.combinat.crystals.affine_factorization.AffineFactorizationCrystal>`
* :class:`AffinizationOf <sage.combinat.crystals.affinization.AffinizationOfCrystal>`
* :class:`AlcovePaths <sage.combinat.crystals.alcove_path.CrystalOfAlcovePaths>`
* :class:`FastRankTwo <sage.combinat.crystals.fast_crystals.FastCrystal>`
* :class:`FullyCommutativeStableGrothendieck
  <sage.combinat.crystals.fully_commutative_stable_grothendieck.FullyCommutativeStableGrothendieckCrystal>`
* :class:`GeneralizedYoungWalls
  <sage.combinat.crystals.generalized_young_walls.CrystalOfGeneralizedYoungWalls>`
* :func:`HighestWeight <sage.combinat.crystals.highest_weight_crystals.HighestWeightCrystal>`
* :class:`Induced <sage.combinat.crystals.induced_structure.InducedCrystal>`
* :class:`KacModule <sage.combinat.crystals.kac_modules.CrystalOfKacModule>`
* :func:`KirillovReshetikhin <sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystal>`
* :class:`KleshchevPartitions <sage.combinat.partition_kleshchev.KleshchevPartitions_all>`
* :class:`KyotoPathModel <sage.combinat.crystals.kyoto_path_model.KyotoPathModel>`
* :class:`Letters <sage.combinat.crystals.letters.CrystalOfLetters>`
* :class:`LSPaths <sage.combinat.crystals.littelmann_path.CrystalOfLSPaths>`
* :class:`Minimaj <sage.combinat.multiset_partition_into_sets_ordered.MinimajCrystal>`
* :class:`NakajimaMonomials <sage.combinat.crystals.monomial_crystals.CrystalOfNakajimaMonomials>`
* :class:`OddNegativeRoots <sage.combinat.crystals.kac_modules.CrystalOfOddNegativeRoots>`
* :class:`ProjectedLevelZeroLSPaths
  <sage.combinat.crystals.littelmann_path.CrystalOfProjectedLevelZeroLSPaths>`
* :class:`RiggedConfigurations
  <sage.combinat.rigged_configurations.rc_crystal.CrystalOfRiggedConfigurations>`
* :class:`ShiftedPrimedTableaux
  <sage.combinat.shifted_primed_tableau.ShiftedPrimedTableaux_shape>`
* :class:`Spins <sage.combinat.crystals.spins.CrystalOfSpins>`
* :class:`SpinsPlus <sage.combinat.crystals.spins.CrystalOfSpinsPlus>`
* :class:`SpinsMinus <sage.combinat.crystals.spins.CrystalOfSpinsMinus>`
* :class:`Tableaux <sage.combinat.crystals.tensor_product.CrystalOfTableaux>`

Subcatalogs:

* :ref:`sage.combinat.crystals.catalog_infinity_crystals`
* :ref:`sage.combinat.crystals.catalog_elementary_crystals`
* :ref:`sage.combinat.crystals.catalog_kirillov_reshetikhin`

Functorial constructions:

* :class:`DirectSum <sage.combinat.crystals.direct_sum.DirectSumOfCrystals>`
* :class:`TensorProduct <sage.combinat.crystals.tensor_product.TensorProductOfCrystals>`
"""

from .letters import CrystalOfLetters as Letters
from .spins import CrystalOfSpins as Spins
from .spins import CrystalOfSpinsPlus as SpinsPlus
from .spins import CrystalOfSpinsMinus as SpinsMinus
from .tensor_product import CrystalOfTableaux as Tableaux
from .fast_crystals import FastCrystal as FastRankTwo
from .affine import AffineCrystalFromClassical as AffineFromClassical
from .affine import AffineCrystalFromClassicalAndPromotion as AffineFromClassicalAndPromotion
from .affine_factorization import AffineFactorizationCrystal as AffineFactorization
from .fully_commutative_stable_grothendieck import FullyCommutativeStableGrothendieckCrystal as FullyCommutativeStableGrothendieck
from sage.combinat.crystals.affinization import AffinizationOfCrystal as AffinizationOf
from .highest_weight_crystals import HighestWeightCrystal as HighestWeight
from .alcove_path import CrystalOfAlcovePaths as AlcovePaths
from .littelmann_path import CrystalOfLSPaths as LSPaths
from .littelmann_path import CrystalOfProjectedLevelZeroLSPaths as ProjectedLevelZeroLSPaths
from .kyoto_path_model import KyotoPathModel
from .generalized_young_walls import CrystalOfGeneralizedYoungWalls as GeneralizedYoungWalls
from .monomial_crystals import CrystalOfNakajimaMonomials as NakajimaMonomials
from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal as KirillovReshetikhin
from sage.combinat.rigged_configurations.rc_crystal import CrystalOfRiggedConfigurations as RiggedConfigurations
from sage.combinat.shifted_primed_tableau import ShiftedPrimedTableaux_shape as ShiftedPrimedTableaux
from sage.combinat.partition_kleshchev import KleshchevPartitions
from sage.combinat.multiset_partition_into_sets_ordered import MinimajCrystal as Minimaj

from sage.combinat.crystals.induced_structure import InducedCrystal as Induced

from sage.combinat.crystals.kac_modules import CrystalOfKacModule as KacModule
from sage.combinat.crystals.kac_modules import CrystalOfOddNegativeRoots as OddNegativeRoots

from .tensor_product import TensorProductOfCrystals as TensorProduct
from .direct_sum import DirectSumOfCrystals as DirectSum

from . import catalog_kirillov_reshetikhin as kirillov_reshetikhin
from . import catalog_infinity_crystals as infinity
from . import catalog_elementary_crystals as elementary
