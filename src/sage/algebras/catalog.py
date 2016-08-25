r"""
Catalog of Algebras

The ``algebras`` object may be used to access examples of various algebras
currently implemented in Sage. Using tab-completion on this object is an
easy way to discover and quickly create the algebras that are available
(as listed here).

Let ``<tab>`` indicate pressing the tab key.  So begin by typing
``algebras.<tab>`` to the see the currently implemented named algebras.

- :class:`algebras.Brauer <sage.combinat.diagram_algebras.BrauerAlgebra>`
- :class:`algebras.Clifford <sage.algebras.clifford_algebra.CliffordAlgebra>`
- :class:`algebras.Descent <sage.combinat.descent_algebra.DescentAlgebra>`
- :class:`algebras.DifferentialWeyl
  <sage.algebras.weyl_algebra.DifferentialWeylAlgebra>`
- :class:`algebras.Exterior <sage.algebras.clifford_algebra.ExteriorAlgebra>`
- :class:`algebras.FiniteDimensional
  <sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra.FiniteDimensionalAlgebra>`
- :class:`algebras.Free <sage.algebras.free_algebra.FreeAlgebraFactory>`
- :class:`algebras.FreeZinbiel <sage.algebras.free_zinbiel_algebra.FreeZinbielAlgebra>`
- :class:`algebras.PreLieAlgebra <sage.combinat.free_prelie_algebra.FreePreLieAlgebra>`
- :func:`algebras.GradedCommutative
  <sage.algebras.commutative_dga.GradedCommutativeAlgebra>`
- :class:`algebras.Group <sage.algebras.group_algebra.GroupAlgebra>`
- :class:`algebras.Hall <sage.algebras.hall_algebra.HallAlgebra>`
- :class:`algebras.Incidence <sage.combinat.posets.incidence_algebras.IncidenceAlgebra>`
- :class:`algebras.IwahoriHecke
  <sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra>`
- :class:`algebras.Moebius <sage.combinat.posets.moebius_algebra.MoebiusAlgebra>`
- :class:`algebras.Jordan
  <sage.algebras.jordan_algebra.JordanAlgebra>`
- :class:`algebras.NilCoxeter
  <sage.algebras.nil_coxeter_algebra.NilCoxeterAlgebra>`
- :class:`algebras.OrlikSolomon
  <sage.algebras.orlik_solomon.OrlikSolomonAlgebra>`
- :class:`algebras.Partition <sage.combinat.diagram_algebras.PartitionAlgebra>`
- :class:`algebras.PlanarPartition <sage.combinat.diagram_algebras.PlanarAlgebra>`
- :func:`algebras.Quaternion
  <sage.algebras.quatalg.quaternion_algebra.QuaternionAlgebraFactory>`
- :class:`algebras.Schur <sage.algebras.schur_algebra.SchurAlgebra>`
- :class:`algebras.Shuffle <sage.algebras.shuffle_algebra.ShuffleAlgebra>`
- :class:`algebras.Steenrod
  <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
- :class:`algebras.TemperleyLieb <sage.combinat.diagram_algebras.TemperleyLiebAlgebra>`
- :class:`algebras.YokonumaHecke
  <sage.algebras.yokonuma_hecke_algebra.YokonumaHeckeAlgebra>`
"""

from sage.algebras.free_algebra import FreeAlgebra as Free
from sage.algebras.iwahori_hecke_algebra import IwahoriHeckeAlgebra as IwahoriHecke
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra as Quaternion
from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra as Steenrod
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra as FiniteDimensional
from sage.algebras.group_algebra import GroupAlgebra as Group
from sage.algebras.clifford_algebra import CliffordAlgebra as Clifford
from sage.algebras.clifford_algebra import ExteriorAlgebra as Exterior
from sage.algebras.weyl_algebra import DifferentialWeylAlgebra as DifferentialWeyl

from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra', 'NilCoxeter')
lazy_import('sage.algebras.free_zinbiel_algebra', 'FreeZinbielAlgebra', 'FreeZinbiel')
lazy_import('sage.algebras.hall_algebra', 'HallAlgebra', 'Hall')
lazy_import('sage.algebras.jordan_algebra', 'JordanAlgebra', 'Jordan')
lazy_import('sage.algebras.orlik_solomon', 'OrlikSolomonAlgebra', 'OrlikSolomon')
lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra', 'Shuffle')
lazy_import('sage.algebras.schur_algebra', 'SchurAlgebra', 'Schur')
lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra', 'GradedCommutative')
lazy_import('sage.algebras.yokonuma_hecke_algebra', 'YokonumaHeckeAlgebra', 'YokonumaHecke')
lazy_import('sage.combinat.posets.incidence_algebras', 'IncidenceAlgebra', 'Incidence')
lazy_import('sage.combinat.descent_algebra', 'DescentAlgebra', 'Descent')
lazy_import('sage.combinat.diagram_algebras', 'BrauerAlgebra', 'Brauer')
lazy_import('sage.combinat.diagram_algebras', 'PartitionAlgebra', 'Partition')
lazy_import('sage.combinat.diagram_algebras', 'PlanarAlgebra', 'PlanarPartition')
lazy_import('sage.combinat.diagram_algebras', 'TemperleyLiebAlgebra', 'TemperleyLieb')
lazy_import('sage.combinat.posets.moebius_algebra', 'MoebiusAlgebra', 'Moebius')
lazy_import('sage.combinat.free_prelie_algebra', 'FreePreLieAlgebra', 'FreePreLie')
del lazy_import # We remove the object from here so it doesn't appear under tab completion

