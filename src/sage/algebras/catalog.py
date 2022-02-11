r"""
Catalog of Algebras

The ``algebras`` object may be used to access examples of various algebras
currently implemented in Sage. Using tab-completion on this object is an
easy way to discover and quickly create the algebras that are available
(as listed here).

Let ``<tab>`` indicate pressing the tab key.  So begin by typing
``algebras.<tab>`` to the see the currently implemented named algebras.

- :class:`algebras.AlternatingCentralExtensionQuantumOnsager
  <sage.algebras.quantum_groups.alt_central_ext_quantum_onsager.ACEQuantumOnsagerAlgebra>`
- :class:`algebras.ArikiKoike
  <sage.algebras.hecke_algebras.ariki_koike_algebra.ArikiKoikeAlgebra>`
- :class:`algebras.AskeyWilson <sage.algebras.askey_wilson.AskeyWilsonAlgebra>`
- :class:`algebras.Blob <sage.combinat.blob_algebra.BlobAlgebra>`
- :class:`algebras.Brauer <sage.combinat.diagram_algebras.BrauerAlgebra>`
- :class:`algebras.Clifford <sage.algebras.clifford_algebra.CliffordAlgebra>`
- :class:`algebras.ClusterAlgebra <sage.algebras.cluster_algebra.ClusterAlgebra>`
- :class:`algebras.Descent <sage.combinat.descent_algebra.DescentAlgebra>`
- :class:`algebras.DifferentialWeyl
  <sage.algebras.weyl_algebra.DifferentialWeylAlgebra>`
- :class:`algebras.Exterior <sage.algebras.clifford_algebra.ExteriorAlgebra>`
- :class:`algebras.FiniteDimensional
  <sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra.FiniteDimensionalAlgebra>`
- :class:`algebras.FQSym <sage.combinat.fqsym.FreeQuasisymmetricFunctions>`
- :class:`algebras.Free <sage.algebras.free_algebra.FreeAlgebraFactory>`
- :class:`algebras.FreeZinbiel <sage.algebras.free_zinbiel_algebra.FreeZinbielAlgebra>`
- :class:`algebras.FreePreLie <sage.combinat.free_prelie_algebra.FreePreLieAlgebra>`
- :class:`algebras.FreeDendriform <sage.combinat.free_dendriform_algebra.FreeDendriformAlgebra>`
- :class:`algebras.FSym <sage.combinat.chas.fsym.FreeSymmetricFunctions>`
- :func:`algebras.GradedCommutative
  <sage.algebras.commutative_dga.GradedCommutativeAlgebra>`
- :class:`algebras.Group <sage.algebras.group_algebra.GroupAlgebra>`
- :class:`algebras.GrossmanLarson <sage.combinat.grossman_larson_algebras.GrossmanLarsonAlgebra>`
- :class:`algebras.Hall <sage.algebras.hall_algebra.HallAlgebra>`
- :class:`algebras.Incidence <sage.combinat.posets.incidence_algebras.IncidenceAlgebra>`
- :class:`algebras.IwahoriHecke
  <sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra>`
- :class:`algebras.Moebius <sage.combinat.posets.moebius_algebra.MoebiusAlgebra>`
- :class:`algebras.Jordan
  <sage.algebras.jordan_algebra.JordanAlgebra>`
- :class:`algebras.Lie <sage.algebras.lie_algebras.lie_algebra.LieAlgebra>`
- :class:`algebras.MalvenutoReutenauer <sage.combinat.fqsym.FreeQuasisymmetricFunctions>`
- :class:`algebras.NilCoxeter
  <sage.algebras.nil_coxeter_algebra.NilCoxeterAlgebra>`
- :class:`algebras.OrlikTerao
  <sage.algebras.orlik_terao.OrlikTeraoAlgebra>`
- :class:`algebras.OrlikSolomon
  <sage.algebras.orlik_solomon.OrlikSolomonAlgebra>`
- :class:`algebras.QuantumClifford
  <sage.algebras.quantum_clifford.QuantumCliffordAlgebra'>`
- :class:`algebras.QuantumGL
  <sage.algebras.quantum_matrix_coordinate_algebra.QuantumGL>`
- :class:`algebras.QuantumMatrixCoordinate
  <sage.algebras.quantum_matrix_coordinate_algebra.QuantumMatrixCoordinateAlgebra>`
- :class:`algebras.QSym <sage.combinat.ncsf_qsym.qsym.QuasiSymmetricFunctions>`
- :class:`algebras.Partition <sage.combinat.diagram_algebras.PartitionAlgebra>`
- :class:`algebras.PlanarPartition <sage.combinat.diagram_algebras.PlanarAlgebra>`
- :class:`algebras.QuantumGroup
  <sage.algebras.quantum_groups.quantum_group_gap.QuantumGroup>`
- :func:`algebras.Quaternion
  <sage.algebras.quatalg.quaternion_algebra.QuaternionAlgebraFactory>`
- :class:`algebras.RationalCherednik
  <sage.algebras.rational_cherednik_algebra.RationalCherednikAlgebra>`
- :class:`algebras.Schur <sage.algebras.schur_algebra.SchurAlgebra>`
- :class:`algebras.Shuffle <sage.algebras.shuffle_algebra.ShuffleAlgebra>`
- :class:`algebras.Steenrod
  <sage.algebras.steenrod.steenrod_algebra.SteenrodAlgebra>`
- :class:`algebras.TemperleyLieb <sage.combinat.diagram_algebras.TemperleyLiebAlgebra>`
- :class:`algebras.WQSym <sage.combinat.chas.wqsym.WordQuasiSymmetricFunctions>`
- :class:`algebras.Yangian <sage.algebras.yangian.Yangian>`
- :class:`algebras.YokonumaHecke
  <sage.algebras.yokonuma_hecke_algebra.YokonumaHeckeAlgebra>`
- :class:`algebras.Tensor <sage.algebras.tensor_algebra.TensorAlgebra>`
"""

from sage.algebras.free_algebra import FreeAlgebra as Free
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra as Quaternion
from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra as Steenrod
from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra as FiniteDimensional
from sage.algebras.group_algebra import GroupAlgebra as Group
from sage.algebras.clifford_algebra import CliffordAlgebra as Clifford
from sage.algebras.clifford_algebra import ExteriorAlgebra as Exterior
from sage.algebras.weyl_algebra import DifferentialWeylAlgebra as DifferentialWeyl
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra as Lie

from sage.misc.lazy_import import lazy_import
lazy_import('sage.algebras.iwahori_hecke_algebra', 'IwahoriHeckeAlgebra', 'IwahoriHecke')
lazy_import('sage.algebras.nil_coxeter_algebra', 'NilCoxeterAlgebra', 'NilCoxeter')
lazy_import('sage.algebras.free_zinbiel_algebra', 'FreeZinbielAlgebra', 'FreeZinbiel')
lazy_import('sage.algebras.askey_wilson', 'AskeyWilsonAlgebra', 'AskeyWilson')
lazy_import('sage.algebras.hall_algebra', 'HallAlgebra', 'Hall')
lazy_import('sage.algebras.jordan_algebra', 'JordanAlgebra', 'Jordan')
lazy_import('sage.algebras.orlik_solomon', 'OrlikSolomonAlgebra', 'OrlikSolomon')
lazy_import('sage.algebras.orlik_terao', 'OrlikTeraoAlgebra', 'OrlikTerao')
lazy_import('sage.algebras.shuffle_algebra', 'ShuffleAlgebra', 'Shuffle')
lazy_import('sage.algebras.schur_algebra', 'SchurAlgebra', 'Schur')
lazy_import('sage.algebras.commutative_dga', 'GradedCommutativeAlgebra', 'GradedCommutative')
lazy_import('sage.algebras.hecke_algebras.ariki_koike_algebra', 'ArikiKoikeAlgebra', 'ArikiKoike')
lazy_import('sage.algebras.rational_cherednik_algebra', 'RationalCherednikAlgebra', 'RationalCherednik')
lazy_import('sage.algebras.yokonuma_hecke_algebra', 'YokonumaHeckeAlgebra', 'YokonumaHecke')
lazy_import('sage.combinat.posets.incidence_algebras', 'IncidenceAlgebra', 'Incidence')
lazy_import('sage.combinat.descent_algebra', 'DescentAlgebra', 'Descent')
lazy_import('sage.combinat.diagram_algebras', 'BrauerAlgebra', 'Brauer')
lazy_import('sage.combinat.diagram_algebras', 'PartitionAlgebra', 'Partition')
lazy_import('sage.combinat.diagram_algebras', 'PlanarAlgebra', 'PlanarPartition')
lazy_import('sage.combinat.diagram_algebras', 'TemperleyLiebAlgebra', 'TemperleyLieb')
lazy_import('sage.combinat.blob_algebra', 'BlobAlgebra', 'Blob')
lazy_import('sage.combinat.posets.moebius_algebra', 'MoebiusAlgebra', 'Moebius')
lazy_import('sage.combinat.free_prelie_algebra', 'FreePreLieAlgebra', 'FreePreLie')
lazy_import('sage.combinat.free_dendriform_algebra', 'FreeDendriformAlgebra', 'FreeDendriform')
lazy_import('sage.combinat.fqsym', 'FreeQuasisymmetricFunctions', 'FQSym')
lazy_import('sage.combinat.fqsym', 'FreeQuasisymmetricFunctions', 'MalvenutoReutenauer')
lazy_import('sage.combinat.chas.wqsym', 'WordQuasiSymmetricFunctions', 'WQSym')
lazy_import('sage.combinat.chas.fsym', 'FreeSymmetricFunctions', 'FSym')
lazy_import('sage.combinat.ncsf_qsym.qsym', 'QuasiSymmetricFunctions', 'QSym')
lazy_import('sage.combinat.grossman_larson_algebras', 'GrossmanLarsonAlgebra', 'GrossmanLarson')
lazy_import('sage.algebras.quantum_clifford', 'QuantumCliffordAlgebra', 'QuantumClifford')
lazy_import('sage.algebras.quantum_matrix_coordinate_algebra',
            'QuantumMatrixCoordinateAlgebra', 'QuantumMatrixCoordinate')
lazy_import('sage.algebras.quantum_matrix_coordinate_algebra', 'QuantumGL')
lazy_import('sage.algebras.tensor_algebra', 'TensorAlgebra', 'Tensor')
lazy_import('sage.algebras.quantum_groups.quantum_group_gap', 'QuantumGroup')
lazy_import('sage.algebras.quantum_groups.ace_quantum_onsager',
           'ACEQuantumOnsagerAlgebra', 'AlternatingCentralExtensionQuantumOnsager')
lazy_import('sage.algebras.yangian', 'Yangian')
del lazy_import # We remove the object from here so it doesn't appear under tab completion

