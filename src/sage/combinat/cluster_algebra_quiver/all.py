r"""
Cluster algebras and quivers

- A compendium on the cluster algebra and quiver package in Sage [MS2011]_

- :ref:`sage.combinat.cluster_algebra_quiver.quiver_mutation_type`
- :ref:`sage.combinat.cluster_algebra_quiver.quiver`
- :ref:`sage.combinat.cluster_algebra_quiver.cluster_seed`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.misc.lazy_import import lazy_import
lazy_import("sage.combinat.cluster_algebra_quiver.quiver_mutation_type",
            "QuiverMutationType")
lazy_import("sage.combinat.cluster_algebra_quiver.quiver", "ClusterQuiver")
lazy_import("sage.combinat.cluster_algebra_quiver.cluster_seed", "ClusterSeed")

del install_doc
del lazy_import
