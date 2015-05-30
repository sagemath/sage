"""
Cluster algebra and quivers features that are imported by default in the interpreter namespace
"""
from sage.misc.lazy_import import lazy_import
lazy_import("sage.combinat.cluster_algebra_quiver.quiver_mutation_type", "QuiverMutationType")
lazy_import("sage.combinat.cluster_algebra_quiver.quiver", "ClusterQuiver")
lazy_import("sage.combinat.cluster_algebra_quiver.cluster_seed", "ClusterSeed")
