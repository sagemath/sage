from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .sandpile import Sandpile, SandpileDivisor, SandpileConfig, firing_graph, parallel_firing_graph, wilmes_algorithm, random_tree, random_digraph, random_DAG, triangle_sandpile

lazy_import('sage.sandpiles.examples', 'sandpiles')

lazy_import('sage.sandpiles.sandpile', 'sandlib', deprecation=(18618,'sandlib() will soon be removed.  Use sandpile() instead.'))
lazy_import('sage.sandpiles.sandpile', 'grid_sandpile', deprecation=(18618,'grid_sandpile() will soon be removed.  Use sandpile.Grid() instead.'))
lazy_import('sage.sandpiles.sandpile', 'complete_sandpile', deprecation=(18618,'complete_sandpile() will soon be removed.  Use sandpile.Complete() instead.'))
lazy_import('sage.sandpiles.sandpile', 'firing_vector', deprecation=(18618,'firing_vector() will soon be removed.  Use SandpileDivisor.is_linearly_equivalent() instead.'))

lazy_import('sage.sandpiles.sandpile', ['admissible_partitions','partition_sandpile','min_cycles','glue_graphs','aztec_sandpile','triangle_sandpile'], deprecation=18618)
