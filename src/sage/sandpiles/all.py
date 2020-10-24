from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

from .sandpile import (Sandpile,
                       SandpileDivisor,
                       SandpileConfig,
                       random_DAG,
                       firing_graph,
                       parallel_firing_graph,
                       wilmes_algorithm,
                       triangle_sandpile)

lazy_import('sage.sandpiles.examples', 'sandpiles')
