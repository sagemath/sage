# sage_setup: distribution = sagemath-graphs
from sage.combinat.all__sagemath_categories import *

from sage.misc.lazy_import import lazy_import

# Posets
from sage.combinat.posets.all import *

# Trees and Tamari interval posets
from sage.combinat.ordered_tree import (OrderedTree, OrderedTrees,
                                        LabelledOrderedTree, LabelledOrderedTrees)
from sage.combinat.binary_tree import (BinaryTree, BinaryTrees,
                                       LabelledBinaryTree, LabelledBinaryTrees)
lazy_import('sage.combinat.rooted_tree', ('RootedTree', 'RootedTrees',
                                          'LabelledRootedTree', 'LabelledRootedTrees'))
lazy_import('sage.combinat.interval_posets', [
            'TamariIntervalPoset', 'TamariIntervalPosets'])

from sage.combinat.graph_path import GraphPaths

from sage.combinat.yang_baxter_graph import YangBaxterGraph

# block designs, etc
from sage.combinat.designs.all import *

# Cluster Algebras and Quivers
from sage.combinat.cluster_algebra_quiver.all import *

# Finite State Machines (Automaton, Transducer)
lazy_import('sage.combinat.finite_state_machine',
            ['Automaton', 'Transducer', 'FiniteStateMachine'])
lazy_import('sage.combinat.finite_state_machine_generators',
            ['automata', 'transducers'])
del lazy_import
