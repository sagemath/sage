dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

exclude *.m4
include requirements.txt

graft sage/combinat/posets
exclude sage/combinat/posets/hasse_cython_flint.p*      # needs flint
include sage/combinat/abstract_tree.p*
include sage/combinat/binary_tree.p*
include sage/combinat/ordered_tree.p*
include sage/combinat/rooted_tree.p*
include sage/combinat/graph_path.p*
include sage/combinat/shard_order.p*
include sage/combinat/tamari_lattices.p*
include sage/combinat/nu_tamari_lattice.p*
include sage/combinat/interval_posets.p*                # check if 'element in DyckWords()' can be tested better
include sage/combinat/yang_baxter_graph.p*
include sage/combinat/cluster_algebra_quiver/mutation_class.p*  # more from there?

graft sage/combinat/designs
graft sage/combinat/cluster_algebra_quiver

#include sage/combinat/root_system/dynkin_diagram.p*     # want?
#include sage/combinat/root_system/cartan_type.p*        # dep of dynkin_diagram

include sage/combinat/finite_state_machine*.p*

#include src/sage/combinat/rigged_configurations/kleber_tree.py # want?

graft sage/graphs
exclude sage/graphs/bliss.p*
prune sage/graphs/bliss_cpp
exclude sage/graphs/mcqd.p*
exclude sage/graphs/graph_decompositions/*tdlib*.*

exclude sage/graphs/chrompoly.p*                        # needs flint
exclude sage/graphs/matchpoly.p*                        # needs flint

include sage/groups/perm_gps/partn_ref/refinement_graphs.p*

graft sage/knots

graft sage/sandpiles

# quivers use bounded_integer_sequences, which depends on flint (for a bad reason)
#include sage/quivers/paths.p*
#include sage/quivers/path_semigroup.p*                  # but not their representations, algebras

# simplicial complexes
graft sage/topology            # depends on sage.combinat.subset (now in sagemath-categories)

# Could also try to add:
# sage/geometry/polyhedron/combinatorial_polyhedron

include sage/databases/knotinfo_db.p*


global-exclude all__sagemath_*.py
global-include all__sagemath_graphs.py

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
