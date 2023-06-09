dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_combinat.py

graft sage/algebras
graft sage/combinat
graft sage/monoids
graft sage/games

include sage/data_structures/stream.p*
include sage/rings/lazy_series*.p*

# included in sagemath-categories
prune sage/combinat/integer_lists
exclude sage/combinat/integer_vector.p*
exclude sage/combinat/backtrack.p*
exclude sage/combinat/combinat.p*
exclude sage/combinat/combinat_cython.p*
exclude sage/combinat/combination.p*
exclude sage/combinat/combinatorial_map.p*
exclude sage/combinat/composition.p*
exclude sage/combinat/permutation.p*
exclude sage/combinat/permutation_cython.p*
exclude sage/combinat/ranker.py                                         # for sage.sets.finite_set_map_cy
exclude sage/combinat/subset.p*
exclude sage/combinat/tools.p*
exclude sage/combinat/tuple.p*
exclude sage/combinat/subsets_hereditary.p*
exclude sage/combinat/subsets_pairwise.p*
exclude sage/combinat/dlx.p*
exclude sage/combinat/matrices/dancing_links.p*
exclude sage/combinat/matrices/dancing_links_c.h
exclude sage/combinat/matrices/dlxcpp.p*

# included in sagemath-graphs
prune sage/combinat/designs
prune sage/combinat/posets
prune sage/combinat/cluster_algebra_quiver
exclude sage/combinat/abstract_tree.p*
exclude sage/combinat/binary_tree.p*
exclude sage/combinat/ordered_tree.p*
exclude sage/combinat/rooted_tree.p*
exclude sage/combinat/graph_path.p*
exclude sage/combinat/shard_order.p*
exclude sage/combinat/tamari_lattices.p*
exclude sage/combinat/nu_tamari_lattice.p*
exclude sage/combinat/interval_posets.p*
exclude sage/combinat/yang_baxter_graph.p*
exclude sage/combinat/finite_state_machine*.p*

# included in sagemath-modules
prune sage/combinat/root_system
exclude sage/combinat/free_module.py
exclude sage/combinat/cartesian_product.py
exclude sage/combinat/family.py
exclude sage/algebras/algebra.py
prune sage/algebras/finite_dimensional_algebras
exclude sage/algebras/group_algebra.py
exclude sage/algebras/orlik_solomon.p*
exclude sage/algebras/orlik_terao.p*

exclude sage/algebras/clifford_algebra*.p*
exclude sage/algebras/exterior_algebra*.p*
exclude sage/algebras/octonion_algebra.p*
exclude sage/algebras/weyl_algebra.p*

# included in sagemath-groups
exclude sage/combinat/enumeration_mod_permgroup.p*

# compile-time library dependencies
prune sage/algebras/quatalg                                 # flint, singular
prune sage/algebras/letterplace                             # singular
prune sage/algebras/fusion_rings                            # number_field (ntl), singular
prune sage/algebras/lie_algebras                # needs modules


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
