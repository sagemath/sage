dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_graphs.py

graft sage/combinat/posets
exclude sage/combinat/posets/hasse_cython_flint.p*      # needs flint
graft sage/graphs
exclude sage/graphs/chrompoly.p*                        # needs flint
exclude sage/graphs/matchpoly.p*                        # needs flint
exclude sage/graphs/convexity_properties.p*             # cimports sage.numerical.backends.generic
include src/sage/groups/perm_gps/partn_ref/refinement_graphs.p*
# simplicial complexes
graft sage/topology            # depends on sage.combinat.subset (now in sagemath-categories)

# Could also try to add:
# sage/geometry/polyhedron/combinatorial_polyhedron


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
