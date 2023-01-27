dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_polyhedra.py

graft sage/misc
# see sage.misc.all__sagemath_polyhedra -- exclude dev tools, session management
exclude sage/misc/dev_tools.p*
exclude sage/misc/sage_timeit_class.p*
exclude sage/misc/edit_module.p*
exclude sage/misc/session.p*
exclude sage/misc/remote_file.p*
exclude sage/misc/profiler.p*
exclude sage/misc/dist.p*
exclude sage/misc/pager.p*
#??? sage/misc/sagedoc.p*
exclude sage/misc/classgraph.p*
exclude sage/misc/reset.p*
exclude sage/misc/cython.p*
exclude sage/misc/trace.p*
exclude sage/misc/explain_pickle.p*
exclude sage/misc/inline_fortran.p*
exclude sage/misc/benchmark.p*
exclude sage/misc/citation.p*
exclude sage/misc/copying.p*
exclude sage/misc/gperftools.p*
exclude sage/misc/messaging.p*
exclude sage/misc/python.p*
exclude sage/misc/sh.p*
exclude sage/misc/viewer.p*
# lost and abandoned
prune sage/misc/notes
# belongs to sagemath-symbolics
exclude sage/misc/derivative.p*
# can hopefully be eliminated
exclude sage/misc/sage_ostools.p*
# still needed for the doctester
##exclude sage/misc/package.p*


graft sage/parallel

include sage/interfaces/process.p*
include sage/interfaces/latte.p*
include sage/interfaces/four_ti_2.p*

include sage/sets/set.py

graft sage/modules
exclude sage/modules/vector_*double*.*  # depends on numpy
exclude sage/modules/vector_numpy*.*    # depends on numpy
exclude sage/modules/vector_mod2*.*     # depends on m4ri
exclude sage/modules/vector_*symbol*.*  # --> sagemath-symbolics

# Need sage.combinat.free_module for polyhedral modules
include sage/combinat/free_module.py
include sage/combinat/ranker.py
include sage/combinat/cartesian_product.py
include sage/combinat/family.py       # until https://trac.sagemath.org/ticket/32624 is done
# could easily include all of sage/sets except disjoint_set (cimports from sage.groups.perm_gps.partn_ref.data_structures)
include sage/sets/family.p*
include sage/sets/finite_enumerated_set.py
include sage/sets/disjoint_union_enumerated_sets.py
include sage/sets/set_from_iterator.py
include sage/sets/non_negative_integers.p*
include sage/misc/lazy_list.p*
include sage/misc/mrange.p*
include sage/misc/callable_dict.p*

graft sage/matrix
exclude sage/matrix/misc.*  # until refactored
exclude sage/matrix/matrix_gap.*
exclude sage/matrix/matrix_*ball*.*     # depends on arb
exclude sage/matrix/matrix_*double_dense.*   # depends on numpy
exclude sage/matrix/matrix_numpy*.*  # depends on numpy
exclude sage/matrix/matrix_*cyclo*.*    # depends on ntl
exclude sage/matrix/matrix_*gap*.*      # depends on gap
exclude sage/matrix/matrix_gf2*.*       # depends on m4ri, m4rie
exclude sage/matrix/matrix_gfpn*.*      # depends on meataxe
exclude sage/matrix/matrix_integer_*.*  # depends on flint, pari, iml, linbox
exclude sage/matrix/matrix_mod2*.*      # depends on m4ri
exclude sage/matrix/matrix_modn*.*      # depends on linbox or flint
exclude sage/matrix/matrix_mpolynom*.*  # depends on singular
exclude sage/matrix/matrix_rational_*.* # depends on flint, pari
exclude sage/matrix/matrix_symbolic_*.* # --> sagemath-symbolics
exclude sage/matrix/change_ring.*       # depends on matrix_integer_*

graft sage/data_structures
exclude sage/data_structures/bounded_integer_sequences.*   # depends on flint

graft sage/geometry
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
exclude sage/geometry/integral_points.pyx  # depends on matrix_integer_dense


include src/sage/numerical/mip.p*
include src/sage/numerical/interactive_simplex_method.p*
include src/sage/numerical/linear_*.p*
include src/sage/numerical/sdp.p*
include src/sage/numerical/backends/cvxopt_*.p*
include src/sage/numerical/backends/cvxpy_*.p*
include src/sage/numerical/backends/generic_*.p*
include src/sage/numerical/backends/interactivelp_*.p*
include src/sage/numerical/backends/logging_*.p*
include src/sage/numerical/backends/matrix_*.p*
include src/sage/numerical/backends/ppl_*.p*
include src/sage/numerical/backends/scip_*.p*


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
