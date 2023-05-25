dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_polyhedra.py

graft sage/misc
exclude sage/misc/all.py
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


include sage/interfaces/latte.p*
include sage/interfaces/four_ti_2.p*

graft sage/modules
exclude sage/modules/vector_*double*.*  # depends on numpy
exclude sage/modules/vector_numpy*.*    # depends on numpy
exclude sage/modules/vector_mod2*.*     # depends on m4ri
exclude sage/modules/vector_*symbol*.*  # --> sagemath-symbolics
prune sage/modules/fp_graded

# Also just modules
graft sage/groups/additive_abelian

# Need sage.combinat.free_module for polyhedral modules
include sage/combinat/free_module.py
include sage/combinat/ranker.py
include sage/combinat/cartesian_product.py
include sage/combinat/family.py       # until https://trac.sagemath.org/ticket/32624 is done
include sage/sets/finite_set_map*.p*
include sage/misc/callable_dict.p*

include sage/algebras/algebra.py
graft sage/algebras/finite_dimensional_algebras  # for hyperplane arrangements
include sage/algebras/group_algebra.py

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

# Can add sage/calculus/functions.p* (jacobian, wronskian) -- excluded from sagemath-categories because it needs matrices

graft sage/geometry
exclude sage/geometry/all.py
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
exclude sage/geometry/ribbon_graph.p*                    # depends on sage.groups.perm_gps
exclude sage/geometry/integral_points_integer_dense.pyx  # depends on matrix_integer_dense

graft sage/quadratic_forms
prune sage/quadratic_forms/genera       # this and below are lazy-imported and can only be tested with pari present
exclude sage/quadratic_forms/quadratic_form__automorphisms.p*
exclude sage/quadratic_forms/quadratic_form__genus.p*
exclude sage/quadratic_forms/quadratic_form__local_density_interfaces.p*
exclude sage/quadratic_forms/quadratic_form__local_normal_form.p*
exclude sage/quadratic_forms/quadratic_form__local_representation_conditions.p*
exclude sage/quadratic_forms/quadratic_form__mass*.p*
exclude sage/quadratic_forms/quadratic_form__siegel_product.p*
exclude sage/quadratic_forms/qfsolve.p*
exclude sage/quadratic_forms/special_values.p*

graft sage/groups/affine_gps
include sage/groups/matrix_gps/all.p*
include sage/groups/matrix_gps/finitely_generated.p*
include sage/groups/matrix_gps/group_element.p*
include sage/groups/matrix_gps/linear.p*
include sage/groups/matrix_gps/matrix_group.p*
include sage/groups/matrix_gps/named_group.p*
include sage/groups/matrix_gps/orthogonal.p*
include sage/groups/matrix_gps/symplectic.p*
include sage/groups/matrix_gps/unitary.p*

#graft sage/tensor              # depends on permutations (now in sagemath-categories),
                                # but also on SymmetricGroup (could be easily fixed)
                                # and sage.sets.disjoint_set

#graft sage/matroids            # many doctests use graphs, finite fields

include sage/numerical/mip.p*
include sage/numerical/interactive_simplex_method.p*
include sage/numerical/linear_*.p*
include sage/numerical/sdp.p*
include sage/numerical/backends/cvxopt_*.p*
include sage/numerical/backends/cvxpy_*.p*
include sage/numerical/backends/generic_*.p*
include sage/numerical/backends/interactivelp_*.p*
include sage/numerical/backends/logging_*.p*
include sage/numerical/backends/matrix_*.p*
include sage/numerical/backends/ppl_*.p*
include sage/numerical/backends/scip_*.p*

# just modules
graft sage/homology
# exclude stuff moved to sage/topology
exclude sage/homology/cell_*.p*
exclude sage/homology/cubical_*.p*
exclude sage/homology/delta_*.p*
exclude sage/homology/examples.p*
exclude sage/homology/simplicial_*.p*
exclude sage/homology/tests.p*

include sage/rings/derivation.p*
include sage/rings/function_field/differential.p*
include sage/rings/function_field/derivations.p*                         # module elements
include sage/rings/function_field/derivations_rational.p*                # module elements
include sage/rings/function_field/differential.p*                        # module elements
include sage/rings/function_field/divisor.p*                             # module elements
include sage/rings/function_field/hermite_form_polynomial.p*             # cimports Matrix
include sage/rings/function_field/valuation.p*                           # ??

include sage/rings/polynomial/laurent_polynomial_mpair.p*                # cimports Matrix
include sage/rings/valuation.p*


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
