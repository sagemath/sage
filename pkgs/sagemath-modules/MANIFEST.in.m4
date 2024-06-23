dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

graft sage/misc
exclude sage/misc/all.py
# exclude what's in sagemath-objects
exclude sage/misc/classcall_metaclass.*
exclude sage/misc/inherit_comparison*.*
exclude sage/misc/weak_dict.*
exclude sage/misc/nested_class*.*
exclude sage/misc/test_nested_class*.p*
exclude sage/misc/abstract_method.*
exclude sage/misc/cachefunc.*
exclude sage/misc/decorators.*
exclude sage/misc/c3_controlled.*
exclude sage/misc/lazy_attribute.*
exclude sage/misc/function_mangling.*
exclude sage/misc/lazy_string.*
exclude sage/misc/lazy_format.*
exclude sage/misc/unknown.*
exclude sage/misc/fast_methods.*
exclude sage/misc/constant_function.*
exclude sage/misc/call.*
exclude sage/misc/bindable_class.*
exclude sage/misc/namespace_package.p*
exclude sage/misc/package_dir.py
exclude sage/misc/verbose.*
exclude sage/misc/repr.*
exclude sage/misc/superseded.*
exclude sage/misc/misc_c.*
exclude sage/misc/flatten.*
exclude sage/misc/lazy_list.p*
exclude sage/misc/lazy_import*.*
exclude sage/misc/sageinspect.*
exclude sage/misc/instancedoc.*
exclude sage/misc/persist.*
exclude sage/misc/sage_unittest.*
exclude sage/misc/fpickle.p*
exclude sage/misc/randstate.*
exclude sage/misc/prandom.*
exclude sage/misc/misc.*
exclude sage/misc/timing.p*
exclude sage/misc/globals.p*
exclude sage/misc/sage_timeit*.p*
exclude sage/misc/session.p*
exclude sage/misc/reset.p*
exclude sage/misc/sage_ostools.p*

# exclude what's in sagemath-categories
exclude sage/misc/object_multiplexer.p*
exclude sage/misc/multireplace.p*
exclude sage/misc/allocator.*
exclude sage/misc/latex*.*
exclude sage/misc/html.p*
exclude sage/misc/mathml.p*
exclude sage/misc/table.p*
exclude sage/misc/map_threaded.p*
exclude sage/misc/mrange.p*
exclude sage/misc/defaults.p*
exclude sage/misc/converting_dict.p*
exclude sage/misc/parser.p*
exclude sage/misc/method_decorator.p*
exclude sage/misc/random_testing.p*
exclude sage/misc/rest_index_of_methods.p*
exclude sage/misc/callable_dict.p*
exclude sage/misc/search.p*
exclude sage/misc/derivative.p*
exclude sage/misc/functional.p*
exclude sage/misc/binary_tree.p*
exclude sage/misc/stopgap.p*

# Exclude what's in sagemath-environment
exclude sage/misc/package.py
exclude sage/misc/package_dir.py
exclude sage/misc/temporary_file.py
exclude sage/misc/viewer.py


# Exclude what's in sagemath-repl
exclude sage/misc/banner.py
exclude sage/misc/sagedoc.py
exclude sage/misc/sage_input.py
exclude sage/misc/sage_eval.py
exclude sage/misc/explain_pickle.py
exclude sage/misc/trace.py
exclude sage/misc/profiler.py
exclude sage/misc/dev_tools.py
exclude sage/misc/edit_module.py
exclude sage/misc/pager.py
exclude sage/misc/cython.py
exclude sage/misc/inline_fortran.py

# see sage.misc.all__sagemath_modules -- exclude dev tools, session management not already excluded above.
exclude sage/misc/sage_timeit_class.p*
exclude sage/misc/edit_module.p*
exclude sage/misc/remote_file.p*
exclude sage/misc/dist.p*
#??? sage/misc/sagedoc.p*
exclude sage/misc/classgraph.p*
exclude sage/misc/benchmark.p*
exclude sage/misc/citation.p*
exclude sage/misc/copying.p*
exclude sage/misc/gperftools.p*
exclude sage/misc/messaging.p*
exclude sage/misc/python.p*
exclude sage/misc/sh.p*
# lost and abandoned
prune sage/misc/notes

graft sage/modules
exclude sage/modules/module.*           # in sagemath-objects
exclude sage/modules/vector_mod2*.*     # depends on m4ri
exclude sage/modules/vector_*symbol*.*  # --> sagemath-symbolics
prune sage/modules/fp_graded

include sage/geometry/toric_lattice*.p*

# Also just modules
graft sage/groups/additive_abelian
graft sage/groups/abelian_gps
include sage/groups/galois_group.p*
exclude sage/groups/abelian_gps/abelian_group_morphism.p*       # gap
exclude sage/groups/abelian_gps/abelian_aut.p*                  # gap
exclude sage/groups/abelian_gps/*gap*.p*                        # gap
exclude sage/groups/abelian_gps/all.p*

# Need sage.combinat.free_module for polyhedral modules
include sage/combinat/free_module.py
include sage/combinat/cartesian_product.py
include sage/combinat/family.py       # until https://trac.sagemath.org/ticket/32624 is done

# root_system; could also instead go to sagemath-polyhedra (which has hyperplane arrangements)
graft sage/combinat/root_system
exclude sage/combinat/root_system/reflection_group*.p*                 # cimports PermutationGroupElement, depends on gap3
exclude sage/combinat/root_system/weyl_group*.p*  # gap

include sage/algebras/algebra.py
include sage/algebras/catalog.py
graft sage/algebras/finite_dimensional_algebras  # for hyperplane arrangements
include sage/algebras/group_algebra.py

graft sage/matrix
exclude sage/matrix/misc.p*             # Matrix_integer_sparse
exclude sage/matrix/misc_flint.p*
exclude sage/matrix/matrix_gap.*
exclude sage/matrix/matrix_*ball*.*     # depends on arb
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
include sage/groups/matrix_gps/catalog.p*
include sage/groups/matrix_gps/finitely_generated.p*
include sage/groups/matrix_gps/group_element.p*
include sage/groups/matrix_gps/linear.p*
include sage/groups/matrix_gps/matrix_group.p*
include sage/groups/matrix_gps/named_group.p*
include sage/groups/matrix_gps/orthogonal.p*
include sage/groups/matrix_gps/symplectic.p*
include sage/groups/matrix_gps/unitary.p*
include sage/groups/matrix_gps/coxeter_group.p*

include sage/groups/perm_gps/partn_ref/refinement_matrices.p*
include sage/groups/perm_gps/partn_ref/refinement_binary.p*

graft sage/tensor
graft sage/matroids            # though many doctests use graphs, finite fields
include sage/algebras/orlik_solomon.p*
include sage/algebras/orlik_terao.p*


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
include sage/rings/finite_rings/maps_finite_field.p*                     # vector space morphisms
include sage/rings/function_field/differential.p*
include sage/rings/function_field/derivations.p*                         # module elements
include sage/rings/function_field/derivations_rational.p*                # module elements
include sage/rings/function_field/differential.p*                        # module elements
include sage/rings/function_field/divisor.p*                             # module elements
include sage/rings/function_field/hermite_form_polynomial.p*             # cimports Matrix
#include sage/rings/function_field/valuation.p*                  -> sagemath-pari

include sage/rings/polynomial/laurent_polynomial_mpair.p*                # cimports Matrix
include sage/rings/polynomial/ore_*.p*
include sage/rings/polynomial/skew_*.p*

include sage/rings/ring_extension*.p*

graft sage/rings/invariants

# mpfr, gsl
graft sage/libs/mpfr
graft sage/libs/gsl
include sage/rings/real_mpfr.p*
include sage/rings/real_field.p*
include sage/rings/polynomial/*_mpfr_*.p*
include sage/rings/cc.p*
include sage/rings/complex_double.p*
include sage/rings/complex_field.p*
include sage/rings/complex_mpfr.p*
include sage/rings/complex_conversion.p*
include sage/rings/real_double_element_gsl.p*
# uses gsl
graft sage/calculus
# exclude what is included in sagemath-categories
exclude sage/calculus/functional.p*
exclude sage/calculus/functions.p*
# exclude symbolics
exclude sage/calculus/all.*
exclude sage/calculus/calculus.*
exclude sage/calculus/desolvers.*
exclude sage/calculus/predefined.*
exclude sage/calculus/tests.*
exclude sage/calculus/var.*

# More from algebras with compile-time dependencies on sagemath-modules
include sage/algebras/clifford_algebra*.p*
include sage/algebras/exterior_algebra*.p*
include sage/algebras/octonion_algebra.p*
include sage/algebras/weyl_algebra.p*
include sage/algebras/lie_algebras/lie_algebra_element.p*


graft sage/coding
prune sage/coding/codecan               # needs sage.groups

graft sage/crypto
## # parts of sage.crypto using linear algebra
## include sage/crypto/lattice.p*
## include sage/crypto/sbox*.p*            # cimports Matrix
## include sage/crypto/boolean_function.p* # dep of sbox
## graft sage/crypto/mq
## # parts of sage.crypto using StringMonoid
## # they could as well go to sagemath-combinat, but sage.crypto.__init__ is nonempty
## include sage/crypto/block_cipher/sdes.p*
## include sage/crypto/stream*.p*
## include sage/crypto/util*.p*
## graft sage/crypto/public_key

# sage.stats.distributions.discrete_gaussian is needed for sage.crypto.lwe.
# We include other parts of sage.stats for simplicity,
# although they cimports numpy directly or via sage.modules.vector_real_double_dense
graft sage/stats
graft sage/probability                                  # uses gsl

# only needed by interpreter
graft sage/libs/mpc
include sage/rings/complex_mpc.p*

# because it needs real_mpfr
include sage/numerical/gauss_legendre.p*
# uses vector
include sage/numerical/optimize.p*

global-exclude all__sagemath_*.*
global-include all__sagemath_modules.py

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist


# TODO:
# include sage/geometry/linear_expression.py
# include sage/geometry/hyperplane_arrangement/affine_subspace.py
# include sage/geometry/hyperplane_arrangement/hyperplane.py
# include sage/numerical/linear_functions.p*
# include sage/numerical/linear_tensor*.p*
