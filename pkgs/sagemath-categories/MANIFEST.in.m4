dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

graft sage/categories
# Exclude what is already shipped in sagemath-objects
exclude sage/categories/action.*
exclude sage/categories/algebra_functor.*
exclude sage/categories/basic.*
exclude sage/categories/cartesian_product.*
exclude sage/categories/category*.*
exclude sage/categories/covariant_functorial_construction.*
exclude sage/categories/facade_sets.*
exclude sage/categories/functor.*
exclude sage/categories/homset.*
exclude sage/categories/homsets.*
exclude sage/categories/map.*
exclude sage/categories/morphism.*
exclude sage/categories/isomorphic_objects.*
exclude sage/categories/objects.*
exclude sage/categories/primer.*
exclude sage/categories/pushout.*
exclude sage/categories/quotients.*
exclude sage/categories/realizations.*
exclude sage/categories/sets_cat.*
exclude sage/categories/sets_with_partial_maps.*
exclude sage/categories/subobjects.*
exclude sage/categories/subquotients.*
exclude sage/categories/with_realizations.*

include sage/geometry/abc.p*

# Interfaces
include sage/interfaces/all.p*
include sage/interfaces/abc.p*
include sage/interfaces/process.p*                      # needed for sage.parallel
include sage/interfaces/tab_completion.p*
include sage/misc/object_multiplexer.p*
include sage/misc/multireplace.p*
include sage/interfaces/sagespawn.p*
include sage/interfaces/quit.p*
include sage/interfaces/cleaner.p*
include sage/interfaces/expect.p*
include sage/interfaces/interface.p*
include sage/interfaces/sage0.p*

graft sage/parallel

graft sage/typeset                     # dep of sage.categories.tensor

include sage/groups/generic.p*
include sage/groups/groups_catalog.p*

include sage/monoids/monoid.p*

include sage/rings/ring.*
include sage/rings/quotient_ring*.p*
include sage/rings/homset.p*
include sage/rings/ideal*.p*
include sage/rings/noncommutative_ideals.p*
include sage/rings/localization.p*
include sage/rings/morphism.p*

include sage/rings/abc.*
include sage/rings/integer*.p*
exclude sage/rings/integer_fake.pxd     # in sagemath-objects
include sage/rings/rational*.*
include sage/rings/infinity.*
include sage/rings/factorint.p*
include sage/rings/sum_of_squares.p*
include sage/rings/generic.p*

include sage/misc/allocator.*
include sage/misc/latex*.*
include sage/misc/html.p*
include sage/misc/mathml.p*
include sage/misc/table.p*
include sage/misc/map_threaded.p*
include sage/misc/mrange.p*
include sage/misc/defaults.p*
include sage/misc/converting_dict.p*
include sage/misc/parser.p*
include sage/misc/method_decorator.p*
include sage/misc/random_testing.p*
include sage/misc/rest_index_of_methods.p*
include sage/misc/callable_dict.p*
include sage/misc/search.p*
include sage/misc/stopgap.p*

## Data structures
include sage/misc/binary_tree.p*
graft sage/data_structures               # bitset needed by sage.graphs and sage.geometry.polyhedron
exclude sage/data_structures/bounded_integer_sequences.*   # depends on flint
exclude sage/data_structures/stream.*
graft sage/groups/perm_gps/partn_ref  # but not partn_ref2, which depends on GAP
exclude sage/groups/perm_gps/partn_ref/refinement_graphs.p*    # sagemath-graphs
exclude sage/groups/perm_gps/partn_ref/refinement_matrices.p*  # sagemath-modules
exclude sage/groups/perm_gps/partn_ref/refinement_binary.p*    # sagemath-modules


# These might later go to a separate distribution sagemath-functions (> sagemath-objects);
# but sage.functions currently depends on basic rings (QQ etc)
graft sage/arith
# Exclude what is included in sagemath-objects already
exclude sage/arith/long.p*
exclude sage/arith/numerical_approx.p*
exclude sage/arith/power.p*

include sage/calculus/functional.p*
include sage/calculus/functions.p*
include sage/misc/derivative.p*
include sage/misc/functional.p*
include sage/symbolic/symbols.p*
include sage/symbolic/function.p*
graft sage/functions


include sage/rings/finite_rings/element_base.*
include sage/rings/finite_rings/stdint.*
include sage/rings/finite_rings/finite_field_base.p*
include sage/rings/finite_rings/finite_field_constructor.p*
include sage/rings/finite_rings/finite_field_prime_modn.p*
include sage/rings/finite_rings/residue_field.p*
include sage/rings/finite_rings/hom_finite_field.p*
include sage/rings/finite_rings/hom_prime_finite_field.p*
include sage/rings/finite_rings/homset.p*
include sage/rings/fast_arith.*
include sage/rings/finite_rings/integer_mod_limits.h
include sage/rings/finite_rings/integer_mod.p*
include sage/rings/finite_rings/integer_mod_ring.p*
include sage/rings/finite_rings/conway_polynomials.p*
include sage/rings/finite_rings/galois_group.p*
include sage/rings/algebraic_closure_finite_field.p*

include sage/rings/number_field/number_field_base.p*
include sage/rings/number_field/number_field_element_base.p*
include sage/rings/number_field/number_field_ideal.p*
include sage/rings/real_double.p*

include sage/rings/real_lazy.p*

include sage/rings/fraction_field.p*
include sage/rings/fraction_field_element.p*
include sage/rings/qqbar_decorators.p*

include sage/rings/padics/padic_generic.p*
include sage/rings/padics/local_generic.p*
include sage/rings/padics/local_generic_element.p*
include sage/rings/padics/precision_error.p*
include sage/rings/padics/misc.p*

include sage/rings/polynomial/polynomial_ring.p*
include sage/rings/polynomial/polynomial_ring_constructor.p*
include sage/rings/polynomial/polynomial_ring_homomorphism.p*
include sage/rings/polynomial/polynomial_quotient_ring.p*
include sage/rings/polynomial/polynomial_singular_interface.p*
include sage/rings/polynomial/multi_polynomial_ring.p*
include sage/rings/polynomial/multi_polynomial_ring_base.p*
include sage/rings/polynomial/multi_polynomial_sequence.p*
include sage/rings/polynomial/multi_polynomial_ideal.p*
include sage/rings/polynomial/infinite_polynomial_*.p*
include sage/rings/polynomial/commutative_polynomial.p*
include sage/rings/polynomial/polynomial_compiled.p*
include sage/rings/polynomial/polynomial_element.p*
include sage/rings/polynomial/polynomial_element_generic.p*
include sage/rings/polynomial/polynomial_fateman.p*
include sage/rings/polynomial/polynomial_quotient_ring_element.p*
include sage/rings/polynomial/multi_polynomial.p*
include sage/rings/polynomial/multi_polynomial_element.p*
include sage/rings/polynomial/polydict.p*
include sage/rings/polynomial/term_order.p*
include sage/rings/polynomial/flatten.p*
include sage/rings/polynomial/cyclotomic.p*
include sage/rings/polynomial/laurent_polynomial_ring*.p*
include sage/rings/polynomial/laurent_polynomial_ideal.p*
include sage/rings/polynomial/laurent_polynomial.p*                     # but not laurent_polynomial_mpair, which needs Matrix
include sage/rings/polynomial/ideal.p*
include sage/rings/polynomial/toy*.p*
include sage/rings/polynomial/symmetric_*.p*
## include sage/rings/polynomial/ore*.p*                                   # need for sage.categories.drinfeld_modules

include sage/rings/continued_fraction*.p*

graft sage/rings/function_field
exclude sage/rings/function_field/*_polymod.*                           # needs Singular
exclude sage/rings/function_field/derivations*.*                        # module elements
exclude sage/rings/function_field/differential.*                        # module elements
exclude sage/rings/function_field/divisor.*                             # module elements
exclude sage/rings/function_field/hermite_form_polynomial.*             # cimports Matrix
exclude sage/rings/function_field/valuation*.*                          # ??
exclude sage/rings/function_field/khuri_makdisi.p*
prune sage/rings/function_field/drinfeld_modules                        # needs ore_polynomial etc.

include sage/rings/power_series_mpoly.p*
include sage/rings/power_series_poly.p*
include sage/rings/power_series_ring_element.p*
include sage/rings/power_series_ring.p*
include sage/rings/multi_power_series_ring.py
include sage/rings/multi_power_series_ring_element.py
include sage/rings/laurent_series_ring*.p*
include sage/rings/puiseux_series_ring*.p*

graft sage/rings/semirings

include sage/rings/tests.p*
include sage/rings/big_oh.p*

include sage/ext/fast_*.p*

include sage/combinat/integer_vector.p*
graft sage/combinat/integer_lists
include sage/combinat/backtrack.p*
include sage/combinat/combinat.p*
include sage/combinat/combinat_cython.p*
include sage/combinat/combination.p*
include sage/combinat/combinatorial_map.p*
include sage/combinat/composition.p*
include sage/combinat/permutation.p*
include sage/combinat/permutation_cython.p*
include sage/combinat/ranker.py                                         # for sage.sets.finite_set_map_cy
include sage/combinat/subset.p*
include sage/combinat/tools.p*
include sage/combinat/tuple.p*
# leave out partition - has complicated deps
# leave out integer_vector_weighted - needs combinat.words.word
include sage/combinat/subsets_hereditary.p*
include sage/combinat/subsets_pairwise.p*
include sage/combinat/dlx.p*
include sage/combinat/matrices/dancing_links.p*
include sage/combinat/matrices/dancing_links_c.h
include sage/combinat/matrices/dlxcpp.p*

# see src/sage/schemes/generic/notes/imports.txt
graft sage/schemes/generic
graft sage/schemes/affine
graft sage/schemes/projective
graft sage/schemes/product_projective

graft sage/sets
exclude sage/sets/pythonclass.*                                 # sagemath-objects

include sage/numerical/backends/generic*backend.p*

## These two should probably go to sagemath-combinat instead.
# include sage/data_structures/stream.p*
# include sage/rings/lazy_series*.p*

global-exclude *.c
global-exclude *.cpp

global-exclude all__sagemath_*.*
global-include all__sagemath_categories.py

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
