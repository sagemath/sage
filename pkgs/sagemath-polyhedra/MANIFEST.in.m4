dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-categories (via m4 include)
include(`../sagemath_categories/src/MANIFEST.in')

prune .tox
exclude *.m4
include requirements.txt

# Extra in sagemath-polyhedra:

global-include all__sagemath_polyhedra.py

include sage/rings/quotient_ring*.p*
include sage/rings/homset.p*
include sage/rings/ideal*.p*
include sage/rings/localization.p*
include sage/rings/morphism.p*

include sage/rings/abc.*
include sage/rings/integer*.*
include sage/rings/rational*.*
include sage/rings/infinity.*
include sage/arith/*.*
include sage/misc/allocator.*
include sage/misc/functional.p*
include sage/ext/mod_int.*

include sage/misc/latex*.*

graft sage/parallel

include sage/interfaces/__init__.py
include sage/interfaces/process.p*
include sage/interfaces/latte.p*
include sage/interfaces/four_ti_2.p*

include sage/sets/__init__.py
include sage/sets/set.py

include sage/rings/number_field/__init__.py
include sage/rings/number_field/number_field_base.*

include sage/rings/finite_rings/__init__.py
include sage/rings/finite_rings/element_base.*
include sage/rings/finite_rings/stdint.*
include sage/rings/finite_rings/finite_field_base.p*
include sage/rings/finite_rings/finite_field_constructor.py
include sage/rings/fast_arith.*
include sage/rings/finite_rings/integer_mod_limits.h
include sage/rings/finite_rings/integer_mod.pxd   # .pyx depends on pari

graft sage/modules
exclude sage/modules/vector_*double*.*  # depends on numpy
exclude sage/modules/vector_mod2*.*     # depends on m4ri
exclude sage/modules/vector_*symbol*.*  # --> sagemath-symbolics

# Need sage.combinat.free_module for polyhedral modules
include sage/combinat/__init__.py
include sage/combinat/quickref.py     # pulled in by __init__
include sage/combinat/tutorial.py     # pulled in by __init__
include sage/combinat/free_module.py
include sage/combinat/ranker.py
include sage/combinat/cartesian_product.py
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
exclude sage/matrix/matrix_*double*.*   # depends on numpy
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





global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
