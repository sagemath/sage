dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

global-include all__sagemath_categories.py
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
# Exclude to make it a namespace package
exclude sage/categories/__init__.py

include sage/rings/ideal.*
include sage/rings/ring.*
graft sage/typeset                     # dep of sage.categories.tensor

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
##################################### NEED TO ADD a file sage/misc/all__sagemath_categories
include sage/misc/allocator.*
include sage/misc/functional.p*
include sage/misc/latex*.*
include sage/misc/html.p*
include sage/misc/table.p*
include sage/misc/map_threaded.p*

include sage/rings/finite_rings/__init__.py
include sage/rings/finite_rings/element_base.*
include sage/rings/finite_rings/stdint.*
include sage/rings/finite_rings/finite_field_base.p*
include sage/rings/finite_rings/finite_field_constructor.py
include sage/rings/fast_arith.*
include sage/rings/finite_rings/integer_mod_limits.h
include sage/rings/finite_rings/integer_mod.pxd   # .pyx depends on pari

include sage/rings/number_field/__init__.py
include sage/rings/number_field/number_field_base.*

include sage/rings/real_double.p*

include sage/rings/fraction_field.p*
include sage/rings/fraction_field_element.p*

include sage/rings/polynomial/__init__.py
include sage/rings/polynomial/polynomial_ring.p*
include sage/rings/polynomial/polynomial_ring_constructor.p*
include sage/rings/polynomial/polynomial_singular_interface.p*
# include sage/rings/integer*.*     # depends on cypari, flint - https://trac.sagemath.org/ticket/30022
# include sage/rings/rational*.*
# include sage/rings/infinity.*

global-exclude *.c
global-exclude *.cpp

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
