dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-objects (via m4 include)
include(`../sagemath_objects/src/MANIFEST.in')

# Extra in sagemath-categories:
global-include all__sagemath_categories.py
graft sage/categories
include sage/misc/prandom.*              # dep of sage/rings/ring
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
include sage/misc/allocator.*
include sage/misc/functional.p*
include sage/misc/latex*.*
include sage/misc/html.p*
include sage/misc/table.p*
include sage/misc/map_threaded.p*
include sage/ext/mod_int.*

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
