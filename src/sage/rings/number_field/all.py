
from sage.rings.number_field.number_field import (NumberField, NumberFieldTower, CyclotomicField, QuadraticField,
                           is_fundamental_discriminant, is_real_place)
from sage.rings.number_field.number_field_element import NumberFieldElement

from sage.rings.number_field.order import EquationOrder, GaussianIntegers, EisensteinIntegers

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import('sage.rings.number_field.totallyreal', 'enumerate_totallyreal_fields_prim')
_lazy_import('sage.rings.number_field.totallyreal_data', 'hermite_constant')
_lazy_import('sage.rings.number_field.totallyreal_rel', 'enumerate_totallyreal_fields_all')
_lazy_import('sage.rings.number_field.totallyreal_rel', 'enumerate_totallyreal_fields_rel')

from sage.rings.number_field.unit_group import UnitGroup
