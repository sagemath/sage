from typing import assert_type

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.infinity import PlusInfinity
from sage.schemes.elliptic_curves.constructor import EllipticCurve, EllipticCurve_from_j
from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.structure.sequence import Sequence_generic


# The factory call is unannotated upstream, so pin the receiver classes the
# same way a statically-typed caller would:
E: EllipticCurve_finite_field = EllipticCurve(GF(11), [1, 1])
P = E(0, 1)

assert_type(EllipticCurve_from_j(1728), EllipticCurve_generic)
assert_type(E.quadratic_twist(), EllipticCurve_field)
assert_type(E(0, 1), EllipticCurvePoint_field)
assert_type(E.lift_x(0), EllipticCurvePoint_field | list[EllipticCurvePoint_field])
assert_type(E.points(), Sequence_generic)
assert_type(E.count_points(), Integer | list[Integer])
assert_type(E.cardinality(), Integer)
assert_type(E.trace_of_frobenius(), Integer)
assert_type(E.is_supersingular(), bool)
assert_type(P.xy(), tuple)
assert_type(P._add_(P), EllipticCurvePoint_field)
assert_type(P.division_points(5), list)
assert_type(P.order(), Integer | PlusInfinity)

E2: EllipticCurve_field = EllipticCurve([1, 1])  # over QQ
P2 = E2.lift_x(0)
assert_type(P2, EllipticCurvePoint_field | list[EllipticCurvePoint_field])
assert_type(E2(0).order(), Integer | PlusInfinity)
