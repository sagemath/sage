from typing import assert_type

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.schemes.elliptic_curves.constructor import EllipticCurve, EllipticCurve_from_j
from sage.schemes.elliptic_curves.ell_field import EllipticCurve_field
from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
from sage.structure.sequence import Sequence_generic


assert_type(EllipticCurve([1, 1]), EllipticCurve_rational_field)
assert_type(EllipticCurve([0, 0, 1, -7, 6]), EllipticCurve_rational_field)
assert_type(EllipticCurve('37a1'), EllipticCurve_rational_field)
assert_type(EllipticCurve(j=1728), EllipticCurve_rational_field)
assert_type(EllipticCurve(QQ, [1, 1]), EllipticCurve_rational_field)
assert_type(EllipticCurve(GF(11), [1, 1]), EllipticCurve_finite_field)
assert_type(EllipticCurve_from_j(1728), EllipticCurve_generic)

E = EllipticCurve(GF(11), [1, 1])
P = E(0, 1)
assert_type(E, EllipticCurve_finite_field)
assert_type(P, EllipticCurvePoint_field)
assert_type(E.quadratic_twist(), EllipticCurve_field)
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

E2 = EllipticCurve([1, 1])
assert_type(E2.lift_x(0), EllipticCurvePoint_field | list[EllipticCurvePoint_field])
assert_type(E2(0).order(), Integer | PlusInfinity)
