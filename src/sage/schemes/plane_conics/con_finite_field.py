# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.finite_rings
r"""
Projective plane conics over finite fields

AUTHORS:

- Marco Streng (2010-07-20)
"""
# ****************************************************************************
#       Copyright (C) 2009/2010 Marco Streng <marco.streng@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.curves.projective_curve import ProjectivePlaneCurve_finite_field
from .con_field import ProjectiveConic_field


class ProjectiveConic_finite_field(ProjectiveConic_field, ProjectivePlaneCurve_finite_field):
    r"""
    Create a projective plane conic curve over a finite field.

    See ``Conic`` for full documentation.

    EXAMPLES::

        sage: K.<a> = FiniteField(9, 'a')
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - a*Z^2)
        Projective Conic Curve over Finite Field in a of size 3^2
         defined by X^2 + Y^2 + (-a)*Z^2

    ::

        sage: P.<X, Y, Z> = FiniteField(5)[]
        sage: Conic(X^2 + Y^2 - 2*Z^2)
        Projective Conic Curve over Finite Field of size 5 defined by X^2 + Y^2 - 2*Z^2
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        TESTS::

            sage: conic = Conic([GF(7)(1), 1, -1]); conic
            Projective Conic Curve over Finite Field of size 7 defined by x^2 + y^2 - z^2
            sage: conic._test_pickling()
        """
        ProjectiveConic_field.__init__(self, A, f)

    def count_points(self, n):
        r"""
        If the base field `B` of ``self`` is finite of order `q`,
        then returns the number of points over `\GF{q}, ..., \GF{q^n}`.

        EXAMPLES::

            sage: P.<x,y,z> = GF(3)[]
            sage: c = Curve(x^2+y^2+z^2); c
            Projective Conic Curve over Finite Field of size 3 defined by x^2 + y^2 + z^2
            sage: c.count_points(4)
            [4, 10, 28, 82]
        """
        F = self.base_ring()
        q = F.cardinality()
        return [q**i + 1 for i in range(1, n + 1)]

    def has_rational_point(self, point=False, read_cache=True,
                           algorithm='default'):
        r"""
        Always returns ``True`` because ``self`` has a point defined over
        its finite base field `B`.

        If ``point`` is True, then returns a second output `S`, which is a
        rational point if one exists.

        Points are cached. If ``read_cache`` is True, then cached information
        is used for the output if available. If no cached point is available
        or ``read_cache`` is False, then random `y`-coordinates are tried
        if ``self`` is smooth and a singular point is returned otherwise.

        EXAMPLES::

            sage: Conic(FiniteField(37), [1, 2, 3, 4, 5, 6]).has_rational_point()
            True

            sage: C = Conic(FiniteField(2), [1, 1, 1, 1, 1, 0]); C
            Projective Conic Curve over Finite Field of size 2
             defined by x^2 + x*y + y^2 + x*z + y*z
            sage: C.has_rational_point(point = True)  # output is random
            (True, (0 : 0 : 1))

            sage: p = next_prime(10^50)
            sage: F = FiniteField(p)
            sage: C = Conic(F, [1, 2, 3]); C
            Projective Conic Curve over Finite Field
             of size 100000000000000000000000000000000000000000000000151
             defined by x^2 + 2*y^2 + 3*z^2
            sage: C.has_rational_point(point = True)  # output is random
            (True,
             (14971942941468509742682168602989039212496867586852
              : 75235465708017792892762202088174741054630437326388 : 1)

            sage: F.<a> = FiniteField(7^20)
            sage: C = Conic([1, a, -5]); C
            Projective Conic Curve over Finite Field in a of size 7^20
             defined by x^2 + a*y^2 + 2*z^2
            sage: C.has_rational_point(point = True)  # output is random
            (True,
             (a^18 + 2*a^17 + 4*a^16 + 6*a^13 + a^12 + 6*a^11 + 3*a^10 + 4*a^9 + 2*a^8
                + 4*a^7 + a^6 + 4*a^4 + 6*a^2 + 3*a + 6
              : 5*a^19 + 5*a^18 + 5*a^17 + a^16 + 2*a^15 + 3*a^14 + 4*a^13 + 5*a^12
                + a^11 + 3*a^10 + 2*a^8 + 3*a^7 + 4*a^6 + 4*a^5 + 6*a^3 + 5*a^2 + 2*a + 4
              : 1))

        TESTS::

            sage: l = Sequence(cartesian_product_iterator([[0, 1] for i in range(6)]))
            sage: bigF = GF(next_prime(2^100))
            sage: bigF2 = GF(next_prime(2^50)^2, 'b')
            sage: m = [[F(b) for b in a]
            ....:      for a in l
            ....:      for F in [GF(2), GF(4, 'a'), GF(5), GF(9, 'a'), bigF, bigF2]]
            sage: m += [[F.random_element() for i in range(6)]
            ....:       for j in range(20) for F in [GF(5), bigF]]
            sage: c = [Conic(a) for a in m if a != [0,0,0,0,0,0]]
            sage: assert all(C.has_rational_point() for C in c)
            sage: r = randrange(0, 5)
            sage: assert all(C.defining_polynomial()(  # long time (1.4s on sage.math, 2013)
            ....:                Sequence(C.has_rational_point(point=True)[1])) == 0
            ....:            for C in c[r::5])
        """
        if not point:
            return True
        if read_cache:
            if self._rational_point is not None:
                return True, self._rational_point
        B = self.base_ring()
        s, pt = self.has_singular_point(point=True)
        if s:
            return True, pt
        while True:
            x = B.random_element()
            Y = PolynomialRing(B, 'Y').gen()
            r = self.defining_polynomial()([x, Y, 1]).roots()
            if r:
                return True, self.point([x, r[0][0], B(1)])
