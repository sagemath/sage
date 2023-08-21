"""
Elliptic curve constructor

AUTHORS:

- William Stein (2005): Initial version

- John Cremona (2008-01): EllipticCurve(j) fixed for all cases
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.rational_field import RationalField

import sage.rings.abc
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.ring import is_Ring

from sage.categories.fields import Fields
_Fields = Fields()

from sage.structure.sequence import Sequence
from sage.structure.element import parent, Expression
from sage.structure.factory import UniqueFactory


class EllipticCurveFactory(UniqueFactory):
    r"""
    Construct an elliptic curve.

    In Sage, an elliptic curve is always specified by
    (the coefficients of) a long Weierstrass equation

    .. MATH::

        y^2 + a_1 xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6.

    INPUT:

    There are several ways to construct an elliptic curve:

    - ``EllipticCurve([a1,a2,a3,a4,a6])``: Elliptic curve with given
      `a`-invariants. The invariants are coerced into a common parent.
      If all are integers, they are coerced into the rational numbers.

    - ``EllipticCurve([a4,a6])``: Same as above, but `a_1=a_2=a_3=0`.

    - ``EllipticCurve(label)``: Returns the elliptic curve over `\QQ`
      from the Cremona database with the given label. The label is a
      string, such as ``"11a"`` or ``"37b2"``. The letters in the
      label *must* be lower case (Cremona's new labeling).

    - ``EllipticCurve(R, [a1,a2,a3,a4,a6])``: Create the elliptic
      curve over `R` with given `a`-invariants. Here `R` can be an
      arbitrary commutative ring, although most functionality is only
      implemented over fields.

    - ``EllipticCurve(j=j0)`` or ``EllipticCurve_from_j(j0)``: Return
      an elliptic curve with `j`-invariant ``j0``.

    - ``EllipticCurve(polynomial)``: Read off the `a`-invariants from
      the polynomial coefficients, see
      :func:`EllipticCurve_from_Weierstrass_polynomial`.

    - ``EllipticCurve(cubic, point)``: The elliptic curve defined by a
      plane cubic (homogeneous polynomial in three variables), with a
      rational point.

    Instead of giving the coefficients as a *list* of length 2 or 5,
    one can also give a *tuple*.

    EXAMPLES:

    We illustrate creating elliptic curves::

        sage: EllipticCurve([0,0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    We create a curve from a Cremona label::

        sage: EllipticCurve('37b2')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational Field
        sage: EllipticCurve('5077a')
        Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        sage: EllipticCurve('389a')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    Old Cremona labels are allowed::

        sage: EllipticCurve('2400FF')
        Elliptic Curve defined by y^2 = x^3 + x^2 + 2*x + 8 over Rational Field

    Unicode labels are allowed::

        sage: EllipticCurve(u'389a')
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    We create curves over a finite field as follows::

        sage: EllipticCurve([GF(5)(0),0,1,-1,0])                                        # optional - sage.rings.finite_rings
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5
        sage: EllipticCurve(GF(5), [0, 0,1,-1,0])                                       # optional - sage.rings.finite_rings
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    Elliptic curves over `\ZZ/N\ZZ` with `N` prime are of type
    "elliptic curve over a finite field"::

        sage: F = Zmod(101)
        sage: EllipticCurve(F, [2, 3])
        Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Ring of integers modulo 101
        sage: E = EllipticCurve([F(2), F(3)])
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field_with_category'>
        sage: E.category()
        Category of schemes over Ring of integers modulo 101

    In contrast, elliptic curves over `\ZZ/N\ZZ` with `N` composite
    are of type "generic elliptic curve"::

        sage: F = Zmod(95)
        sage: EllipticCurve(F, [2, 3])
        Elliptic Curve defined by y^2 = x^3 + 2*x + 3 over Ring of integers modulo 95
        sage: E = EllipticCurve([F(2), F(3)])
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_generic.EllipticCurve_generic_with_category'>
        sage: E.category()
        Category of schemes over Ring of integers modulo 95

    The following is a curve over the complex numbers::

        sage: E = EllipticCurve(CC, [0,0,1,-1,0])
        sage: E
        Elliptic Curve defined by y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x
         over Complex Field with 53 bits of precision
        sage: E.j_invariant()
        2988.97297297297

    We can also create elliptic curves by giving the Weierstrass equation::

        sage: R2.<x,y> = PolynomialRing(QQ,2)
        sage: EllipticCurve(y^2 + y - ( x^3 + x - 9 ))
        Elliptic Curve defined by y^2 + y = x^3 + x - 9 over Rational Field

        sage: R.<x,y> = GF(5)[]                                                         # optional - sage.rings.finite_rings
        sage: EllipticCurve(x^3 + x^2 + 2 - y^2 - y*x)                                  # optional - sage.rings.finite_rings
        Elliptic Curve defined by y^2 + x*y  = x^3 + x^2 + 2 over Finite Field of size 5

    We can also create elliptic curves by giving a smooth plane cubic with a rational point::

        sage: R3.<x,y,z> = PolynomialRing(QQ,3)
        sage: F = x^3 + y^3 + 30*z^3
        sage: P = [1,-1,0]
        sage: EllipticCurve(F,P)
        Elliptic Curve defined by y^2 - 270*y = x^3 - 24300 over Rational Field

    We can explicitly specify the `j`-invariant::

        sage: E = EllipticCurve(j=1728); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        1728
        '32a2'

        sage: E = EllipticCurve(j=GF(5)(2)); E; E.j_invariant()                         # optional - sage.rings.finite_rings
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
        2

    See :trac:`6657` ::

        sage: EllipticCurve(GF(144169), j=1728)                                         # optional - sage.rings.finite_rings
        Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 144169

    Elliptic curves over the same ring with the same Weierstrass
    coefficients are identical, even when they are constructed in
    different ways (see :trac:`11474`)::

        sage: EllipticCurve('11a3') is EllipticCurve(QQ, [0, -1, 1, 0, 0])
        True

    By default, when a rational value of `j` is given, the constructed
    curve is a minimal twist (minimal conductor for curves with that
    `j`-invariant).  This can be changed by setting the optional
    parameter ``minimal_twist``, which is True by default, to False::

        sage: EllipticCurve(j=100)
        Elliptic Curve defined by y^2 = x^3 + x^2 + 3392*x + 307888 over Rational Field
        sage: E =EllipticCurve(j=100); E
        Elliptic Curve defined by y^2 = x^3 + x^2 + 3392*x + 307888 over Rational Field
        sage: E.conductor()
        33129800
        sage: E.j_invariant()
        100
        sage: E =EllipticCurve(j=100, minimal_twist=False); E
        Elliptic Curve defined by y^2 = x^3 + 488400*x - 530076800 over Rational Field
        sage: E.conductor()
        298168200
        sage: E.j_invariant()
        100

    Without this option, constructing the curve could take a long time
    since both `j` and `j-1728` have to be factored to compute the
    minimal twist (see :trac:`13100`)::

       sage: E = EllipticCurve_from_j(2^256+1, minimal_twist=False)
       sage: E.j_invariant() == 2^256+1
       True

    TESTS::

        sage: R = ZZ['u', 'v']
        sage: EllipticCurve(R, [1,1])
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Multivariate Polynomial Ring in u, v
        over Integer Ring

    We create a curve and a point over ``QQbar`` (see :trac:`6879`)::

        sage: E = EllipticCurve(QQbar, [0,1])                                           # optional - sage.rings.number_field
        sage: E(0)                                                                      # optional - sage.rings.number_field
        (0 : 1 : 0)
        sage: E.base_field()                                                            # optional - sage.rings.number_field
        Algebraic Field

        sage: E = EllipticCurve(RR, [1,2]); E; E.base_field()
        Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 2.00000000000000
         over Real Field with 53 bits of precision
        Real Field with 53 bits of precision
        sage: E = EllipticCurve(CC,[3,4]); E; E.base_field()
        Elliptic Curve defined by y^2 = x^3 + 3.00000000000000*x + 4.00000000000000
         over Complex Field with 53 bits of precision
        Complex Field with 53 bits of precision
        sage: E = EllipticCurve(QQbar, [5,6]); E; E.base_field()                        # optional - sage.rings.number_field
        Elliptic Curve defined by y^2 = x^3 + 5*x + 6 over Algebraic Field
        Algebraic Field

    See :trac:`6657` ::

        sage: EllipticCurve(3, j=1728)
        Traceback (most recent call last):
        ...
        ValueError: First parameter (if present) must be a ring when j is specified

        sage: EllipticCurve(GF(5), j=3/5)                                               # optional - sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: First parameter must be a ring containing 3/5

    If the universe of the coefficients is a general field, the object
    constructed has type :class:`EllipticCurve_field`.  Otherwise it is
    :class:`EllipticCurve_generic`.  See :trac:`9816` ::

        sage: E = EllipticCurve([QQbar(1), 3]); E                                       # optional - sage.rings.number_field
        Elliptic Curve defined by y^2 = x^3 + x + 3 over Algebraic Field
        sage: type(E)                                                                   # optional - sage.rings.number_field
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>

        sage: E = EllipticCurve([RR(1), 3]); E
        Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 3.00000000000000
        over Real Field with 53 bits of precision
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>

        sage: E = EllipticCurve([SR(i),i]); E                                           # optional - sage.symbolic
        Elliptic Curve defined by y^2 = x^3 + I*x + I over Symbolic Ring
        sage: type(E)                                                                   # optional - sage.symbolic
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>
        sage: E.category()                                                              # optional - sage.symbolic
        Category of schemes over Symbolic Ring
        sage: SR in Fields()                                                            # optional - sage.symbolic
        True

        sage: F = FractionField(PolynomialRing(QQ,'t'))
        sage: t = F.gen()
        sage: E = EllipticCurve([t,0]); E
        Elliptic Curve defined by y^2 = x^3 + t*x
         over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>
        sage: E.category()
        Category of schemes over
         Fraction Field of Univariate Polynomial Ring in t over Rational Field

    See :trac:`12517`::

        sage: E = EllipticCurve([1..5])
        sage: EllipticCurve(E.a_invariants())
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

    See :trac:`11773`::

        sage: E = EllipticCurve()
        Traceback (most recent call last):
        ...
        TypeError: invalid input to EllipticCurve constructor
    """
    def create_key_and_extra_args(self, x=None, y=None, j=None, minimal_twist=True, **kwds):
        r"""
        Return a ``UniqueFactory`` key and possibly extra parameters.

        INPUT: See the documentation for :class:`EllipticCurveFactory`.

        OUTPUT:

        A pair ``(key, extra_args)``:

        - ``key`` has the form `(R, (a_1, a_2, a_3, a_4, a_6))`,
          representing a ring and the Weierstrass coefficients of an
          elliptic curve over that ring;

        - ``extra_args`` is a dictionary containing additional data to
          be inserted into the elliptic curve structure.

        EXAMPLES::

            sage: EllipticCurve.create_key_and_extra_args(j=8000)
            ((Rational Field, (0, 1, 0, -3, 1)), {})

        When constructing a curve over `\QQ` from a Cremona or LMFDB
        label, the invariants from the database are returned as
        ``extra_args``::

            sage: key, data = EllipticCurve.create_key_and_extra_args('389.a1')
            sage: key
            (Rational Field, (0, 1, 1, -2, 0))
            sage: data['conductor']
            389
            sage: data['cremona_label']
            '389a1'
            sage: data['lmfdb_label']
            '389.a1'
            sage: data['rank']
            2
            sage: data['torsion_order']
            1

        User-specified keywords are also included in ``extra_args``::

            sage: key, data = EllipticCurve.create_key_and_extra_args((0, 0, 1, -23737, 960366), rank=4)
            sage: data['rank']
            4

        Furthermore, keywords takes precedence over data from the
        database, which can be used to specify an alternative set of
        generators for the Mordell-Weil group::

            sage: key, data = EllipticCurve.create_key_and_extra_args('5077a1', gens=[[1, -1], [-2, 3], [4, -7]])
            sage: data['gens']
            [[1, -1], [-2, 3], [4, -7]]
            sage: E = EllipticCurve.create_object(0, key, **data)
            sage: E.gens()
            [(-2 : 3 : 1), (1 : -1 : 1), (4 : -7 : 1)]

        Note that elliptic curves are equal if and only they have the
        same base ring and Weierstrass equation; the data in
        ``extra_args`` do not influence comparison of elliptic curves.
        A consequence of this is that passing keyword arguments only
        works when constructing an elliptic curve the first time::

            sage: E = EllipticCurve('433a1', gens=[[-1, 1], [3, 4]])
            sage: E.gens()
            [(-1 : 1 : 1), (3 : 4 : 1)]
            sage: E = EllipticCurve('433a1', gens=[[-1, 0], [0, 1]])
            sage: E.gens()
            [(-1 : 1 : 1), (3 : 4 : 1)]

        .. WARNING::

            Manually specifying extra data is almost never necessary
            and is not guaranteed to have any effect, as the above
            example shows.  Almost no checking is done, so specifying
            incorrect data may lead to wrong results of computations
            instead of errors or warnings.

        TESTS::

            sage: var('x', 'y', 'v', 'w')
            (x, y, v, w)
            sage: EllipticCurve(y^2 + y > x^3 + x - 9)
            Traceback (most recent call last):
            ...
            ValueError: no symbolic relations other than equalities are allowed
            sage: E = EllipticCurve(y^2 + y == x^3 + x - 9)
            sage: E is EllipticCurve(y^2 + y - ( x^3 + x - 9 ))
            True
            sage: R.<x,y> = QQ[]
            sage: E is EllipticCurve(y^2 + y - ( x^3 + x - 9 ))
            True
        """
        R = None
        if is_Ring(x):
            (R, x) = (x, y)

        if j is not None:
            if R is not None:
                try:
                    j = R(j)
                except (ZeroDivisionError, ValueError, TypeError):
                    raise ValueError("First parameter must be a ring containing %s" % j)
            elif x is not None:
                raise ValueError("First parameter (if present) must be a ring when j is specified")
            x = coefficients_from_j(j, minimal_twist)

        if isinstance(x, Expression) and x.is_relational():
            import operator
            if x.operator() != operator.eq:
                raise ValueError("no symbolic relations other than equalities are allowed")
            x = x.lhs() - x.rhs()

        if isinstance(parent(x), sage.rings.abc.SymbolicRing):
            x = x._polynomial_(QQ['x', 'y'])

        if isinstance(x, MPolynomial):
            if y is None:
                x = coefficients_from_Weierstrass_polynomial(x)
            else:
                # x is a cubic, y a rational point
                x = EllipticCurve_from_cubic(x, y, morphism=False).ainvs()

        if isinstance(x, str):
            # Interpret x as a Cremona or LMFDB label.
            from sage.databases.cremona import CremonaDatabase
            x, data = CremonaDatabase().coefficients_and_data(x)
            # User-provided keywords may override database entries.
            data.update(kwds)
            kwds = data

        if not isinstance(x, (list, tuple)):
            raise TypeError("invalid input to EllipticCurve constructor")

        if len(x) == 2:
            x = (0, 0, 0, x[0], x[1])
        elif len(x) != 5:
            raise ValueError("sequence of coefficients must have length 2 or 5")

        if R is None:
            R = Sequence(x).universe()
            if R in (ZZ, int):
                R = QQ

        return (R, tuple(R(a) for a in x)), kwds

    def create_object(self, version, key, **kwds):
        r"""
        Create an object from a ``UniqueFactory`` key.

        EXAMPLES::

            sage: E = EllipticCurve.create_object(0, (GF(3), (1, 2, 0, 1, 2)))          # optional - sage.rings.finite_rings
            sage: type(E)                                                               # optional - sage.rings.finite_rings
            <class 'sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field_with_category'>

        .. NOTE::

            Keyword arguments are currently only passed to the
            constructor for elliptic curves over `\QQ`; elliptic
            curves over other fields do not support them.
        """
        R, x = key

        if R is QQ:
            from .ell_rational_field import EllipticCurve_rational_field
            return EllipticCurve_rational_field(x, **kwds)
        elif isinstance(R, NumberField):
            from .ell_number_field import EllipticCurve_number_field
            return EllipticCurve_number_field(R, x)
        elif isinstance(R, sage.rings.abc.pAdicField):
            from .ell_padic_field import EllipticCurve_padic_field
            return EllipticCurve_padic_field(R, x)
        elif isinstance(R, FiniteField) or (isinstance(R, sage.rings.abc.IntegerModRing) and R.characteristic().is_prime()):
            from .ell_finite_field import EllipticCurve_finite_field
            return EllipticCurve_finite_field(R, x)
        elif R in _Fields:
            from .ell_field import EllipticCurve_field
            return EllipticCurve_field(R, x)
        from .ell_generic import EllipticCurve_generic
        return EllipticCurve_generic(R, x)


EllipticCurve = EllipticCurveFactory('sage.schemes.elliptic_curves.constructor.EllipticCurve')


def EllipticCurve_from_Weierstrass_polynomial(f):
    r"""
    Return the elliptic curve defined by a cubic in (long) Weierstrass
    form.

    INPUT:

    - ``f`` -- a inhomogeneous cubic polynomial in long Weierstrass
      form.

    OUTPUT: The elliptic curve defined by it.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: f = y^2 + 1*x*y + 3*y - (x^3 + 2*x^2 + 4*x + 6)
        sage: EllipticCurve(f)
        Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 6 over Rational Field
        sage: EllipticCurve(f).a_invariants()
        (1, 2, 3, 4, 6)

    The polynomial ring may have extra variables as long as they
    do not occur in the polynomial itself::

        sage: R.<x,y,z,w> = QQ[]
        sage: EllipticCurve(-y^2 + x^3 + 1)
        Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
        sage: EllipticCurve(-x^2 + y^3 + 1)
        Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
        sage: EllipticCurve(-w^2 + z^3 + 1)
        Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field

    TESTS::

        sage: from sage.schemes.elliptic_curves.constructor import EllipticCurve_from_Weierstrass_polynomial
        sage: EllipticCurve_from_Weierstrass_polynomial(-w^2 + z^3 + 1)
        Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
    """
    return EllipticCurve(coefficients_from_Weierstrass_polynomial(f))

def coefficients_from_Weierstrass_polynomial(f):
    r"""
    Return the coefficients `[a_1, a_2, a_3, a_4, a_6]` of a cubic in
    Weierstrass form.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import coefficients_from_Weierstrass_polynomial
        sage: R.<w,z> = QQ[]
        sage: coefficients_from_Weierstrass_polynomial(-w^2 + z^3 + 1)
        [0, 0, 0, 0, 1]
    """
    R = f.parent()
    cubic_variables = [ x for x in R.gens() if f.degree(x) == 3 ]
    quadratic_variables = [ y for y in R.gens() if f.degree(y) == 2 ]
    try:
        x = cubic_variables[0]
        y = quadratic_variables[0]
    except IndexError:
        raise ValueError('polynomial is not in long Weierstrass form')

    a1 = a2 = a3 = a4 = a6 = 0
    x3 = y2 = None
    for coeff, mon in f:
        if mon == x**3:
            x3 = coeff
        elif mon == x**2:
            a2 = coeff
        elif mon == x:
            a4 = coeff
        elif mon == 1:
            a6 = coeff
        elif mon == y**2:
            y2 = -coeff
        elif mon == x*y:
            a1 = -coeff
        elif mon == y:
            a3 = -coeff
        else:
            raise ValueError('polynomial is not in long Weierstrass form')

    if x3 != y2:
        raise ValueError('the coefficient of x^3 and -y^2 must be the same')
    elif x3 != 1:
        a1, a2, a3, a4, a6 = a1/x3, a2/x3, a3/x3, a4/x3, a6/x3
    return [a1, a2, a3, a4, a6]


def EllipticCurve_from_c4c6(c4, c6):
    r"""
    Return an elliptic curve with given `c_4` and
    `c_6` invariants.

    EXAMPLES::

        sage: E = EllipticCurve_from_c4c6(17, -2005)
        sage: E
        Elliptic Curve defined by y^2  = x^3 - 17/48*x + 2005/864 over Rational Field
        sage: E.c_invariants()
        (17, -2005)
    """
    try:
        K = c4.parent()
    except AttributeError:
        K = RationalField()
    if K not in _Fields:
        K = K.fraction_field()
    return EllipticCurve([-K(c4)/K(48), -K(c6)/K(864)])


def EllipticCurve_from_j(j, minimal_twist=True):
    r"""
    Return an elliptic curve with given `j`-invariant.

    INPUT:

    - ``j`` -- an element of some field.

    - ``minimal_twist`` (boolean, default True) -- If True and ``j``
      is in `\QQ`, the curve returned is a minimal twist, i.e. has
      minimal conductor; when there is more than one curve with
      minimal conductor, the curve returned is the one whose label
      comes first if the curves are in the CremonaDatabase, otherwise
      the one whose minimal a-invariants are first lexicographically.
      If `j` is not in `\QQ` this parameter is ignored.

    OUTPUT:

    An elliptic curve with `j`-invariant `j`.

    EXAMPLES::

        sage: E = EllipticCurve_from_j(0); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 + y = x^3 over Rational Field
        0
        '27a3'

        sage: E = EllipticCurve_from_j(1728); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        1728
        '32a2'

        sage: E = EllipticCurve_from_j(1); E; E.j_invariant()
        Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field
        1

    The ``minimal_twist`` parameter (ignored except over `\QQ` and
    True by default) controls whether or not a minimal twist is
    computed::

        sage: EllipticCurve_from_j(100)
        Elliptic Curve defined by y^2 = x^3 + x^2 + 3392*x + 307888 over Rational Field
        sage: _.conductor()
        33129800
        sage: EllipticCurve_from_j(100, minimal_twist=False)
        Elliptic Curve defined by y^2 = x^3 + 488400*x - 530076800 over Rational Field
        sage: _.conductor()
        298168200

    Since computing the minimal twist requires factoring both `j` and
    `j-1728` the following example would take a long time without
    setting ``minimal_twist`` to False::

       sage: E = EllipticCurve_from_j(2^256+1, minimal_twist=False)
       sage: E.j_invariant() == 2^256+1
       True
    """
    return EllipticCurve(coefficients_from_j(j, minimal_twist))


def coefficients_from_j(j, minimal_twist=True):
    r"""
    Return Weierstrass coefficients `(a_1, a_2, a_3, a_4, a_6)` for an
    elliptic curve with given `j`-invariant.

    INPUT: See :func:`EllipticCurve_from_j`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import coefficients_from_j
        sage: coefficients_from_j(0)
        [0, 0, 1, 0, 0]
        sage: coefficients_from_j(1728)
        [0, 0, 0, -1, 0]
        sage: coefficients_from_j(1)
        [1, 0, 0, 36, 3455]

    The ``minimal_twist`` parameter (ignored except over `\QQ` and
    True by default) controls whether or not a minimal twist is
    computed::

        sage: coefficients_from_j(100)
        [0, 1, 0, 3392, 307888]
        sage: coefficients_from_j(100, minimal_twist=False)
        [0, 0, 0, 488400, -530076800]
    """
    try:
        K = j.parent()
    except AttributeError:
        K = RationalField()
    if K not in _Fields:
        K = K.fraction_field()

    char = K.characteristic()
    if char == 2:
        if j == 0:
            return Sequence([0, 0, 1, 0, 0], universe=K)
        else:
            return Sequence([1, 0, 0, 0, 1/j], universe=K)
    if char == 3:
        if j == 0:
            return Sequence([0, 0, 0, 1, 0], universe=K)
        else:
            return Sequence([0, j, 0, 0, -j**2], universe=K)

    if K is RationalField():
        # we construct the minimal twist, i.e. the curve with minimal
        # conductor with this j_invariant:
        if j == 0:
            return Sequence([0, 0, 1, 0, 0], universe=K) # 27a3
        if j == 1728:
            return Sequence([0, 0, 0, -1, 0], universe=K) # 32a2

        if not minimal_twist:
            k = j-1728
            return Sequence([0, 0, 0, -3*j*k, -2*j*k**2], universe=K)

        n = j.numerator()
        m = n-1728*j.denominator()
        a4 = -3*n*m
        a6 = -2*n*m**2

        # Now E=[0,0,0,a4,a6] has j-invariant j=n/d
        from sage.sets.set import Set
        for p in Set(n.prime_divisors()+m.prime_divisors()):
            e = min(a4.valuation(p)//2,a6.valuation(p)//3)
            if e > 0:
                p = p**e
                a4 /= p**2
                a6 /= p**3

        # Now E=[0,0,0,a4,a6] is minimal at all p != 2,3
        tw = [-1,2,-2,3,-3,6,-6]
        E1 = EllipticCurve([0,0,0,a4,a6])
        Elist = [E1] + [E1.quadratic_twist(t) for t in tw]
        min_cond = min(E.conductor() for E in Elist)
        Elist = [E for E in Elist if E.conductor() == min_cond]
        if len(Elist) > 1:
            from sage.databases.cremona import CremonaDatabase, parse_cremona_label
            if min_cond <= CremonaDatabase().largest_conductor():
                sorter = lambda E: parse_cremona_label(E.label(), numerical_class_code=True)
            else:
                sorter = lambda E: E.ainvs()
            Elist.sort(key=sorter)
        return Sequence(Elist[0].ainvs())

    # defaults for all other fields:
    if j == 0:
        return Sequence([0, 0, 0, 0, 1], universe=K)
    if j == 1728:
        return Sequence([0, 0, 0, 1, 0], universe=K)
    k = j-1728
    return Sequence([0, 0, 0, -3*j*k, -2*j*k**2], universe=K)


def EllipticCurve_from_cubic(F, P=None, morphism=True):
    r"""
    Construct an elliptic curve from a ternary cubic with a rational point.

    If you just want the Weierstrass form and are not interested in
    the morphism then it is easier to use the function
    :func:`~sage.schemes.elliptic_curves.jacobian.Jacobian`
    instead. If there is a rational point on the given cubic, this
    function will construct the same elliptic curve but you do not have to
    supply the point ``P``.

    INPUT:

    - ``F`` -- a homogeneous cubic in three variables with rational
      coefficients, as a polynomial ring element, defining a smooth
      plane cubic curve `C`.

    - ``P`` -- a 3-tuple `(x,y,z)` defining a projective point on `C`,
      or ``None``.  If ``None`` then a rational flex will be used as a
      base point if one exists, otherwise an error will be raised.

    - ``morphism`` -- boolean (default: ``True``).  If ``True``
      returns a birational isomorphism from `C` to a Weierstrass
      elliptic curve `E`, otherwise just returns `E`.

    OUTPUT:

    Either (when ``morphism``=``False``) an elliptic curve `E` in long
    Weierstrass form isomorphic to the plane cubic curve `C` defined
    by the equation `F=0`.

    Or (when ``morphism=True``), a birational isomorphism from `C` to
    the elliptic curve `E`. If the given point is a flex, this is a
    linear isomorphism.

    .. NOTE::

        The function
        :func:`~sage.schemes.elliptic_curves.jacobian.Jacobian` may be
        used instead.  It constructs the same elliptic curve (which is in
        all cases the Jacobian of `(F=0)`) and needs no base point to be
        provided, but also returns no isomorphism since in general there
        is none: the plane cubic is only isomorphic to its Jacobian when
        it has a rational point.

    .. NOTE::

        When ``morphism=True``, a birational isomorphism between the
        curve `F=0` and the Weierstrass curve is returned. If the point
        happens to be a flex, then this is a linear isomorphism.  The
        morphism does not necessarily take the given point `P` to the
        point at infinity on `E`, since we always use a rational flex
        on `C` as base-point when one exists.

    EXAMPLES:

    First we find that the Fermat cubic is isomorphic to the curve
    with Cremona label 27a1::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^3 + y^3 + z^3
        sage: P = [1,-1,0]
        sage: E = EllipticCurve_from_cubic(cubic, P, morphism=False); E
        Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
        sage: E.cremona_label()
        '27a1'
        sage: EllipticCurve_from_cubic(cubic, [0,1,-1], morphism=False).cremona_label()
        '27a1'
        sage: EllipticCurve_from_cubic(cubic, [1,0,-1], morphism=False).cremona_label()
        '27a1'

    Next we find the minimal model and conductor of the Jacobian of the
    Selmer curve::

        sage: R.<a,b,c> = QQ[]
        sage: cubic = a^3 + b^3 + 60*c^3
        sage: P = [1,-1,0]
        sage: E = EllipticCurve_from_cubic(cubic, P, morphism=False);  E
        Elliptic Curve defined by y^2 - 540*y = x^3 - 97200 over Rational Field
        sage: E.minimal_model()
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
        sage: E.conductor()
        24300

    We can also get the birational isomorphism to and from the
    Weierstrass form. We start with an example where ``P`` is a flex
    and the equivalence is a linear isomorphism::

        sage: f = EllipticCurve_from_cubic(cubic, P, morphism=True)
        sage: f
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by a^3 + b^3 + 60*c^3
          To:   Elliptic Curve defined by y^2 - 540*y = x^3 - 97200 over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (-c : 3*a : 1/180*a + 1/180*b)

        sage: finv = f.inverse();  finv
        Scheme morphism:
          From: Elliptic Curve defined by y^2 - 540*y = x^3 - 97200 over Rational Field
          To:   Projective Plane Curve over Rational Field defined by a^3 + b^3 + 60*c^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (1/3*y : -1/3*y + 180*z : -x)

        Scheme morphism:
          From: Elliptic Curve defined by y^2 + 2*x*y + 20*y = x^3 - x^2 - 20*x - 400/3
                over Rational Field
          To:   Closed subscheme of Projective Space of dimension 2 over Rational Field
                defined by: a^3 + b^3 + 60*c^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (x + y + 20*z : -x - y : -x)

    We verify that `f` maps the chosen point `P=(1,-1,0)` on the cubic
    to the origin of the elliptic curve::

        sage: f([1,-1,0])
        (0 : 1 : 0)
        sage: finv([0,1,0])
        (-1 : 1 : 0)

    To verify the output, we plug in the polynomials to check that
    this indeed transforms the cubic into Weierstrass form::

        sage: cubic(finv.defining_polynomials()) * finv.post_rescaling()
        -x^3 + y^2*z - 540*y*z^2 + 97200*z^3

        sage: E.defining_polynomial()(f.defining_polynomials()) * f.post_rescaling()
        a^3 + b^3 + 60*c^3

    If the given point is not a flex and the cubic has no rational
    flexes, then the cubic can not be transformed to a Weierstrass
    equation by a linear transformation. The general birational
    transformation is still a birational isomorphism, but is
    quadratic::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^2*y + 4*x*y^2 + x^2*z + 8*x*y*z + 4*y^2*z + 9*x*z^2 + 9*y*z^2
        sage: f = EllipticCurve_from_cubic(cubic, [1,-1,1], morphism=True); f
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined
                by x^2*y + 4*x*y^2 + x^2*z + 8*x*y*z + 4*y^2*z + 9*x*z^2 + 9*y*z^2
          To:   Elliptic Curve defined
                by y^2 + 7560/19*x*y + 552960000000/2352637*y = x^3 - 3445200/133*x^2
                over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
                (2527/17280*x^2 + 133/2160*x*y + 133/108000*y^2 + 133/2880*x*z
                   + 931/18000*y*z - 3857/48000*z^2
                 : -6859/288*x^2 + 323/36*x*y + 359/1800*y^2 + 551/48*x*z
                   + 2813/300*y*z + 24389/800*z^2
                 : -2352637/99532800000*x^2 - 2352637/124416000000*x*y
                   - 2352637/622080000000*y^2 + 2352637/82944000000*x*z
                   + 2352637/207360000000*y*z - 2352637/276480000000*z^2)

    Note that the morphism returned cannot be evaluated directly at
    the given point ``P=(1:-1:1)`` since the polynomials defining it
    all vanish there::

        sage: f([1,-1,1])
        Traceback (most recent call last):
        ...
        ValueError: [0, 0, 0] does not define a valid point since all entries are 0

    Using the group law on the codomain elliptic curve, which has rank
    1 and full 2-torsion, and the inverse morphism, we can find many
    points on the cubic.  First we find the preimages of multiples of
    the generator::

        sage: E = f.codomain()
        sage: E.label()
        '720e2'
        sage: E.rank()
        1
        sage: R = E.gens()[0]; R
        (-17280000/2527 : 9331200000/6859 : 1)
        sage: finv = f.inverse()
        sage: [finv(k*R) for k in range(1,10)]
        [(-4 : 1 : 0),
        (-1 : 4 : 1),
        (-20 : -55/76 : 1),
        (319/399 : -11339/7539 : 1),
        (159919/14360 : -4078139/1327840 : 1),
        (-27809119/63578639 : 1856146436/3425378659 : 1),
        (-510646582340/56909753439 : 424000923715/30153806197284 : 1),
        (-56686114363679/4050436059492161 : -2433034816977728281/1072927821085503881 : 1),
        (650589589099815846721/72056273157352822480 : -347376189546061993109881/194127383495944026752320 : 1)]

    The elliptic curve also has torsion, which we can map back::

        sage: E.torsion_points()
        [(-144000000/17689 : 3533760000000/2352637 : 1),
        (-92160000/17689 : 2162073600000/2352637 : 1),
        (-5760000/17689 : -124070400000/2352637 : 1),
        (0 : 1 : 0)]
        sage: [finv(Q) for Q in E.torsion_points() if Q]
        [(9 : -9/4 : 1), (-9 : 0 : 1), (0 : 1 : 0)]

    In this example, the given point ``P`` is not a flex but the cubic
    does have a rational flex, ``(-4:0:1)``.  We return a linear
    isomorphism which maps this flex to the point at infinity on the
    Weierstrass model::

        sage: R.<a,b,c> = QQ[]
        sage: cubic =  a^3 + 7*b^3 + 64*c^3
        sage: P = [2,2,-1]
        sage: f = EllipticCurve_from_cubic(cubic, P, morphism=True)
        sage: E = f.codomain();  E
        Elliptic Curve defined by y^2 - 258048*y = x^3 - 22196256768 over Rational Field
        sage: E.minimal_model()
        Elliptic Curve defined by y^2 + y = x^3 - 331 over Rational Field

        sage: f
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by a^3 + 7*b^3 + 64*c^3
          To:   Elliptic Curve defined by y^2 - 258048*y = x^3 - 22196256768 over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (b : -48*a : -1/5376*a - 1/1344*c)

        sage: finv = f.inverse();  finv
        Scheme morphism:
          From: Elliptic Curve defined by y^2 - 258048*y = x^3 - 22196256768 over Rational Field
          To:   Projective Plane Curve over Rational Field defined by a^3 + 7*b^3 + 64*c^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (-1/48*y : x : 1/192*y - 1344*z)

        sage: cubic(finv.defining_polynomials()) * finv.post_rescaling()
        -x^3 + y^2*z - 258048*y*z^2 + 22196256768*z^3

        sage: E.defining_polynomial()(f.defining_polynomials()) * f.post_rescaling()
        a^3 + 7*b^3 + 64*c^3

        sage: f(P)
        (5376 : -258048 : 1)
        sage: f([-4,0,1])
        (0 : 1 : 0)

    It is possible to not provide a base point ``P`` provided that the
    cubic has a rational flex.  In this case the flexes will be found
    and one will be used as a base point::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^3 + y^3 + z^3
        sage: f = EllipticCurve_from_cubic(cubic, morphism=True)
        sage: f
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
          To:   Elliptic Curve defined by y^2 - 9*y = x^3 - 27 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
                (y : -3*x : -1/3*x - 1/3*z)

    An error will be raised if no point is given and there are no rational flexes::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = 3*x^3 + 4*y^3 + 5*z^3
        sage: EllipticCurve_from_cubic(cubic)
        Traceback (most recent call last):
        ...
        ValueError: A point must be given when the cubic has no rational flexes

    An example over a finite field, using a flex::

        sage: K = GF(17)                                                                # optional - sage.rings.finite_rings
        sage: R.<x,y,z> = K[]                                                           # optional - sage.rings.finite_rings
        sage: cubic = 2*x^3 + 3*y^3 + 4*z^3                                             # optional - sage.rings.finite_rings
        sage: EllipticCurve_from_cubic(cubic, [0,3,1])                                  # optional - sage.rings.finite_rings
        Scheme morphism:
          From: Projective Plane Curve over Finite Field of size 17
                defined by 2*x^3 + 3*y^3 + 4*z^3
          To:   Elliptic Curve defined by y^2 + 16*y = x^3 + 11 over Finite Field of size 17
          Defn: Defined on coordinates by sending (x : y : z) to
                (-x : 4*y : 4*y + 5*z)

    An example in characteristic 3::

        sage: K = GF(3)                                                                 # optional - sage.rings.finite_rings
        sage: R.<x,y,z> = K[]                                                           # optional - sage.rings.finite_rings
        sage: cubic = x^3 + y^3 + z^3 + x*y*z                                           # optional - sage.rings.finite_rings
        sage: EllipticCurve_from_cubic(cubic, [0,1,-1])                                 # optional - sage.rings.finite_rings
        Scheme morphism:
          From: Projective Plane Curve over Finite Field of size 3
                defined by x^3 + y^3 + x*y*z + z^3
          To:   Elliptic Curve defined by y^2 + x*y = x^3 + 1 over Finite Field of size 3
          Defn: Defined on coordinates by sending (x : y : z) to
                (y + z : -y : x)

    An example over a number field, using a non-flex and where there are no rational flexes::

        sage: K.<a> = QuadraticField(-3)                                                # optional - sage.rings.number_field
        sage: R.<x,y,z> = K[]                                                           # optional - sage.rings.number_field
        sage: cubic = 2*x^3 + 3*y^3 + 5*z^3                                             # optional - sage.rings.number_field
        sage: EllipticCurve_from_cubic(cubic, [1,1,-1])                                 # optional - sage.rings.number_field
        Scheme morphism:
          From: Projective Plane Curve over Number Field in a
                with defining polynomial x^2 + 3 with a = 1.732050807568878?*I
                defined by 2*x^3 + 3*y^3 + 5*z^3
          To:   Elliptic Curve defined by
                y^2 + 1754460/2053*x*y + 5226454388736000/8653002877*y
                = x^3 + (-652253285700/4214809)*x^2
                over Number Field in a with defining polynomial x^2 + 3
                with a = 1.732050807568878?*I
          Defn: Defined on coordinates by sending (x : y : z) to
                (-16424/127575*x^2 - 231989/680400*x*y - 14371/64800*y^2 - 26689/81648*x*z - 10265/27216*y*z - 2053/163296*z^2
                 : 24496/315*x^2 + 119243/840*x*y + 4837/80*y^2 + 67259/504*x*z + 25507/168*y*z + 5135/1008*z^2
                 : 8653002877/2099914709760000*x^2 + 8653002877/699971569920000*x*y + 8653002877/933295426560000*y^2 + 8653002877/419982941952000*x*z + 8653002877/279988627968000*y*z + 8653002877/335986353561600*z^2)

    An example over a function field, using a non-flex::

        sage: K.<t> = FunctionField(QQ)
        sage: R.<x,y,z> = K[]
        sage: cubic = x^3 + t*y^3 + (1+t)*z^3
        sage: EllipticCurve_from_cubic(cubic, [1,1,-1], morphism=False)
        Elliptic Curve defined by y^2 + ((162*t^6+486*t^5+810*t^4+810*t^3+486*t^2+162*t)/(t^6+12*t^5-3*t^4-20*t^3-3*t^2+12*t+1))*x*y + ((314928*t^14+4094064*t^13+23462136*t^12+78102144*t^11+167561379*t^10+243026001*t^9+243026001*t^8+167561379*t^7+78102144*t^6+23462136*t^5+4094064*t^4+314928*t^3)/(t^14+40*t^13+577*t^12+3524*t^11+8075*t^10+5288*t^9-8661*t^8-17688*t^7-8661*t^6+5288*t^5+8075*t^4+3524*t^3+577*t^2+40*t+1))*y = x^3 + ((2187*t^12+13122*t^11-17496*t^10-207765*t^9-516132*t^8-673596*t^7-516132*t^6-207765*t^5-17496*t^4+13122*t^3+2187*t^2)/(t^12+24*t^11+138*t^10-112*t^9-477*t^8+72*t^7+708*t^6+72*t^5-477*t^4-112*t^3+138*t^2+24*t+1))*x^2
         over Rational function field in t over Rational Field

    TESTS:

    Here is a test for :trac:`21092`::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = -3*x^2*y + 3*x*y^2 + 4*x^2*z + 4*y^2*z - 3*x*z^2 + 3*y*z^2 - 8*z^3
        sage: EllipticCurve_from_cubic(cubic, (-4/5, 4/5, 3/5) )
        Scheme morphism:
          From: Projective Plane Curve over Rational Field defined
                by -3*x^2*y + 3*x*y^2 + 4*x^2*z + 4*y^2*z - 3*x*z^2 + 3*y*z^2 - 8*z^3
          To:   Elliptic Curve defined by y^2 + 24*x*y + 3024*y = x^3 + 495*x^2 + 36288*x
                over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
                (-1/3*z : 3*x : -1/1008*x + 1/1008*y + 1/378*z)
    """
    from sage.schemes.curves.constructor import Curve
    from sage.matrix.constructor import Matrix
    from sage.schemes.elliptic_curves.weierstrass_transform import \
        WeierstrassTransformationWithInverse

    # check the input
    R = F.parent()
    K = R.base_ring()
    if not is_MPolynomialRing(R):
        raise TypeError('equation must be a polynomial')
    if R.ngens() != 3 or F.nvariables() != 3:
        raise TypeError('equation must be a polynomial in three variables')
    if not F.is_homogeneous():
        raise TypeError('equation must be a homogeneous polynomial')

    C = Curve(F)
    if P:
        try:
            CP = C(P)
        except (TypeError, ValueError):
            raise TypeError('{} does not define a point on a projective curve over {} defined by {}'.format(P,K,F))

    x, y, z = R.gens()

    # Test whether P is a flex; if not test whether there are any rational flexes:

    hessian = Matrix([[F.derivative(v1, v2) for v1 in R.gens()] for v2 in R.gens()]).det()
    if P and hessian(P) == 0:
        flex_point = P
    else:
        flexes = C.intersection(Curve(hessian)).rational_points()
        if flexes:
            flex_point = list(flexes[0])
            if not P:
                P = flex_point
                CP = C(P)
        else:
            flex_point = None

    if flex_point is not None: # first case: base point is a flex
        P = flex_point
        L = tangent_at_smooth_point(C,P)
        dx, dy, dz = [L.coefficient(v) for v in R.gens()]

        # find an invertible matrix M such that (0,1,0)M=P and
        # ML'=(0,0,1)' where L=[dx,dy,dx].  Then the linear transform
        # by M takes P to [0,1,0] and L to Z=0:

        if P[0]:
            Q1 = [0,-dz,dy]
            Q2 = [0,1,0] if dy else [0,0,1]
        elif P[1]:
            Q1 = [dz,0,-dx]
            Q2 = [1,0,0] if dx else [0,0,1]
        else:
            Q1 = [-dy,dx,0]
            Q2 = [1,0,0] if dx else [0,1,0]

        M = Matrix(K,[Q1,P,Q2])
        # assert M.is_invertible()
        # assert list(vector([0,1,0])*M) == P
        # assert list(M*vector([dx,dy,dz]))[:2] == [0,0]

        M = M.transpose()
        F2 = R(M.act_on_polynomial(F))

        # scale and dehomogenise
        a = K(F2.coefficient(x**3))
        b = K(F2.coefficient(y*y*z))

        F3 = F2([-x, y/b, z*a*b]) / a
        # assert F3.coefficient(x**3) == -1
        # assert F3.coefficient(y*y*z) == 1
        E = EllipticCurve(F3([x,y,1]))
        if not morphism:
            return E

        # Construct the (linear) morphism
        M = M * Matrix(K,[[-1,0,0],[0,1/b,0],[0,0,a*b]])
        inv_defining_poly = [ M[i,0]*x + M[i,1]*y + M[i,2]*z for i in range(3) ]
        inv_post = 1/a
        M = M.inverse()
        fwd_defining_poly = [ M[i,0]*x + M[i,1]*y + M[i,2]*z for i in range(3) ]
        fwd_post = a

    else: # Second case: no flexes
        if not P:
            raise ValueError('A point must be given when the cubic has no rational flexes')
        L = tangent_at_smooth_point(C,P)
        Qlist = [Q for Q in C.intersection(Curve(L)).rational_points() if C(Q) != CP]
        # assert Qlist
        P2 = C(Qlist[0])
        L2 = tangent_at_smooth_point(C,P2)
        Qlist = [Q for Q in C.intersection(Curve(L2)).rational_points() if C(Q) != P2]
        # assert Qlist
        P3 = C(Qlist[0])

        # NB This construction of P3 relies on P2 not being a flex.
        # If we want to use a non-flex as P when there are rational
        # flexes this would be a problem.  However, the only condition
        # which P3 must satisfy is that it is on the tangent at P2, it
        # need not lie on the cubic.

        # send P, P2, P3 to (1:0:0), (0:1:0), (0:0:1) respectively
        M = Matrix(K, [P, list(P2), list(P3)]).transpose()
        F2 = M.act_on_polynomial(F)
        xyzM = [ M[i,0]*x + M[i,1]*y + M[i,2]*z for i in range(3) ]
        # assert F(xyzM)==F2

        # substitute x = U^2, y = V*W, z = U*W, and rename (x,y,z)=(U,V,W)
        T1 = [x*x,y*z,x*z]
        S1 = x**2*z
        F3 = F2(T1) // S1
        xyzC = [ t(T1) for t in xyzM ]
        # assert F3 == F(xyzC) // S1

        # scale and dehomogenise
        a = K(F3.coefficient(x**3))
        b = K(F3.coefficient(y*y*z))
        ab = a*b

        T2 = [-x, y/b, ab*z]
        F4 = F3(T2) / a
        # assert F4.coefficient(x**3) == -1
        # assert F4.coefficient(y*y*z) == 1
        xyzW = [ t(T2) for t in xyzC ]
        S2 = a*S1(T2)
        # assert F4 == F(xyzW) // S2

        E = EllipticCurve(F4([x,y,1]))
        if not morphism:
            return E

        inv_defining_poly = xyzW
        inv_post = 1/S2
        # assert F4==F(inv_defining_poly)*inv_post
        MI = M.inverse()
        xyzI = [ (MI[i,0]*x + MI[i,1]*y + MI[i,2]*z) for i in range(3) ]
        T1I = [x*z,x*y,z*z] # inverse of T1
        xyzIC = [ t(xyzI) for t in T1I ]
        T2I = [-x, b*y, z/ab] # inverse of T2
        xyzIW = [ t(xyzIC) for t in T2I ]
        fwd_defining_poly = xyzIW
        fwd_post = a/(x*z*z)(xyzI)
        # assert F4(fwd_defining_poly)*fwd_post == F

    # Construct the morphism

    return WeierstrassTransformationWithInverse(
        C, E, fwd_defining_poly, fwd_post, inv_defining_poly, inv_post)


def tangent_at_smooth_point(C,P):
    r"""Return the tangent at the smooth point `P` of projective curve `C`.

    INPUT:

    - ``C`` -- a projective plane curve.

    - ``P`` -- a 3-tuple `(x,y,z)` defining a projective point on `C`.

    OUTPUT:

    The linear form defining the tangent at `P` to `C`.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: from sage.schemes.elliptic_curves.constructor import tangent_at_smooth_point
        sage: C = Curve(x^3 + y^3 + 60*z^3)
        sage: tangent_at_smooth_point(C, [1,-1,0])
        x + y

        sage: K.<t> = FunctionField(QQ)
        sage: R.<x,y,z> = K[]
        sage: C = Curve(x^3 + 2*y^3 + 3*z^3)
        sage: from sage.schemes.elliptic_curves.constructor import tangent_at_smooth_point
        sage: tangent_at_smooth_point(C,[1,1,-1])
        3*x + 6*y + 9*z
    """
    # Over function fields such as QQ(t) an error is raised with the
    # default (factor=True).  Note that factor=False returns the
    # product of the tangents in case of a multiple point, while here
    # `P` is assumed smooth so factorization is unnecessary, but over
    # QQ (for example) including the factorization gives better
    # results, for example returning x+y instead of 3x+3y in the
    # doctest.
    try:
        return C.tangents(P)[0]
    except NotImplementedError:
        return C.tangents(P,factor=False)[0]

def chord_and_tangent(F, P):
    r"""Return the third point of intersection of a cubic with the tangent at one point.

    INPUT:

    - ``F`` -- a homogeneous cubic in three variables with rational
      coefficients, as a polynomial ring element, defining a smooth
      plane cubic curve.

    - ``P`` -- a 3-tuple `(x,y,z)` defining a projective point on the
      curve `F=0`.

    OUTPUT:

    A point ``Q`` such that ``F(Q)=0``, namely the third point of
    intersection of the tangent at ``P`` with the curve ``F=0``, so
    ``Q=P`` if and only if ``P`` is a flex.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: from sage.schemes.elliptic_curves.constructor import chord_and_tangent
        sage: F = x^3 + y^3 + 60*z^3
        sage: chord_and_tangent(F, [1,-1,0])
        (-1 : 1 : 0)

        sage: F = x^3 + 7*y^3 + 64*z^3
        sage: p0 = [2,2,-1]
        sage: p1 = chord_and_tangent(F, p0);  p1
        (5 : -3 : 1)
        sage: p2 = chord_and_tangent(F, p1);  p2
        (-1265/314 : 183/314 : 1)

    TESTS::

        sage: F(list(p2))
        0
        sage: list(map(type, p2))
        [<... 'sage.rings.rational.Rational'>,
         <... 'sage.rings.rational.Rational'>,
         <... 'sage.rings.rational.Rational'>]

    See :trac:`16068`::

        sage: F = x**3 - 4*x**2*y - 65*x*y**2 + 3*x*y*z - 76*y*z**2
        sage: chord_and_tangent(F, [0, 1, 0])
        (0 : 0 : 1)
    """
    from sage.schemes.curves.constructor import Curve
    # check the input
    R = F.parent()
    if not is_MPolynomialRing(R):
        raise TypeError('equation must be a polynomial')
    if R.ngens() != 3:
        raise TypeError('{} is not a polynomial in three variables'.format(F))
    if not F.is_homogeneous():
        raise TypeError('{} is not a homogeneous polynomial'.format(F))
    x, y, z = R.gens()
    if len(P) != 3:
        raise TypeError('{} is not a projective point'.format(P))
    K = R.base_ring()
    try:
        C = Curve(F)
        P = C(P)
    except (TypeError, ValueError):
        raise TypeError('{} does not define a point on a projective curve over {} defined by {}'.format(P,K,F))

    L = Curve(tangent_at_smooth_point(C,P))
    Qlist = [Q for Q in C.intersection(L).rational_points() if Q != P]
    if Qlist:
        return Qlist[0]
    return P


def projective_point(p):
    r"""
    Return equivalent point with denominators removed

    INPUT:

    - ``P``, ``Q`` -- list/tuple of projective coordinates.

    OUTPUT:

    List of projective coordinates.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import projective_point
        sage: projective_point([4/5, 6/5, 8/5])
        [2, 3, 4]
        sage: F = GF(11)                                                                # optional - sage.rings.finite_rings
        sage: projective_point([F(4), F(8), F(2)])                                      # optional - sage.rings.finite_rings
        [4, 8, 2]
    """
    from sage.rings.integer import GCD_list
    from sage.arith.functions import LCM_list
    try:
        p_gcd = GCD_list([x.numerator() for x in p])
        p_lcm = LCM_list(x.denominator() for x in p)
    except AttributeError:
        return p
    scale = p_lcm / p_gcd
    return [scale * x for x in p]


def are_projectively_equivalent(P, Q, base_ring):
    r"""
    Test whether ``P`` and ``Q`` are projectively equivalent.

    INPUT:

    - ``P``, ``Q`` -- list/tuple of projective coordinates.

    - ``base_ring`` -- the base ring.

    OUTPUT: A boolean.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import are_projectively_equivalent
        sage: are_projectively_equivalent([0,1,2,3], [0,1,2,2], base_ring=QQ)
        False
        sage: are_projectively_equivalent([0,1,2,3], [0,2,4,6], base_ring=QQ)
        True
    """
    from sage.matrix.constructor import matrix
    return matrix(base_ring, [P, Q]).rank() < 2


def EllipticCurves_with_good_reduction_outside_S(S=[], proof=None, verbose=False):
    r"""
    Return a sorted list of all elliptic curves defined over `\QQ`
    with good reduction outside the set `S` of primes.

    INPUT:

    - ``S`` -- list of primes (default: empty list)

    - ``proof`` -- boolean (default ``True``): the MW basis for
      auxiliary curves will be computed with this proof flag

    - ``verbose`` -- boolean (default ``False``): if ``True``, some details
      of the computation will be output

    .. NOTE::

        Proof flag: The algorithm used requires determining all
        S-integral points on several auxiliary curves, which in turn
        requires the computation of their generators.  This is not
        always possible (even in theory) using current knowledge.

        The value of this flag is passed to the function which
        computes generators of various auxiliary elliptic curves, in
        order to find their S-integral points.  Set to ``False`` if the
        default (``True``) causes warning messages, but note that you can
        then not rely on the set of curves returned being
        complete.

    EXAMPLES::

        sage: EllipticCurves_with_good_reduction_outside_S([])
        []
        sage: elist = EllipticCurves_with_good_reduction_outside_S([2])
        sage: elist
        [Elliptic Curve defined by y^2 = x^3 + 4*x over Rational Field,
         Elliptic Curve defined by y^2 = x^3 - x over Rational Field,
         ...
         Elliptic Curve defined by y^2 = x^3 - x^2 - 13*x + 21 over Rational Field]
        sage: len(elist)
        24
        sage: ', '.join(e.label() for e in elist)
        '32a1, 32a2, 32a3, 32a4, 64a1, 64a2, 64a3, 64a4, 128a1, 128a2, 128b1, 128b2, 128c1, 128c2, 128d1, 128d2, 256a1, 256a2, 256b1, 256b2, 256c1, 256c2, 256d1, 256d2'

    Without ``Proof=False``, this example gives two warnings::

        sage: elist = EllipticCurves_with_good_reduction_outside_S([11], proof=False)   # long time (14s on sage.math, 2011)
        sage: len(elist)                                                                # long time
        12
        sage: ', '.join(e.label() for e in elist)                                       # long time
        '11a1, 11a2, 11a3, 121a1, 121a2, 121b1, 121b2, 121c1, 121c2, 121d1, 121d2, 121d3'

        sage: elist = EllipticCurves_with_good_reduction_outside_S([2,3])               # long time (26s on sage.math, 2011)
        sage: len(elist)                                                                # long time
        752
        sage: conds = sorted(set([e.conductor() for e in elist]))                       # long time
        sage: max(conds)                                                                # long time
        62208
        sage: [N.factor() for N in conds]                                               # long time
        [2^3 * 3,
         3^3,
         2^5,
         2^2 * 3^2,
         2^4 * 3,
         2 * 3^3,
         2^6,
         2^3 * 3^2,
         2^5 * 3,
         2^2 * 3^3,
         2^7,
         2^4 * 3^2,
         2 * 3^4,
         2^6 * 3,
         2^3 * 3^3,
         3^5,
         2^8,
         2^5 * 3^2,
         2^2 * 3^4,
         2^7 * 3,
         2^4 * 3^3,
         2 * 3^5,
         2^6 * 3^2,
         2^3 * 3^4,
         2^8 * 3,
         2^5 * 3^3,
         2^2 * 3^5,
         2^7 * 3^2,
         2^4 * 3^4,
         2^6 * 3^3,
         2^3 * 3^5,
         2^8 * 3^2,
         2^5 * 3^4,
         2^7 * 3^3,
         2^4 * 3^5,
         2^6 * 3^4,
         2^8 * 3^3,
         2^5 * 3^5,
         2^7 * 3^4,
         2^6 * 3^5,
         2^8 * 3^4,
         2^7 * 3^5,
         2^8 * 3^5]
    """
    from .ell_egros import egros_from_jlist, egros_get_j
    return egros_from_jlist(egros_get_j(S, proof=proof, verbose=verbose), S)
