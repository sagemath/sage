"""
Elliptic curve constructor

AUTHORS:

- William Stein (2005): Initial version

- John Cremona (2008-01): EllipticCurve(j) fixed for all cases
"""
from __future__ import absolute_import

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

import six
import sage.rings.all as rings

from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
from sage.rings.rational_field import is_RationalField
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.number_field.number_field import is_NumberField
from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
from sage.rings.ring import is_Ring

from sage.categories.fields import Fields
_Fields = Fields()

from sage.structure.sequence import Sequence
from sage.structure.element import parent
from sage.structure.factory import UniqueFactory
from sage.symbolic.ring import SR
from sage.symbolic.expression import is_SymbolicEquation


class EllipticCurveFactory(UniqueFactory):
    r"""
    Construct an elliptic curve.

    In Sage, an elliptic curve is always specified by
    (the coefficients of) a long Weierstrass equation

    .. math::

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

        sage: EllipticCurve([GF(5)(0),0,1,-1,0])
        Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5
        sage: EllipticCurve(GF(5), [0, 0,1,-1,0])
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
        Elliptic Curve defined by y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x over Complex Field with 53 bits of precision
        sage: E.j_invariant()
        2988.97297297297

    We can also create elliptic curves by giving the Weierstrass equation::

        sage: x, y = var('x,y')
        sage: EllipticCurve(y^2 + y ==  x^3 + x - 9)
        Elliptic Curve defined by y^2 + y = x^3 + x - 9 over Rational Field

        sage: R.<x,y> = GF(5)[]
        sage: EllipticCurve(x^3 + x^2 + 2 - y^2 - y*x)
        Elliptic Curve defined by y^2 + x*y  = x^3 + x^2 + 2 over Finite Field of size 5

    We can explicitly specify the `j`-invariant::

        sage: E = EllipticCurve(j=1728); E; E.j_invariant(); E.label()
        Elliptic Curve defined by y^2 = x^3 - x over Rational Field
        1728
        '32a2'

        sage: E = EllipticCurve(j=GF(5)(2)); E; E.j_invariant()
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 5
        2

    See :trac:`6657` ::

        sage: EllipticCurve(GF(144169),j=1728)
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

       sage: E = EllipticCurve_from_j(2^256+1,minimal_twist=False)
       sage: E.j_invariant() == 2^256+1
       True

    TESTS::

        sage: R = ZZ['u', 'v']
        sage: EllipticCurve(R, [1,1])
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Multivariate Polynomial Ring in u, v
        over Integer Ring

    We create a curve and a point over ``QQbar`` (see :trac:`6879`)::

        sage: E = EllipticCurve(QQbar,[0,1])
        sage: E(0)
        (0 : 1 : 0)
        sage: E.base_field()
        Algebraic Field

        sage: E = EllipticCurve(RR,[1,2]); E; E.base_field()
        Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 2.00000000000000 over Real Field with 53 bits of precision
        Real Field with 53 bits of precision
        sage: EllipticCurve(CC,[3,4]); E; E.base_field()
        Elliptic Curve defined by y^2 = x^3 + 3.00000000000000*x + 4.00000000000000 over Complex Field with 53 bits of precision
        Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 2.00000000000000 over Real Field with 53 bits of precision
        Real Field with 53 bits of precision
        sage: E = EllipticCurve(QQbar,[5,6]); E; E.base_field()
        Elliptic Curve defined by y^2 = x^3 + 5*x + 6 over Algebraic Field
        Algebraic Field

    See :trac:`6657` ::

        sage: EllipticCurve(3,j=1728)
        Traceback (most recent call last):
        ...
        ValueError: First parameter (if present) must be a ring when j is specified

        sage: EllipticCurve(GF(5),j=3/5)
        Traceback (most recent call last):
        ...
        ValueError: First parameter must be a ring containing 3/5

    If the universe of the coefficients is a general field, the object
    constructed has type EllipticCurve_field.  Otherwise it is
    EllipticCurve_generic.  See :trac:`9816` ::

        sage: E = EllipticCurve([QQbar(1),3]); E
        Elliptic Curve defined by y^2 = x^3 + x + 3 over Algebraic Field
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>

        sage: E = EllipticCurve([RR(1),3]); E
        Elliptic Curve defined by y^2 = x^3 + 1.00000000000000*x + 3.00000000000000 over Real Field with 53 bits of precision
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>

        sage: E = EllipticCurve([i,i]); E
        Elliptic Curve defined by y^2 = x^3 + I*x + I over Symbolic Ring
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>
        sage: E.category()
        Category of schemes over Symbolic Ring
        sage: SR in Fields()
        True

        sage: F = FractionField(PolynomialRing(QQ,'t'))
        sage: t = F.gen()
        sage: E = EllipticCurve([t,0]); E
        Elliptic Curve defined by y^2 = x^3 + t*x over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: type(E)
        <class 'sage.schemes.elliptic_curves.ell_field.EllipticCurve_field_with_category'>
        sage: E.category()
        Category of schemes over Fraction Field of Univariate Polynomial Ring in t over Rational Field

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
        """
        Return a ``UniqueFactory`` key and possibly extra parameters.

        INPUT:

        See the documentation for :class:`EllipticCurveFactory`.

        OUTPUT:

        A pair ``(key, extra_args)``:

        - ``key`` has the form `(R, (a_1, a_2, a_3, a_4, a_6))`,
          representing a ring and the Weierstrass coefficients of an
          elliptic curve over that ring;

        - ``extra_args`` is a dictionary containing additional data to
          be inserted into the elliptic curve structure.

        EXAMPLES::

            sage: EllipticCurve.create_key_and_extra_args(j=8000)
            ((Rational Field, (0, -1, 0, -3, -1)), {})

        When constructing a curve over `\\QQ` from a Cremona or LMFDB
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

        if is_SymbolicEquation(x):
            x = x.lhs() - x.rhs()

        if parent(x) is SR:
            x = x._polynomial_(rings.QQ['x', 'y'])

        if is_MPolynomial(x):
            if y is None:
                x = coefficients_from_Weierstrass_polynomial(x)
            else:
                x = coefficients_from_cubic(x, y, morphism=False)

        if isinstance(x, six.string_types):
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
            if R in (rings.ZZ, int, long):
                R = rings.QQ

        return (R, tuple(R(a) for a in x)), kwds

    def create_object(self, version, key, **kwds):
        """
        Create an object from a ``UniqueFactory`` key.

        EXAMPLES::

            sage: E = EllipticCurve.create_object(0, (GF(3), (1, 2, 0, 1, 2)))
            sage: type(E)
            <class 'sage.schemes.elliptic_curves.ell_finite_field.EllipticCurve_finite_field_with_category'>

        .. NOTE::

            Keyword arguments are currently only passed to the
            constructor for elliptic curves over `\\QQ`; elliptic
            curves over other fields do not support them.

        """
        R, x = key

        if R is rings.QQ:
            from .ell_rational_field import EllipticCurve_rational_field
            return EllipticCurve_rational_field(x, **kwds)
        elif is_NumberField(R):
            from .ell_number_field import EllipticCurve_number_field
            return EllipticCurve_number_field(R, x)
        elif rings.is_pAdicField(R):
            from .ell_padic_field import EllipticCurve_padic_field
            return EllipticCurve_padic_field(R, x)
        elif is_FiniteField(R) or (is_IntegerModRing(R) and R.characteristic().is_prime()):
            from .ell_finite_field import EllipticCurve_finite_field
            return EllipticCurve_finite_field(R, x)
        elif R in _Fields:
            from .ell_field import EllipticCurve_field
            return EllipticCurve_field(R, x)
        from .ell_generic import EllipticCurve_generic
        return EllipticCurve_generic(R, x)


EllipticCurve = EllipticCurveFactory('sage.schemes.elliptic_curves.constructor.EllipticCurve')


def EllipticCurve_from_Weierstrass_polynomial(f):
    """
    Return the elliptic curve defined by a cubic in (long) Weierstrass
    form.

    INPUT:

    - ``f`` -- a inhomogeneous cubic polynomial in long Weierstrass
      form.

    OUTPUT:

    The elliptic curve defined by it.

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
    """
    Return the coefficients `(a_1, a_2, a_3, a_4, a_5)` for a cubic in
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
    """
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
        K = rings.RationalField()
    if K not in _Fields:
        K = K.fraction_field()
    return EllipticCurve([-K(c4)/K(48), -K(c6)/K(864)])


def EllipticCurve_from_j(j, minimal_twist=True):
    """
    Return an elliptic curve with given `j`-invariant.

    INPUT:

    - ``j`` -- an element of some field.

    - ``minimal_twist`` (boolean, default True) -- If True and ``j`` is in `\QQ`, the curve returned is a
      minimal twist, i.e. has minimal conductor.  If `j` is not in `\QQ` this parameter is ignored.

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

       sage: E = EllipticCurve_from_j(2^256+1,minimal_twist=False)
       sage: E.j_invariant() == 2^256+1
       True

    """
    return EllipticCurve(coefficients_from_j(j, minimal_twist))

def coefficients_from_j(j, minimal_twist=True):
    """
    Return Weierstrass coefficients `(a_1, a_2, a_3, a_4, a_6)` for an
    elliptic curve with given `j`-invariant.

    INPUT:

    See :func:`EllipticCurve_from_j`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import coefficients_from_j
        sage: coefficients_from_j(0)
        [0, 0, 1, 0, 0]
        sage: coefficients_from_j(1728)
        [0, 0, 0, -1, 0]
        sage: coefficients_from_j(1)
        [1, 0, 0, 36, 3455]

    The ``minimal_twist`` parameter (ignored except over `\\QQ` and
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
        K = rings.RationalField()
    if K not in _Fields:
        K = K.fraction_field()

    char=K.characteristic()
    if char==2:
        if j == 0:
            return Sequence([0, 0, 1, 0, 0], universe=K)
        else:
            return Sequence([1, 0, 0, 0, 1/j], universe=K)
    if char == 3:
        if j==0:
            return Sequence([0, 0, 0, 1, 0], universe=K)
        else:
            return Sequence([0, j, 0, 0, -j**2], universe=K)

    if K is rings.RationalField():
        # we construct the minimal twist, i.e. the curve with minimal
        # conductor with this j_invariant:
        if j == 0:
            return Sequence([0, 0, 1, 0, 0], universe=K) # 27a3
        if j == 1728:
            return Sequence([0, 0, 0, -1, 0], universe=K) # 32a2

        if not minimal_twist:
            k=j-1728
            return Sequence([0, 0, 0, -3*j*k, -2*j*k**2], universe=K)

        n = j.numerator()
        m = n-1728*j.denominator()
        a4 = -3*n*m
        a6 = -2*n*m**2

        # Now E=[0,0,0,a4,a6] has j-invariant j=n/d
        from sage.sets.set import Set
        for p in Set(n.prime_divisors()+m.prime_divisors()):
            e = min(a4.valuation(p)//2,a6.valuation(p)//3)
            if e>0:
                p  = p**e
                a4 /= p**2
                a6 /= p**3

        # Now E=[0,0,0,a4,a6] is minimal at all p != 2,3
        tw = [-1,2,-2,3,-3,6,-6]
        E1 = EllipticCurve([0,0,0,a4,a6])
        Elist = [E1] + [E1.quadratic_twist(t) for t in tw]
        Elist.sort(key=lambda E: E.conductor())
        return Sequence(Elist[0].ainvs())

    # defaults for all other fields:
    if j == 0:
        return Sequence([0, 0, 0, 0, 1], universe=K)
    if j == 1728:
        return Sequence([0, 0, 0, 1, 0], universe=K)
    k=j-1728
    return Sequence([0, 0, 0, -3*j*k, -2*j*k**2], universe=K)


def EllipticCurve_from_cubic(F, P, morphism=True):
    r"""
    Construct an elliptic curve from a ternary cubic with a rational point.

    If you just want the Weierstrass form and are not interested in
    the morphism then it is easier to use
    :func:`~sage.schemes.elliptic_curves.jacobian.Jacobian`
    instead. This will construct the same elliptic curve but you don't
    have to supply the point ``P``.

    INPUT:

    - ``F`` -- a homogeneous cubic in three variables with rational
      coefficients, as a polynomial ring element, defining a smooth
      plane cubic curve.

    - ``P`` -- a 3-tuple `(x,y,z)` defining a projective point on the
      curve `F=0`. Need not be a flex, but see caveat on output.

    - ``morphism`` -- boolean (default: ``True``). Whether to return
      the morphism or just the elliptic curve.

    OUTPUT:

    An elliptic curve in long Weierstrass form isomorphic to the curve
    `F=0`.

    If ``morphism=True`` is passed, then a birational equivalence
    between F and the Weierstrass curve is returned. If the point
    happens to be a flex, then this is an isomorphism.

    EXAMPLES:

    First we find that the Fermat cubic is isomorphic to the curve
    with Cremona label 27a1::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^3+y^3+z^3
        sage: P = [1,-1,0]
        sage: E = EllipticCurve_from_cubic(cubic, P, morphism=False); E
        Elliptic Curve defined by y^2 + 2*x*y + 1/3*y = x^3 - x^2 - 1/3*x - 1/27 over Rational Field
        sage: E.cremona_label()
        '27a1'
        sage: EllipticCurve_from_cubic(cubic, [0,1,-1], morphism=False).cremona_label()
        '27a1'
        sage: EllipticCurve_from_cubic(cubic, [1,0,-1], morphism=False).cremona_label()
        '27a1'

    Next we find the minimal model and conductor of the Jacobian of the
    Selmer curve::

        sage: R.<a,b,c> = QQ[]
        sage: cubic = a^3+b^3+60*c^3
        sage: P = [1,-1,0]
        sage: E = EllipticCurve_from_cubic(cubic, P, morphism=False);  E
        Elliptic Curve defined by y^2 + 2*x*y + 20*y = x^3 - x^2 - 20*x - 400/3 over Rational Field
        sage: E.minimal_model()
        Elliptic Curve defined by y^2 = x^3 - 24300 over Rational Field
        sage: E.conductor()
        24300

    We can also get the birational equivalence to and from the
    Weierstrass form. We start with an example where ``P`` is a flex
    and the equivalence is an isomorphism::

        sage: f = EllipticCurve_from_cubic(cubic, P, morphism=True)
        sage: f
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
                a^3 + b^3 + 60*c^3
          To:   Elliptic Curve defined by y^2 + 2*x*y + 20*y = x^3 - x^2 - 20*x - 400/3
                over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (-c : -b + c : 1/20*a + 1/20*b)

        sage: finv = f.inverse();  finv
        Scheme morphism:
          From: Elliptic Curve defined by y^2 + 2*x*y + 20*y = x^3 - x^2 - 20*x - 400/3
                over Rational Field
          To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
          a^3 + b^3 + 60*c^3
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
        -x^3 + x^2*z + 2*x*y*z + y^2*z + 20*x*z^2 + 20*y*z^2 + 400/3*z^3

        sage: E.defining_polynomial()(f.defining_polynomials()) * f.post_rescaling()
        a^3 + b^3 + 60*c^3

    If the point is not a flex then the cubic can not be transformed
    to a Weierstrass equation by a linear transformation. The general
    birational transformation is quadratic::

        sage: cubic =  a^3+7*b^3+64*c^3
        sage: P = [2,2,-1]
        sage: f = EllipticCurve_from_cubic(cubic, P, morphism=True)
        sage: E = f.codomain();  E
        Elliptic Curve defined by y^2 - 722*x*y - 21870000*y = x^3
        + 23579*x^2 over Rational Field
        sage: E.minimal_model()
        Elliptic Curve defined by y^2 + y = x^3 - 331 over Rational Field

        sage: f
        Scheme morphism:
          From: Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
                a^3 + 7*b^3 + 64*c^3
          To:   Elliptic Curve defined by y^2 - 722*x*y - 21870000*y =
                x^3 + 23579*x^2 over Rational Field
          Defn: Defined on coordinates by sending (a : b : c) to
                (-5/112896*a^2 - 17/40320*a*b - 1/1280*b^2 - 29/35280*a*c
                 - 13/5040*b*c - 4/2205*c^2 :
                 -4055/112896*a^2 - 4787/40320*a*b - 91/1280*b^2 - 7769/35280*a*c
                 - 1993/5040*b*c - 724/2205*c^2 :
                 1/4572288000*a^2 + 1/326592000*a*b + 1/93312000*b^2 + 1/142884000*a*c
                 + 1/20412000*b*c + 1/17860500*c^2)

        sage: finv = f.inverse();  finv
        Scheme morphism:
          From: Elliptic Curve defined by y^2 - 722*x*y - 21870000*y =
                x^3 + 23579*x^2 over Rational Field
          To:   Closed subscheme of Projective Space of dimension 2 over Rational Field defined by:
                a^3 + 7*b^3 + 64*c^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (2*x^2 + 227700*x*z - 900*y*z :
                 2*x^2 - 32940*x*z + 540*y*z :
                 -x^2 - 56520*x*z - 180*y*z)

        sage: cubic(finv.defining_polynomials()) * finv.post_rescaling()
        -x^3 - 23579*x^2*z - 722*x*y*z + y^2*z - 21870000*y*z^2

        sage: E.defining_polynomial()(f.defining_polynomials()) * f.post_rescaling()
        a^3 + 7*b^3 + 64*c^3

    TESTS::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^2*y + 4*x*y^2 + x^2*z + 8*x*y*z + 4*y^2*z + 9*x*z^2 + 9*y*z^2
        sage: EllipticCurve_from_cubic(cubic, [1,-1,1], morphism=False)
        Elliptic Curve defined by y^2 - 882*x*y - 2560000*y = x^3 - 127281*x^2 over Rational Field
    """
    import sage.matrix.all as matrix

    # check the input
    R = F.parent()
    if not is_MPolynomialRing(R):
        raise TypeError('equation must be a polynomial')
    if R.ngens() != 3:
        raise TypeError('equation must be a polynomial in three variables')
    if not F.is_homogeneous():
        raise TypeError('equation must be a homogeneous polynomial')
    K = F.parent().base_ring()
    try:
        P = [K(c) for c in P]
    except TypeError:
        raise TypeError('cannot convert %s into %s'%(P,K))
    if F(P) != 0:
        raise ValueError('%s is not a point on %s'%(P,F))
    if len(P) != 3:
        raise TypeError('%s is not a projective point'%P)
    x, y, z = R.gens()

    # First case: if P = P2 then P is a flex
    P2 = chord_and_tangent(F, P)
    if are_projectively_equivalent(P, P2, base_ring=K):
        # find the tangent to F in P
        dx = K(F.derivative(x)(P))
        dy = K(F.derivative(y)(P))
        dz = K(F.derivative(z)(P))
        # find a second point Q on the tangent line but not on the cubic
        for tangent in [[dy, -dx, K.zero()], [dz, K.zero(), -dx], [K.zero(), -dz, dx]]:
            tangent = projective_point(tangent)
            Q = [tangent[0]+P[0], tangent[1]+P[1], tangent[2]+P[2]]
            F_Q = F(Q)
            if F_Q != 0:  # At most one further point may accidentally be on the cubic
                break
        assert F_Q != 0
        # pick linearly independent third point
        for third_point in [(1,0,0), (0,1,0), (0,0,1)]:
            M = matrix.matrix(K, [Q, P, third_point]).transpose()
            if M.is_invertible():
                break
        F2 = R(M.act_on_polynomial(F))
        # scale and dehomogenise
        a = K(F2.coefficient(x**3))
        F3 = F2/a
        b = K(F3.coefficient(y*y*z))
        S = rings.PolynomialRing(K, 'x,y,z')
        # elliptic curve coordinates
        X, Y, Z = S.gen(0), S.gen(1), S(-1/b)*S.gen(2)
        F4 = F3(X, Y, Z)
        E = EllipticCurve(F4.subs(z=1))
        if not morphism:
            return E
        inv_defining_poly = [ M[i,0]*X + M[i,1]*Y + M[i,2]*Z for i in range(3) ]
        inv_post = -1/a
        M = M.inverse()
        trans_x, trans_y, trans_z = [ M[i,0]*x + M[i,1]*y + M[i,2]*z for i in range(3) ]
        fwd_defining_poly = [trans_x, trans_y, -b*trans_z]
        fwd_post = -a

    # Second case: P is not a flex, then P, P2, P3 are different
    else:
        P3 = chord_and_tangent(F, P2)
        # send P, P2, P3 to (1:0:0), (0:1:0), (0:0:1) respectively
        M = matrix.matrix(K, [P, P2, P3]).transpose()
        F2 = M.act_on_polynomial(F)
        # substitute x = U^2, y = V*W, z = U*W, and rename (x,y,z)=(U,V,W)
        F3 = F2.substitute({x:x**2, y:y*z, z:x*z}) // (x**2*z)
        # scale and dehomogenise
        a = K(F3.coefficient(x**3))
        F4 = F3/a
        b = K(F4.coefficient(y*y*z))
        # change to a polynomial in only two variables
        S = rings.PolynomialRing(K, 'x,y,z')
        # elliptic curve coordinates
        X, Y, Z = S.gen(0), S.gen(1), S(-1/b)*S.gen(2)
        F5 = F4(X, Y, Z)
        E = EllipticCurve(F5.subs(z=1))
        if not morphism:
            return E
        inv_defining_poly = [ M[i,0]*X*X + M[i,1]*Y*Z + M[i,2]*X*Z for i in range(3) ]
        inv_post = -1/a/(X**2)/Z
        M = M.inverse()
        trans_x, trans_y, trans_z = [
            (M[i,0]*x + M[i,1]*y + M[i,2]*z) for i in range(3) ]
        fwd_defining_poly = [ trans_x*trans_z, trans_x*trans_y, -b*trans_z*trans_z ]
        fwd_post = -a/(trans_x*trans_z*trans_z)

    # Construct the morphism
    from sage.schemes.projective.projective_space import ProjectiveSpace
    P2 = ProjectiveSpace(2, K, names=[str(_) for _ in R.gens()])
    cubic = P2.subscheme(F)
    from sage.schemes.elliptic_curves.weierstrass_transform import \
        WeierstrassTransformationWithInverse
    return WeierstrassTransformationWithInverse(
        cubic, E, fwd_defining_poly, fwd_post, inv_defining_poly, inv_post)


def chord_and_tangent(F, P):
    """
    Use the chord and tangent method to get another point on a cubic.

    INPUT:

    - ``F`` -- a homogeneous cubic in three variables with rational
      coefficients, as a polynomial ring element, defining a smooth
      plane cubic curve.

    - ``P`` -- a 3-tuple `(x,y,z)` defining a projective point on the
      curve `F=0`.

    OUTPUT:

    Another point satisfying the equation ``F``.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: from sage.schemes.elliptic_curves.constructor import chord_and_tangent
        sage: F = x^3+y^3+60*z^3
        sage: chord_and_tangent(F, [1,-1,0])
        [1, -1, 0]

        sage: F = x^3+7*y^3+64*z^3
        sage: p0 = [2,2,-1]
        sage: p1 = chord_and_tangent(F, p0);  p1
        [-5, 3, -1]
        sage: p2 = chord_and_tangent(F, p1);  p2
        [1265, -183, -314]

    TESTS::

        sage: F(p2)
        0
        sage: map(type, p2)
        [<type 'sage.rings.rational.Rational'>,
         <type 'sage.rings.rational.Rational'>,
         <type 'sage.rings.rational.Rational'>]

    See :trac:`16068`::

        sage: F = x**3 - 4*x**2*y - 65*x*y**2 + 3*x*y*z - 76*y*z**2
        sage: chord_and_tangent(F, [0, 1, 0])
        [0, 0, -1]
    """
    # check the input
    R = F.parent()
    if not is_MPolynomialRing(R):
        raise TypeError('equation must be a polynomial')
    if R.ngens() != 3:
        raise TypeError('%s is not a polynomial in three variables'%F)
    if not F.is_homogeneous():
        raise TypeError('%s is not a homogeneous polynomial'%F)
    x, y, z = R.gens()
    if len(P) != 3:
        raise TypeError('%s is not a projective point'%P)
    K = R.base_ring()
    try:
        P = [K(c) for c in P]
    except TypeError:
        raise TypeError('cannot coerce %s into %s'%(P,K))
    if F(P) != 0:
        raise ValueError('%s is not a point on %s'%(P,F))

    # find the tangent to F in P
    dx = K(F.derivative(x)(P))
    dy = K(F.derivative(y)(P))
    dz = K(F.derivative(z)(P))
    # if dF/dy(P) = 0, change variables so that dF/dy != 0
    if dy == 0:
        if dx != 0:
            g = F.substitute({x:y, y:x})
            Q = [P[1], P[0], P[2]]
            R = chord_and_tangent(g, Q)
            return [R[1], R[0], R[2]]
        elif dz != 0:
            g = F.substitute({y:z, z:y})
            Q = [P[0], P[2], P[1]]
            R = chord_and_tangent(g, Q)
            return [R[0], R[2], R[1]]
        else:
            raise ValueError('%s is singular at %s'%(F, P))

    # t will be our choice of parmeter of the tangent plane
    #     dx*(x-P[0]) + dy*(y-P[1]) + dz*(z-P[2])
    # through the point P
    t = rings.PolynomialRing(K, 't').gen(0)
    Ft = F(dy*t+P[0], -dx*t+P[1], P[2])
    if Ft == 0:   # (dy, -dx, 0) is projectively equivalent to P
        # then (0, -dz, dy) is not projectively equivalent to P
        g = F.substitute({x:z, z:x})
        Q = [P[2], P[1], P[0]]
        R = chord_and_tangent(g, Q)
        return [R[2], R[1], R[0]]
    # Ft has a double zero at t=0 by construction, which we now remove
    Ft = Ft // t**2

    # first case: the third point is at t=infinity
    if Ft.is_constant():
        return projective_point([dy, -dx, K(0)])
    # second case: the third point is at finite t
    else:
        assert Ft.degree() == 1
        t0 = Ft.roots()[0][0]
        return projective_point([dy*t0+P[0], -dx*t0+P[1], P[2]])


def projective_point(p):
    """
    Return equivalent point with denominators removed

    INPUT:

    - ``P``, ``Q`` -- list/tuple of projective coordinates.

    OUTPUT:

    List of projective coordinates.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import projective_point
        sage: projective_point([4/5, 6/5, 8/5])
        [2, 3, 4]
        sage: F = GF(11)
        sage: projective_point([F(4), F(8), F(2)])
        [4, 8, 2]
    """
    from sage.rings.integer import GCD_list, LCM_list
    try:
        p_gcd = GCD_list([x.numerator() for x in p])
        p_lcm = LCM_list([x.denominator() for x in p])
    except AttributeError:
        return p
    scale = p_lcm / p_gcd
    return [scale * x for x in p]


def are_projectively_equivalent(P, Q, base_ring):
    """
    Test whether ``P`` and ``Q`` are projectively equivalent.

    INPUT:

    - ``P``, ``Q`` -- list/tuple of projective coordinates.

    - ``base_ring`` -- the base ring.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.constructor import are_projectively_equivalent
        sage: are_projectively_equivalent([0,1,2,3], [0,1,2,2], base_ring=QQ)
        False
        sage: are_projectively_equivalent([0,1,2,3], [0,2,4,6], base_ring=QQ)
        True
    """
    from sage.matrix.constructor import matrix
    return matrix(base_ring, [P, Q]).rank() < 2


def EllipticCurve_from_plane_curve(C, P):
    """
    Deprecated way to construct an elliptic curve.

    Use :meth:`~sage.schemes.elliptic_curves.jacobian.Jacobian` instead.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: C = Curve(x^3+y^3+z^3)
        sage: P = C(1,-1,0)
        sage: E = EllipticCurve_from_plane_curve(C,P); E  # long time (3s on sage.math, 2013)
        doctest:...: DeprecationWarning: use Jacobian(C) instead
        See http://trac.sagemath.org/3416 for details.
        Elliptic Curve defined by y^2 = x^3 - 27/4 over Rational Field
    """
    from sage.misc.superseded import deprecation
    deprecation(3416, 'use Jacobian(C) instead')
    # Note: this function never used the rational point
    from sage.schemes.elliptic_curves.jacobian import Jacobian
    return Jacobian(C)


def EllipticCurves_with_good_reduction_outside_S(S=[], proof=None, verbose=False):
    r"""
    Returns a sorted list of all elliptic curves defined over `Q`
    with good reduction outside the set `S` of primes.

    INPUT:

    -  ``S`` - list of primes (default: empty list).

    - ``proof`` - True/False (default True): the MW basis for
      auxiliary curves will be computed with this proof flag.

    - ``verbose`` - True/False (default False): if True, some details
      of the computation will be output.

    .. note::

        Proof flag: The algorithm used requires determining all
        S-integral points on several auxiliary curves, which in turn
        requires the computation of their generators.  This is not
        always possible (even in theory) using current knowledge.

        The value of this flag is passed to the function which
        computes generators of various auxiliary elliptic curves, in
        order to find their S-integral points.  Set to False if the
        default (True) causes warning messages, but note that you can
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
        sage: ', '.join([e.label() for e in elist])
        '32a1, 32a2, 32a3, 32a4, 64a1, 64a2, 64a3, 64a4, 128a1, 128a2, 128b1, 128b2, 128c1, 128c2, 128d1, 128d2, 256a1, 256a2, 256b1, 256b2, 256c1, 256c2, 256d1, 256d2'

    Without ``Proof=False``, this example gives two warnings::

        sage: elist = EllipticCurves_with_good_reduction_outside_S([11],proof=False)  # long time (14s on sage.math, 2011)
        sage: len(elist)  # long time
        12
        sage: ', '.join([e.label() for e in elist])  # long time
        '11a1, 11a2, 11a3, 121a1, 121a2, 121b1, 121b2, 121c1, 121c2, 121d1, 121d2, 121d3'

        sage: elist = EllipticCurves_with_good_reduction_outside_S([2,3]) # long time (26s on sage.math, 2011)
        sage: len(elist) # long time
        752
        sage: max([e.conductor() for e in elist]) # long time
        62208
        sage: [N.factor() for N in Set([e.conductor() for e in elist])] # long time
        [2^7,
        2^8,
        2^3 * 3^4,
        2^2 * 3^3,
        2^8 * 3^4,
        2^4 * 3^4,
        2^3 * 3,
        2^7 * 3,
        2^3 * 3^5,
        3^3,
        2^8 * 3,
        2^5 * 3^4,
        2^4 * 3,
        2 * 3^4,
        2^2 * 3^2,
        2^6 * 3^4,
        2^6,
        2^7 * 3^2,
        2^4 * 3^5,
        2^4 * 3^3,
        2 * 3^3,
        2^6 * 3^3,
        2^6 * 3,
        2^5,
        2^2 * 3^4,
        2^3 * 3^2,
        2^5 * 3,
        2^7 * 3^4,
        2^2 * 3^5,
        2^8 * 3^2,
        2^5 * 3^2,
        2^7 * 3^5,
        2^8 * 3^5,
        2^3 * 3^3,
        2^8 * 3^3,
        2^5 * 3^5,
        2^4 * 3^2,
        2 * 3^5,
        2^5 * 3^3,
        2^6 * 3^5,
        2^7 * 3^3,
        3^5,
        2^6 * 3^2]
    """
    from .ell_egros import (egros_from_jlist, egros_get_j)
    return egros_from_jlist(egros_get_j(S, proof=proof, verbose=verbose), S)
