# sage_setup: distribution = sagemath-schemes
r"""
Projective plane conics over a rational function field

The class :class:`ProjectiveConic_rational_function_field` represents a
projective plane conic over a rational function field `F(t)`, where `F`
is any field. Instances can be created using :func:`Conic`.

AUTHORS:

- Lennart Ackermans (2016-02-07): initial version

EXAMPLES:

Create a conic::

    sage: K = FractionField(PolynomialRing(QQ, 't'))
    sage: P.<X, Y, Z> = K[]
    sage: Conic(X^2 + Y^2 - Z^2)
    Projective Conic Curve over Fraction Field of Univariate
    Polynomial Ring in t over Rational Field defined by
    X^2 + Y^2 - Z^2

Points can be found using :meth:`has_rational_point`::

    sage: K.<t> = FractionField(QQ['t'])
    sage: C = Conic([1, -t, t])
    sage: C.has_rational_point(point=True)                                              # needs sage.libs.singular
    (True, (0 : 1 : 1))
"""

# ****************************************************************************
#       Copyright (C) 2016 Lennart Ackermans
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import diagonal_matrix, matrix, block_matrix
from sage.schemes.plane_conics.con_field import ProjectiveConic_field
from sage.arith.functions import lcm
from sage.arith.misc import GCD as gcd
from sage.modules.free_module_element import vector
from sage.rings.fraction_field import FractionField_generic


class ProjectiveConic_rational_function_field(ProjectiveConic_field):
    r"""
    Create a projective plane conic curve over a rational function field
    `F(t)`, where `F` is any field.

    The algorithms used in this class come mostly from [HC2006]_.

    EXAMPLES::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: P.<X, Y, Z> = K[]
        sage: Conic(X^2 + Y^2 - Z^2)
        Projective Conic Curve over Fraction Field of Univariate
        Polynomial Ring in t over Rational Field defined by
        X^2 + Y^2 - Z^2

    TESTS::

        sage: K = FractionField(PolynomialRing(QQ, 't'))
        sage: Conic([K(1), 1, -1])._test_pickling()

    REFERENCES:

    - [HC2006]_
    - [Ack2016]_
    """
    def __init__(self, A, f):
        r"""
        See ``Conic`` for full documentation.

        EXAMPLES::

            sage: c = Conic([1, 1, 1]); c
            Projective Conic Curve over Rational Field defined by
            x^2 + y^2 + z^2
        """
        ProjectiveConic_field.__init__(self, A, f)

    def has_rational_point(self, point=False, algorithm='default',
                           read_cache=True):
        r"""
        Return ``True`` if and only if the conic ``self``
        has a point over its base field `F(t)`, which is a field of rational
        functions.

        If ``point`` is True, then returns a second output, which is
        a rational point if one exists.

        Points are cached whenever they are found. Cached information
        is used if and only if ``read_cache`` is True.

        The default algorithm does not (yet) work for all base fields `F`.
        In particular, sage is required to have:

        * an algorithm for finding the square root of elements in finite
          extensions of `F`;

        * a factorization and gcd algorithm for `F[t]`;

        * an algorithm for solving conics over `F`.

        ALGORITHM:

        The parameter ``algorithm`` specifies the algorithm
        to be used:

        * ``'default'`` -- use a native Sage implementation, based on the
          algorithm Conic in [HC2006]_.

        * ``'magma'`` (requires Magma to be installed) --
          delegates the task to the Magma computer algebra
          system.

        EXAMPLES:

        We can find points for function fields over (extensions of) `\QQ`
        and finite fields::

            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: C = Conic(K, [t^2 - 2, 2*t^3, -2*t^3 - 13*t^2 - 2*t + 18])
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (True, (-3 : (t + 1)/t : 1))

            sage: R.<t> = FiniteField(23)[]
            sage: C = Conic([2, t^2 + 1, t^2 + 5])
            sage: C.has_rational_point()                                                # needs sage.libs.singular
            True
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (True, (5*t : 8 : 1))

            sage: # needs sage.rings.number_field
            sage: F.<i> = QuadraticField(-1)
            sage: R.<t> = F[]
            sage: C = Conic([1, i*t, -t^2 + 4])
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (True, (-t - 2*i : -2*i : 1))

        It works on non-diagonal conics as well::

            sage: K.<t> = QQ[]
            sage: C = Conic([4, -4, 8, 1, -4, t + 4])
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (True, (1/2 : 1 : 0))

        If no point exists output still depends on the argument ``point``::

            sage: K.<t> = QQ[]
            sage: C = Conic(K, [t^2, (t-1), -2*(t-1)])
            sage: C.has_rational_point()                                                # needs sage.libs.singular
            False
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (False, None)

        Due to limitations in Sage of algorithms we depend on, it is not
        yet possible to find points on conics over multivariate function fields
        (see the requirements above)::

            sage: F.<t1> = FractionField(QQ['t1'])
            sage: K.<t2> = FractionField(F['t2'])
            sage: a = K(1)
            sage: b = 2*t2^2 + 2*t1*t2 - t1^2
            sage: c = -3*t2^4 - 4*t1*t2^3 + 8*t1^2*t2^2 + 16*t1^3 - t2 - 48*t1^4
            sage: C = Conic([a,b,c])
            sage: C.has_rational_point()                                                # needs sage.libs.singular
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square() not implemented for elements of
            Univariate Quotient Polynomial Ring in tbar over Fraction Field
            of Univariate Polynomial Ring in t1 over Rational Field with
            modulus tbar^2 + t1*tbar - 1/2*t1^2

        In some cases, the algorithm requires us to be
        able to solve conics over `F`. In particular, the following does not
        work::

            sage: P.<u> = QQ[]
            sage: E = P.fraction_field()
            sage: Q.<Y> = E[]
            sage: F.<v> = E.extension(Y^2 - u^3 - 1)
            sage: R.<t> = F[]
            sage: K = R.fraction_field()                                                # needs sage.rings.function_field
            sage: C = Conic(K, [u, v, 1])                                               # needs sage.rings.function_field
            sage: C.has_rational_point()                                                # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: has_rational_point not implemented for conics
            over base field Univariate Quotient Polynomial Ring in v over
            Fraction Field of Univariate Polynomial Ring in u over Rational
            Field with modulus v^2 - u^3 - 1

        TESTS::

            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: a = (2*t^2 - 3/2*t + 1)/(37/3*t^2 + t - 1/4)
            sage: b = (1/2*t^2 + 1/3)/(-73*t^2 - 2*t + 11/4)
            sage: c = (6934/3*t^6 + 8798/3*t^5 - 947/18*t^4 + 3949/9*t^3 + 20983/18*t^2 + 28/3*t - 131/3)/(-2701/3*t^4 - 293/3*t^3 + 301/6*t^2 + 13/4*t - 11/16)
            sage: C = Conic([a,b,c])
            sage: C.has_rational_point(point=True)                                      # needs sage.libs.singular
            (True, (4*t + 4 : 2*t + 2 : 1))

        A long time test::

            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: a = (-1/3*t^6 - 14*t^5 - 1/4*t^4 + 7/2*t^2 - 1/2*t - 1)/(24/5*t^6 - t^5 - 1/4*t^4 + t^3 - 3*t^2 + 8/5*t + 5)
            sage: b = (-3*t^3 + 8*t + 1/2)/(-1/3*t^3 + 3/2*t^2 + 1/12*t + 1/2)
            sage: c = (1232009/225*t^25 - 1015925057/8100*t^24 + 1035477411553/1458000*t^23 + 7901338091/30375*t^22 - 1421379260447/729000*t^21 + 266121260843/972000*t^20 + 80808723191/486000*t^19 - 516656082523/972000*t^18 + 21521589529/40500*t^17 + 4654758997/21600*t^16 - 20064038625227/9720000*t^15 - 173054270347/324000*t^14 + 536200870559/540000*t^13 - 12710739349/50625*t^12 - 197968226971/135000*t^11 - 134122025657/810000*t^10 + 22685316301/120000*t^9 - 2230847689/21600*t^8 - 70624099679/270000*t^7 - 4298763061/270000*t^6 - 41239/216000*t^5 - 13523/36000*t^4 + 493/36000*t^3 + 83/2400*t^2 + 1/300*t + 1/200)/(-27378/125*t^17 + 504387/500*t^16 - 97911/2000*t^15 + 1023531/4000*t^14 + 1874841/8000*t^13 + 865381/12000*t^12 + 15287/375*t^11 + 6039821/6000*t^10 + 599437/1500*t^9 + 18659/250*t^8 + 1218059/6000*t^7 + 2025127/3000*t^6 + 1222759/6000*t^5 + 38573/200*t^4 + 8323/125*t^3 + 15453/125*t^2 + 17031/500*t + 441/10)
            sage: C = Conic([a,b,c])
            sage: C.has_rational_point(point=True)      # long time (4 seconds)         # needs sage.libs.singular
            (True,
             ((-2/117*t^8 + 304/1053*t^7 + 40/117*t^6 - 1/27*t^5 - 110/351*t^4 - 2/195*t^3 + 11/351*t^2 + 1/117)/(t^4 + 2/39*t^3 + 4/117*t^2 + 2/39*t + 14/39) : -5/3*t^4 + 19*t^3 : 1))

        ``has_rational_point`` used to fail for some conics over function fields
        over finite fields, due to :issue:`20003`::

            sage: K.<t> = PolynomialRing(GF(7))
            sage: C = Conic([5*t^2 + 4, t^2 + 3*t + 3, 6*t^2 + 3*t + 2,
            ....:            5*t^2 + 5, 4*t + 3, 4*t^2 + t + 5])
            sage: C.has_rational_point()
            True
        """
        from .constructor import Conic

        if read_cache:
            if self._rational_point is not None:
                return (True, self._rational_point) if point else True

        if algorithm != 'default':
            return ProjectiveConic_field.has_rational_point(self, point,
                algorithm, read_cache)

        # Default algorithm
        if self.base_ring().characteristic() == 2:
            raise NotImplementedError("has_rational_point not implemented \
for function field of characteristic 2.")
        new_conic, transformation, inverse = self.diagonalization()
        coeff = new_conic.coefficients()
        if coeff[0] == 0:
            return (True, transformation([1, 0, 0])) if point else True
        elif coeff[3] == 0:
            return (True, transformation([0, 1, 0])) if point else True
        elif coeff[5] == 0:
            return (True, transformation([0, 0, 1])) if point else True

        # We save the coefficients of the reduced form in coeff
        # A zero of the reduced conic can be multiplied by multipliers
        # to get a zero of the old conic
        (coeff, multipliers) = new_conic._reduce_conic()
        new_conic = Conic(coeff)
        transformation = transformation \
            * new_conic.hom(diagonal_matrix(multipliers))
        if coeff[0].degree() % 2 == coeff[1].degree() % 2 and \
                coeff[1].degree() % 2 == coeff[2].degree() % 2:
            case = 0
        else:
            case = 1

        t, = self.base_ring().base().gens()  # t in F[t]
        supp = []
        roots = [[], [], []]
        remove = None
        # loop through the coefficients and find a root of f_i (as in
        # [HC2006]) modulo each element in the coefficients' support
        for i in (0, 1, 2):
            supp.append(list(coeff[i].factor()))
            for p in supp[i]:
                if p[1] != 1:
                    raise ValueError("Expected factor of exponent 1.")
                # Convert to monic factor
                x = p[0] / list(p[0])[-1]
                N = p[0].base_ring().extension(x, 'tbar')
                R = PolynomialRing(N, 'u')
                u, = R.gens()
                # If p[0] has degree 1, sage might forget the "defining
                # polynomial" of N, so we define our own modulo operation
                if p[0].degree() == 1:
                    mod = t.parent().hom([-x[0]])
                else:
                    mod = N
                if i == 0:
                    x = -mod(coeff[2]) / mod(coeff[1])
                elif i == 1:
                    x = -mod(coeff[0]) / mod(coeff[2])
                else:
                    x = -mod(coeff[1]) / mod(coeff[0])
                if x.is_square():
                    root = N(x.sqrt())
                else:
                    return (False, None) if point else False
                # if case == 0 and p[0] has degree 1, we switch to case
                # 1 and remove this factor out of the support. In [HC2006]
                # this is done later, in FindPoint.
                if case == 0 and p[0].degree() == 1:
                    case = 1
                    # remove later so the loop iterator stays in place.
                    remove = (i, p)
                else:
                    roots[i].append(root)
        if remove:
            supp[remove[0]].remove(remove[1])
        supp = [[p[0] for p in supp[i]] for i in (0, 1, 2)]

        if case == 0:
            # Find a solution of (5) in [HC2006]
            leading_conic = Conic(self.base_ring().base_ring(),
                        [coeff[0].leading_coefficient(),
                        coeff[1].leading_coefficient(),
                        coeff[2].leading_coefficient()])
            has_point = leading_conic.has_rational_point(True)
            if has_point[0]:
                if point:
                    pt = new_conic.find_point(supp, roots, case,
                        has_point[1])
                else:
                    pt = True
                return (True, transformation(pt)) if point else True
            else:
                return (False, None) if point else False
        # case == 1:
        if point:
            pt = new_conic.find_point(supp, roots, case)
        else:
            pt = True
        return (True, transformation(pt)) if point else True

    def _reduce_conic(self):
        r"""
        Return the reduced form of the conic, i.e. a conic with base field
        `K=F(t)` and coefficients `a,b,c` such that `a,b,c \in F[t]`,
        `\gcd(a,b)=\gcd(b,c)=\gcd(c,a)=1` and `abc` is square-free.

        Assumes ``self`` is in diagonal form.

        OUTPUT:

        A tuple (coefficients, multipliers), the coefficients of the conic
        in reduced form and multipliers `\lambda, \mu, \nu \in F(t)^*` such
        that `(x,y,z) \in F(t)` is a solution of the reduced conic if and only
        if `(\lambda x, \mu y, \nu z)` is a solution of ``self``.

        ALGORITHM:

        The algorithm used is the algorithm ReduceConic in [HC2006]_.

        EXAMPLES::

            sage: K.<t> = FractionField(PolynomialRing(QQ, 't'))
            sage: C = Conic(K, [t^2 - 2, 2*t^3, -2*t^3 - 13*t^2 - 2*t + 18])
            sage: C._reduce_conic()
            ([t^2 - 2, 2*t, -2*t^3 - 13*t^2 - 2*t + 18], [t, 1, t])
        """

        # start with removing fractions
        coeff = [self.coefficients()[0], self.coefficients()[3],
                self.coefficients()[5]]
        coeff = lcm(lcm(coeff[0].denominator(), coeff[1].denominator()),
                coeff[2].denominator()) * vector(coeff)
        # go to base ring of fraction field
        coeff = [self.base().base()(x) for x in coeff]
        coeff = vector(coeff) / gcd(coeff)
        # remove common divisors
        labda = mu = nu = 1
        g1 = g2 = g3 = 0
        ca, cb, cc = coeff
        while g1 != 1 or g2 != 1 or g3 != 1:
            g1 = gcd(ca, cb)
            ca = ca / g1
            cb = cb / g1
            cc = cc * g1
            nu = g1 * nu
            g2 = gcd(ca, cc)
            ca = ca / g2
            cc = cc / g2
            cb = cb * g2
            mu = g2 * mu
            g3 = gcd(cb, cc)
            cb = cb / g3
            cc = cc / g3
            ca = ca * g3
            labda = g3 * labda
        coeff = [ca, cb, cc]
        multipliers = [labda, mu, nu]

        # remove squares
        for i, x in enumerate(coeff):
            if isinstance(x.parent(), FractionField_generic):
                # go to base ring of fraction field
                x = self.base().base()(x)

            try:
                decom = x.squarefree_decomposition()
            except (NotImplementedError, AttributeError):
                decom = x.factor()
            x = decom.unit()
            x2 = 1
            for factor in decom:
                if factor[1] > 1:
                    if factor[1] % 2 == 0:
                        x2 *= factor[0] ** (factor[1] // 2)
                    else:
                        x *= factor[0]
                        x2 *= factor[0] ** ((factor[1] - 1) // 2)
                else:
                    x *= factor[0]
            for j, y in enumerate(multipliers):
                if j != i:
                    multipliers[j] = y * x2
            coeff[i] = self.base_ring().base().coerce(x)

        return (coeff, multipliers)

    def find_point(self, supports, roots, case, solution=0):
        r"""
        Given a solubility certificate like in [HC2006]_, find a point on
        ``self``. Assumes ``self`` is in reduced form (see [HC2006]_ for a
        definition).

        If you don't have a solubility certificate and just want to find a
        point, use the function :meth:`has_rational_point` instead.

        INPUT:

        - ``self`` -- conic in reduced form
        - ``supports`` -- 3-tuple where ``supports[i]`` is a list of all monic
          irreducible `p \in F[t]` that divide the `i`-th of the 3 coefficients
        - ``roots`` -- 3-tuple containing lists of roots of all elements of
          ``supports[i]``, in the same order
        - ``case`` -- 1 or 0, as in [HC2006]_
        - ``solution`` -- (default: 0) a solution of (5) in [HC2006]_, if
          ``case`` = 0, 0 otherwise

        OUTPUT:

        A point `(x,y,z) \in F(t)` of ``self``. Output is undefined when the
        input solubility certificate is incorrect.

        ALGORITHM:

        The algorithm used is the algorithm FindPoint in [HC2006]_, with
        a simplification from [Ack2016]_.

        EXAMPLES::

            sage: K.<t> = FractionField(QQ['t'])
            sage: C = Conic(K, [t^2 - 2, 2*t^3, -2*t^3 - 13*t^2 - 2*t + 18])
            sage: C.has_rational_point(point=True)  # indirect test                     # needs sage.libs.singular
            (True, (-3 : (t + 1)/t : 1))

        Different solubility certificates give different points::

            sage: # needs sage.rings.number_field
            sage: K.<t> = PolynomialRing(QQ, 't')
            sage: C = Conic(K, [t^2 - 2, 2*t, -2*t^3 - 13*t^2 - 2*t + 18])
            sage: supp = [[t^2 - 2], [t], [t^3 + 13/2*t^2 + t - 9]]
            sage: tbar1 = QQ.extension(supp[0][0], 'tbar').gens()[0]
            sage: tbar2 = QQ.extension(supp[1][0], 'tbar').gens()[0]
            sage: tbar3 = QQ.extension(supp[2][0], 'tbar').gens()[0]
            sage: roots = [[tbar1 + 1], [1/3*tbar2^0], [2/3*tbar3^2 + 11/3*tbar3 - 3]]
            sage: C.find_point(supp, roots, 1)
            (3 : t + 1 : 1)
            sage: roots = [[-tbar1 - 1], [-1/3*tbar2^0], [-2/3*tbar3^2 - 11/3*tbar3 + 3]]
            sage: C.find_point(supp, roots, 1)
            (3 : -t - 1 : 1)
        """
        Ft = self.base().base()
        F = Ft.base()
        t, = Ft.gens()
        coefficients = [Ft(self.coefficients()[0]), Ft(self.coefficients()[3]),
            Ft(self.coefficients()[5])]
        deg = [coefficients[0].degree(), coefficients[1].degree(),
               coefficients[2].degree()]
        # definitions as in [HC2006] and [Ack2016]
        A = ((deg[1] + deg[2]) / 2).ceil() - case
        B = ((deg[2] + deg[0]) / 2).ceil() - case
        C = ((deg[0] + deg[1]) / 2).ceil() - case

        # For all roots as calculated by has_rational_point(), we create
        # a system of linear equations. As in [Ack2016], we do this
        # by calculating the matrices for all phi_p, with basis consisting
        # of monomials of x, y and z in the space V of potential solutions:
        # t^0, ..., t^A, t^0, ..., t^B and t^0, ..., t^C.
        phi = []
        for (i, p) in enumerate(supports[0]):
            # lift to F[t] and map to R, with R as defined above
            if roots[0][i].parent().is_finite():
                root = roots[0][i].polynomial()
            else:
                root = roots[0][i].lift()
            alpha = root.parent().hom([t])(root)
            d = p.degree()
            # Calculate y - alpha*z mod p for all basis vectors
            phi_p = [[] for i in range(A + B + C + 4)]
            phi_p[0:A + 1] = [vector(F, d)] * (A + 1)
            phi_p[A + 1] = vector(F, d, {0: F(1)})
            lastpoly = F(1)
            for n in range(B):
                lastpoly = (lastpoly * t) % p
                phi_p[A + 2 + n] = vector(F, d, lastpoly.monomial_coefficients())
            lastpoly = -alpha % p
            phi_p[A + B + 2] = vector(F, d, lastpoly.monomial_coefficients())
            for n in range(C):
                lastpoly = (lastpoly * t) % p
                phi_p[A + B + 3 + n] = vector(F, d, lastpoly.monomial_coefficients())
            phi_p[A + B + C + 3] = vector(F, d)
            phi.append(matrix(phi_p).transpose())
        for (i, p) in enumerate(supports[1]):
            if roots[1][i].parent().is_finite():
                root = roots[1][i].polynomial()
            else:
                root = roots[1][i].lift()
            alpha = root.parent().hom([t])(root)
            d = p.degree()
            # Calculate z - alpha*x mod p for all basis vectors
            phi_p = [[] for i in range(A + B + C + 4)]
            phi_p[A + 1:A + B + 2] = [vector(F, d)] * (B + 1)
            phi_p[A + B + 2] = vector(F, d, {0: F(1)})
            lastpoly = F(1)
            for n in range(C):
                lastpoly = (lastpoly * t) % p
                phi_p[A + B + 3 + n] = vector(F, d, lastpoly.monomial_coefficients())
            lastpoly = -alpha % p
            phi_p[0] = vector(F, d, lastpoly.monomial_coefficients())
            for n in range(A):
                lastpoly = (lastpoly * t) % p
                phi_p[1 + n] = vector(F, d, lastpoly.monomial_coefficients())
            phi_p[A + B + C + 3] = vector(F, d)
            phi.append(matrix(phi_p).transpose())
        for (i, p) in enumerate(supports[2]):
            if roots[2][i].parent().is_finite():
                root = roots[2][i].polynomial()
            else:
                root = roots[2][i].lift()
            alpha = root.parent().hom([t])(root)
            d = p.degree()
            # Calculate x - alpha*y mod p for all basis vectors
            phi_p = [[] for i in range(A + B + C + 4)]
            phi_p[A + B + 2:A + B + C + 3] = [vector(F, d)] * (C + 1)
            phi_p[0] = vector(F, d, {0: F(1)})
            lastpoly = F(1)
            for n in range(A):
                lastpoly = (lastpoly * t) % p
                phi_p[1 + n] = vector(F, d, lastpoly.monomial_coefficients())
            lastpoly = -alpha % p
            phi_p[A + 1] = vector(F, d, lastpoly.monomial_coefficients())
            for n in range(B):
                lastpoly = (lastpoly * t) % p
                phi_p[A + 2 + n] = vector(F, d, lastpoly.monomial_coefficients())
            phi_p[A + B + C + 3] = vector(F, d)
            phi.append(matrix(phi_p).transpose())
        if case == 0:
            # We need three more equations
            lx = Ft(solution[0]).leading_coefficient()
            ly = Ft(solution[1]).leading_coefficient()
            lz = Ft(solution[2]).leading_coefficient()
            ABC = A + B + C
            phi.append(matrix([vector(F, ABC + 4, {A: 1, ABC + 3: -lx}),
                vector(F, ABC + 4, {A + B + 1: 1, ABC + 3: -ly}),
                vector(F, ABC + 4, {ABC + 2: 1, ABC + 3: -lz})]))
        # Create the final matrix which we will solve
        M = block_matrix(phi, ncols=1, subdivide=False)
        solution_space = M.right_kernel()
        for v in solution_space.basis():
            if v[:A + B + C + 3] != 0:
                # we do not want to return a trivial solution
                X = Ft(list(v[:A + 1]))
                Y = Ft(list(v[A + 1:A + B + 2]))
                Z = Ft(list(v[A + B + 2:A + B + C + 3]))
                return self.point([X, Y, Z])

        raise RuntimeError("No solution has been found: possibly incorrect "
                           "solubility certificate.")
